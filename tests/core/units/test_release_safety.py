"""Release blockers: reject unsafe semantics and verify numeric/backend contracts."""

from dataclasses import replace
from decimal import Decimal
import math
import subprocess
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import pint

import core.units as units
from core.units import agricultural_specs as agricultural
from core.units import conversions as conversions
from core.units import heat_content as heat
from core.units import registry


class ReleaseSafetyTests(unittest.TestCase):
    def test_carbon_and_certificates_are_not_physical_or_fungible(self):
        for source, target in (
            ("EUA", "UKA"), ("EUA", "CCA"), ("EUA", "tonne_CO2e"),
            ("VCU", "GS_VER"), ("REC", "MWh"), ("REC", "GO"),
            ("USD/EUA", "USD/UKA"), ("USD/REC", "USD/MWh"),
        ):
            with self.subTest(source=source, target=target), self.assertRaises(pint.DimensionalityError):
                units.Q_(1, source).to(target)
        self.assertEqual(units.Q_(1, "eua").to("EU_allowance").magnitude, 1)

    def test_direct_pint_mass_and_pressure_bypasses_fail(self):
        for source, target in (
            ("wmt", "dmt"), ("dmt", "dmtu"), ("wmt/day", "kg/hour"),
            ("USD/wmt", "USD/dmt"), ("admt", "bdmt"),
            ("psig", "psia"), ("Pa", "barg"), ("kpag", "kpaa"),
        ):
            with self.subTest(source=source, target=target):
                with self.assertRaises(pint.DimensionalityError):
                    units.Q_(1, source).to(target)
                with self.assertRaises(pint.DimensionalityError):
                    units.ureg.convert(1, source, target)
        self.assertEqual(units.convert(units.Q_(Decimal(100), "wmt"), "dmt", moisture_fraction=Decimal(".08")).magnitude, 92)
        self.assertAlmostEqual(float(units.convert(units.Q_(0, "psig"), "psia").magnitude), conversions.STANDARD_ATMOSPHERE_PSI)
        with self.assertRaises(pint.DimensionalityError):
            units.Q_(1, "wmt").to_base_units()

    def test_unknown_moisture_does_not_mean_zero(self):
        with self.assertRaisesRegex(ValueError, "moisture"):
            conversions.convert_bulk_ore(units.Q_(100, "wmt"), "dmt")
        with self.assertRaisesRegex(ValueError, "moisture"):
            conversions.BulkOreAssay().wmt_to_dmt(100)
        known_dry = conversions.BulkOreAssay(moisture_fraction=0)
        self.assertEqual(known_dry.wmt_to_dmt(Decimal(100)), 100)
        # Contained-metal extraction on an already dry basis does not need moisture.
        self.assertEqual(conversions.convert_bulk_ore(units.Q_(Decimal(1), "dmt"), "dmtu", grade=Decimal(".62")).magnitude, 62)

    def test_gas_thermal_dispatch_preserves_every_reference_state(self):
        gas = conversions.GasAssay(gcv_btu_per_scf=1000, relative_density=.6)
        for source in ("scf", "scf_tx", "scf_60F_14696psi_dry", "Sm3", "Nm3", "Sm3_20C", "m3_Groningen"):
            with self.subTest(source=source):
                original = units.Q_(Decimal("2.5"), source)
                expected = conversions.convert_gas_energy(original, "MJ", assay=gas, basis="gcv")
                actual = conversions.convert(original, "MJ", assay=gas, basis="gcv")
                restored = conversions.convert(actual, source, assay=gas, basis="gcv")
                self.assertAlmostEqual(float(actual.magnitude), float(expected.magnitude), places=8)
                self.assertAlmostEqual(float(restored.magnitude), 2.5, places=8)
                self.assertIsInstance(actual.magnitude, Decimal)

    def test_dak_must_establish_convergence(self):
        state = (units.Q_(5000, "kPa"), units.Q_(300, "K"))
        result = conversions.calculate_z_factor(*state, method="dak")
        self.assertTrue(math.isfinite(result) and result > 0)
        with self.assertRaises(ArithmeticError):
            conversions.calculate_z_factor(*state, method="dak", max_iter=1, tol=1e-16)
        for kwargs in ({"max_iter": 0}, {"max_iter": True}, {"tol": float("nan")}, {"tol": 0}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                conversions.calculate_z_factor(*state, method="dak", **kwargs)

    def test_unknown_components_fail_before_optional_backend_calls(self):
        fake = SimpleNamespace(Detail=lambda: self.fail("Backend must not be invoked"))
        with patch.dict(sys.modules, {"pyaga8": fake}):
            for composition in ({"CH4": .9, "UNKNOWN": .1}, {"CH4": 1, "UNKNOWN": 0}, {"CH4": float("nan")}):
                with self.subTest(composition=composition), self.assertRaises(ValueError):
                    conversions.calculate_z_factor(units.Q_(1, "MPa"), units.Q_(300, "K"), method="aga8", molar_composition=composition)

    def test_backend_identity_is_explicit_and_invalid_results_fail(self):
        backend = SimpleNamespace(PropsSI=lambda *args: .95)
        with patch.dict(sys.modules, {"CoolProp": SimpleNamespace(), "CoolProp.CoolProp": backend}):
            with patch.object(backend, "PropsSI", return_value=.95) as calculate:
                for method, prefix in (("heos", "HEOS::"), ("gerg2008", "GERG2008::")):
                    result = conversions.calculate_z_factor(units.Q_(1, "MPa"), units.Q_(300, "K"), method=method, molar_composition={"CH4": 1})
                    self.assertEqual(result, .95)
                    self.assertTrue(calculate.call_args.args[-1].startswith(prefix))
            with patch.object(backend, "PropsSI", return_value=float("nan")):
                with self.assertRaises(ArithmeticError):
                    conversions.calculate_z_factor(units.Q_(1, "MPa"), units.Q_(300, "K"), method="heos", molar_composition={"CH4": 1})
        with patch.dict(sys.modules, {"pyaga8": None}):
            with self.assertRaises(ImportError):
                conversions.calculate_z_factor(units.Q_(1, "MPa"), units.Q_(300, "K"), method="aga8", molar_composition={"CH4": 1})

    def test_assay_typeerrors_are_not_retried_or_suppressed(self):
        class BrokenAssay:
            calls = 0
            def as_context(self, *, basis="lhv"):
                self.calls += 1
                raise TypeError("internal assay defect")
        assay = BrokenAssay()
        with self.assertRaisesRegex(TypeError, "internal assay defect"):
            conversions.convert_commodity(units.Q_(1, "kg"), "m^3", assay=assay, basis="lhv")
        self.assertEqual(assay.calls, 1)

    def test_config_failure_rolls_back_units_symbols_metadata_and_freeze(self):
        quantity = units.Q_(Decimal("1.25"), "kg")
        payload = {
            "currencies": [{"code": "ROLLBACKCCY", "symbols": ["RB$"], "decimals": 7}],
            "commodities": [
                {"definition": "rollback_parcel = 123 * kilogram"},
                {"definition": "kilogram = second"},
            ],
        }
        with self.assertRaises((ValueError, pint.errors.RedefinitionError)):
            units.load_desk_config(payload)
        self.assertNotIn("ROLLBACKCCY", units.ureg)
        self.assertNotIn("rollback_parcel", units.ureg)
        self.assertNotIn("RB$", units.CURRENCY_SYMBOL_TO_CODE)
        self.assertNotIn("ROLLBACKCCY", units.CURRENCY_MINOR_DIGITS)
        self.assertFalse(units.is_registry_frozen())
        self.assertEqual(quantity.to("g").magnitude, 1250)

    def test_failed_individual_registration_also_rolls_back(self):
        with self.assertRaises(pint.errors.RedefinitionError):
            units.register_currency("PARTIALCCY", decimals=4, sub_units={"kilogram": "PARTIALCCY"})
        self.assertNotIn("PARTIALCCY", units.ureg)
        self.assertNotIn("PARTIALCCY", units.CURRENCY_MINOR_DIGITS)

    def test_config_readers_cannot_observe_uncommitted_definitions(self):
        from concurrent.futures import ThreadPoolExecutor
        from threading import Event
        from core.units import config
        written, release = Event(), Event()
        original_register = registry.register_commodity_unit
        def staged_register(definition, **kwargs):
            original_register(definition, **kwargs)
            written.set()
            if not release.wait(5):
                raise RuntimeError("Test writer timed out")
            raise ValueError("Simulated late configuration failure")
        def writer():
            with patch.object(config, "register_commodity_unit", staged_register):
                with self.assertRaisesRegex(ValueError, "Simulated"):
                    units.load_desk_config({"commodities": [{"definition": "uncommitted_parcel = 5 * kilogram"}]})
        def reader():
            with self.assertRaises(pint.errors.UndefinedUnitError):
                units.parse_unit("uncommitted_parcel")
        with ThreadPoolExecutor(max_workers=2) as pool:
            write = pool.submit(writer)
            self.assertTrue(written.wait(5))
            read = pool.submit(reader)
            release.set()
            write.result(timeout=10)
            read.result(timeout=10)

    def test_freezing_blocks_direct_pint_mutators(self):
        script = '''
import core.units as u
u.freeze_registry()
for action in (lambda: u.ureg.define("forbidden_direct_unit = meter"),
               lambda: u.ureg.load_definitions(["forbidden_loaded_unit = meter"]),
               lambda: setattr(u.ureg, "default_system", "cgs"),
               lambda: u.ureg.add_context(__import__("pint").Context("forbidden_context"))):
    try: action()
    except RuntimeError: pass
    else: raise AssertionError("Direct Pint mutation bypassed freeze")
assert u.parse_quantity("1 km").magnitude == 1
with u.sandbox_registry() as sandbox:
    sandbox.define("sandbox_mutable_unit = meter")
assert "sandbox_mutable_unit" not in u.ureg
'''
        result = subprocess.run([sys.executable, "-B", "-c", script], capture_output=True, text=True, timeout=30)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_scalar_and_batch_numeric_policies_match(self):
        for decimal in (True, False):
            for value in (2, 2.5, Decimal("2.50"), "2.50"):
                scalar = units.parse_quantity(value, "kg", as_decimal=decimal)
                batch = units.parse_quantities([value], "kg", as_decimal=decimal)[0]
                self.assertIsInstance(scalar.magnitude, Decimal if decimal else float)
                self.assertEqual(scalar, batch)
            for expression in ("2.50 kg", "kg", "2.50"):
                self.assertIsInstance(units.parse_quantity(expression, as_decimal=decimal).magnitude, Decimal if decimal else float)
        exact = "1.234567890123456789012345678"
        self.assertEqual(units.parse_quantity(exact).magnitude, Decimal(exact))
        for text in ("1,500.25", "(1,500.25)"):
            self.assertEqual(units.parse_quantity(text, "kg"), units.parse_quantities([text], "kg")[0])
        for invalid in (float("nan"), Decimal("Infinity"), True):
            with self.assertRaises((ValueError, TypeError)):
                units.parse_quantity(invalid, "kg")

    def test_finite_properties_and_mass_balance_fields(self):
        for assay, field in (
            (conversions.CrudeAssay.brent(), "sulfur_percentage"),
            (conversions.BiofuelAssay.ethanol(), "density_kg_per_m3"),
            (agricultural.GrainAssay.corn_56lb(), "actual_moisture_pct"),
        ):
            for invalid in (float("nan"), float("inf")):
                with self.subTest(assay=type(assay).__name__, field=field), self.assertRaises((ValueError, TypeError)):
                    replace(assay, **{field: invalid})
        assay = agricultural.OilseedCrushAssay.canola_50lb()
        with self.assertRaises(ValueError):
            replace(assay, oil_yield_pct=90, meal_yield_pct=90)

    def test_reference_states_and_precision_are_not_silently_collapsed(self):
        with self.assertRaisesRegex(ValueError, "reference states"):
            heat.HeatContentAssay(gcv=units.Q_(40, "MJ/scf"), ncv=units.Q_(35, "MJ/Nm3"))
        crude = conversions.CrudeAssay(specific_gravity=.812345678901234)
        self.assertEqual(crude.specific_gravity, .812345678901234)
        self.assertEqual(crude.density_kg_per_m3, crude.specific_gravity * conversions.WATER_DENSITY_60F_KG_M3)

    def test_measured_reference_metadata_must_agree(self):
        with self.assertRaisesRegex(ValueError, "reference states"):
            heat.HeatContentAssay(
                gcv=units.Q_(20, "MJ/kg"), ncv=units.Q_(18, "MJ/kg"),
                gcv_reference_temperature=units.Q_(288.15, "K"),
                ncv_reference_temperature=units.Q_(273.15, "K"),
            )
        with self.assertRaisesRegex(ValueError, "absolute"):
            heat.HeatContentAssay(ncv=units.Q_(18, "MJ/kg"), ncv_reference_pressure=units.Q_(1, "mpag"))

    def test_bootstrap_collisions_are_explicitly_allowlisted(self):
        definitions = registry.COMMODITY_UNIT_DEFINITIONS + ("unauthorized_alias = meter = second",)
        with patch.object(registry, "COMMODITY_UNIT_DEFINITIONS", definitions):
            with self.assertRaises(pint.errors.RedefinitionError):
                with units.sandbox_registry():
                    self.fail("Unexpected collision was accepted")

    def test_shared_symbol_contexts_and_semantic_compatibility(self):
        for expression, code in (("kr25/MWh", "NOK"), ("Rs25/mt", "PKR"), ("Fr25/mt", "XOF"), ("\u00a525/MWh", "CNY")):
            with self.subTest(expression=expression):
                parsed = units.parse_quantity(expression, default_currency=code)
                self.assertIn(f"[currency_{code}]", parsed.dimensionality)
        for source, target in (("wmt", "dmt"), ("psig", "psia"), ("REC", "MWh"), ("EUA", "UKA")):
            self.assertFalse(units.are_convertible(source, target))
            self.assertFalse(units.Q_(1, source).is_compatible_with(target))
            self.assertFalse(units.parse_unit(source).is_compatible_with(target))

    def test_semantic_authorization_does_not_leak_between_threads(self):
        from concurrent.futures import ThreadPoolExecutor
        from core.units._safety import semantic_conversion
        def unauthorized():
            with self.assertRaises(pint.DimensionalityError):
                units.Q_(1, "wmt").to("dmt")
        with semantic_conversion():
            with ThreadPoolExecutor(max_workers=1) as pool:
                pool.submit(unauthorized).result(timeout=5)

    def test_yield_and_negative_mass_guards_apply_to_direct_assay_methods(self):
        with self.assertRaises(ValueError):
            replace(agricultural.CottonSeedLintAssay.us_upland_480lb(), lint_turnout_pct=80, cottonseed_turnout_pct=80)
        assay = agricultural.SeafoodAssay.atlantic_salmon()
        with self.assertRaises(ValueError):
            assay.whole_to_edible(units.Q_(float("nan"), "kg"))
        with self.assertRaises(ValueError):
            assay.whole_to_edible(units.Q_(Decimal(-1), "kg"))
