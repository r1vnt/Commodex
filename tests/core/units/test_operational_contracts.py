"""Public imports, assay round trips, quote arbitrage and startup lifecycle."""

from decimal import Decimal
import inspect
import subprocess
import sys
import unittest
from unittest.mock import patch

import pint

import core.units as units
from core.units import agricultural_specs as agricultural
from core.units import config, conversions, converters, currencies, heat_content, quote_parsing, registry
from tests.core.units.test_agricultural_specs import agricultural_presets
from tests.core.units.test_registry_compatibility import assay_presets


class OperationalContracts(unittest.TestCase):
    def test_public_surface_is_complete_and_unique(self):
        from core.units import _spec_engine, material_specs, chemical_specs, energy_environment_specs, materials, agricultural_quality, moisture, exceptions, energy_systems, precision, audit
        from core.units import reference
        modules = (registry, conversions, agricultural, heat_content, converters, config, currencies, quote_parsing, reference,
                   _spec_engine, material_specs, chemical_specs, energy_environment_specs, materials, agricultural_quality, moisture, exceptions, energy_systems, precision, audit)
        expected = {name for module in modules for name in module.__all__}
        self.assertEqual(set(units.__all__), expected)
        self.assertEqual(len(units.__all__), len(expected))
        namespace = {}
        exec("from core.units import *", namespace)
        self.assertEqual(set(namespace) - {"__builtins__"}, expected)
        for module in modules:
            for name in module.__all__:
                self.assertIs(namespace[name], getattr(module, name))

    def test_new_units_and_case_sensitive_nuclear_labels(self):
        for source, target, expected in (
            ("bushel_canola", "pound", 50), ("bushel_rice", "pound", 45),
            ("kgU", "kilogram", 1), ("mtU", "kilogram", 1000),
            ("tU", "kilogram", 1000), ("mtu", "kilogram", 10),
            ("kVA", "watt", 1000), ("MVA", "watt", 1e6),
            ("kVAR", "watt", 1000), ("MVAR", "watt", 1e6),
        ):
            with self.subTest(source=source):
                self.assertEqual(units.Q_(1, source).to(target).magnitude, expected)
        self.assertEqual(units.parse_unit("SWU").dimensionality, {"[separative_work]": 1})
        for minor in ("zac", "zarc", "ZAc", "ZARC"):
            self.assertEqual(units.parse_quantity(f"25 {minor}/mt").units, units.parse_unit("ZAR_cent/mt"))
        for absolute, gauge, neutral in (("kpaa", "kpag", "kPa"), ("mpaa", "mpag", "MPa"), ("paa", "pag", "Pa")):
            converted = units.convert(units.Q_(Decimal(0), gauge), absolute)
            self.assertAlmostEqual(float(converted.to("Pa").magnitude), 101325)
            restored = units.convert(converted, gauge)
            self.assertAlmostEqual(float(restored.magnitude), 0)
            self.assertEqual(units.Q_(1, absolute).to(neutral).magnitude, 1)
        units.verify_registry_integrity()

    def test_every_assay_class_has_a_decimal_round_trip(self):
        # One representative per class; the existing preset suites exhaust the
        # individual products. Explicit processing models use explicit inverses.
        examples = {}
        for _, assay in list(assay_presets()) + list(agricultural_presets()):
            examples.setdefault(type(assay), assay)
        examples[conversions.BulkOreAssay] = conversions.BulkOreAssay(moisture_fraction=.08)
        examples[conversions.RefineryCrudeAssay] = conversions.RefineryCrudeAssay(
            conversions.CrudeAssay.brent(),
            (conversions.DistillationCut("whole", 0, 800, 100, 100, .84),),
        )
        examples[heat_content.HeatContentAssay] = heat_content.HeatContentAssay(ncv=units.Q_(18, "MJ/kg"))
        expected_classes = {
            value for module in (conversions, agricultural, heat_content)
            for value in vars(module).values()
            if inspect.isclass(value) and value.__name__ in getattr(module, "_LEGACY_EXPORTS", module.__all__) and value.__name__.endswith("Assay")
        }
        self.assertEqual(set(examples), expected_classes)
        special_methods = {
            conversions.BatteryMaterialAssay: ("wmt_to_dmt", "dmt_to_wmt"),
            conversions.BullionAssay: ("gross_to_fine", "fine_to_gross"),
            agricultural.LivestockDairyAssay: ("live_to_carcass_mass", "carcass_to_live_mass"),
            agricultural.SeafoodAssay: ("whole_to_edible", "edible_to_whole"),
            agricultural.PoultryEggAssay: ("live_to_rtc_mass", "rtc_to_live_mass"),
        }
        representatives = {}
        for name in tuple(units.ureg._units):
            unit = units.ureg.Unit(name)
            representatives.setdefault(unit.dimensionality, unit)
        for name in ("kg", "m^3", "mol", "MJ", "scf", "Nm3"):
            unit = units.parse_unit(name)
            representatives[unit.dimensionality] = unit
        for cls, assay in examples.items():
            with self.subTest(assay=cls.__name__):
                original = units.Q_(Decimal("2.5"), "kg")
                if cls in special_methods:
                    forward, backward = special_methods[cls]
                    if forward == "wmt_to_dmt":
                        original = units.Q_(Decimal("2.5"), "wmt")
                    converted = getattr(assay, forward)(original)
                    restored = getattr(assay, backward)(converted)
                elif cls is conversions.RefineryCrudeAssay:
                    converted = units.convert(original, "oil_barrel", assay=assay.bulk)
                    restored = units.convert(converted, original.units, assay=assay.bulk)
                elif cls is conversions.BulkOreAssay:
                    original = units.Q_(Decimal("2.5"), "wmt")
                    converted = units.convert(original, "dmt", assay=assay)
                    restored = units.convert(converted, original.units, assay=assay)
                else:
                    context = assay.as_context()
                    edges = [(a, b) for a, b in context.funcs if (b, a) in context.funcs and not any("currency_" in str(dim) for dim in (a, b))]
                    self.assertTrue(edges, cls.__name__)
                    source, target = edges[0]
                    original = units.Q_(Decimal("2.5"), representatives[units.ureg.get_dimensionality(source)])
                    with units.ureg.context(context):
                        converted = original.to(representatives[units.ureg.get_dimensionality(target)])
                        restored = converted.to(original.units)
                self.assertIsInstance(converted.magnitude, Decimal)
                self.assertIsInstance(restored.magnitude, Decimal)
                self.assertAlmostEqual(float(restored.to(original.units).magnitude), 2.5, places=7)

    def test_compound_quotes_through_master_router(self):
        rates = {"EUR": Decimal(".9"), "GBP": Decimal(".8")}
        cases = [
            ("EUR/MWh", "USD/MMBtu", None),
            ("p/th", "$/MMBtu", None),
            ("c/bu", "$/mt", agricultural.GrainAssay.corn_56lb()),
            ("$/gal", "$/mt", conversions.BiofuelAssay.ethanol()),
            ("$/gal", "$/mt", conversions.NGLAssay.propane()),
            ("$/gal", "$/mt", agricultural.OliveOilAssay.extra_virgin_evoo()),
        ]
        for source, target, assay in cases:
            for amount in (Decimal("125.123456789"), Decimal(0), Decimal(-10)):
                with self.subTest(source=source, target=target, assay=type(assay).__name__, amount=amount):
                    original = units.Q_(amount, units.parse_unit(source))
                    converted = units.convert(original, target, assay=assay, rates=rates)
                    restored = units.convert(converted, original.units, assay=assay, rates=rates)
                    self.assertIsInstance(converted.magnitude, Decimal)
                    self.assertIsInstance(restored.magnitude, Decimal)
                    self.assertAlmostEqual(float(restored.magnitude), float(amount), places=7)
        energy = units.convert(units.Q_(Decimal(1), "EUR/MWh"), "USD/MMBtu", rates=rates)
        expected = Decimal(1) / rates["EUR"] * Decimal("1055.05585262") / Decimal(3600)
        self.assertAlmostEqual(float(energy.magnitude), float(expected), places=12)

    def test_auto_wrappers_use_master_but_processing_modes_remain_explicit(self):
        cases = (
            (agricultural.convert_softs, agricultural.SoftsAssay.raw_sugar_no11(), "raw_to_processed"),
            (agricultural.convert_livestock_dairy, agricultural.LivestockDairyAssay.live_cattle_standard(), "live_to_carcass"),
        )
        for converter, assay, mode in cases:
            with patch.object(conversions, "convert", wraps=conversions.convert) as master:
                result = converter(units.Q_(Decimal(10), "kg"), "g", assay=assay)
                master.assert_called_once()
                self.assertEqual(result.magnitude, 10000)
            with patch.object(conversions, "convert", side_effect=AssertionError("processing used master")):
                result = converter(units.Q_(Decimal(10), "kg"), "kg", assay=assay, mode=mode)
                self.assertIsInstance(result.magnitude, Decimal)

    def test_custom_latent_heat_round_trip_and_explicit_crude_calorific_data(self):
        for magnitude in (18., Decimal(18)):
            original = units.Q_(magnitude, "MJ/kg")
            dry = heat_content.as_received_to_dry(original, 20, kind="ncv", latent_heat_mj_per_kg_water=Decimal("2.5"), method="latent_water")
            restored = heat_content.dry_to_as_received(dry, 20, kind="ncv", latent_heat_mj_per_kg_water=Decimal("2.5"), method="latent_water")
            self.assertAlmostEqual(float(restored.magnitude), float(magnitude), places=12)
            self.assertIsInstance(restored.magnitude, type(magnitude))
        for invalid in (-1, 0, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                heat_content.dry_to_as_received(units.Q_(18, "MJ/kg"), 20, latent_heat_mj_per_kg_water=invalid)
        crude = conversions.CrudeAssay.brent()
        with self.assertRaisesRegex(ValueError, "Supply at least one"):
            units.convert(units.Q_(1, "oil_barrel"), "MMBtu", assay=crude)
        explicit = heat_content.heat_content_from_assay(
            crude, ncv=units.Q_(42, "MJ/kg"), density=units.Q_(crude.density_kg_per_m3, "kg/m^3"),
        )
        original = units.Q_(Decimal(1), "oil_barrel")
        energy = units.convert(original, "MMBtu", assay=crude, heat_content=explicit)
        restored = units.convert(energy, original.units, assay=crude, heat_content=explicit)
        self.assertAlmostEqual(float(restored.magnitude), 1, places=12)

    def test_incompatible_and_nonfinite_inputs_fail_loudly(self):
        with self.assertRaisesRegex(ValueError, "FX rates"):
            units.convert(units.Q_(1, "USD"), "EUR")
        with self.assertRaises(pint.DimensionalityError):
            units.convert(units.Q_(1, "oil_barrel"), "mt")
        with self.assertRaisesRegex(ValueError, "Supply at least one"):
            units.convert(units.Q_(1, "kg"), "MJ", assay=agricultural.GrainAssay.corn_56lb())
        # Tagged pressure dimensionality alone never replaces datum routing.
        with self.assertRaises(ValueError):
            units.convert(units.Q_(-1, "psia"), "psig")
        absolute = units.convert(units.Q_(0, "psig"), "psia")
        self.assertAlmostEqual(float(absolute.magnitude), conversions.STANDARD_ATMOSPHERE_PSI)
        for magnitude in (Decimal(-1), float("nan"), Decimal("Infinity"), float("-inf")):
            with self.subTest(magnitude=magnitude), self.assertRaises(ValueError):
                units.convert(units.Q_(magnitude, "kg"), "g")

    def test_desk_lifecycle_and_lossless_custom_serialisation(self):
        script = '''
from decimal import Decimal
import pint
import core.units as u
application = pint.get_application_registry().get()
u.load_desk_config({
    "currencies": [{"code": "LIFECYCLECCY", "decimals": 8, "symbols": ["LC$"]}],
    "commodities": [{"definition": "lifecycle_parcel = 123 * kilogram"}],
})
assert u.is_registry_frozen()
for action in (lambda: u.register_currency("FROZENCCY", decimals=2),
               lambda: u.register_commodity_unit("frozen_parcel = kilogram")):
    try: action()
    except RuntimeError: pass
    else: raise AssertionError("Registration after freeze succeeded")
original = u.Q_(Decimal("1234567890.000000000123400"), "LIFECYCLECCY/lifecycle_parcel")
restored = u.quantity_from_dict(u.quantity_to_dict(original))
assert restored.units == original.units
assert str(restored.magnitude) == str(original.magnitude)
with u.sandbox_registry() as sandbox:
    sandbox.define("sandbox_lifecycle_parcel = 7 * kilogram")
    assert "sandbox_lifecycle_parcel" not in u.ureg
assert u.is_registry_frozen()
assert pint.get_application_registry().get() is application
u.verify_registry_integrity()
'''
        completed = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)


if __name__ == "__main__":
    unittest.main()
