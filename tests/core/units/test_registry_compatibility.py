"""Broad compatibility coverage against the real commodity unit registry."""

import ast
from decimal import Decimal
import inspect
import math
from pathlib import Path
import subprocess
import sys
import unittest

import pint

from core.units import conversions as c
from core.units import registry as r


def assay_presets():
    for name, cls in vars(c).items():
        if not inspect.isclass(cls) or name not in c.__all__ or not name.endswith("Assay"):
            continue
        for method, descriptor in vars(cls).items():
            if isinstance(descriptor, classmethod):
                factory = getattr(cls, method)
                if all(p.default is not inspect.Parameter.empty or p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD) for p in inspect.signature(factory).parameters.values()):
                    yield name + "." + method, factory()


class RegistryCompatibilityTests(unittest.TestCase):
    coverage_counts = {}

    def test_all_fx_context_currency_pairs(self):
        currencies = sorted({c._normalize_currency(code) for code in r.CURRENCY_CODES})
        fixings = {currency: Decimal(index + 1) / 10 for index, currency in enumerate(currencies) if currency != "USD"}
        fixings["USD"] = Decimal(1)
        context = c.create_fx_context(fixings)
        checked = 0
        with r.ureg.context(context):
            for source in currencies:
                for target in currencies:
                    for magnitude in (1.25, Decimal("1.25")):
                        with self.subTest(source=source, target=target, magnitude=magnitude):
                            quantity = r.ureg.Quantity(magnitude, source)
                            converted = quantity.to(target)
                            expected = Decimal("1.25") * fixings[target] / fixings[source]
                            self.assertAlmostEqual(float(converted.magnitude), float(expected), places=7)
                            if isinstance(magnitude, Decimal):
                                self.assertIsInstance(converted.magnitude, Decimal)
                            checked += 1
        self.coverage_counts["fx_currency_pair_cases"] = checked

    def test_gas_reference_unit_and_rate_aliases(self):
        bases = ("scf", "scf_tx", "scf_60F_14696psi_dry", "Sm3", "Nm3", "Sm3_20C", "m3_Groningen")
        gas_dimensions = {r.ureg.Unit(name).dimensionality for name in bases}
        time_dimension = r.ureg.second.dimensionality
        assay = c.GasAssay(gcv_btu_per_scf=1000)
        checked = 0
        for name in tuple(r.ureg._units):
            unit = r.ureg.Unit(r.ureg.get_name(name))
            if unit.dimensionality in gas_dimensions:
                target = r.ureg.MJ
            elif unit.dimensionality * time_dimension in gas_dimensions:
                target = r.ureg.MJ / r.ureg.second
            else:
                continue
            for magnitude in (1.25, Decimal("1.25")):
                with self.subTest(unit=name, magnitude=magnitude):
                    original = r.ureg.Quantity(magnitude, unit)
                    converted = c.convert_gas_energy(original, target, assay=assay)
                    restored = c.convert_gas_energy(converted, unit, assay=assay)
                    self.assertAlmostEqual(float(restored.magnitude), 1.25, places=7)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(converted.magnitude, Decimal)
                    checked += 1
        self.coverage_counts["gas_unit_and_rate_alias_cases"] = checked
    def test_gas_reference_states_and_gauge_inputs(self):
        temperature = r.ureg.Quantity(300, "K")
        gauge = r.ureg.Quantity(100, "psig")
        absolute = c.convert_pressure(gauge, "psia")
        self.assertEqual(c.calculate_z_factor(gauge, temperature), c.calculate_z_factor(absolute, temperature))
        for standard in ("ferc", "tx", "iso15", "din0", "gost20"):
            for magnitude in (100.0, Decimal("100")):
                with self.subTest(standard=standard, magnitude=magnitude):
                    original = r.ureg.Quantity(magnitude, "m^3")
                    gas = c.convert_actual_to_standard_gas_volume(original, gauge, temperature, target_standard=standard, z_actual=0.95, z_standard=1)
                    equivalent = c.convert_actual_to_standard_gas_volume(original, absolute, temperature, target_standard=standard, z_actual=0.95, z_standard=1)
                    self.assertAlmostEqual(float(gas.magnitude), float(equivalent.magnitude), places=7)
                    restored = c.convert_standard_to_actual_gas_volume(gas, absolute, temperature, source_standard=standard, target_unit="m^3", z_actual=0.95, z_standard=1)
                    self.assertAlmostEqual(float(restored.magnitude), 100)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(restored.magnitude, Decimal)
                    self.assertNotEqual(gas.dimensionality, original.dimensionality)
        for bad in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                c.convert_actual_to_standard_gas_volume(r.ureg.Quantity(1, "m^3"), absolute, temperature, z_actual=bad, z_standard=1)

    def test_runtime_currency_aliases_and_frozen_registry(self):
        # Lifecycle changes stay in a separate interpreter, outside other tests.
        script = '''
from decimal import Decimal
from core.units import registry as r, conversions as c
r.register_currency("TSTFX", decimals=3, aliases=("FXALIAS",), sub_units={"TSTFX_minor": "0.001 * TSTFX"})
r.freeze_registry()
q=r.ureg.Quantity(Decimal("1000"), "TSTFX_minor/oil_barrel")
converted=c.convert_fx(q,"USD",{"FXALIAS/USD":Decimal("2")})
assert converted.magnitude==2 and converted.units==r.ureg.USD/r.ureg.oil_barrel
assert c.triangulate_fx_rate({"USD/FXALIAS":Decimal("0.5")},"USD","TSTFX")==Decimal("0.5")
converted=c.convert_currency(q,"USD",rates={"FXALIAS":Decimal("0.5")})
assert converted.magnitude==2
with r.ureg.context(c.create_fx_context({"FXALIAS":Decimal("0.5")})):
    assert r.ureg.Quantity(Decimal("1"),"USD").to("TSTFX").magnitude==Decimal("0.5")
assert c.convert(r.ureg.Quantity(1,"m^3"),"kg",assay=c.BitumenAssay.pen_60_70()).magnitude>0
custom_quote=r.ureg.Quantity(Decimal("10"),"TSTFX/oil_barrel")
crude=c.CrudeAssay.brent()
converted=c.convert_hydrocarbon(custom_quote,"TSTFX/mt",assay=crude)
restored=c.convert_hydrocarbon(converted,custom_quote.units,assay=crude)
assert abs(restored.magnitude-10)<Decimal("1e-12")
'''
        result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True)
        self.assertEqual(result.returncode, 0, result.stderr)

    def test_flow_rates_and_quotes_on_capacity_units(self):
        for magnitude in (100.0, Decimal("100")):
            for source, target, assay in (
                ("oil_barrel/day", "kg/hour", c.CrudeAssay.brent()),
                ("scfd", "MWh/day", c.GasAssay(gcv_btu_per_scf=1000)),
                ("m^3/day", "kg/hour", c.BitumenAssay.pen_60_70()),
                ("kg/day", "mol/hour", c.ChemicalAssay.caustic_soda_50pct()),
                ("kg/hour", "MWh/day", c.HydrogenAssay.ammonia_carrier()),
            ):
                with self.subTest(source=source, target=target, magnitude=magnitude):
                    original = r.ureg.Quantity(magnitude, source)
                    converted = c.convert(original, target, assay=assay)
                    restored = c.convert(converted, source, assay=assay)
                    self.assertAlmostEqual(float(restored.magnitude), float(magnitude), places=7)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(converted.magnitude, Decimal)
            quote = r.ureg.Quantity(magnitude, "USD/(m^3/day)")
            assay = c.BitumenAssay.pen_60_70()
            converted = c.convert_commodity_quote(quote, "EUR/(kg/hour)", assay=assay, rates={"EUR": Decimal("0.9")})
            restored = c.convert_commodity_quote(converted, quote.units, assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(restored.magnitude), float(magnitude), places=7)
            direct = c.convert_hydrocarbon(r.ureg.Quantity(magnitude, "oil_barrel/day"), "kg/hour", assay=c.CrudeAssay.brent())
            self.assertGreater(direct.magnitude, 0)
            direct = c.convert_gas_energy(r.ureg.Quantity(magnitude, "scfd"), "MWh/day", btu_per_scf=1000)
            self.assertGreater(direct.magnitude, 0)

    def test_registry_singleton_exports_and_referenced_units(self):
        self.assertIs(c.ureg, r.ureg)
        self.assertIs(c.Quantity, r.Quantity)
        self.assertIs(c.Unit, r.Unit)
        r.verify_registry_integrity()
        for module in (r, c):
            self.assertEqual(len(module.__all__), len(set(module.__all__)))
            self.assertFalse(set(module.__all__) - vars(module).keys())
        tree = ast.parse(Path(c.__file__).read_text(encoding="utf-8"))
        references = {node.attr for node in ast.walk(tree) if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "ureg"}
        for attribute in references:
            with self.subTest(attribute=attribute):
                self.assertTrue(hasattr(r.ureg, attribute))

    def test_all_preset_context_edges(self):
        # Index real registry units by dimension; prefer convenient magnitudes.
        representatives = {}
        for name in tuple(r.ureg._units):
            unit = r.ureg.Unit(r.ureg.get_name(name))
            representatives.setdefault(unit.dimensionality, unit)
        for name in ("kg", "m^3", "MJ", "mol", "scf", "scf_60F_14696psi_dry", "scf_tx", "Sm3", "Nm3", "m3_Groningen"):
            unit = r.ureg.Unit(name)
            representatives[unit.dimensionality] = unit
        # Hydrocarbon contexts also register every currency / physical quote.
        for currency in r.CURRENCY_CODES:
            for physical in ("m^3", "kg"):
                unit = r.ureg.Unit(currency) / r.ureg.Unit(physical)
                representatives[unit.dimensionality] = unit
        count = 0
        for label, assay in assay_presets():
            if not hasattr(assay, "as_context"):
                continue
            context = assay.as_context()
            with r.ureg.context(context):
                edges = tuple(context.funcs)
                for source_dim, target_dim in edges:
                    source = representatives[r.ureg.get_dimensionality(source_dim)]
                    target = representatives[r.ureg.get_dimensionality(target_dim)]
                    for magnitude in (2.5, Decimal("2.5")):
                        with self.subTest(preset=label, source=str(source), target=str(target), magnitude=magnitude):
                            quantity = r.ureg.Quantity(magnitude, source)
                            result = quantity.to(target)
                            self.assertTrue(math.isfinite(float(result.magnitude)))
                            self.assertGreater(result.magnitude, 0)
                            if isinstance(magnitude, Decimal):
                                self.assertIsInstance(result.magnitude, Decimal)
                            if (target_dim, source_dim) in context.funcs:
                                self.assertAlmostEqual(float(result.to(source).magnitude), 2.5, places=7)
                            count += 1
        self.coverage_counts["preset_context_edge_cases"] = count
        self.assertGreater(count, 1000)

    def test_all_registry_physical_aliases_through_master_scaling(self):
        dimensions = {r.ureg.Unit(unit).dimensionality: unit for unit in ("kg", "m^3", "MJ", "mol", "Pa")}
        checked = 0
        for name in tuple(r.ureg._units):
            unit = r.ureg.Unit(r.ureg.get_name(name))
            target = dimensions.get(unit.dimensionality)
            if target is None:
                continue
            for magnitude in (1.25, Decimal("1.25")):
                with self.subTest(unit=name, magnitude=magnitude):
                    source = r.ureg.Quantity(magnitude, unit)
                    if r.ureg._semantic_tags(unit._units) and unit.dimensionality != r.ureg.pascal.dimensionality:
                        with self.assertRaises(pint.DimensionalityError):
                            c.convert(source, target)
                        checked += 1
                        continue
                    converted = c.convert(source, target)
                    restored = c.convert(converted, unit)
                    self.assertAlmostEqual(float(restored.magnitude), 1.25, places=7)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(converted.magnitude, Decimal)
                    checked += 1
        self.assertGreater(checked, 100)
        self.coverage_counts["physical_unit_alias_cases"] = checked

    def test_every_currency_and_registered_minor_unit(self):
        currencies = sorted(r.CURRENCY_CODES)
        fixings = {currency: Decimal(index + 1) / 10 for index, currency in enumerate(currencies) if currency != "USD"}
        checked = 0
        for name in tuple(r.ureg._units):
            unit = r.ureg.Unit(r.ureg.get_name(name))
            dimensions = unit.dimensionality
            if len(dimensions) != 1:
                continue
            dim, exponent = next(iter(dimensions.items()))
            if not dim.startswith("[currency_") or exponent != 1:
                continue
            iso = dim[len("[currency_"):-1]
            if iso not in r.CURRENCY_CODES:
                continue
            for magnitude in (1.25, Decimal("1.25")):
                with self.subTest(unit=name, magnitude=magnitude):
                    source = r.ureg.Quantity(magnitude, unit / r.ureg.oil_barrel)
                    usd = c.convert_currency(source, "USD", rates=fixings)
                    restored = c.convert_currency(usd, unit, rates=fixings)
                    self.assertEqual(restored.units, source.units)
                    self.assertAlmostEqual(float(restored.magnitude), 1.25, places=7)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(restored.magnitude, Decimal)
                    checked += 1
        self.coverage_counts["currency_unit_alias_cases"] = checked

    def test_ngl_shrinkage_and_distinct_gas_reference(self):
        for factory in (c.NGLAssay.ethane, c.NGLAssay.propane, c.NGLAssay.n_butane, c.NGLAssay.y_grade):
            assay = factory()
            for magnitude in (100.0, Decimal("100")):
                liquid = r.ureg.Quantity(magnitude, "gallon")
                gas = assay.liquid_to_gas_shrinkage(liquid)
                self.assertEqual(gas.units, r.ureg.scf_60F_14696psi_dry)
                self.assertAlmostEqual(float(assay.gas_shrinkage_to_liquid(gas).magnitude), 100)
        gas = r.ureg.Quantity(1, "scf_60F_14696psi_dry")
        with self.assertRaises(pint.DimensionalityError):
            gas.to("scf")
        with r.ureg.context(c.GasAssay(gcv_btu_per_scf=1000).as_context()):
            self.assertLess(gas.to("scf").magnitude, 1)
            self.assertAlmostEqual(gas.to("scf").to(gas.units).magnitude, 1)


if __name__ == "__main__":
    unittest.main()
