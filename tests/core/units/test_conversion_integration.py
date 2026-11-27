"""Whole-module consistency and float/Decimal commodity round trips."""

import ast
from collections import Counter
from decimal import Decimal
import dis
import inspect
from pathlib import Path
from types import CodeType
import unittest

from core.units import conversions as c
from core.units.registry import ureg


class ConversionIntegrationTests(unittest.TestCase):
    def test_symbols_helpers_and_global_references(self):
        for name in (*c.__all__, "CrudeAssay", "GasAssay", "CoalAssay", "CleanFuelAssay", "_scale_factor", "parse_unit"):
            self.assertIn(name, vars(c))
        self.assertEqual(len(c.__all__), len(set(c.__all__)))

        def audit(code, namespace, builtins, label):
            missing = {instruction.argval for instruction in dis.get_instructions(code) if instruction.opname == "LOAD_GLOBAL" and instruction.argval not in namespace and instruction.argval not in builtins}
            self.assertFalse(missing, f"Undefined helpers in {label}: {sorted(missing)}")
            for constant in code.co_consts:
                if isinstance(constant, CodeType):
                    audit(constant, namespace, builtins, label + "." + constant.co_name)

        def check(function, label):
            if inspect.isfunction(function) and function.__module__ == c.__name__:
                # Generated dataclass methods have their own globals dictionary.
                audit(function.__code__, function.__globals__, function.__builtins__, label)

        for name, value in vars(c).items():
            check(value, name)
            if inspect.isclass(value) and value.__module__ == c.__name__:
                for method, function in vars(value).items():
                    if isinstance(function, (classmethod, staticmethod)):
                        function = function.__func__
                    if isinstance(function, property):
                        function = function.fget
                    check(function, name + "." + method)

    def test_no_duplicate_definitions(self):
        tree = ast.parse(Path(c.__file__).read_text(encoding="utf-8"))

        def check(body, label):
            names = Counter(node.name for node in body if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)))
            duplicates = [name for name, count in names.items() if count > 1]
            self.assertFalse(duplicates, f"Duplicate definitions in {label}: {duplicates}")
        check(tree.body, "module")
        for node in tree.body:
            if isinstance(node, ast.ClassDef):
                check(node.body, node.name)

    def test_float_and_decimal_round_trip_matrix(self):
        bitumen = c.BitumenAssay.pen_60_70()
        petcoke = c.PetcokeAssay.usgc_fuel_grade_65s()
        chemical = c.ChemicalAssay.caustic_soda_50pct()
        hydrogen = c.HydrogenAssay.ammonia_carrier()
        ferro = c.FerroalloyAssay.silicomanganese()
        metal = c.MetalAlloyAssay.aluminium_6061()
        coal = c.CoalAssay(nar_kcal_per_kg=6000)
        clean = c.CleanFuelAssay.hydrogen()

        for magnitude in (100.0, Decimal("100")):
            cases = (
                ("hot bitumen", "m^3", lambda q: c.convert(q, "kg", assay=bitumen, observed_temp=150), lambda q: c.convert(q, "m^3", assay=bitumen, observed_temp=150)),
                ("petcoke moisture", "wmt", petcoke.wmt_to_dmt, petcoke.dmt_to_wmt),
                ("petcoke MMBtu", "wmt", lambda q: c.convert(q, "MMBtu_it", assay=petcoke), lambda q: c.convert(q, "wmt", assay=petcoke)),
                ("chemical active", "mt", chemical.solution_to_active_mass, chemical.active_to_solution_mass),
                ("carrier cracking", "mt", hydrogen.carrier_to_h2_mass, hydrogen.h2_to_carrier_mass),
                ("ferroalloy contained", "mt", lambda q: ferro.alloy_to_contained_metal(q, element="Si"), lambda q: ferro.contained_metal_to_alloy(q, element="Si")),
                ("metal contained", "mt", lambda q: metal.alloy_to_contained_element(q, "Al"), lambda q: metal.contained_element_to_alloy(q, "Al")),
                ("metal remelt", "mt", metal.gross_to_recovered_melt_mass, metal.recovered_melt_to_gross_mass),
                ("pressure", "psig", lambda q: c.convert_pressure(q, "bara"), lambda q: c.convert_pressure(q, "psig")),
                ("FX", "US_cent/oil_barrel", lambda q: c.convert_fx(q, "GBP_pence", {"USD/GBP": Decimal("0.8")}), lambda q: c.convert_fx(q, "US_cent", {"USD/GBP": Decimal("0.8")})),
                ("quote", "USD/m^3", lambda q: c.convert_commodity_quote(q, "EUR/kg", assay=bitumen, rates={"EUR": Decimal("0.9")}), lambda q: c.convert_commodity_quote(q, "USD/m^3", assay=bitumen, rates={"EUR": Decimal("0.9")})),
                ("coal context", "mt", lambda q: c.convert(q, "GJ", assay=coal), lambda q: c.convert(q, "mt", assay=coal)),
                ("clean fuel context", "Nm3_H2", lambda q: c.convert(q, "MWh", assay=clean), lambda q: c.convert(q, "Nm3_H2", assay=clean)),
            )
            for name, unit, forward, reverse in cases:
                with self.subTest(case=name, magnitude=magnitude):
                    original = ureg.Quantity(magnitude, unit)
                    converted = forward(original)
                    restored = reverse(converted).to(original.units)
                    self.assertAlmostEqual(float(restored.magnitude), float(magnitude), places=8)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(converted.magnitude, Decimal)
                        self.assertIsInstance(restored.magnitude, Decimal)
            feed = ureg.Quantity(magnitude, "mt")
            water = chemical.water_required_to_dilute(feed, 25)
            diluted = chemical.dilute_to(25)
            active = diluted.solution_to_active_mass(feed + water)
            self.assertAlmostEqual(float(active.magnitude), float(chemical.solution_to_active_mass(feed).magnitude))

    def test_earlier_context_reference_values_and_options(self):
        coal = c.CoalAssay(nar_kcal_per_kg=6000)
        self.assertAlmostEqual(ureg.Quantity(1, "mt").to("GJ", coal.as_context()).magnitude, coal.gj_per_mt_nar)
        self.assertGreater(c.convert(ureg.Quantity(1, "mt"), "GJ", assay=coal, basis="gar").magnitude, coal.gj_per_mt_nar)
        for fuel in (c.CleanFuelAssay.hydrogen(), c.CleanFuelAssay.ammonia(), c.CleanFuelAssay.methanol()):
            self.assertAlmostEqual(ureg.Quantity(1, "kg").to("MJ", fuel.as_context()).magnitude, fuel.lhv_mj_per_kg)
            with self.assertRaises(ValueError):
                fuel.as_context(basis="invalid")
        with self.assertRaises(ValueError):
            coal.as_context(basis="invalid")


if __name__ == "__main__":
    unittest.main()
