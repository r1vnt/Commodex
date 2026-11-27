"""Module-wide agricultural API, preset, registry and conversion checks."""

import ast
from collections import Counter
from decimal import Decimal
import dis
import inspect
import math
from pathlib import Path
from types import CodeType
import unittest

import pint

import core.units as package
from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units import heat_content as h
from core.units import registry as registry
from core.units import config as config
from core.units import currencies
from core.units import quote_parsing
from core.units import converters
from core.units.registry import ureg


def agricultural_presets():
    for name, cls in vars(a).items():
        if not inspect.isclass(cls) or name not in a._LEGACY_EXPORTS or not name.endswith("Assay"):
            continue
        for method, descriptor in vars(cls).items():
            if isinstance(descriptor, classmethod):
                factory = getattr(cls, method)
                if all(p.default is not inspect.Parameter.empty or p.kind in (p.VAR_POSITIONAL, p.VAR_KEYWORD) for p in inspect.signature(factory).parameters.values()):
                    yield name + "." + method, factory()


class AgriculturalIntegrationTests(unittest.TestCase):
    coverage_counts = {}

    def test_complete_public_api_and_package_reexports(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        from core.units import _spec_engine, material_specs, chemical_specs, energy_environment_specs, materials, agricultural_quality, moisture, exceptions, energy_systems, precision, audit
        from core.units import reference
        modules = (a, h, c, registry, config, currencies, quote_parsing, converters, reference,
                   _spec_engine, material_specs, chemical_specs, energy_environment_specs, materials, agricultural_quality, moisture, exceptions, energy_systems, precision, audit)
        self.assertEqual(set(package.__all__), {name for module in modules for name in module.__all__})
        classes = []
        for name, value in vars(a).items():
            if name not in a._LEGACY_EXPORTS:
                continue
            if inspect.isclass(value) and name.endswith("Assay"):
                classes.append(name)
            if inspect.isclass(value) or name.startswith(("create_", "convert_")):
                self.assertIn(name, a.__all__)
        self.assertEqual(len(classes), 15)
        namespace = {}
        exec("from core.units import *", namespace)
        for name in a.__all__:
            self.assertIs(namespace[name], getattr(a, name))
            self.assertIs(getattr(package, name), getattr(a, name))

    def test_no_duplicate_definitions_missing_helpers_or_registry_units(self):
        tree = ast.parse(Path(a.__file__).read_text(encoding="utf-8"))
        for body in [tree.body] + [node.body for node in tree.body if isinstance(node, ast.ClassDef)]:
            names = Counter(node.name for node in body if isinstance(node, (ast.FunctionDef, ast.ClassDef)))
            self.assertFalse([name for name, count in names.items() if count > 1])
        for node in ast.walk(tree):
            if isinstance(node, ast.Attribute) and isinstance(node.value, ast.Name) and node.value.id == "ureg":
                self.assertTrue(hasattr(ureg, node.attr), node.attr)

        def audit_code(code, namespace, builtins, label):
            missing = {item.argval for item in dis.get_instructions(code) if item.opname == "LOAD_GLOBAL" and item.argval not in namespace and item.argval not in builtins}
            self.assertFalse(missing, f"{label}: {missing}")
            for constant in code.co_consts:
                if isinstance(constant, CodeType):
                    audit_code(constant, namespace, builtins, label)

        def audit(value, label):
            if inspect.isfunction(value) and value.__module__ == a.__name__:
                audit_code(value.__code__, value.__globals__, value.__builtins__, label)

        for name, value in vars(a).items():
            audit(value, name)
            if inspect.isclass(value) and value.__module__ == a.__name__:
                for method, function in vars(value).items():
                    if isinstance(function, (classmethod, staticmethod)):
                        function = function.__func__
                    if isinstance(function, property):
                        function = function.fget
                    audit(function, name + "." + method)

    def test_all_presets_properties_and_bidirectional_contexts(self):
        presets = list(agricultural_presets())
        self.assertGreater(len(presets), 30)
        checked = 0
        for label, assay in presets:
            for name, descriptor in vars(type(assay)).items():
                if isinstance(descriptor, property):
                    value = getattr(assay, name)
                    if isinstance(value, (int, float, Decimal)):
                        self.assertTrue(math.isfinite(float(value)), label + "." + name)
            context = assay.as_context()
            if not context.funcs:
                with self.assertRaises(pint.DimensionalityError):
                    c.convert(ureg.Quantity(1, "m^3"), "kg", assay=assay)
                continue
            for volume_unit in ("m^3", "liter", "gallon", "foot^3"):
                for mass_unit in ("kg", "lb", "metric_ton", "short_ton"):
                    for magnitude in (2.5, Decimal("2.5")):
                        with self.subTest(preset=label, volume=volume_unit, mass=mass_unit, magnitude=magnitude):
                            original = ureg.Quantity(magnitude, volume_unit)
                            mass = c.convert(original, mass_unit, assay=assay)
                            restored = c.convert(mass, volume_unit, assay=assay)
                            self.assertAlmostEqual(float(restored.magnitude), 2.5, places=7)
                            self.assertGreater(mass.magnitude, 0)
                            if isinstance(magnitude, Decimal):
                                self.assertIsInstance(mass.magnitude, Decimal)
                            checked += 1
            rate = c.convert(ureg.Quantity(Decimal("24"), "m^3/day"), "kg/hour", assay=assay)
            restored = c.convert(rate, "m^3/day", assay=assay)
            self.assertAlmostEqual(float(restored.magnitude), 24, places=7)
            quote = ureg.Quantity(Decimal("100"), "USD/m^3")
            converted = c.convert_commodity_quote(quote, "EUR/kg", assay=assay, rates={"EUR": Decimal("0.9")})
            restored = c.convert_commodity_quote(converted, quote.units, assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(restored.magnitude), 100, places=7)
            with self.assertRaises(pint.DimensionalityError):
                ureg.Quantity(1, "m^3").to("kg")
            checked += 2
        self.coverage_counts.update(presets=len(presets), preset_conversion_cases=checked)

    def test_legacy_softs_oilseed_and_livestock_workflows(self):
        sugar = a.SoftsAssay.raw_sugar_no11()
        coffee = a.SoftsAssay.green_coffee_60kg_bag()
        soy = a.OilseedCrushAssay.soybeans_cb_60lb()
        cattle = a.LivestockDairyAssay.live_cattle_standard()
        milk = a.LivestockDairyAssay.whole_milk_standard()
        for magnitude in (100.0, Decimal("100")):
            mass = ureg.Quantity(magnitude, "mt")
            refined = a.convert_softs(mass, "mt", assay=sugar, mode="raw_to_processed")
            self.assertAlmostEqual(float(refined.magnitude), 92.5)
            restored = a.convert_softs(refined, "mt", assay=sugar, mode="processed_to_raw")
            self.assertAlmostEqual(float(restored.magnitude), 100)
            bag_mass = coffee.packages_to_mass(magnitude)
            self.assertAlmostEqual(float(bag_mass.to("kg").magnitude), 6000)
            self.assertAlmostEqual(coffee.mass_to_packages(bag_mass), 100)
            oil = a.convert_oilseed_crush(mass, "mt", assay=soy, mode="seed_to_oil")
            meal = a.convert_oilseed_crush(mass, "mt", assay=soy, mode="seed_to_meal")
            self.assertAlmostEqual(float(oil.magnitude), 100 * 11 / 60)
            self.assertAlmostEqual(float(meal.magnitude), 100 * 44 / 60)
            carcass = a.convert_livestock_dairy(mass, "mt", assay=cattle, mode="live_to_carcass")
            self.assertAlmostEqual(float(carcass.magnitude), 62)
            restored = a.convert_livestock_dairy(carcass, "mt", assay=cattle, mode="carcass_to_live")
            self.assertAlmostEqual(float(restored.magnitude), 100)
            milk_volume = ureg.Quantity(magnitude, "liter")
            solids = a.convert_livestock_dairy(milk_volume, "kg", assay=milk, mode="milk_to_solids")
            self.assertAlmostEqual(float(solids.magnitude), 7.519)
        self.assertAlmostEqual(soy.gross_crush_margin_usd_per_mt(400, 1000, 400), 76.6667, places=4)

    def test_legacy_count_precision_and_oilseed_thermal_routes(self):
        count = Decimal("10.123456789123456789")
        coffee = a.SoftsAssay.green_coffee_60kg_bag()
        mass = coffee.packages_to_mass(count)
        recovered = coffee.mass_to_packages(mass)
        self.assertIsInstance(recovered, Decimal)
        self.assertLess(abs(recovered - count), Decimal("1e-20"))
        soy = a.OilseedCrushAssay.soybeans_cb_60lb()
        mass = soy.bushels_to_seed_mass(count)
        recovered = soy.seed_mass_to_bushels(mass)
        self.assertIsInstance(recovered, Decimal)
        # Pint's pound/metric conversion factors limit the final few digits.
        self.assertLess(abs(recovered - count), Decimal("1e-14"))
        for magnitude in (2.5, Decimal("2.5")):
            for source, target in (("kg", "lb"), ("m^3", "liter")):
                original = ureg.Quantity(magnitude, source)
                converted = a.convert_oilseed_crush(original, target, assay=soy, temp_c=60)
                self.assertAlmostEqual(float(converted.magnitude), float(original.to(target).magnitude))
            volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_oilseed_crush(volume, "kg", assay=soy, temp_c=60)
            expected_density = 920 * (1 - 0.0008 * (60 - 25))
            self.assertAlmostEqual(float(mass.magnitude), float(magnitude) * expected_density)
            recovered = a.convert_oilseed_crush(mass, "m^3", assay=soy, temp_c=60)
            self.assertAlmostEqual(float(recovered.magnitude), float(magnitude))
        with self.assertRaises(ValueError):
            soy.as_context(temp_c=float("nan"))

    def test_product_temperature_variants_and_downstream_cottonseed_crush(self):
        checked = 0
        for label, assay in agricultural_presets():
            parameters = inspect.signature(assay.as_context).parameters
            options = []
            if "product" in parameters:
                options.extend(({"product": "meal"}, {"product": "oil"}))
            if "temp_c" in parameters:
                options.extend(({"temp_c": 20}, {"temp_c": 40}))
            for option in options:
                with self.subTest(preset=label, option=option):
                    context = assay.as_context(**option)
                    original = ureg.Quantity(Decimal("1"), "m^3")
                    with ureg.context(context):
                        direct = original.to("kg")
                        restored = direct.to("m^3")
                    kwargs = {"assay": assay}
                    if "product" in option:
                        kwargs["product"] = option["product"]
                    if "temp_c" in option:
                        kwargs["observed_temp"] = option["temp_c"]
                    routed = c.convert(original, "kg", **kwargs)
                    self.assertAlmostEqual(float(routed.magnitude), float(direct.magnitude))
                    self.assertAlmostEqual(float(restored.magnitude), 1)
                    checked += 1
        cotton = a.CottonSeedLintAssay.us_upland_480lb()
        # Explicit test-lot extraction assumptions, rather than a market preset.
        crush = a.OilseedCrushAssay("Test cottonseed lot", "cottonseed", 50, 18, 45)
        raw = ureg.Quantity(Decimal("100"), "mt")
        seed = cotton.gin_seed_cotton(raw)["cottonseed"]
        oil = a.convert_oilseed_crush(seed, "mt", assay=crush, mode="seed_to_oil")
        meal = a.convert_oilseed_crush(seed, "mt", assay=crush, mode="seed_to_meal")
        self.assertAlmostEqual(float(oil.magnitude), 9.72)
        self.assertAlmostEqual(float(meal.magnitude), 24.3)
        self.coverage_counts["product_temperature_variants"] = checked


if __name__ == "__main__":
    unittest.main()
