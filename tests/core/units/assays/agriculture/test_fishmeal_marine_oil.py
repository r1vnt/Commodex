"""Fishmeal/marine-oil reduction and explicit product-density regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class FishmealMarineOilAssayTests(unittest.TestCase):
    def test_exports_presets_and_reduction_round_trips(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("FishmealMarineOilAssay", "create_fishmeal_marine_oil_context", "convert_fishmeal_marine_oil"):
            self.assertIn(name, a.__all__)
        for assay in (a.FishmealMarineOilAssay.peruvian_anchoveta_super_prime(), a.FishmealMarineOilAssay.menhaden_standard()):
            for magnitude in (100.0, Decimal("100")):
                whole = ureg.Quantity(magnitude, "mt")
                products = assay.whole_to_meal_and_oil(whole)
                self.assertAlmostEqual(float(products["meal"].magnitude), assay.meal_yield_pct)
                self.assertAlmostEqual(float(products["oil"].magnitude), assay.oil_yield_pct)
                self.assertAlmostEqual(float(sum(q.magnitude for q in products.values())), assay.combined_yield_pct)
                self.assertLessEqual(sum(q.magnitude for q in products.values()), magnitude)
                for forward, reverse in (("whole_to_meal", "meal_to_whole"), ("whole_to_oil", "oil_to_whole")):
                    with self.subTest(species=assay.species, mode=forward, magnitude=magnitude):
                        output = a.convert_fishmeal_marine_oil(whole, "kg", assay=assay, mode=forward)
                        restored = a.convert_fishmeal_marine_oil(output, "mt", assay=assay, mode=reverse)
                        self.assertAlmostEqual(float(restored.magnitude), 100)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(output.magnitude, Decimal)

    def test_meal_and_oil_density_modes_and_product_selection(self):
        assay = a.FishmealMarineOilAssay.peruvian_anchoveta_super_prime()
        for product, density in (("meal", 640), ("oil", 925)):
            for magnitude in (2.5, Decimal("2.5")):
                with self.subTest(product=product, magnitude=magnitude):
                    volume = ureg.Quantity(magnitude, "m^3")
                    mass = a.convert_fishmeal_marine_oil(volume, "kg", assay=assay, product=product)
                    self.assertAlmostEqual(float(mass.magnitude), 2.5 * density)
                    restored = a.convert_fishmeal_marine_oil(mass, "m^3", assay=assay, product=product)
                    self.assertAlmostEqual(float(restored.magnitude), 2.5)
                    explicit = a.convert_fishmeal_marine_oil(volume, "kg", assay=assay, mode=product + "_vol_to_mass")
                    self.assertEqual(explicit, mass)
                    explicit = a.convert_fishmeal_marine_oil(mass, "m^3", assay=assay, mode=product + "_mass_to_vol")
                    self.assertAlmostEqual(float(explicit.magnitude), 2.5)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(mass.magnitude, Decimal)
        # Pint's litre -> cubic metre factor carries floating-definition rounding.
        self.assertAlmostEqual(float(assay.oil_volume_to_mass(ureg.Quantity(Decimal("1000"), "liter")).magnitude), 0.925, places=12)

    def test_product_selection_through_rates_quotes_and_master(self):
        assay = a.FishmealMarineOilAssay.peruvian_anchoveta_super_prime()
        for product, density in (("meal", 640), ("oil", 925)):
            rate = a.convert_fishmeal_marine_oil(ureg.Quantity(Decimal("24"), "m^3/day"), "kg/hour", assay=assay, product=product)
            self.assertAlmostEqual(float(rate.magnitude), density)
            result = c.convert(ureg.Quantity(Decimal("1"), "m^3"), "kg", assay=assay, product=product)
            self.assertEqual(result.magnitude, density)
            price = ureg.Quantity(Decimal("1000"), "USD/m^3")
            converted = c.convert_commodity_quote(price, "EUR/kg", assay=assay, product=product, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(converted.magnitude), 900 / density)
            master = c.convert(price, "EUR/kg", assay=assay, product=product, rates={"USD/EUR": Decimal("0.9")})
            self.assertEqual(master, converted)

    def test_standard_and_jumbo_bag_counts(self):
        standard = a.FishmealMarineOilAssay.peruvian_anchoveta_super_prime()
        for assay in (standard, replace(standard, bag_mass_kg=1000)):
            for count in (20.0, Decimal("20.123456789")):
                mass = assay.bags_to_meal_mass(count)
                self.assertAlmostEqual(float(mass.to("kg").magnitude), float(count) * assay.bag_mass_kg)
                self.assertAlmostEqual(float(assay.meal_mass_to_bags(mass)), float(count))
                if isinstance(count, Decimal):
                    self.assertIsInstance(mass.magnitude, Decimal)
                    self.assertIsInstance(assay.meal_mass_to_bags(mass), Decimal)
        self.assertEqual(standard.bags_to_meal_mass(Decimal("20")).magnitude, 1)

    def test_gross_product_value_and_zero_oil(self):
        assay = a.FishmealMarineOilAssay.peruvian_anchoveta_super_prime()
        self.assertEqual(assay.gross_product_value_per_mt_whole(Decimal("1500"), Decimal("2000")), 445)
        self.assertEqual(assay.gross_product_value_per_mt_whole(1500.0, 2000.0), 445)
        self.assertIsInstance(assay.gross_product_value_per_mt_whole(Decimal("1500"), 2000), Decimal)
        zero = replace(assay, oil_yield_pct=0)
        self.assertEqual(zero.whole_to_oil_mass(ureg.Quantity(Decimal("1"), "mt")).magnitude, 0)
        self.assertEqual(zero.gross_product_value_per_mt_whole(Decimal("1500"), Decimal("2000")), 345)
        with self.assertRaises(ValueError):
            zero.oil_to_whole_mass(ureg.Quantity(1, "mt"))

    def test_validation_and_context_isolation(self):
        assay = a.FishmealMarineOilAssay.peruvian_anchoveta_super_prime()
        invalid = {"species": ("", None), "meal_yield_pct": (0, 101, float("nan")), "oil_yield_pct": (-1, 78), "protein_pct": (39, 81), "oil_density_kg_per_m3": (749, 1101), "meal_bulk_density_kg_per_m3": (299, 1201), "bag_mass_kg": (0, -1)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.whole_to_meal_mass(ureg.Quantity(bad, "kg"))
            with self.assertRaises(ValueError):
                assay.oil_volume_to_mass(ureg.Quantity(bad, "m^3"))
            with self.assertRaises(ValueError):
                assay.bags_to_meal_mass(bad)
        with self.assertRaises(ValueError):
            assay.gross_product_value_per_mt_whole(Decimal("NaN"), 2000)
        for product in ("invalid", None):
            with self.assertRaises(ValueError):
                assay.as_context(product=product)
        with self.assertRaises(ValueError):
            c.convert(ureg.Quantity(1, "m^3"), "kg", assay=assay, product="invalid")
        with self.assertRaises(pint.DimensionalityError):
            ureg.Quantity(1, "m^3").to("kg")
        with self.assertRaises(TypeError):
            a.convert_fishmeal_marine_oil(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_fishmeal_marine_oil(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()
