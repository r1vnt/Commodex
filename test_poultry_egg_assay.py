"""Poultry carcass and egg count/processing integration regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class PoultryEggAssayTests(unittest.TestCase):
    def test_exports_and_poultry_round_trips(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("PoultryEggAssay", "create_poultry_egg_context", "convert_poultry_egg"):
            self.assertIn(name, a.__all__)
        for assay in (a.PoultryEggAssay.broiler_chicken_standard(), a.PoultryEggAssay.turkey_standard()):
            for magnitude in (100.0, Decimal("100")):
                for forward, reverse, percentage in (("live_to_rtc", "rtc_to_live", assay.rtc_dressing_pct), ("live_to_boneless", "boneless_to_live", assay.boneless_meat_pct_of_live)):
                    with self.subTest(category=assay.category, mode=forward, magnitude=magnitude):
                        live = ureg.Quantity(magnitude, "kg")
                        changed = a.convert_poultry_egg(live, "kg", assay=assay, mode=forward)
                        self.assertAlmostEqual(float(changed.magnitude), percentage)
                        restored = a.convert_poultry_egg(changed, "kg", assay=assay, mode=reverse)
                        self.assertAlmostEqual(float(restored.magnitude), 100)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(changed.magnitude, Decimal)

    def test_egg_dozen_and_case_counts(self):
        assay = a.PoultryEggAssay.us_large_shell_eggs()
        self.assertEqual(assay.eggs_to_shell_mass(Decimal("1")).magnitude, Decimal("0.0567"))
        self.assertEqual(assay.dozens_to_shell_mass(Decimal("1")).magnitude, Decimal("0.6804"))
        self.assertEqual(assay.cases_to_shell_mass(Decimal("1")).magnitude, Decimal("20.412"))
        self.assertEqual(assay.dozen_mass_kg, 0.6804)
        self.assertEqual(assay.case_mass_kg, 20.412)
        for count in (10.0, Decimal("10.123456789")):
            for forward, reverse in ((assay.eggs_to_shell_mass, assay.shell_mass_to_eggs), (assay.dozens_to_shell_mass, assay.shell_mass_to_dozens), (assay.cases_to_shell_mass, assay.shell_mass_to_cases)):
                mass = forward(count)
                restored = reverse(mass)
                self.assertAlmostEqual(float(restored), float(count))
                if isinstance(count, Decimal):
                    self.assertIsInstance(mass.magnitude, Decimal)
                    self.assertIsInstance(restored, Decimal)

    def test_egg_recoveries_are_relative_to_shell_mass(self):
        assay = a.PoultryEggAssay.us_large_shell_eggs()
        for magnitude in (100.0, Decimal("100")):
            shell = ureg.Quantity(magnitude, "kg")
            for forward, reverse, expected in (("shell_to_liquid", "liquid_to_shell", 88), ("shell_to_dried", "dried_to_shell", 24)):
                result = a.convert_poultry_egg(shell, "kg", assay=assay, mode=forward)
                self.assertAlmostEqual(float(result.magnitude), expected)
                restored = a.convert_poultry_egg(result, "kg", assay=assay, mode=reverse)
                self.assertAlmostEqual(float(restored.magnitude), 100)
                if isinstance(magnitude, Decimal):
                    self.assertIsInstance(result.magnitude, Decimal)
        self.assertEqual(a.convert_poultry_egg(ureg.Quantity(100, "kg"), "kg", assay=assay).magnitude, 100)

    def test_poultry_and_egg_price_parity(self):
        broiler = a.PoultryEggAssay.broiler_chicken_standard()
        self.assertEqual(broiler.live_price_to_rtc_price(Decimal("0.74")), 1)
        self.assertEqual(broiler.rtc_price_to_live_price(Decimal("1")), Decimal("0.74"))
        eggs = a.PoultryEggAssay.us_large_shell_eggs()
        for price in (2.0, Decimal("2")):
            liquid_price = eggs.usd_per_dozen_to_liquid_usd_per_lb(price)
            pounds_per_dozen = eggs.shell_to_liquid_egg_mass(eggs.dozens_to_shell_mass(Decimal("1"))).to("lb").magnitude
            self.assertAlmostEqual(float(liquid_price), 2 / float(pounds_per_dozen), places=4)
            restored = eggs.liquid_usd_per_lb_to_usd_per_dozen(liquid_price)
            self.assertAlmostEqual(float(restored), 2, delta=0.0002)
            if isinstance(price, Decimal):
                self.assertIsInstance(liquid_price, Decimal)

    def test_liquid_density_context_rates_and_quote_router(self):
        eggs = a.PoultryEggAssay.us_large_shell_eggs()
        for magnitude in (2.5, Decimal("2.5")):
            liquid_volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_poultry_egg(liquid_volume, "kg", assay=eggs)
            self.assertAlmostEqual(float(mass.magnitude), 2587.5)
            self.assertAlmostEqual(float(a.convert_poultry_egg(mass, "m^3", assay=eggs).magnitude), 2.5)
            rate = a.convert_poultry_egg(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=eggs)
            self.assertAlmostEqual(float(rate.magnitude), 2587.5 / 24)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/m^3"), "EUR/kg", assay=eggs, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 * 0.9 / 1035)
        broiler = a.PoultryEggAssay.broiler_chicken_standard()
        self.assertFalse(broiler.as_context().funcs)
        with self.assertRaises(ValueError):
            broiler.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(pint.DimensionalityError):
            a.convert_poultry_egg(ureg.Quantity(1, "m^3"), "kg", assay=broiler)

    def test_invalid_specs_masses_counts_categories_and_modes(self):
        eggs = a.PoultryEggAssay.us_large_shell_eggs()
        invalid = {"category": ("invalid", None), "rtc_dressing_pct": (0, 101, float("nan")), "boneless_meat_pct_of_live": (0, 75), "egg_mass_grams": (19, 151), "liquid_egg_yield_pct": (0, 101), "dried_egg_yield_pct": (0, 89), "bulk_density_kg_per_m3": (99, 2001)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(eggs, **{field: value})
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                eggs.eggs_to_shell_mass(bad)
            with self.assertRaises(ValueError):
                eggs.cases_to_shell_mass(bad)
            with self.assertRaises(ValueError):
                eggs.shell_to_liquid_egg_mass(ureg.Quantity(bad, "kg"))
        with self.assertRaises(ValueError):
            eggs.live_to_rtc_mass(ureg.Quantity(1, "kg"))
        with self.assertRaises(ValueError):
            a.PoultryEggAssay.broiler_chicken_standard().dozens_to_shell_mass(1)
        with self.assertRaises(ValueError):
            eggs.usd_per_dozen_to_liquid_usd_per_lb(Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_poultry_egg(1, "kg", assay=eggs)
        with self.assertRaises(ValueError):
            a.convert_poultry_egg(ureg.Quantity(1, "kg"), "kg", assay=eggs, mode="invalid")


if __name__ == "__main__":
    unittest.main()
