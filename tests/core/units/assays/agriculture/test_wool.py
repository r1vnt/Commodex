"""Wool clean/top yields, nominal bales and same-currency price parity tests."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class WoolAssayTests(unittest.TestCase):
    def test_exports_presets_and_processing_round_trips(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("WoolAssay", "create_wool_context", "convert_wool"):
            self.assertIn(name, a.__all__)
        for preset in ("australian_merino_19micron", "medium_merino_21micron", "crossbred_28micron"):
            assay = getattr(a.WoolAssay, preset)()
            for magnitude in (100.0, Decimal("100")):
                for forward, reverse in (("greasy_to_clean", "clean_to_greasy"), ("clean_to_top", "top_to_clean"), ("greasy_to_top", "top_to_greasy")):
                    with self.subTest(preset=preset, mode=forward, magnitude=magnitude):
                        original = ureg.Quantity(magnitude, "kg")
                        changed = a.convert_wool(original, "lb", assay=assay, mode=forward)
                        restored = a.convert_wool(changed, "kg", assay=assay, mode=reverse)
                        self.assertAlmostEqual(float(restored.magnitude), 100)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(changed.magnitude, Decimal)
                clean = assay.greasy_to_clean_mass(ureg.Quantity(magnitude, "kg"))
                self.assertAlmostEqual(float(clean.magnitude), assay.clean_yield_pct)
                top = assay.greasy_to_top_mass(ureg.Quantity(magnitude, "kg"))
                self.assertAlmostEqual(float(top.magnitude), 100 * assay.greasy_to_top_fraction)

    def test_lot_override_and_exact_decimal_yield(self):
        assay = a.WoolAssay.australian_merino_19micron()
        original = ureg.Quantity(Decimal("100"), "kg")
        percentage = Decimal("65.123456789123456789")
        result = assay.greasy_to_clean_mass(original, clean_yield_pct=percentage)
        self.assertEqual(result.magnitude, percentage)
        restored = assay.clean_to_greasy_mass(result, clean_yield_pct=percentage)
        self.assertAlmostEqual(float(restored.magnitude), 100)
        converted = a.convert_wool(original, "kg", assay=assay, mode="greasy_to_top", clean_yield_pct=50)
        self.assertEqual(converted.magnitude, Decimal("45.5"))
        self.assertEqual(a.convert_wool(original, "kg", assay=assay).magnitude, 100)

    def test_bale_counts_and_custom_dumped_bales(self):
        for assay in (a.WoolAssay.australian_merino_19micron(), replace(a.WoolAssay.australian_merino_19micron(), bale_mass_kg=250)):
            for count in (10.0, Decimal("10.123456789")):
                mass = assay.bales_to_mass(count)
                self.assertAlmostEqual(float(mass.magnitude), float(count) * assay.bale_mass_kg)
                self.assertAlmostEqual(float(assay.mass_to_bales(mass)), float(count))
                if isinstance(count, Decimal):
                    self.assertIsInstance(mass.magnitude, Decimal)
                    self.assertIsInstance(assay.mass_to_bales(mass), Decimal)
        self.assertEqual(a.WoolAssay.australian_merino_19micron().bales_to_mass(Decimal("10")).magnitude, 1750)

    def test_clean_greasy_bale_and_cent_price_parity(self):
        assay = a.WoolAssay.australian_merino_19micron()
        for price in (10.0, Decimal("10")):
            greasy = assay.clean_price_to_greasy_price(price)
            self.assertAlmostEqual(float(greasy), 6.6)
            self.assertAlmostEqual(float(assay.greasy_price_to_clean_price(greasy)), 10)
            bale = assay.greasy_price_per_kg_to_bale_price(greasy)
            self.assertAlmostEqual(float(bale), 1155)
            self.assertAlmostEqual(float(assay.bale_price_to_greasy_price_per_kg(bale)), 6.6)
            if isinstance(price, Decimal):
                self.assertIsInstance(greasy, Decimal)
                self.assertIsInstance(bale, Decimal)
        self.assertEqual(assay.clean_cents_per_kg_to_greasy_price(Decimal("1000")), Decimal("6.6"))
        self.assertEqual(assay.greasy_price_to_clean_cents_per_kg(Decimal("6.6")), 1000)
        expected = ureg.Quantity(Decimal("100"), "US_cent/lb").to("USD/kg").magnitude
        converted = assay.cents_per_lb_to_price_per_kg(Decimal("100"))
        self.assertAlmostEqual(float(converted), float(expected), places=4)
        reverse = assay.price_per_kg_to_cents_per_lb(converted)
        self.assertAlmostEqual(float(reverse), 100, delta=0.002)

    def test_bulk_context_master_rates_and_missing_density(self):
        assay = a.WoolAssay.australian_merino_19micron()
        for magnitude in (2.5, Decimal("2.5")):
            volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_wool(volume, "kg", assay=assay)
            self.assertEqual(mass.magnitude, 800)
            self.assertAlmostEqual(float(a.convert_wool(mass, "m^3", assay=assay).magnitude), 2.5)
            rate = a.convert_wool(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay)
            self.assertAlmostEqual(float(rate.magnitude), 800 / 24)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "AUD/m^3"), "USD/kg", assay=assay, rates={"AUD": Decimal("1.5")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 / 1.5 / 320)
        missing = replace(assay, bulk_density_kg_per_m3=None)
        self.assertFalse(missing.as_context().funcs)
        with self.assertRaises(ValueError):
            missing.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(pint.DimensionalityError):
            a.convert_wool(ureg.Quantity(1, "m^3"), "kg", assay=missing)

    def test_invalid_specs_yields_values_and_modes(self):
        assay = a.WoolAssay.australian_merino_19micron()
        invalid = {"wool_type": ("", None), "clean_yield_pct": (0, 101, Decimal("NaN")), "combing_top_yield_pct": (0, 101), "micron": (9, 61, float("inf")), "bale_mass_kg": (0, -1), "bulk_density_kg_per_m3": (99, 2001)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for bad in (0, 101, Decimal("NaN"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.greasy_to_clean_mass(ureg.Quantity(1, "kg"), clean_yield_pct=bad)
            with self.assertRaises(ValueError):
                assay.clean_price_to_greasy_price(10, clean_yield_pct=bad)
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.bales_to_mass(bad)
            with self.assertRaises(ValueError):
                assay.clean_to_top_mass(ureg.Quantity(bad, "kg"))
        with self.assertRaises(ValueError):
            assay.clean_price_to_greasy_price(Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_wool(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_wool(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()
