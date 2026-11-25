"""Rice milling recovery, bag/bushel units and pricing integration tests."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class RiceAssayTests(unittest.TestCase):
    def test_exports_and_presets(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        namespace = {}
        exec("from core.units.agricultural_specs import *", namespace)
        for name in ("RiceAssay", "create_rice_context", "convert_rice"):
            self.assertIn(name, a.__all__)
            self.assertIs(namespace[name], getattr(a, name))
        for preset in ("thai_white_5pct_broken", "us_cme_rough_rice"):
            self.assertTrue(getattr(a.RiceAssay, preset)().as_context().funcs)

    def test_milling_round_trips_and_output_split(self):
        for assay in (a.RiceAssay.thai_white_5pct_broken(), a.RiceAssay.us_cme_rough_rice()):
            for magnitude in (100.0, Decimal("100")):
                paddy = ureg.Quantity(magnitude, "mt")
                for forward, reverse, expected in (("paddy_to_brown", "brown_to_paddy", assay.brown_rice_yield_pct), ("paddy_to_milled", "milled_to_paddy", assay.milling_yield_pct)):
                    with self.subTest(variety=assay.name, mode=forward, magnitude=magnitude):
                        output = a.convert_rice(paddy, "mt", assay=assay, mode=forward)
                        self.assertAlmostEqual(float(output.magnitude), expected)
                        restored = a.convert_rice(output, "mt", assay=assay, mode=reverse)
                        self.assertAlmostEqual(float(restored.magnitude), 100)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(output.magnitude, Decimal)
                brown = assay.paddy_to_brown_mass(paddy)
                milled = assay.brown_to_milled_mass(brown)
                self.assertAlmostEqual(float(milled.magnitude), assay.milling_yield_pct)
                self.assertAlmostEqual(float(assay.milled_to_brown_mass(milled).magnitude), assay.brown_rice_yield_pct)
                split = assay.milled_to_head_and_broken(milled)
                self.assertAlmostEqual(float(sum(q.magnitude for q in split.values())), float(milled.magnitude))
                self.assertAlmostEqual(float(split["head_rice"].magnitude), 100 * assay.head_rice_fraction_of_paddy)
                if assay.head_rice_pct == 95:
                    self.assertAlmostEqual(float(split["broken_rice"].magnitude), 3.5)
        for head in (0, 100):
            assay = replace(a.RiceAssay.us_cme_rough_rice(), head_rice_pct=head)
            split = assay.milled_to_head_and_broken(ureg.Quantity(Decimal("100"), "mt"))
            self.assertEqual(split["head_rice"].magnitude, head)
            self.assertEqual(split["broken_rice"].magnitude, 100 - head)

    def test_us_bushel_and_export_bag_counts(self):
        assay = a.RiceAssay.us_cme_rough_rice()
        self.assertAlmostEqual(float(assay.mass_to_bushels(ureg.Quantity(45, "lb"))), 1)
        self.assertAlmostEqual(float(assay.bushels_to_mass(1).to("lb").magnitude), 45)
        for count in (100.0, Decimal("100.123456789")):
            for forward, reverse, unit_weight in ((assay.bushels_to_mass, assay.mass_to_bushels, None), (assay.bags_to_mass, assay.mass_to_bags, 50)):
                result = forward(count)
                self.assertAlmostEqual(float(reverse(result)), float(count))
                if unit_weight is not None:
                    self.assertAlmostEqual(float(result.to("kg").magnitude), float(count) * 50)
                if isinstance(count, Decimal):
                    self.assertIsInstance(result.magnitude, Decimal)
                    self.assertIsInstance(reverse(result), Decimal)
        self.assertEqual(assay.bags_to_mass(Decimal("20")).magnitude, 1)

    def test_price_parity_matches_registry(self):
        assay = a.RiceAssay.thai_white_5pct_broken()
        for price in (10.0, Decimal("10")):
            converted = assay.usd_per_cwt_to_usd_per_mt(price)
            expected = ureg.Quantity(price, "USD/commodity_short_hundredweight").to("USD/mt")
            self.assertAlmostEqual(float(converted), float(expected.magnitude), places=4)
            self.assertAlmostEqual(float(assay.usd_per_mt_to_usd_per_cwt(converted)), 10, places=4)
            if isinstance(price, Decimal):
                self.assertIsInstance(converted, Decimal)
        self.assertEqual(assay.paddy_price_to_milled_parity_per_mt(Decimal("350")), 500)

    def test_bulk_context_master_quotes_and_missing_density(self):
        assay = a.RiceAssay.us_cme_rough_rice()
        for magnitude in (2.5, Decimal("2.5")):
            volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_rice(volume, "kg", assay=assay)
            self.assertEqual(mass.magnitude, 1450)
            self.assertAlmostEqual(float(a.convert_rice(mass, "m^3", assay=assay).magnitude), 2.5)
            rate = a.convert_rice(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay)
            self.assertAlmostEqual(float(rate.magnitude), 1450 / 24)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/m^3"), "EUR/kg", assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 * 0.9 / 580)
        missing = replace(assay, bulk_density_kg_per_m3=None)
        self.assertFalse(missing.as_context().funcs)
        with self.assertRaises(ValueError):
            missing.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(pint.DimensionalityError):
            a.convert_rice(ureg.Quantity(1, "m^3"), "kg", assay=missing)

    def test_invalid_specs_inputs_and_modes(self):
        assay = a.RiceAssay.us_cme_rough_rice()
        invalid = {"variety": ("", " ", None), "milling_yield_pct": (0, 101, Decimal("NaN")), "head_rice_pct": (-1, 101), "brown_rice_yield_pct": (0, 67, 101), "bushel_weight_lb": (0, float("inf")), "bag_mass_kg": (0, -1), "bulk_density_kg_per_m3": (99, 2001)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.paddy_to_milled_mass(ureg.Quantity(bad, "kg"))
            with self.assertRaises(ValueError):
                assay.bushels_to_mass(bad)
            with self.assertRaises(ValueError):
                assay.bags_to_mass(bad)
        with self.assertRaises(ValueError):
            assay.usd_per_cwt_to_usd_per_mt(Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_rice(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_rice(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()
