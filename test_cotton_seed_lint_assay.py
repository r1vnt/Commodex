"""Cotton ginning, lint bales, seed density and netback regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class CottonSeedLintAssayTests(unittest.TestCase):
    def test_ginning_mass_balance_and_inverse_yields(self):
        for assay in (a.CottonSeedLintAssay.us_upland_480lb(), a.CottonSeedLintAssay.indian_shankar6_170kg(), a.CottonSeedLintAssay.brazil_australian_227kg()):
            for magnitude in (100.0, Decimal("100")):
                raw = ureg.Quantity(magnitude, "mt")
                split = assay.gin_seed_cotton(raw)
                self.assertAlmostEqual(float(split["lint"].magnitude), assay.lint_turnout_pct)
                self.assertAlmostEqual(float(split["cottonseed"].magnitude), assay.cottonseed_turnout_pct)
                self.assertAlmostEqual(float(split["trash"].magnitude), assay.gin_trash_loss_pct)
                self.assertAlmostEqual(float(sum(q.magnitude for q in split.values())), 100)
                for forward, reverse in (("seed_cotton_to_lint", "lint_to_seed_cotton"), ("seed_cotton_to_seed", "seed_to_seed_cotton")):
                    with self.subTest(origin=assay.origin, magnitude=magnitude, mode=forward):
                        processed = a.convert_cotton_seed_lint(raw, "kg", assay=assay, mode=forward)
                        restored = a.convert_cotton_seed_lint(processed, "mt", assay=assay, mode=reverse)
                        self.assertAlmostEqual(float(restored.magnitude), 100)
                if isinstance(magnitude, Decimal):
                    for output in split.values():
                        self.assertIsInstance(output.magnitude, Decimal)
        no_trash = replace(a.CottonSeedLintAssay.us_upland_480lb(), cottonseed_turnout_pct=60)
        self.assertEqual(no_trash.gin_seed_cotton(ureg.Quantity(Decimal("1"), "mt"))["trash"].magnitude, 0)

    def test_bale_standards_and_reverse_raw_cotton_orders(self):
        us = a.CottonSeedLintAssay.us_upland_480lb()
        self.assertAlmostEqual(float(us.bales_to_lint_mass(Decimal("1")).to("lb").magnitude), 480)
        self.assertAlmostEqual(float(us.bales_to_lint_mass(Decimal("1")).to("kg").magnitude), ureg.Quantity(1, "cotton_bale").to("kg").magnitude)
        for assay in (us, a.CottonSeedLintAssay.indian_shankar6_170kg(), a.CottonSeedLintAssay.brazil_australian_227kg()):
            for count in (10.0, Decimal("10.123456789")):
                lint = assay.bales_to_lint_mass(count)
                self.assertAlmostEqual(float(assay.lint_mass_to_bales(lint)), float(count))
                raw = assay.bales_to_seed_cotton_mass(count)
                self.assertAlmostEqual(float(assay.seed_cotton_to_lint_mass(raw).magnitude), float(lint.magnitude))
                if isinstance(count, Decimal):
                    self.assertIsInstance(lint.magnitude, Decimal)
                    self.assertIsInstance(assay.lint_mass_to_bales(lint), Decimal)

    def test_price_parity_and_gross_ginning_value(self):
        assay = a.CottonSeedLintAssay.us_upland_480lb()
        for price in (80.0, Decimal("80")):
            per_tonne = assay.lint_cents_per_lb_to_usd_per_mt(price)
            expected = ureg.Quantity(price, "US_cent/lb").to("USD/mt")
            self.assertAlmostEqual(float(per_tonne), float(expected.magnitude), places=4)
            self.assertAlmostEqual(float(assay.lint_usd_per_mt_to_cents_per_lb(per_tonne)), 80, places=4)
            per_bale = assay.lint_cents_per_lb_to_usd_per_bale(price)
            self.assertEqual(per_bale, 384)
            self.assertAlmostEqual(float(assay.lint_usd_per_bale_to_cents_per_lb(per_bale)), 80)
            if isinstance(price, Decimal):
                self.assertIsInstance(per_tonne, Decimal)
        self.assertEqual(assay.gross_ginning_value_usd_per_mt_seed_cotton(Decimal("2000"), Decimal("300")), 962)

    def test_bulk_seed_context_rates_and_quotes(self):
        assay = a.CottonSeedLintAssay.us_upland_480lb()
        for magnitude in (2.5, Decimal("2.5")):
            mass = a.convert_cotton_seed_lint(ureg.Quantity(magnitude, "m^3"), "kg", assay=assay)
            self.assertAlmostEqual(float(mass.magnitude), 1050)
            back = a.convert_cotton_seed_lint(mass, "m^3", assay=assay)
            self.assertAlmostEqual(float(back.magnitude), 2.5)
            rate = a.convert_cotton_seed_lint(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay)
            self.assertAlmostEqual(float(rate.magnitude), 1050 / 24)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/m^3"), "EUR/kg", assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 * 0.9 / 420)

    def test_validation_and_zero_seed_turnout(self):
        assay = a.CottonSeedLintAssay.us_upland_480lb()
        invalid = {"origin": ("", None), "lint_turnout_pct": (0, 101, float("nan")), "cottonseed_turnout_pct": (-1, 61), "bale_mass_kg": (0, -1), "seed_bulk_density_kg_per_m3": (199, 1001)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.gin_seed_cotton(ureg.Quantity(bad, "kg"))
            with self.assertRaises(ValueError):
                assay.bales_to_lint_mass(bad)
        zero = replace(assay, cottonseed_turnout_pct=0)
        self.assertEqual(zero.seed_cotton_to_cottonseed_mass(ureg.Quantity(1, "mt")).magnitude, 0)
        with self.assertRaises(ValueError):
            zero.cottonseed_to_seed_cotton_mass(ureg.Quantity(1, "mt"))
        with self.assertRaises(ValueError):
            assay.lint_cents_per_lb_to_usd_per_mt(Decimal("NaN"))
        with self.assertRaises(ValueError):
            a.convert_cotton_seed_lint(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()
