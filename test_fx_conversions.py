"""FX minor-unit scaling, quote preservation and sparse-rate regressions."""

from decimal import Decimal
import unittest

from core.units import conversions as c
from core.units.registry import ureg


class FxConversionTests(unittest.TestCase):
    def test_context_binds_all_pairwise_rates(self):
        rates = {"USD": Decimal("1"), "EUR": Decimal("0.8"), "GBP": Decimal("0.6"), "JPY": Decimal("150")}
        context = c.create_fx_context(rates)
        # Exercise every callback after all registration loops have completed.
        for source, source_rate in rates.items():
            for target, target_rate in rates.items():
                for magnitude in (10.0, Decimal("10")):
                    with self.subTest(source=source, target=target, magnitude=magnitude):
                        result = ureg.Quantity(magnitude, source).to(target, context)
                        expected = Decimal("10") * target_rate / source_rate
                        self.assertAlmostEqual(float(result.magnitude), float(expected))
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(result.magnitude, Decimal)

    def test_settlement_float_multiplication_and_half_even_boundaries(self):
        # Binary float multiplication yields 0.5250000000000001 here.
        result = c.convert_currency(ureg.Quantity(0.07, "USD"), "EUR", rates={"EUR": 7.5}, quantize_settlement=True)
        self.assertEqual(result.magnitude, Decimal("0.52"))
        for amount, expected in ((0.125, "0.12"), (0.135, "0.14"), (-0.125, "-0.12"), (0.12500000000000002, "0.13")):
            with self.subTest(amount=amount):
                result = c.convert_currency(ureg.Quantity(amount, "GBP"), "GBP", rates={}, base_currency="GBP", quantize_settlement=True)
                self.assertEqual(result.magnitude, Decimal(expected))
        ureg.define("section6_two_pence = 0.02 * GBP")
        result = c.convert_currency(ureg.Quantity(Decimal("1.255"), "GBP"), "section6_two_pence", rates={}, base_currency="GBP", quantize_settlement=True)
        self.assertEqual(result.magnitude, Decimal("63"))
        result = c.convert_currency(ureg.Quantity(Decimal("1.245"), "GBP"), "section6_two_pence", rates={}, base_currency="GBP", quantize_settlement=True)
        self.assertEqual(result.magnitude, Decimal("62"))

    def test_nonfinite_amounts_and_currency_powers(self):
        for amount in (float("nan"), float("inf"), -float("inf"), Decimal("NaN"), Decimal("Infinity")):
            with self.subTest(amount=amount):
                with self.assertRaisesRegex(ValueError, "FX amount must be finite"):
                    c.convert_fx(amount, "EUR", {"USD/EUR": 0.9}, source_currency="USD")
                for convert in (
                    lambda: c.convert_fx(ureg.Quantity(amount, "USD / oil_barrel"), "EUR", {"USD/EUR": 0.9}),
                    lambda: c.convert_currency(ureg.Quantity(amount, "USD"), "EUR", rates={"EUR": 0.9}),
                ):
                    with self.assertRaisesRegex(ValueError, "FX amount must be finite"):
                        convert()
        for units in ("USD / USD", "1 / USD", "USD ** 2"):
            with self.subTest(units=units), self.assertRaises(ValueError):
                c.convert_fx(ureg.Quantity(100, units), "EUR", {"USD/EUR": 0.9})
        quantity = ureg.Quantity(Decimal("2"), "kUSD / oil_barrel")
        result = c.convert_fx(quantity, "EUR", {"USD/EUR": Decimal("0.9")})
        self.assertEqual(result.magnitude, Decimal("1800"))
        self.assertEqual(result.units, ureg.EUR / ureg.oil_barrel)

    def test_minor_source_units_preserve_denominators(self):
        for magnitude in (125.0, Decimal("125")):
            for denominator in ("", " / therm_it", " / oil_barrel", " / metric_ton"):
                with self.subTest(magnitude=magnitude, denominator=denominator):
                    quantity = ureg.Quantity(magnitude, "GBP_pence" + denominator)
                    expected_unit = ureg.Unit("USD" + denominator)
                    anchor = c.convert_currency(quantity, "USD", rates={"GBP": Decimal("0.8")})
                    sparse = c.convert_fx(quantity, "USD", {"GBP/USD": Decimal("1.25")})
                    for result in (anchor, sparse):
                        self.assertEqual(result.units, expected_unit)
                        self.assertAlmostEqual(float(result.magnitude), 1.5625)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(result.magnitude, Decimal)
        cents = c.convert_fx(ureg.Quantity(Decimal("100"), "US_cent / gallon"), "EUR", {"USD/EUR": Decimal("0.9")})
        self.assertEqual(cents.units, ureg.EUR / ureg.gallon)
        self.assertEqual(cents.magnitude, Decimal("0.9"))

    def test_minor_targets_and_same_currency_scaling(self):
        for source, target, expected in (("GBP", "GBP_pence", 100), ("GBP_pence", "GBP", 0.01), ("GBP_pence", "GBP_pence", 1)):
            quantity = ureg.Quantity(Decimal("1"), source + " / therm_it")
            for result in (c.convert_currency(quantity, target, rates={}, base_currency="GBP"), c.convert_fx(quantity, target, {})):
                self.assertEqual(result.units, ureg.Unit(target + " / therm_it"))
                self.assertEqual(result.magnitude, Decimal(str(expected)))
        result = c.convert_fx(Decimal("1"), "GBP_pence", {"USD/GBP": Decimal("0.8")}, source_currency="USD")
        self.assertEqual(result.magnitude, 80)
        self.assertEqual(result.units, ureg.GBP_pence)

    def test_price_quotes_and_multihop_rates(self):
        quote = ureg.Quantity(Decimal("75.50"), "USD / oil_barrel")
        result = c.convert_fx(quote, "GBP", {"USD/EUR": Decimal("0.9"), ("EUR", "GBP"): Decimal("0.8")})
        self.assertEqual(result.units, ureg.GBP / ureg.oil_barrel)
        self.assertEqual(result.magnitude, Decimal("54.3600"))
        self.assertEqual(c.triangulate_fx_rate({"USD/EUR": Decimal("0.8")}, "EUR", "USD"), Decimal("1.25"))
        self.assertEqual(c.convert_fx(10.0, "EUR", {"USD/EUR": 0.9}, source_currency="USD").magnitude, 9)

    def test_settlement_rounding_uses_target_scale(self):
        for target, expected in (("GBP", Decimal("1.24")), ("GBP_pence", Decimal("124"))):
            result = c.convert_currency(ureg.Quantity(Decimal("1.245"), "GBP"), target, rates={}, base_currency="GBP", quantize_settlement=True)
            self.assertEqual(result.magnitude, expected)
        result = c.convert_currency(ureg.Quantity(Decimal("1.255"), "GBP"), "GBP_pence", rates={}, base_currency="GBP", quantize_settlement=True)
        self.assertEqual(result.magnitude, 126)

    def test_invalid_pairs_rates_and_currency_dimensions(self):
        for key in ("EURUSD", "EU", "EUR/USD/GBP", "EUR/", "/USD", ("EUR",), ("EUR", "USD", "GBP"), 42):
            with self.subTest(key=key), self.assertRaises(ValueError):
                c.triangulate_fx_rate({key: 1}, "EUR", "USD")
        for rate in (0, -1, float("nan"), float("inf")):
            with self.subTest(rate=rate), self.assertRaises(ValueError):
                c.triangulate_fx_rate({"USD/EUR": rate}, "USD", "EUR")
        with self.assertRaises(ValueError):
            c.convert_fx(ureg.Quantity(1, "USD"), "EUR", {"USD/EUR": 0.9}, source_currency="GBP")
        for unit in ("kg", "1 / USD", "USD * EUR", "USD / EUR"):
            with self.subTest(unit=unit), self.assertRaises(ValueError):
                c.convert_fx(ureg.Quantity(1, unit), "EUR", {"USD/EUR": 0.9})
        with self.assertRaises(ValueError):
            c.convert_fx(1, "EUR", {"USD/EUR": 0.9})
        with self.assertRaises(ValueError):
            c.convert_fx(ureg.Quantity(1, "USD"), "EUR / kg", {"USD/EUR": 0.9})


if __name__ == "__main__":
    unittest.main()
