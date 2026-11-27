"""Source-aware parsing must resolve original tokens before compatibility defaults."""

from concurrent.futures import ThreadPoolExecutor
from dataclasses import FrozenInstanceError, replace
from datetime import date
from decimal import Decimal
import unittest

import pint

import core.units as u
from core.units import currencies, registry


class StrictQuoteParsingTests(unittest.TestCase):
    def test_every_requested_ambiguous_physical_token_is_rejected(self):
        for token in ("MT", "M", "MM", "m", "MB", "gal", "oz", "bu", "GT", "lot"):
            with self.subTest(token=token), self.assertRaisesRegex(ValueError, "Ambiguous quote token"):
                u.parse_quote(f"25 USD/{token}")
        for token in ("c", "p", "$", "¥"):
            with self.subTest(token=token), self.assertRaisesRegex(ValueError, "Ambiguous"):
                u.parse_quote(f"25 {token}/kilogram")

    def test_source_labels_alone_never_guess_a_physical_convention(self):
        context = u.QuoteParsingContext(commodity_id="corn", currency="USD", exchange="CBOT", vendor="feed", region="US")
        for token in ("gal", "oz", "bu", "MT", "lot"):
            with self.subTest(token=token), self.assertRaises(ValueError):
                u.parse_quote(f"1 USD/{token}", context=context)

    def test_explicit_source_rule_resolves_each_ambiguous_unit(self):
        cases = (("MT", "tonne"), ("M", "million_metric_tonnes"), ("MM", "MMbbl"),
                 ("m", "meter"), ("MB", "kbbl"), ("gal", "imperial_gallon"),
                 ("oz", "toz"), ("bu", "bushel_corn"), ("GT", "gross_tonnage"),
                 ("lot", "kbbl"))
        for token, target in cases:
            with self.subTest(token=token):
                rule = u.QuoteTokenRule(token, target, vendor="example-feed", commodity_id="CONTRACT-1")
                context = u.QuoteParsingContext(vendor="example-feed", commodity_id="CONTRACT-1", token_rules=(rule,))
                result = u.parse_quote(f"10 USD/{token}", context=context)
                expected = u.Q_(Decimal(10), u.ureg.USD / u.ureg.Unit(target))
                self.assertEqual(result, expected)
        # Roman multipliers are explicit named units, not inferred SI prefixes.
        for multiplier, expected in (("kilocount", ".1"), ("megacount", ".0001")):
            context = u.QuoteParsingContext(token_rules=(u.QuoteTokenRule("M", multiplier),))
            result = u.parse_quote("100 USD/(M*oil_barrel)", context=context)
            self.assertEqual(result.to("USD/oil_barrel").magnitude, Decimal(expected))

    def test_currency_context_resolves_shared_symbols_and_minor_denominations(self):
        for token, code, target in (("$", "CAD", "CAD"), ("¥", "CNY", "CNY"),
                                    ("¥", "JPY", "JPY"), ("£", "GBP", "GBP"),
                                    ("c", "CAD", "CAD_cent"), ("p", "GBP", "GBP_penny")):
            with self.subTest(token=token, code=code):
                context = u.QuoteParsingContext(currency=code)
                result = u.parse_quote(f"125 {token}/kilogram", context=context)
                self.assertEqual(result.units, u.ureg.Unit(target) / u.ureg.kilogram)
        with self.assertRaises(ValueError):
            u.parse_quote("1 p/kilogram", context=u.QuoteParsingContext(currency="USD"))
        rules = (u.QuoteTokenRule("$", "CAD", vendor="feed"),)
        result = u.parse_quote("$25/kilogram", context=u.QuoteParsingContext(vendor="feed", token_rules=rules))
        self.assertEqual(result.units, u.ureg.CAD / u.ureg.kilogram)

    def test_scope_requires_all_source_selectors_to_match(self):
        rule = u.QuoteTokenRule("bu", "bushel_corn", commodity_id="corn", currency="USD",
                               exchange="CBOT", vendor="feed-a", region="US")
        context = u.QuoteParsingContext(commodity_id="corn", currency="USD", exchange="CBOT",
                                        vendor="feed-a", region="US", token_rules=(rule,))
        self.assertEqual(u.parse_quote("5 c/bu", context=context).units, u.ureg.USD_cent / u.ureg.bushel_corn)
        for name, value in (("commodity_id", "wheat"), ("currency", "CAD"), ("exchange", "OTHER"),
                            ("vendor", "feed-b"), ("region", "UK"), ("vendor", None)):
            with self.subTest(name=name), self.assertRaisesRegex(ValueError, "Ambiguous quote token 'bu'"):
                u.parse_quote("5 c/bu", context=replace(context, **{name: value}))

    def test_dated_rules_use_half_open_intervals_and_require_a_date(self):
        old = u.QuoteTokenRule("gal", "gallon_petroleum", vendor="feed", valid_from="2020-01-01", valid_to="2021-01-01")
        new = u.QuoteTokenRule("gal", "imperial_gallon", vendor="feed", valid_from="2021-01-01")
        context = u.QuoteParsingContext(vendor="feed", token_rules=(old, new))
        with self.assertRaisesRegex(ValueError, "reference_date"):
            u.parse_quote("1 EUR/gal", context=context)
        for when, target in (("2020-01-01", "gallon_petroleum"), ("2020-12-31", "gallon_petroleum"), ("2021-01-01", "imperial_gallon")):
            with self.subTest(when=when):
                result = u.parse_quote("1 EUR/gal", context=replace(context, reference_date=when))
                self.assertEqual(result.units, u.ureg.EUR / u.ureg.Unit(target))
        with self.assertRaisesRegex(ValueError, "Ambiguous"):
            u.parse_quote("1 EUR/gal", context=replace(context, reference_date="2019-12-31"))

    def test_overlapping_conflicts_fail_without_rule_order_priority(self):
        first = u.QuoteTokenRule("oz", "toz", vendor="feed")
        second = u.QuoteTokenRule("oz", "ounce", commodity_id="gold")
        for rules in ((first, second), (second, first)):
            with self.subTest(rules=rules), self.assertRaisesRegex(ValueError, "Conflicting source"):
                u.parse_quote("1 USD/oz", context=u.QuoteParsingContext(vendor="feed", commodity_id="gold", token_rules=rules))
        compatible = replace(second, unit="troy_ounce")
        result = u.parse_quote("1 USD/oz", context=u.QuoteParsingContext(vendor="feed", commodity_id="gold", token_rules=(first, compatible)))
        self.assertEqual(result.units, u.ureg.USD / u.ureg.toz)

    def test_strict_path_bypasses_primed_legacy_quote_cache(self):
        for expression in ("$/MT", "$/oz", "c/bu", "p/th"):
            u.parse_unit(expression)
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                u.parse_quote_unit(expression)
        self.assertEqual(u.parse_unit("$/oz"), u.ureg.USD / u.ureg.toz)
        self.assertEqual(u.parse_unit("c/bu"), u.ureg.USD_cent / u.ureg.bushel)

    def test_all_parser_entry_points_share_context_and_decimal_policy(self):
        context = u.QuoteParsingContext(currency="USD", vendor="feed", token_rules=(u.QuoteTokenRule("gal", "imperial_gallon", vendor="feed"),))
        expected = u.ureg.USD / u.ureg.imperial_gallon
        self.assertEqual(u.parse_unit("$/gal", context=context), expected)
        self.assertEqual(u.parse_quote_unit("$/gal", context=context), expected)
        self.assertEqual(u.parse_quantity("$1,234.50/gal", context=context).magnitude, Decimal("1234.50"))
        self.assertEqual(u.parse_quote("($1,234.50/gal)", context=context).magnitude, Decimal("-1234.50"))
        self.assertEqual(u.parse_quote("-$.50/gal", context=context).magnitude, Decimal("-.50"))
        self.assertEqual(u.parse_quote("1.5e3 USD per gal", context=context).magnitude, Decimal(1500))
        for as_decimal, kind in ((True, Decimal), (False, float)):
            batches = (u.parse_quantities(["1 USD/gal", "2 USD/gal"], context=context, as_decimal=as_decimal),
                       u.parse_quantities([1, 2], "USD/gal", context=context, as_decimal=as_decimal),
                       u.parse_quantities([1, 2], ["USD/gal", "USD/gal"], context=context, as_decimal=as_decimal))
            for quantities in batches:
                self.assertTrue(all(isinstance(quantity.magnitude, kind) and quantity.units == expected for quantity in quantities))
        components = u.decompose_quote("$/gal", context=context)
        self.assertEqual(components.commodity_unit, "imperial_gallon")

    def test_context_currency_and_default_currency_cannot_disagree(self):
        context = u.QuoteParsingContext(currency="CAD")
        for expression in ("USD/kilogram", "$1 USD/kilogram", "US$/kilogram"):
            with self.subTest(expression=expression), self.assertRaises(ValueError):
                u.parse_quote(expression, context=context)
        with self.assertRaisesRegex(ValueError, "default_currency conflicts"):
            u.parse_unit("$/kilogram", context=context, default_currency="USD")
        self.assertEqual(u.parse_unit("$/kilogram", context=u.QuoteParsingContext(), default_currency="CAD"), u.ureg.CAD / u.ureg.kilogram)

    def test_currency_metadata_is_resolved_at_the_reference_date(self):
        context = u.QuoteParsingContext(currency="BGN", reference_date="2025-12-31")
        self.assertEqual(u.parse_quote("5 BGN/kilogram", context=context).magnitude, 5)
        for when in (None, "2026-01-01", "2027-01-01"):
            with self.subTest(when=when), self.assertRaises(ValueError):
                u.parse_quote("5 BGN/kilogram", context=replace(context, reference_date=when))
        with self.assertRaisesRegex(ValueError, "Insufficient exact validity"):
            u.parse_quote("5 USD/kilogram", context=u.QuoteParsingContext(reference_date="2020-01-01"))
        usd = replace(u.resolve_currency("USD"), valid_from=date(2000, 1, 1))
        catalog = u.CurrencyCatalog((usd,), version="test overlay", observed_on=currencies.DATA_VERSION)
        self.assertEqual(u.parse_quote("5 USD/kilogram", context=u.QuoteParsingContext(reference_date="2020-01-01", currency_catalog=catalog)).magnitude, 5)

    def test_vendor_minor_alias_retains_scale_and_requires_its_namespace(self):
        alias = u.CurrencyAlias("GBp", "GBP", "GBP_penny", Decimal(".01"), vendor="feed")
        catalog = u.CurrencyCatalog(currencies.CURRENCY_CATALOG.records, version="feed-1", observed_on=currencies.DATA_VERSION, aliases=(alias,))
        context = u.QuoteParsingContext(vendor="feed", currency_catalog=catalog)
        result = u.parse_quote("125 GBp/kilogram", context=context)
        self.assertEqual(result.to("GBP/kilogram").magnitude, Decimal("1.25"))
        for vendor in (None, "other"):
            with self.subTest(vendor=vendor), self.assertRaises(ValueError):
                u.parse_quote("125 GBp/kilogram", context=replace(context, vendor=vendor))

    def test_token_boundaries_and_algebra_preserve_canonical_units(self):
        context = u.QuoteParsingContext(token_rules=(u.QuoteTokenRule("M", "meter"),))
        for expression in ("USD/MMBtu_it", "USD/MWh", "USD/meter**3", "USD/kilogram/day", "USD/(kilogram/day)"):
            with self.subTest(expression=expression):
                self.assertEqual(u.parse_quote_unit(expression, context=context), u.ureg.Unit(expression))
        self.assertEqual(u.parse_quote_unit("USD/m^3", context=u.QuoteParsingContext(token_rules=(u.QuoteTokenRule("m", "meter"),))), u.ureg.USD / u.ureg.meter**3)
        rule = u.QuoteTokenRule("bu", "bushel_corn")
        self.assertEqual(u.parse_quote_unit("USD/(bu/day)", context=u.QuoteParsingContext(token_rules=(rule,))), u.ureg.USD * u.ureg.day / u.ureg.bushel_corn)
        inverse = u.QuoteParsingContext(token_rules=(u.QuoteTokenRule("m", "meter**-1"),))
        self.assertEqual(u.parse_quote_unit("USD/m", context=inverse), u.ureg.USD * u.ureg.meter)

    def test_parallel_source_profiles_cannot_contaminate_each_other_or_legacy(self):
        contexts = tuple(u.QuoteParsingContext(vendor=name, token_rules=(u.QuoteTokenRule("gal", unit, vendor=name),))
                         for name, unit in (("us-feed", "gallon_petroleum"), ("uk-feed", "imperial_gallon")))
        original = dict(registry._QUOTE_REPLACEMENTS)
        def parse(index):
            return u.parse_quote_unit("USD/gal", context=contexts[index % 2])
        with ThreadPoolExecutor(max_workers=4) as pool:
            parsed = list(pool.map(parse, range(24)))
        for index, result in enumerate(parsed):
            self.assertEqual(result, u.ureg.USD / u.ureg.Unit("gallon_petroleum" if index % 2 == 0 else "imperial_gallon"))
        self.assertEqual(dict(registry._QUOTE_REPLACEMENTS), original)

    def test_strict_flag_explicitly_enables_compatibility_fallback(self):
        context = u.QuoteParsingContext(strict=False)
        for expression in ("$/oz", "$/gal", "p/th", "c/bu"):
            with self.subTest(expression=expression):
                self.assertEqual(u.parse_quote_unit(expression, context=context), u.parse_unit(expression))
        explicit = replace(context, token_rules=(u.QuoteTokenRule("oz", "ounce"),))
        self.assertEqual(u.parse_quote_unit("USD/oz", context=explicit), u.ureg.USD / u.ureg.ounce)

    def test_unusable_quote_shapes_and_unsafe_syntax_fail(self):
        for expression in ("USD", "kilogram", "USD/EUR", "USD/dimensionless", "USD/kilogram + meter",
                           "USD/__import__('os')", "USD/meter.real", "USD/meter**1.5", "USD/(meter", "USD/[meter]"):
            with self.subTest(expression=expression), self.assertRaises((ValueError, pint.PintError)):
                u.parse_quote_unit(expression)
        with self.assertRaises(TypeError):
            u.parse_unit(u.ureg.USD / u.ureg.kilogram, context=u.QuoteParsingContext())
        for value in (1, "1", True, float("nan")):
            with self.subTest(value=value), self.assertRaises((ValueError, TypeError)):
                u.parse_quote(value)

    def test_context_and_rules_are_validated_and_immutable(self):
        context = u.QuoteParsingContext(token_rules=[u.QuoteTokenRule("MT", "tonne")])
        self.assertIsInstance(context.token_rules, tuple)
        with self.assertRaises(FrozenInstanceError):
            context.strict = False
        with self.assertRaises(FrozenInstanceError):
            context.token_rules[0].unit = "million_metric_tonnes"
        for values in (dict(strict="true"), dict(vendor=" feed "), dict(token_rules=("MT",)), dict(currency="$")):
            with self.subTest(values=values), self.assertRaises((TypeError, ValueError)):
                u.QuoteParsingContext(**values)
        with self.assertRaisesRegex(ValueError, "contract-qualified"):
            u.QuoteTokenRule("lot", "tonne")
        with self.assertRaises(ValueError):
            u.QuoteTokenRule("MT", "tonne", valid_from="2021-01-01", valid_to="2020-01-01")


if __name__ == "__main__":
    unittest.main()
