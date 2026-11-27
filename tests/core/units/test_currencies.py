"""Reference provenance, historical denominations, settlement and FX contracts."""

from dataclasses import FrozenInstanceError, replace
from datetime import datetime, timedelta, timezone
from decimal import Decimal
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import unittest

import pint

import core.units as units
from core.units import currencies
from tests.core.units.reference_data.currencies._test_build import build


class CurrencyReferenceTests(unittest.TestCase):
    def test_raw_and_derived_data_rebuild_reproducibly(self):
        directory = currencies.DATA_DIRECTORY
        with tempfile.TemporaryDirectory() as temporary:
            destination = Path(temporary)
            for filename in ("list-one.xml", "list-three.xml", "overlay.json", "conventions.json"):
                shutil.copyfile(directory / filename, destination / filename)
            build(destination)
            for filename in ("reference.json", "manifest.json"):
                self.assertEqual((destination / filename).read_bytes(), (directory / filename).read_bytes())
            catalog = currencies.load_currency_catalog(destination)
            self.assertEqual(catalog.records, currencies.CURRENCY_CATALOG.records)
            # Changes to either raw or generated facts cannot pass unnoticed.
            for filename in ("list-one.xml", "reference.json", "overlay.json", "conventions.json"):
                original = (destination / filename).read_bytes()
                (destination / filename).write_bytes(original + b" ")
                with self.assertRaisesRegex(ValueError, "checksum"):
                    currencies.load_currency_catalog(destination)
                (destination / filename).write_bytes(original)

    def test_current_historical_and_numeric_codes_are_separate(self):
        self.assertGreater(len(units.CURRENCY_REFERENCE_CODES), len(units.CURRENCY_CODES))
        self.assertNotIn("BGN", units.CURRENCY_CODES)
        self.assertIn("BGN", units.CURRENCY_REFERENCE_CODES)
        self.assertEqual(units.resolve_currency("008").code, "ALL")
        self.assertEqual(units.resolve_currency("ALL").numeric_code, "008")
        self.assertEqual(units.resolve_currency("975", as_of="2025-12-31").code, "BGN")
        self.assertEqual(units.resolve_currency("BGN", allow_historical=True).status, "historical")
        self.assertEqual(units.parse_quantity("100 BGN/mt").units, units.ureg.Unit("BGN/mt"))
        with self.assertRaises(ValueError):
            units.resolve_currency("BGN")
        for identifier in ("8", "000", "NOT_A_CURRENCY"):
            with self.subTest(identifier=identifier), self.assertRaises(ValueError):
                units.resolve_currency(identifier)

    def test_point_in_time_boundaries_and_unknown_dates_fail_closed(self):
        self.assertEqual(units.resolve_currency("BGN", as_of="1999-07-05").code, "BGN")
        for when in ("1999-07-04", "2026-01-01"):
            with self.subTest(when=when), self.assertRaises(ValueError):
                units.resolve_currency("BGN", as_of=when)
        self.assertEqual(units.resolve_currency("EUR", as_of="1999-01-01").code, "EUR")
        with self.assertRaisesRegex(ValueError, "Insufficient exact validity"):
            units.resolve_currency("USD", as_of="2020-01-01")
        self.assertEqual(units.resolve_currency("USD", as_of=currencies.DATA_VERSION).code, "USD")
        with self.assertRaisesRegex(ValueError, "later than"):
            units.resolve_currency("USD", as_of="2027-01-01")
        historic = units.resolve_currency("BGN", allow_historical=True)
        self.assertEqual(historic.withdrawal_dates, ("2026-01",))
        self.assertEqual(units.resolve_currency("EUR").withdrawal_dates, ())
        self.assertTrue(any(entry["code"] == "EUR" and "SERBIA AND MONTENEGRO" in entry["entities"]
                            for entry in currencies.CURRENCY_CATALOG.entity_history))

    def test_proprietary_account_units_use_explicit_catalogue_and_registry(self):
        # Isolate startup registration so the shared registry remains unchanged.
        script = '''
from dataclasses import replace
from decimal import Decimal
import core.units as u
u.register_currency("DESKPOINT", decimals=6)
record = replace(u.resolve_currency("EUR"), code="DESKPOINT", numeric_code=None,
                 name="Example settlement points", minor_digits=6, entities=(),
                 kind="proprietary", category="unit_of_account", source_file="desk.json")
catalog = u.CurrencyCatalog((record,), version="example-1", observed_on="2026-11-25")
amount = u.parse_currency_quantity("1.1234567", "DESKPOINT", catalog=catalog, as_of="2026-11-25")
assert u.round_currency(amount, catalog=catalog).magnitude == Decimal("1.123457")
assert "DESKPOINT" not in u.CURRENCY_REFERENCE_CODES
'''
        completed = subprocess.run([sys.executable, "-B", "-c", script], capture_output=True, text=True, timeout=30)
        self.assertEqual(completed.returncode, 0, completed.stderr)

    def test_shared_symbol_families_require_context(self):
        for symbol, code in (("$", "CAD"), ("£", "GBP"), ("¥", "CNY"), ("Rs", "PKR"),
                             ("Fr", "XOF"), ("kr", "NOK"), ("$", "MXN")):
            with self.subTest(symbol=symbol, code=code):
                with self.assertRaisesRegex(ValueError, "currency_hint"):
                    units.resolve_currency(symbol)
                self.assertEqual(units.resolve_currency(symbol, currency_hint=code).code, code)
                self.assertEqual(units.parse_currency_quantity("1.25", symbol, currency_hint=code).units, units.ureg.Unit(code))
        with self.assertRaises(ValueError):
            units.resolve_currency("$", currency_hint="JPY")
        with self.assertRaises(ValueError):
            units.resolve_currency("USD", currency_hint="CAD")
        self.assertEqual(units.resolve_currency("US$").code, "USD")

    def test_cnh_requires_fx_and_xdr_is_iso_account_unit(self):
        self.assertEqual(units.resolve_currency("CNH").kind, "quotation_convention")
        self.assertIsNone(units.resolve_currency("CNH").numeric_code)
        self.assertNotIn("CNH", units.CURRENCY_CODES)
        self.assertIn("CNH", units.CURRENCY_REFERENCE_CODES)
        with self.assertRaises(pint.DimensionalityError):
            units.Q_(100, "CNY").to("CNH")
        result = units.convert_fx(units.Q_(Decimal(100), "CNY"), "CNH", {"CNY/CNH": Decimal(".99")})
        self.assertEqual(result.magnitude, 99)
        xdr = units.resolve_currency("SDR")
        self.assertEqual((xdr.code, xdr.numeric_code, xdr.kind, xdr.category), ("XDR", "960", "iso4217", "unit_of_account"))
        self.assertEqual(units.parse_unit("SDR/mt"), units.parse_unit("XDR/mt"))
        for code in ("XXX", "XTS"):
            with self.assertRaises(ValueError):
                units.parse_currency_quantity(100, code)

    def test_vendor_namespaces_preserve_minor_denomination_case_and_scale(self):
        alias = units.CurrencyAlias("GBp", "GBP", "GBP_penny", Decimal(".01"), vendor="example-feed")
        catalog = units.CurrencyCatalog(currencies.CURRENCY_CATALOG.records, version="desk-1",
                                       observed_on=currencies.DATA_VERSION, aliases=(alias,))
        amount = units.parse_currency_quantity("123.4", "GBp", vendor="example-feed", catalog=catalog)
        self.assertEqual(amount.to("GBP").magnitude, Decimal("1.234"))
        with self.assertRaisesRegex(ValueError, "preserve scaling"):
            catalog.resolve("GBp", vendor="example-feed")
        for vendor in (None, "different-feed"):
            with self.subTest(vendor=vendor), self.assertRaises(ValueError):
                units.parse_currency_quantity(123, "GBp", vendor=vendor, catalog=catalog)
        wrong_scale = replace(alias, scale=Decimal(1))
        invalid = units.CurrencyCatalog(catalog.records, version="bad-scale", observed_on=catalog.observed_on, aliases=(wrong_scale,))
        with self.assertRaisesRegex(ValueError, "scale conflicts"):
            units.parse_currency_quantity(123, "GBp", vendor="example-feed", catalog=invalid)

    def test_dated_vendor_alias_does_not_escape_its_interval(self):
        alias = units.CurrencyAlias("EURO", "EUR", "EUR", vendor="example-feed", valid_from="2020-01-01", valid_to="2021-01-01")
        catalog = units.CurrencyCatalog(currencies.CURRENCY_CATALOG.records, version="desk-1", observed_on=currencies.DATA_VERSION, aliases=(alias,))
        self.assertEqual(units.parse_currency_quantity(5, "EURO", vendor="example-feed", catalog=catalog, as_of="2020-01-01").magnitude, 5)
        for when in (None, "2019-12-31", "2021-01-01"):
            with self.subTest(when=when), self.assertRaises(ValueError):
                units.parse_currency_quantity(5, "EURO", vendor="example-feed", catalog=catalog, as_of=when)

    def test_exact_explicit_redenominations_do_not_round_or_use_fx(self):
        original = units.Q_(Decimal("195.583"), "BGN")
        self.assertEqual(units.redenominate_currency(original, "BGN-EUR-2026", as_of="2026-01-01"), units.Q_(Decimal(100), "EUR"))
        self.assertEqual(units.redenominate_currency(units.Q_(Decimal(100000), "BGL"), "BGL-BGN-1999", as_of="1999-07-05"), units.Q_(Decimal(100), "BGN"))
        with self.assertRaises(pint.DimensionalityError):
            original.to("EUR")
        with self.assertRaises(ValueError):
            units.redenominate_currency(original, "BGN-EUR-2026", as_of="2025-12-31")
        with self.assertRaises(ValueError):
            units.redenominate_currency(units.Q_(100, "USD"), "BGN-EUR-2026", as_of="2026-01-01")
        with self.assertRaises(pint.DimensionalityError):
            units.redenominate_currency(units.Q_(100, "BGN/mt"), "BGN-EUR-2026", as_of="2026-01-01")

    def test_iso_precision_exceptions_are_independent_of_coin_subdivisions(self):
        for code, digits in (("JPY", 0), ("KWD", 3), ("CLF", 4), ("UYW", 4), ("MGA", 2), ("MRU", 2)):
            with self.subTest(code=code):
                self.assertEqual(units.CURRENCY_MINOR_DIGITS[code], digits)
                self.assertEqual(units.resolve_currency(code).minor_digits, digits)
        self.assertEqual(units.Q_(Decimal(1), "MGA_iraimbilanja").to("MGA").magnitude, Decimal(".2"))
        self.assertEqual(units.round_currency(units.Q_(Decimal("1.23456"), "CLF")).magnitude, Decimal("1.2346"))
        self.assertEqual(units.round_currency(units.Q_(Decimal("1.235"), "MGA")).magnitude, Decimal("1.24"))

    def test_cash_increment_and_dated_policy_are_explicit(self):
        amount = units.Q_(Decimal("1.03"), "CAD")
        self.assertEqual(units.round_currency(amount).magnitude, Decimal("1.03"))
        # An explicit example contract policy, not a universal national cash rule.
        rule = units.CurrencySettlementRule("CAD", Decimal(".05"), method="cash", source="example contract",
                                            valid_from="2026-01-01", valid_to="2027-01-01")
        rounded = units.round_currency(amount.to("CAD_cent"), rule=rule, as_of="2026-11-25")
        self.assertEqual(rounded.magnitude, Decimal(105))
        self.assertEqual(rounded.units, units.ureg.Unit("CAD_cent"))
        for when in (None, "2025-12-31", "2027-01-01"):
            with self.subTest(when=when), self.assertRaises(ValueError):
                units.round_currency(amount, rule=rule, as_of=when)
        with self.assertRaises(ValueError):
            units.round_currency(units.Q_(1, "USD"), rule=rule, as_of="2026-11-25")
        with self.assertRaises(ValueError):
            units.round_currency(units.Q_(1, "CAD/mt"), rule=rule, as_of="2026-11-25")

    def test_undefined_precision_requires_policy_but_does_not_block_quotes(self):
        amount = units.Q_(Decimal("1.2346"), "XDR")
        self.assertIsNone(units.decompose_quote("XDR/mt").settlement_decimals)
        with self.assertRaisesRegex(ValueError, "explicit settlement rule"):
            units.round_currency(amount)
        with self.assertRaisesRegex(ValueError, "explicit settlement rule"):
            units.convert_currency(units.Q_(1, "USD"), "XDR", rates={"XDR": 1}, quantize_settlement=True)
        policy = units.CurrencySettlementRule("XDR", Decimal(".001"), method="contract")
        self.assertEqual(units.round_currency(amount, rule=policy).magnitude, Decimal("1.235"))

    def test_snapshots_and_facts_cannot_be_mutated(self):
        with self.assertRaises(AttributeError):
            currencies.CURRENCY_CATALOG.version = "changed"
        with self.assertRaises(AttributeError):
            del currencies.CURRENCY_CATALOG._sealed
        with self.assertRaises(TypeError):
            currencies.CURRENCY_CATALOG.symbols["$"] = "CAD"
        with self.assertRaises(FrozenInstanceError):
            units.resolve_currency("USD").minor_digits = 0
        with self.assertRaises(TypeError):
            currencies.CURRENCY_OVERLAY["redenominations"][0]["old_units_per_new"] = "1"

    def test_invalid_monetary_values_and_metadata_are_rejected(self):
        for value in (True, float("nan"), float("inf"), Decimal("-Infinity")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                units.parse_currency_quantity(value, "USD")
        for increment in (0, -1, float("nan")):
            with self.subTest(increment=increment), self.assertRaises(ValueError):
                units.CurrencySettlementRule("USD", increment)
        with self.assertRaises(ValueError):
            units.CurrencySettlementRule("$", Decimal(".01"))
        with self.assertRaises(ValueError):
            replace(units.resolve_currency("USD"), numeric_code="84")


class FXFixingTests(unittest.TestCase):
    timestamp = datetime(2026, 11, 25, 12, tzinfo=timezone.utc)

    def fixing(self, base="USD", quote="EUR", rate=".8", **kwargs):
        values = dict(base=base, quote=quote, rate=Decimal(rate), source="example fixing", timestamp=self.timestamp,
                      tenor="spot", methodology="published mid")
        values.update(kwargs)
        return units.FXFixing(**values)

    def test_direct_inverse_and_cross_orientation_with_master_router(self):
        fixing = self.fixing()
        self.assertEqual(fixing.inverted().rate, Decimal("1.25"))
        self.assertEqual((fixing.inverted().base, fixing.inverted().quote), ("EUR", "USD"))
        snapshot = units.FXRateSet((fixing, self.fixing("EUR", "GBP", ".9")), as_of=self.timestamp)
        self.assertEqual(units.convert_fx(units.Q_(Decimal(100), "USD"), "GBP", snapshot).magnitude, 72)
        self.assertEqual(units.convert_fx(units.Q_(Decimal(80), "EUR"), "USD", snapshot).magnitude, 100)
        quote = units.convert(units.Q_(Decimal(100), "USD/MWh"), "GBP/MWh", rates=snapshot)
        self.assertEqual(quote.magnitude, 72)
        with self.assertRaises(AttributeError):
            snapshot.as_of = self.timestamp + timedelta(days=1)
        with self.assertRaises(AttributeError):
            del snapshot._sealed
        with self.assertRaises(TypeError):
            snapshot["USD/EUR"] = 1

    def test_fixing_time_validation_prevents_lookahead_and_stale_rates(self):
        fixing = self.fixing()
        with self.assertRaisesRegex(ValueError, "Look-ahead"):
            units.FXRateSet((fixing,), as_of=self.timestamp - timedelta(seconds=1))
        with self.assertRaisesRegex(ValueError, "Stale"):
            units.FXRateSet((fixing,), as_of=self.timestamp + timedelta(hours=2), max_age=timedelta(hours=1))
        expired = replace(fixing, valid_until=self.timestamp + timedelta(hours=1))
        with self.assertRaisesRegex(ValueError, "expired"):
            units.FXRateSet((expired,), as_of=self.timestamp + timedelta(hours=2))
        for values in (dict(timestamp=self.timestamp.replace(tzinfo=None)), dict(rate=0), dict(source=""), dict(tenor=""), dict(methodology=""), dict(base="$")):
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.fixing(**values)

    def test_incoherent_snapshots_cannot_be_triangulated(self):
        first = self.fixing()
        for differences in (dict(source="another source"), dict(timestamp=self.timestamp - timedelta(seconds=1)),
                            dict(tenor="1m"), dict(methodology="bid")):
            with self.subTest(differences=differences), self.assertRaisesRegex(ValueError, "mixed fixing"):
                units.FXRateSet((first, self.fixing("EUR", "GBP", ".9", **differences)), as_of=self.timestamp)
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            units.FXRateSet((first, first), as_of=self.timestamp)
        for extra in ((self.fixing("EUR", "USD", "2"),),
                      (self.fixing("EUR", "GBP", ".9"), self.fixing("USD", "GBP", ".75"))):
            with self.subTest(extra=extra), self.assertRaisesRegex(ValueError, "Inconsistent"):
                units.FXRateSet((first,) + extra, as_of=self.timestamp)


if __name__ == "__main__":
    unittest.main()
