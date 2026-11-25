# Currency reference data and settlement

`currencies.py` loads an offline, versioned snapshot from
`data/currencies/2026-11-25/`. `registry.py` builds independent Pint dimensions
from those records. No import downloads or updates reference data.

## Provenance and scope

- `list-one.xml` and `list-three.xml` are the original current and historical
  ISO 4217 lists from [SIX, the ISO maintenance agency](https://www.six-group.com/en/products-services/financial-information/market-reference-data/data-standards.html).
  Their publication dates and SHA-256 hashes are in `manifest.json`.
- `reference.json` is generated from those files and a reviewed `overlay.json`.
  Numeric codes remain three-character strings, including leading zeros.
  Entity names and the original precision of withdrawal dates are preserved.
- `overlay.json` contains additional facts with source links, legal denomination
  events, dated overrides, and non-ISO quotation conventions.
- `conventions.json` preserves Commodex's existing broker symbols and named
  subunits. These are application conventions, not ISO-provided symbol data.
  Manifest hashes cover the generated data and both overlays as well as the XML.

`CURRENCY_CODES` means current ISO codes. `CURRENCY_REFERENCE_CODES` additionally
includes historical denominations and CNH. The current ISO list includes funds,
precious-metal codes, XDR, and reserved identifiers; category metadata distinguishes
these. Strict monetary ingestion rejects `XXX` and `XTS`.

SIX's historical list also records entity withdrawals of currencies still in
use elsewhere. For example, its EUR entry for Serbia and Montenegro does not
mean EUR itself was withdrawn. Such entries are retained in
`CURRENCY_CATALOG.entity_history`, separate from denomination records.

## Historical resolution

```python
from core.units import resolve_currency, parse_currency_quantity

lev = resolve_currency("975", as_of="2025-12-31")  # BGN, numeric code 975
amount = parse_currency_quantity("100.25", "BGN", as_of="2025-12-31")
historic = resolve_currency("BGN", allow_historical=True)  # undated inspection
```

Validity intervals are `[valid_from, valid_to)`. The reviewed BGN interval ends
on 2026-01-01, when Bulgaria adopted EUR. January's parallel cash circulation is
a separate rule, not an extension of its current ISO status. The overlay also
records EUR's introduction on 1999-01-01.

**This snapshot is not a complete historical validity database.** SIX's XML
does not provide all introduction dates or exact withdrawal days. A withdrawal
month is retained as a month, not expanded into an invented day. Unknown dates
and precision remain `None`. Dated resolution rejects insufficient evidence,
ambiguous code reuse, and dates later than the snapshot. For example, this
snapshot cannot resolve USD for 2020 from an exact validity interval. Current
records can be resolved at the snapshot's observation date. Add sourced intervals
to a new versioned overlay before enabling other historical queries.

Legacy `parse_unit`, `parse_quantity`, and direct `Q_` remain undated dimensional
APIs. They can represent historical currency amounts and keep existing broker
shorthand defaults. Use `resolve_currency` / `parse_currency_quantity` at dated
ingestion boundaries; constructing a Pint quantity does not prove temporal validity.

## Symbols and feed aliases

```python
parse_currency_quantity("25", "$", currency_hint="CAD")
parse_currency_quantity("25", "¥", currency_hint="CNY")
```

Strict ingestion requires a hint for shared dollar/peso, pound, yen/yuan, rupee,
franc and krone symbols. Explicit codes or disambiguated forms such as `US$`
are preferable. Supported candidates are in `conventions.json`; unsupported
symbols or hints fail rather than selecting a currency by geographic guess.

Aliases have an exact, case-sensitive token, vendor namespace, unit, scale,
and optional validity interval. For example, a desk can supply a `CurrencyAlias`
for `GBp`, code `GBP`, unit `GBP_penny`, scale `Decimal(".01")`, vendor
`"example-feed"` in its own `CurrencyCatalog`. Then use
`parse_currency_quantity(123, "GBp", vendor="example-feed", catalog=catalog)`.
The unit's scale must match the alias declaration. A scaled alias cannot pass
through the metadata-only resolver and lose its denomination.

The bundled aliases include SDR and explicitly named desk RMB conventions.
They are **not** a complete Bloomberg, Reuters, exchange or price-reporting
agency alias catalogue. Add licensed, validated feed mappings to reviewed
overlays, with their actual dates and scaling; do not uppercase feed tokens.

CNH has its own quotation dimension and no invented ISO numeric code. CNY/CNH
conversion requires an explicit FX rate. XDR (SDR) already is an ISO-coded unit
of account; it is not classified as a proprietary currency. For proprietary
settlement denominations, register a dedicated currency dimension during startup
and supply a matching `CurrencyDefinition` in a desk `CurrencyCatalog` with
`kind="proprietary"` and `category="unit_of_account"`. A custom catalogue does
not silently mutate the process's unit registry.

## Redenomination is explicit

```python
from decimal import Decimal
from core.units import Q_, redenominate_currency

redenominate_currency(Q_(Decimal("195.583"), "BGN"), "BGN-EUR-2026",
                     as_of="2026-01-01")  # 100 EUR, before settlement rounding
```

The event stores an exact old-units-per-new ratio, effective date and source.
The initial reviewed events are BGL/BGN (1000:1, 1999-07-05) and BGN/EUR
(1.95583:1, 2026-01-01). The latter uses the
[ECB's published irrevocable conversion rate](https://www.ecb.europa.eu/press/pr/date/2026/html/ecb.pr260101~c830245e42.en.html).
Pint never treats these currencies as dimensionally interchangeable. Legal
denomination conversion, market FX, commodity contract rebasing and final
settlement rounding are separate operations. This helper accepts plain amounts.

## ISO precision and settlement increments

`round_currency` defaults to the record's ISO decimal exponent with Decimal
half-even rounding. Examples include JPY 0, KWD 3, and CLF/UYW 4 digits. MGA
and MRU have ISO precision 2; their physical fifth subunits remain separate
unit definitions. Unknown/undefined precision, including XDR, requires an
explicit policy. Quote decomposition reports `settlement_decimals=None` in
that case; FX calculations themselves do not require rounding.

```python
from core.units import CurrencySettlementRule, round_currency

# Example contract policy; choose the actual applicable policy for the payment.
policy = CurrencySettlementRule("CAD", Decimal(".05"), method="cash",
                                source="example contract",
                                valid_from="2026-01-01", valid_to="2027-01-01")
round_currency(Q_(Decimal("1.03"), "CAD"))  # 1.03 CAD, ISO precision
round_currency(Q_(Decimal("1.03"), "CAD"), rule=policy,
               as_of="2026-11-25")  # 1.05 CAD, explicit cash increment
```

Rules retain method, rounding policy, effective interval and optional source.
An interval requires `as_of`. No national cash-rounding table is assumed from
ISO digits, and no jurisdiction, payment method or contract is inferred.
Round final plain amounts, not currency/commodity quotes. Returned quantities
retain their input currency unit, including a named minor denomination.

## FX fixing metadata

`FXFixing(base, quote, rate, source, timestamp, tenor, methodology)` always means
**1 base unit = rate quote units**. It requires explicit major currency codes,
a positive finite Decimal-compatible rate and a timezone-aware timestamp.
Optional `valid_until` bounds the fixing's use. `inverted()` swaps the direction.

`FXRateSet(fixings, as_of=aware_datetime, max_age=timedelta(...))` is an immutable
sparse pair mapping accepted by `convert`, `convert_fx` and `triangulate_fx_rate`.
It preserves the original fixings and rejects look-ahead, expiry, staleness,
duplicate directions and mixtures of source, timestamp, tenor or methodology.
Inverse/cross-rate cycles must agree within relative tolerance `1e-8`. Construct
a separate snapshot for each coherent fixing convention. Legacy unannotated
rate dictionaries remain supported but do not acquire provenance automatically.

## Updating the snapshot

1. Create a new dated directory; keep the previous snapshot for reproducibility.
2. Download the original SIX XML over verified HTTPS and retain the bytes.
3. Review dated overlays, alias scales, precision changes, country history and
   each denomination event against primary sources. Cash rules need their own
   source and jurisdiction/contract policy.
4. Run `python -B tools/build_currency_reference.py <new-directory>` to generate
   `reference.json` and the manifest. This compiler never downloads data.
5. Load with `load_currency_catalog(<new-directory>)` and review/test the diff.
   To change the process default, update `DATA_VERSION` deliberately and package
   the snapshot directory alongside `currencies.py`.

ISO 10962 CFI codes classify financial instruments. They belong in instrument
reference data, not the ISO 4217 currency dimension or currency numeric code.
