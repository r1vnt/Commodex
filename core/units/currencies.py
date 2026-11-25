"""Version-pinned currency facts, settlement policies and auditable FX snapshots.

ISO metadata is distinct from symbols, cash rules, quotation conventions and
redenominations. Nothing in this module downloads data during import.
"""

from __future__ import annotations

from collections.abc import Iterable, Iterator, Mapping
from dataclasses import dataclass, replace
from datetime import date, datetime, timedelta
from decimal import Decimal, ROUND_HALF_EVEN
from hashlib import sha256
import json
from pathlib import Path
from types import MappingProxyType
from typing import Any

import pint


DATA_VERSION = "2026-11-25"
DATA_DIRECTORY = Path(__file__).with_name("data") / "currencies" / DATA_VERSION


def _date(value: date | str) -> date:
    if isinstance(value, datetime):
        raise TypeError("Use a date, not a datetime, for currency validity.")
    return value if isinstance(value, date) else date.fromisoformat(value)


def _number(value: Decimal | float | int | str) -> Decimal:
    if isinstance(value, bool):
        raise ValueError("Boolean is not a monetary value.")
    result = Decimal(str(value))
    if not result.is_finite():
        raise ValueError("Monetary values must be finite.")
    return result


def _currency_code(value: str) -> str:
    """Fixings/policies require a major code, never an implicitly chosen symbol."""
    from .conversions import _normalize_currency
    code = _normalize_currency(value)
    if not isinstance(value, str) or value.strip().upper() != code:
        raise ValueError("Use an explicit major currency code for fixing/settlement metadata.")
    return code


@dataclass(frozen=True, slots=True)
class CurrencyDefinition:
    code: str
    numeric_code: str | None
    name: str
    status: str
    minor_digits: int | None
    valid_from: date | None
    valid_to: date | None
    withdrawal_dates: tuple[str, ...]
    entities: tuple[str, ...]
    kind: str
    source_file: str
    category: str = "currency"
    is_fund: bool = False

    def __post_init__(self) -> None:
        if not isinstance(self.code, str) or not self.code.isidentifier() or not self.code.isupper():
            raise ValueError("Currency codes must be uppercase identifiers.")
        if self.status not in {"current", "historical"}:
            raise ValueError("Currency status must be current or historical.")
        if self.numeric_code is not None and (not isinstance(self.numeric_code, str) or len(self.numeric_code) != 3 or not self.numeric_code.isascii() or not self.numeric_code.isdigit()):
            raise ValueError("ISO numeric codes must retain three ASCII digits.")
        if self.minor_digits is not None and (isinstance(self.minor_digits, bool) or not isinstance(self.minor_digits, int) or self.minor_digits < 0):
            raise ValueError("Invalid minor-unit exponent.")
        for field in ("valid_from", "valid_to"):
            if getattr(self, field) is not None:
                object.__setattr__(self, field, _date(getattr(self, field)))
        if self.valid_from and self.valid_to and self.valid_from >= self.valid_to:
            raise ValueError("Currency validity must be a non-empty half-open interval.")
        for field in ("withdrawal_dates", "entities"):
            object.__setattr__(self, field, tuple(getattr(self, field)))


@dataclass(frozen=True, slots=True)
class CurrencyAlias:
    alias: str
    code: str
    unit: str
    scale: Decimal = Decimal(1)
    vendor: str | None = None
    valid_from: date | None = None
    valid_to: date | None = None

    def __post_init__(self) -> None:
        if any(not isinstance(value, str) or not value for value in (self.alias, self.code, self.unit)):
            raise ValueError("Alias, currency code and unit are required.")
        scale = _number(self.scale)
        if scale <= 0:
            raise ValueError("Alias scale must be positive.")
        object.__setattr__(self, "scale", scale)
        for field in ("valid_from", "valid_to"):
            if getattr(self, field) is not None:
                object.__setattr__(self, field, _date(getattr(self, field)))
        if self.valid_from and self.valid_to and self.valid_from >= self.valid_to:
            raise ValueError("Alias interval must be non-empty.")


class CurrencyCatalog:
    """Immutable reference snapshot; unknown dates fail closed in dated queries."""
    __slots__ = ("records", "version", "observed_on", "aliases", "entity_history", "symbols", "shared_symbols", "_sealed")

    def __init__(self, records: Iterable[CurrencyDefinition], *, version: str, observed_on: date | str,
                 aliases: Iterable[CurrencyAlias] = (), symbols: Mapping[str, str] | None = None,
                 shared_symbols: Mapping[str, Iterable[str]] | None = None,
                 entity_history: Iterable[Mapping[str, Any]] = ()) -> None:
        self.records = tuple(records)
        self.version = version
        self.observed_on = _date(observed_on)
        self.aliases = tuple(aliases)
        self.entity_history = _freeze(tuple(entity_history))
        self.symbols = MappingProxyType(dict(symbols or {}))
        self.shared_symbols = MappingProxyType({key: frozenset(values) for key, values in (shared_symbols or {}).items()})
        codes = {record.code for record in self.records}
        if any(alias.code not in codes for alias in self.aliases):
            raise ValueError("Currency alias refers to an undefined denomination.")
        self._sealed = True

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, "_sealed", False):
            raise AttributeError("Currency catalogs are immutable snapshots.")
        object.__setattr__(self, name, value)

    def __delattr__(self, name: str) -> None:
        raise AttributeError("Currency catalogs are immutable snapshots.")

    def resolve_alias(self, token: str, *, vendor: str | None = None, as_of: date | str | None = None) -> CurrencyAlias | None:
        matches = [alias for alias in self.aliases if alias.alias == token and alias.vendor == vendor]
        if as_of is not None:
            when = _date(as_of)
            matches = [a for a in matches if (a.valid_from is None or when >= a.valid_from) and (a.valid_to is None or when < a.valid_to)]
        elif any(a.valid_from or a.valid_to for a in matches):
            raise ValueError("A dated vendor alias requires as_of.")
        if len(matches) > 1:
            raise ValueError("Ambiguous vendor alias.")
        return matches[0] if matches else None

    def resolve(self, token: str, *, as_of: date | str | None = None, vendor: str | None = None,
                currency_hint: str | None = None, allow_historical: bool = False) -> CurrencyDefinition:
        if not isinstance(token, str) or not token.strip():
            raise ValueError("Currency identifier must be a non-empty string.")
        token = token.strip()
        alias = self.resolve_alias(token, vendor=vendor, as_of=as_of)
        if alias:
            if alias.scale != 1:
                raise ValueError("Minor-unit vendor aliases require parse_currency_quantity to preserve scaling.")
            token = alias.code
        elif vendor is not None:
            # Codes remain valid in vendor feeds, but unknown non-code aliases
            # are never guessed using a different vendor's namespace.
            if not any(r.code == token for r in self.records):
                raise ValueError("Unknown vendor-specific currency alias.")
        choices = self.shared_symbols.get(token)
        if choices:
            if currency_hint is None or currency_hint.upper() not in choices:
                raise ValueError(f"Ambiguous currency symbol {token!r}; supply a valid currency_hint.")
            token = currency_hint.upper()
        else:
            if token not in self.symbols and not alias and token.isalpha() and not (token.isupper() or token.islower()):
                raise ValueError("Mixed-case currency tokens require an explicit vendor alias.")
            token = self.symbols.get(token, token).upper()
            if currency_hint is not None and currency_hint.upper() != token:
                raise ValueError("Currency hint conflicts with the explicit identifier.")
        candidates = [r for r in self.records if r.code == token or r.numeric_code == token]
        if as_of is None:
            current = [r for r in candidates if r.status == "current"]
            candidates = current or (candidates if allow_historical else [])
        else:
            when = _date(as_of)
            if when > self.observed_on:
                raise ValueError("Currency date is later than this reference snapshot.")
            known, uncertain = [], []
            for record in candidates:
                if record.valid_from and when < record.valid_from or record.valid_to and when >= record.valid_to:
                    continue
                if record.valid_from is None and not (record.status == "current" and when == self.observed_on):
                    uncertain.append(record)
                elif record.status == "historical" and record.valid_to is None:
                    uncertain.append(record)
                else:
                    known.append(record)
            if uncertain:
                raise ValueError("Insufficient exact validity data for this point-in-time currency query.")
            candidates = known
        if len(candidates) != 1:
            raise ValueError(f"Unknown, inactive or ambiguous currency identifier: {token!r}")
        return candidates[0]


def _load_reference(directory: str | Path = DATA_DIRECTORY) -> tuple[CurrencyCatalog, Mapping[str, Any], Mapping[str, Any]]:
    directory = Path(directory)
    manifest = json.loads((directory / "manifest.json").read_text(encoding="utf-8"))
    for entry in manifest["sources"]:
        if sha256((directory / entry["file"]).read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError("Currency reference source checksum mismatch.")
    for filename, expected in manifest.get("derived_sha256", {}).items():
        if sha256((directory / filename).read_bytes()).hexdigest() != expected:
            raise ValueError("Currency reference derived-data checksum mismatch.")
    data = json.loads((directory / "reference.json").read_text(encoding="utf-8"))
    conventions = json.loads((directory / "conventions.json").read_text(encoding="utf-8"))
    overlay = json.loads((directory / "overlay.json").read_text(encoding="utf-8"))
    if data["schema_version"] != 1 or len({manifest["version"], data["version"], overlay["version"]}) != 1:
        raise ValueError("Unsupported/mismatched currency snapshot versions.")
    records = []
    for item in data["records"]:
        values = dict(item)
        for key in ("valid_from", "valid_to"):
            values[key] = _date(values[key]) if values[key] else None
        for key in ("entities", "withdrawal_dates"):
            values[key] = tuple(values[key])
        records.append(CurrencyDefinition(**values))
    aliases = [CurrencyAlias(**{**item, "scale": Decimal(item["scale"])}) for item in overlay["aliases"]]
    return CurrencyCatalog(records, version=data["version"], observed_on=data["observed_on"], aliases=aliases,
                           symbols=conventions["symbols"], shared_symbols=conventions["shared_symbols"],
                           entity_history=data["entity_history"]), _freeze(conventions), _freeze(overlay)


def _freeze(value: Any) -> Any:
    if isinstance(value, dict):
        return MappingProxyType({key: _freeze(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    return value


CURRENCY_CATALOG, CURRENCY_CONVENTIONS, CURRENCY_OVERLAY = _load_reference()


def load_currency_catalog(directory: str | Path) -> CurrencyCatalog:
    """Load a checked snapshot without changing the process's unit registry."""
    return _load_reference(directory)[0]


def resolve_currency(identifier: str, *, as_of: date | str | None = None, vendor: str | None = None,
                     currency_hint: str | None = None, allow_historical: bool = False) -> CurrencyDefinition:
    """Resolve currency metadata without guessing symbol context or unknown dates."""
    return CURRENCY_CATALOG.resolve(identifier, as_of=as_of, vendor=vendor,
                                    currency_hint=currency_hint, allow_historical=allow_historical)


def parse_currency_quantity(amount: Decimal | float | int | str, identifier: str, *, catalog: CurrencyCatalog = CURRENCY_CATALOG,
                            vendor: str | None = None, as_of: date | str | None = None,
                            currency_hint: str | None = None, allow_historical: bool = False) -> pint.Quantity:
    """Strict currency-only ingestion; vendor minor denominations retain scale."""
    from .registry import ureg
    alias = catalog.resolve_alias(identifier, vendor=vendor, as_of=as_of)
    record = catalog.resolve(alias.code if alias else identifier, as_of=as_of, vendor=None if alias else vendor,
                             currency_hint=currency_hint, allow_historical=allow_historical)
    if record.category in {"testing", "no_currency"}:
        raise ValueError("ISO testing/no-currency identifiers are not monetary amounts.")
    if alias:
        # Do not trust a vendor unit definition without checking its declared scale.
        measured = ureg.Quantity(Decimal(1), alias.unit).to(record.code).magnitude
        if measured != alias.scale:
            raise ValueError("Vendor alias scale conflicts with the unit registry.")
    return ureg.Quantity(_number(amount), alias.unit if alias else record.code)


@dataclass(frozen=True, slots=True)
class CurrencySettlementRule:
    """A contract/cash increment, independent of ISO decimal places."""
    currency: str
    increment: Decimal
    method: str = "electronic"
    rounding: str = ROUND_HALF_EVEN
    valid_from: date | None = None
    valid_to: date | None = None
    source: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "currency", _currency_code(self.currency))
        for field in ("valid_from", "valid_to"):
            if getattr(self, field) is not None:
                object.__setattr__(self, field, _date(getattr(self, field)))
        number = _number(self.increment)
        if number <= 0 or self.method not in {"cash", "electronic", "contract"}:
            raise ValueError("A positive increment and explicit settlement method are required.")
        object.__setattr__(self, "increment", number)
        # Let Decimal validate the requested rounding policy without rounding an amount.
        Decimal(1).quantize(Decimal(1), rounding=self.rounding)
        if self.valid_from and self.valid_to and self.valid_from >= self.valid_to:
            raise ValueError("Invalid settlement rule interval.")


def round_currency(quantity: pint.Quantity, *, rule: CurrencySettlementRule | None = None,
                   as_of: date | str | None = None, catalog: CurrencyCatalog = CURRENCY_CATALOG) -> pint.Quantity:
    """Round a plain currency amount, never a commodity quote or FX rate."""
    from .conversions import _numerator_currency_iso
    from .registry import ureg
    if not hasattr(quantity, "units"):
        raise TypeError("Expected a currency quantity.")
    code = _numerator_currency_iso(quantity)
    if quantity.dimensionality != {f"[currency_{code}]": 1}:
        raise ValueError("Settlement rounding requires a plain currency amount, not a price.")
    major = quantity.to(code)
    if rule is None:
        record = catalog.resolve(code, as_of=as_of, allow_historical=as_of is None)
        if record.minor_digits is None:
            raise ValueError("Unknown/undefined ISO precision requires an explicit settlement rule.")
        increment, rounding = Decimal(1).scaleb(-record.minor_digits), ROUND_HALF_EVEN
    else:
        if rule.currency != code:
            raise ValueError("Settlement rule currency does not match the amount.")
        if rule.valid_from or rule.valid_to:
            if as_of is None:
                raise ValueError("A dated settlement rule requires as_of.")
            when = _date(as_of)
            if rule.valid_from and when < rule.valid_from or rule.valid_to and when >= rule.valid_to:
                raise ValueError("Settlement rule is not effective on this date.")
        increment, rounding = rule.increment, rule.rounding
    value = (_number(major.magnitude) / increment).quantize(Decimal(1), rounding=rounding) * increment
    return ureg.Quantity(value, code).to(quantity.units)


def redenominate_currency(quantity: pint.Quantity, event_id: str, *, as_of: date | str) -> pint.Quantity:
    """Explicit dated legal denomination conversion, not an implicit FX quote."""
    from .conversions import _numerator_currency_iso
    from .registry import ureg
    events = [e for e in CURRENCY_OVERLAY["redenominations"] if e["id"] == event_id]
    if len(events) != 1:
        raise ValueError("Unknown redenomination event.")
    event = events[0]
    if _date(as_of) < _date(event["effective"]):
        raise ValueError("Redenomination is not yet effective.")
    if _numerator_currency_iso(quantity) != event["old"]:
        raise ValueError("Redenomination source currency does not match the event.")
    # This API intentionally accepts plain amounts; contract rebasing needs policy.
    major = quantity.to(event["old"])
    amount = _number(major.magnitude) / Decimal(event["old_units_per_new"])
    return ureg.Quantity(amount, event["new"])


@dataclass(frozen=True, slots=True)
class FXFixing:
    """One base unit buys rate quote units; provenance is mandatory."""
    base: str
    quote: str
    rate: Decimal
    source: str
    timestamp: datetime
    tenor: str
    methodology: str
    valid_until: datetime | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "base", _currency_code(self.base))
        object.__setattr__(self, "quote", _currency_code(self.quote))
        value = _number(self.rate)
        if value <= 0 or self.base == self.quote:
            raise ValueError("FX fixing requires distinct currencies and a positive rate.")
        object.__setattr__(self, "rate", value)
        if not isinstance(self.timestamp, datetime) or self.timestamp.tzinfo is None or self.timestamp.utcoffset() is None:
            raise ValueError("FX timestamp must be timezone-aware.")
        if any(not isinstance(value, str) or not value.strip() for value in (self.source, self.tenor, self.methodology)):
            raise ValueError("FX source, tenor and methodology are required.")
        if self.valid_until is not None and (not isinstance(self.valid_until, datetime) or self.valid_until.tzinfo is None or self.valid_until.utcoffset() is None or self.valid_until < self.timestamp):
            raise ValueError("Invalid fixing expiry.")

    def inverted(self) -> FXFixing:
        return replace(self, base=self.quote, quote=self.base, rate=Decimal(1) / self.rate)


class FXRateSet(Mapping[str, Decimal]):
    """Coherent, point-in-time sparse FX map accepted by existing convert APIs."""
    __slots__ = ("fixings", "as_of", "_rates", "_sealed")

    def __init__(self, fixings: Iterable[FXFixing], *, as_of: datetime, max_age: timedelta | None = None) -> None:
        if not isinstance(as_of, datetime) or as_of.tzinfo is None or as_of.utcoffset() is None:
            raise ValueError("FX as_of must be timezone-aware.")
        if max_age is not None and max_age < timedelta(0):
            raise ValueError("max_age must be non-negative.")
        self.fixings = tuple(fixings)
        self.as_of = as_of
        groups = {(f.source, f.timestamp, f.tenor, f.methodology) for f in self.fixings}
        if len(groups) > 1:
            raise ValueError("Cannot triangulate mixed fixing sources, times, tenors or methodologies.")
        rates = {}
        for fixing in self.fixings:
            if fixing.timestamp > as_of or fixing.valid_until is not None and as_of > fixing.valid_until:
                raise ValueError("Look-ahead or expired FX fixing.")
            if max_age is not None and as_of - fixing.timestamp > max_age:
                raise ValueError("Stale FX fixing.")
            key = fixing.base + "/" + fixing.quote
            if key in rates:
                raise ValueError("Duplicate directed FX fixing.")
            rates[key] = fixing.rate
        self._rates = MappingProxyType(rates)
        # A coherent fixing cannot contain contradictory inverse/cross rates.
        graph = {}
        for fixing in self.fixings:
            graph.setdefault(fixing.base, []).append((fixing.quote, fixing.rate))
            graph.setdefault(fixing.quote, []).append((fixing.base, Decimal(1) / fixing.rate))
        values = {}
        for start in graph:
            if start in values:
                continue
            values[start] = Decimal(1)
            pending = [start]
            while pending:
                node = pending.pop()
                for destination, rate in graph[node]:
                    expected = values[node] * rate
                    if destination in values:
                        if abs(values[destination] - expected) > max(abs(expected), abs(values[destination])) * Decimal("1e-8"):
                            raise ValueError("Inconsistent inverse/cross rates in fixing snapshot.")
                    else:
                        values[destination] = expected
                        pending.append(destination)
        self._sealed = True

    def __setattr__(self, name: str, value: Any) -> None:
        if getattr(self, "_sealed", False):
            raise AttributeError("FX fixing snapshots are immutable.")
        object.__setattr__(self, name, value)

    def __delattr__(self, name: str) -> None:
        raise AttributeError("FX fixing snapshots are immutable.")

    def __getitem__(self, key: str) -> Decimal:
        return self._rates[key]

    def __iter__(self) -> Iterator[str]:
        return iter(self._rates)

    def __len__(self) -> int:
        return len(self._rates)


__all__ = ("CurrencyDefinition", "CurrencyAlias", "CurrencyCatalog", "CURRENCY_CATALOG",
           "DATA_VERSION", "load_currency_catalog", "resolve_currency", "parse_currency_quantity", "CurrencySettlementRule",
           "round_currency", "redenominate_currency", "FXFixing", "FXRateSet")
