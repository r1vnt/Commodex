"""Commodex quote_parser implementation and public interfaces."""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from datetime import date
from decimal import Decimal
import ast
import pint
import re




# Exact tokens, not substrings: MM must not catch MMBtu or a custom unit name.
_quote_parsing_AMBIGUOUS_TOKENS = frozenset({
    "MT", "M", "MM", "m", "MB", "gal", "Gal", "gallon", "gallons",
    "oz", "bu", "bushel", "c", "cent", "cents", "p", "pence", "penny", "GT", "gt", "lot", "lots",
    "ton", "tons", "cwt", "barrel", "barrels", "ppt", "mu", "ton_mile",
})



_quote_parsing_SELECTORS = ("commodity_id", "currency", "exchange", "vendor", "region")



_quote_parsing_IDENTIFIER = re.compile(r"(?<!\w)[^\W\d]\w*(?!\w)")





def _quote_parsing_selector(value: str | None, name: str) -> str | None:
    if value is not None and (not isinstance(value, str) or not value or value != value.strip()):
        raise ValueError(f"{name} must be a non-empty, trimmed identifier.")
    return value





def _quote_parsing_major_currency(value: str | None) -> str | None:
    if value is None:
        return None
    from core.units.currencies import _currency_code
    return _currency_code(value)





@dataclass(frozen=True, slots=True)
class QuoteTokenRule:
    """An exact token convention, scoped by source/product and optional dates.

    ``unit`` is a registered Pint unit expression, including named scaled units.
    Dates are half-open. Overlapping rules must agree; there is no guessed priority.
    """
    token: str
    unit: str
    commodity_id: str | None = None
    currency: str | None = None
    exchange: str | None = None
    vendor: str | None = None
    region: str | None = None
    valid_from: date | None = None
    valid_to: date | None = None

    def __post_init__(self) -> None:
        from core.units.currencies import _date
        if not isinstance(self.token, str) or not (self.token.isidentifier() or self.token in {"$", "¥", "£", "€", "₨"}):
            raise ValueError("A quote rule must match one complete token.")
        _quote_parsing_selector(self.unit, "unit")
        if self.unit is None:
            raise ValueError("A quote rule requires an explicit unit.")
        for name in _quote_parsing_SELECTORS:
            _quote_parsing_selector(getattr(self, name), name)
        object.__setattr__(self, "currency", _quote_parsing_major_currency(self.currency))
        for name in ("valid_from", "valid_to"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, _date(getattr(self, name)))
        if self.valid_from and self.valid_to and self.valid_from >= self.valid_to:
            raise ValueError("Quote rule validity must be a non-empty half-open interval.")
        if self.token in {"lot", "lots"} and self.commodity_id is None:
            raise ValueError("A lot convention requires a contract-qualified commodity_id.")


from core.units.currencies import CURRENCY_CATALOG


@dataclass(frozen=True, slots=True)
class QuoteParsingContext:
    """Explicit quote origin and immutable conventions for one ingestion stream.

    Source labels have no built-in unit defaults. Supply reviewed token_rules
    for ambiguous physical tokens and a currency_catalog for dated metadata.
    """
    commodity_id: str | None = None
    currency: str | None = None
    exchange: str | None = None
    vendor: str | None = None
    region: str | None = None
    reference_date: date | None = None
    strict: bool = True
    token_rules: tuple[QuoteTokenRule, ...] = ()
    currency_catalog: CurrencyCatalog = field(default=CURRENCY_CATALOG, repr=False)

    def __post_init__(self) -> None:
        from core.units.currencies import CurrencyCatalog
        from core.units.currencies import _date
        for name in _quote_parsing_SELECTORS:
            _quote_parsing_selector(getattr(self, name), name)
        object.__setattr__(self, "currency", _quote_parsing_major_currency(self.currency))
        if self.reference_date is not None:
            object.__setattr__(self, "reference_date", _date(self.reference_date))
        if not isinstance(self.strict, bool):
            raise TypeError("strict must be a boolean.")
        object.__setattr__(self, "token_rules", tuple(self.token_rules))
        if not all(isinstance(rule, QuoteTokenRule) for rule in self.token_rules):
            raise TypeError("token_rules must contain QuoteTokenRule objects.")
        if not isinstance(self.currency_catalog, CurrencyCatalog):
            raise TypeError("currency_catalog must be a CurrencyCatalog.")





def _quote_parsing_source_unit(token: str, context: QuoteParsingContext) -> pint.Unit | None:
    from core.units.exceptions import AmbiguousUnitError
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.registry import ureg
    matches = []
    for rule in context.token_rules:
        if rule.token != token or any(getattr(rule, name) is not None and getattr(rule, name) != getattr(context, name) for name in _quote_parsing_SELECTORS):
            continue
        if rule.valid_from or rule.valid_to:
            if context.reference_date is None:
                raise MissingReferenceConditionsError(f"Dated quote convention for {token!r} requires reference_date.")
            if rule.valid_from and context.reference_date < rule.valid_from or rule.valid_to and context.reference_date >= rule.valid_to:
                continue
        matches.append(ureg.Unit(rule.unit))
    if matches and any(unit != matches[0] for unit in matches[1:]):
        raise AmbiguousUnitError(f"Conflicting source conventions for token {token!r}.")
    return matches[0] if matches else None





def _quote_parsing_unit_grammar(expression: str) -> None:
    """Accept only algebraic unit syntax, without calls, attributes or literals."""
    try:
        tree = ast.parse(expression, mode="eval")
    except (SyntaxError, ValueError) as error:
        raise ValueError("Invalid quote unit syntax.") from error

    def visit(node):
        if isinstance(node, ast.Name):
            return
        if isinstance(node, ast.Constant) and type(node.value) is int and node.value == 1:
            return  # Pint renders reciprocal units with an explicit unity numerator.
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Mult, ast.Div)):
            visit(node.left)
            visit(node.right)
            return
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            visit(node.left)
            power = node.right
            if isinstance(power, ast.UnaryOp) and isinstance(power.op, (ast.UAdd, ast.USub)):
                power = power.operand
            if isinstance(power, ast.Constant) and type(power.value) is int and abs(power.value) <= 64:
                return
        raise ValueError("Quote units allow identifiers, multiplication, division and integer powers only.")

    visit(tree.body)





def _quote_parsing_numerator(token: str, context: QuoteParsingContext) -> str:
    from core.units.exceptions import AmbiguousUnitError
    from core.units.registry import ureg
    from core.units.registry import sanitize_unit_string
    unit = _quote_parsing_source_unit(token, context)
    if unit is None:
        alias = context.currency_catalog.resolve_alias(token, vendor=context.vendor, as_of=context.reference_date)
        if alias:
            from core.units.currencies import parse_currency_quantity
            unit = parse_currency_quantity(1, token, vendor=context.vendor, as_of=context.reference_date,
                                           catalog=context.currency_catalog).units
        elif token in {"c", "cent", "cents", "p", "pence", "penny"} and context.currency is not None:
            candidate = f"{context.currency}_cent" if token in {"c", "cent", "cents"} else f"{context.currency}_penny"
            if candidate and candidate in ureg:
                unit = ureg.Unit(candidate)
        elif token in context.currency_catalog.symbols or token in context.currency_catalog.shared_symbols or token.isdigit():
            if context.strict or context.currency is not None or token not in context.currency_catalog.shared_symbols:
                record = context.currency_catalog.resolve(token, currency_hint=context.currency, as_of=context.reference_date)
                unit = ureg.Unit(record.code)
        elif token not in _quote_parsing_AMBIGUOUS_TOKENS:
            # Direct Pint lookup preserves named minor units; no legacy ticker map.
            if context.strict and token.isalpha() and not (token.isupper() or token.islower()):
                raise AmbiguousUnitError(f"Mixed-case currency token {token!r} requires an explicit source alias or rule.")
            _quote_parsing_unit_grammar(token)
            if token.upper() in {record.code for record in context.currency_catalog.records}:
                token = token.upper()
            if not context.strict and token not in ureg:
                token = sanitize_unit_string(token, context.currency)
            unit = ureg.Unit(token)
    if unit is None and not context.strict:
        unit = ureg.Unit(sanitize_unit_string(f"{token}/kilogram", context.currency).split("/", 1)[0])
    if unit is None:
        raise AmbiguousUnitError(f"Ambiguous quote token {token!r}; supply currency or an applicable source rule.")
    dimensions = unit.dimensionality
    if len(dimensions) != 1:
        raise ValueError("Quote numerator must be an isolated currency denomination.")
    dimension, power = next(iter(dimensions.items()))
    if not dimension.startswith("[currency_") or power != 1:
        raise ValueError("Quote numerator must be an isolated currency denomination.")
    code = dimension[len("[currency_"):-1]
    if context.currency is not None and code != context.currency:
        raise ValueError("Quote currency conflicts with the parsing context.")
    context.currency_catalog.resolve(code, as_of=context.reference_date)
    return str(unit)





def _quote_parsing_normalize_quote_expression(raw: str, context: QuoteParsingContext, default_currency: str | None = None) -> str:
    from core.units.exceptions import AmbiguousUnitError
    from core.units.registry import ureg
    from core.units.registry import sanitize_unit_string
    if not isinstance(context, QuoteParsingContext):
        raise TypeError("context must be a QuoteParsingContext.")
    if default_currency is not None:
        code = _quote_parsing_major_currency(default_currency)
        if context.currency is not None and code != context.currency:
            raise ValueError("default_currency conflicts with the parsing context.")
        context = replace(context, currency=code)
    if not isinstance(raw, str) or not raw.strip():
        raise ValueError("Quote unit expression must be a non-empty string.")
    text = re.sub(r"\s+per\s+", " / ", raw.strip(), flags=re.I).replace("^", "**")
    if "/" not in text:
        raise ValueError("A commodity quote requires currency/commodity units.")
    numerator, denomination = (part.strip() for part in text.split("/", 1))
    numerator = _quote_parsing_numerator(numerator, context)

    def resolve(match):
        from core.units.exceptions import AmbiguousUnitError
        token = match[0]
        unit = _quote_parsing_source_unit(token, context)
        if unit is not None:
            return f"({unit})"
        if token in _quote_parsing_AMBIGUOUS_TOKENS:
            if context.strict:
                raise AmbiguousUnitError(f"Ambiguous quote token {token!r}; supply an applicable source rule or an explicit unit name.")
            return token
        # A directly registered name bypasses prefix-engine collisions, including
        # custom units. This lookup never lowercases a vendor token.
        return ureg._units[token].name if token in ureg._units else token

    denomination = _quote_parsing_IDENTIFIER.sub(resolve, denomination)
    expression = f"{numerator} / {denomination}"
    if not context.strict:
        expression = sanitize_unit_string(expression, context.currency)
    _quote_parsing_unit_grammar(expression)
    parsed = ureg.Unit(expression)
    commodity = ureg.Unit(numerator) / parsed
    if not commodity.dimensionality or any(name.startswith("[currency_") for name in commodity.dimensionality):
        raise ValueError("Quote denomination must represent a commodity, not a currency or a bare multiplier.")
    return expression





def parse_quote(value: int | float | Decimal | str, unit: str | None = None, *,
                context: QuoteParsingContext | None = None, as_decimal: bool = True) -> pint.Quantity:
    """Parse a commodity quote, rejecting unresolved ambiguity by default."""
    from core.units.registry import parse_quantity
    return parse_quantity(value, unit, context=context if context is not None else QuoteParsingContext(), as_decimal=as_decimal)





def parse_quote_unit(expression: str, *, context: QuoteParsingContext | None = None) -> pint.Unit:
    """Parse quote units with the same source conventions as scalar/batch input."""
    from core.units.registry import parse_unit
    return parse_unit(expression, context=context if context is not None else QuoteParsingContext())





parse_market_quote = parse_quote


__all__ = ('QuoteParsingContext', 'QuoteTokenRule', 'parse_market_quote', 'parse_quote', 'parse_quote_unit')
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)

# Static declarations for lazy public interfaces.
from typing import TYPE_CHECKING, Any, Callable
if TYPE_CHECKING:
    from core.units.currencies import CurrencyCatalog
