"""Commodex registry implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Iterable, Mapping
from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import dataclass, replace
from decimal import Decimal, InvalidOperation
from functools import cache, wraps
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
from typing import Any, Final, Iterator
import ast
import copy
import json
import keyword
import math
import pint
import re
import threading
from pint.facets.nonmultiplicative.definitions import OffsetConverter
from core.units.exceptions import CommodityUnitError

class _DecimalOffsetConverter(OffsetConverter):
    """Preserve Decimal affine arithmetic on older supported Pint releases."""

    def to_reference(self, value, inplace=False):
        if isinstance(value, Decimal):
            return value * Decimal(str(self.scale)) + Decimal(str(self.offset))
        return super().to_reference(value, inplace)

    def from_reference(self, value, inplace=False):
        if isinstance(value, Decimal):
            return (value - Decimal(str(self.offset))) / Decimal(str(self.scale))
        return super().from_reference(value, inplace)

class ForeignRegistryError(CommodityUnitError, ValueError):
    """Raised when a Pint Quantity or Unit belongs to a foreign UnitRegistry."""

from core.units.exceptions import RegistryVersionMismatchError

class RegistryCompatibilityError(RegistryVersionMismatchError):
    """Raised when registry schema, version, definition hash, or required units mismatch."""

# Global lock to ensure thread-safe dynamic registration of custom units and currencies
_REGISTRY_LOCK = threading.RLock()
_REGISTRY_FROZEN: bool = False

class _FreezableDict(dict):
    """Dictionary supporting startup registration that seals permanently on freeze_registry()."""

    def _ensure_mutable(self) -> None:
        from core.units.exceptions import RegistryFrozenError
        if _REGISTRY_FROZEN:
            raise RegistryFrozenError("UnitRegistry is frozen; registry metadata mappings cannot be mutated.")

    def _restore(self, snapshot: Mapping[Any, Any]) -> None:
        dict.clear(self)
        dict.update(self, snapshot)

    def __setitem__(self, key: Any, value: Any) -> None:
        self._ensure_mutable()
        super().__setitem__(key, value)

    def __delitem__(self, key: Any) -> None:
        self._ensure_mutable()
        super().__delitem__(key)

    def pop(self, *args: Any) -> Any:
        self._ensure_mutable()
        return super().pop(*args)

    def popitem(self) -> tuple[Any, Any]:
        self._ensure_mutable()
        return super().popitem()

    def clear(self) -> None:
        self._ensure_mutable()
        super().clear()

    def update(self, *args: Any, **kwargs: Any) -> None:
        self._ensure_mutable()
        super().update(*args, **kwargs)

    def setdefault(self, key: Any, default: Any = None) -> Any:
        if key in self:
            return self[key]
        self._ensure_mutable()
        return super().setdefault(key, default)

# ---------------------------------------------------------------------------
# ISO-4217 Currency Codes & Reference Data Integrity
# ---------------------------------------------------------------------------

REGISTRY_SCHEMA_VERSION: Final[int] = 1
from core.units.currencies import CURRENCY_CATALOG
REGISTRY_VERSION: Final[str] = f"commodex-units-v{REGISTRY_SCHEMA_VERSION}+{CURRENCY_CATALOG.version}"
_DATA_ROOT: Final[Path] = (Path(__file__).resolve().parents[0] / "reference_data")

from core.units.currencies import DATA_DIRECTORY
REFERENCE_DATA_HASHES: Final[Mapping[str, str]] = MappingProxyType({
    "currencies_manifest": sha256((DATA_DIRECTORY / "manifest.json").read_bytes()).hexdigest(),
    "currencies_reference": sha256((DATA_DIRECTORY / "reference.json").read_bytes()).hexdigest(),
    "currencies_conventions": sha256((DATA_DIRECTORY / "conventions.json").read_bytes()).hexdigest(),
    "currencies_overlay": sha256((DATA_DIRECTORY / "overlay.json").read_bytes()).hexdigest(),
    "energy_specs": sha256((_DATA_ROOT / "units" / "energy_specs.json").read_bytes()).hexdigest(),
    "assay_reference_manifest": sha256(((Path(__file__).resolve().parents[0] / "reference_data") / "standards" / "2026-11-27" / "manifest.json").read_bytes()).hexdigest(),
    "gas_components": sha256((_DATA_ROOT / "assays" / "gas_components.json").read_bytes()).hexdigest(),
    "astm_d2161": sha256((_DATA_ROOT / "units" / "astm_d2161.json").read_bytes()).hexdigest(),
    "icumsa_brix": sha256((_DATA_ROOT / "units" / "icumsa_brix.json").read_bytes()).hexdigest(),
})

from core.units.currencies import CURRENCY_CATALOG

# Current ISO codes and all parseable reference denominations are separate views.
CURRENCY_CODES: Final[tuple[str, ...]] = tuple(sorted({record.code for record in CURRENCY_CATALOG.records if record.status == "current" and record.kind == "iso4217"}))
from core.units.currencies import CURRENCY_CATALOG
CURRENCY_REFERENCE_CODES: Final[tuple[str, ...]] = tuple(sorted({record.code for record in CURRENCY_CATALOG.records}))
from core.units.currencies import CURRENCY_CONVENTIONS
CURRENCY_SYMBOL_TO_CODE: Final[dict[str, str]] = _FreezableDict(CURRENCY_CONVENTIONS["symbols"])
from core.units.currencies import CURRENCY_CATALOG

# ISO decimal exponents are not physical coin subdivisions or cash increments.
CURRENCY_MINOR_UNIT_EXPONENT: Final[Mapping[str, int]] = MappingProxyType({record.code: record.minor_digits for record in CURRENCY_CATALOG.records if record.minor_digits is not None})
CURRENCY_MINOR_DIGITS: dict[str, int] = _FreezableDict(CURRENCY_MINOR_UNIT_EXPONENT)
_ZERO_MINOR_UNIT = frozenset(code for code, digits in CURRENCY_MINOR_DIGITS.items() if digits == 0)
_THREE_MINOR_UNIT = frozenset(code for code, digits in CURRENCY_MINOR_DIGITS.items() if digits == 3)
_NON_DECIMAL_MINOR_UNIT = frozenset({"MGA", "MRU"})  # physical fifth subdivisions, not ISO precision

# Each monetary denomination has an independent dimension. Only undated,
# unscaled aliases in the default namespace are safe in the generic parser.
CURRENCY_UNIT_DEFINITIONS: Final[tuple[str, ...]] = tuple(
    f"{code} = [currency_{code}]" + "".join(
        f" = {alias.alias}" for alias in CURRENCY_CATALOG.aliases
        if alias.code == code and alias.unit == code and alias.scale == 1
        and alias.vendor is None and alias.valid_from is None and alias.valid_to is None
    ) for code in CURRENCY_REFERENCE_CODES
)

from core.units.currencies import CURRENCY_CONVENTIONS
# Constructs Pint definitions for the minor units of each ISO 4217 currency, allowing for precise representation of fractional currency amounts in calculations.
CURRENCY_MINOR_UNIT_DEFINITIONS: Final[tuple[str, ...]] = tuple(CURRENCY_CONVENTIONS["minor_definitions"])

# ---------------------------------------------------------------------------
# Master Commodity Definitions
# ---------------------------------------------------------------------------

# Mutable hook for desk configurations or application settings prior to registry build
CUSTOM_COMMODITY_UNIT_DEFINITIONS: tuple[str, ...] = ()

from core.units.standards import _reference_reference_document

COMMODITY_UNIT_DEFINITIONS: Final[tuple[str, ...]] = (
    *CURRENCY_UNIT_DEFINITIONS, *CURRENCY_MINOR_UNIT_DEFINITIONS,
    *(record["definition"] for record in _reference_reference_document("units", "registry_units")["records"]),
    *CUSTOM_COMMODITY_UNIT_DEFINITIONS,
)

# ---------------------------------------------------------------------------
# Registry Singleton & Serialisation Bindings
# ---------------------------------------------------------------------------

# Only these observed Pint-default collisions are deliberate. New collisions
# (including desk extensions) must fail instead of silently replacing a unit.
PINT_COLLISION_ALLOWLIST: Final[Mapping[str, frozenset[str]]] = MappingProxyType({
    "USD_mill": frozenset({"mil"}), "basis_point": frozenset({"bp", "bps"}),
    "percent": frozenset({"%", "percent"}), "btu_it": frozenset({"BTU", "Btu"}),
    "quad_it": frozenset({"quad"}), "therm_it": frozenset({"th", "therm"}),
    "toe": frozenset({"toe", "tonne_of_oil_equivalent"}),
    "cal_it": frozenset({"cal", "cal_it"}), "bbl": frozenset({"bbl"}),
    "gallon_petroleum": frozenset({"Gal", "gal"}), "mt": frozenset({"tonne"}),
    "kmt": frozenset({"kt"}), "commodity_short_hundredweight": frozenset({"cwt", "hundredweight"}),
    "commodity_short_ton": frozenset({"short_ton"}), "commodity_long_ton": frozenset({"long_ton"}),
    "toz": frozenset({"ozt", "toz", "troy_ounce"}), "carat": frozenset({"carat", "ct"}),
    "grain": frozenset({"gr", "grain"}), "pennyweight": frozenset({"pennyweight"}),
    "dwt": frozenset({"dwt"}), "commodity_peck": frozenset({"peck", "pk"}),
    "board_foot": frozenset({"BF", "FBM", "board_feet", "board_foot"}),
})

_BOOTSTRAP_COLLISIONS: Final[Mapping[str, frozenset[str]]] = PINT_COLLISION_ALLOWLIST
_FORBIDDEN_DEFINITION_TOKENS: Final[frozenset[str]] = frozenset({"ppt", "mu", "ton_mile"})
_MAX_UNIT_STRING_LENGTH: Final[int] = 256
_MAX_UNIT_EXPONENT: Final[int] = 12

def _compute_definition_hash(definitions: Iterable[str], custom_registrations: Iterable[str] = (),) -> str:
    """Compute a deterministic SHA-256 fingerprint of registry definitions and reference snapshots."""
    from core.units.currencies import CURRENCY_CATALOG
    canonical_payload = {
        "schema_version": REGISTRY_SCHEMA_VERSION,
        "registry_version": REGISTRY_VERSION,
        "currency_data_version": CURRENCY_CATALOG.version,
        "reference_data_hashes": dict(sorted(REFERENCE_DATA_HASHES.items())),
        "collision_allowlist": {
            key: sorted(tokens) for key, tokens in sorted(PINT_COLLISION_ALLOWLIST.items())
        },
        "commodity_definitions": list(definitions),
        "custom_registrations": list(custom_registrations),
    }
    encoded = json.dumps(canonical_payload, sort_keys=True, separators=(",", ":"), ensure_ascii=True).encode("utf-8")
    return sha256(encoded).hexdigest()

BASELINE_DEFINITION_HASH: Final[str] = _compute_definition_hash(COMMODITY_UNIT_DEFINITIONS, ())

def _clone_context_chain(chain: pint.facets.context.objects.ContextChain,) -> pint.facets.context.objects.ContextChain:
    cloned = pint.facets.context.objects.ContextChain()
    cloned.contexts = list(chain.contexts)
    cloned.maps = list(chain.maps)
    cloned._graph = None
    return cloned

class _SemanticCompatibility:
    def is_compatible_with(self, other, *contexts, **kwargs):
        from core.units.validation import _safety_SEMANTIC_SCOPE
        if not super().is_compatible_with(other, *contexts, **kwargs):
            return False
        if _safety_SEMANTIC_SCOPE.get():
            return True
        registry = self._REGISTRY
        source = self.units if isinstance(self, pint.Quantity) else self
        if isinstance(other, str):
            target = registry.parse_expression(other).units
        elif isinstance(other, pint.Quantity):
            target = other.units
        elif isinstance(other, pint.Unit):
            target = other
        else:
            target = registry.dimensionless
        return registry._semantic_tags(source._units) == registry._semantic_tags(target._units)

class _SemanticUnit(_SemanticCompatibility, pint.Unit):
    pass

class _SemanticQuantity(_SemanticCompatibility, pint.Quantity):
    pass

class _GuardedUnitRegistry(pint.UnitRegistry):
    """Guard ordinary Pint APIs, isolate active contexts per thread/task, and seal on freeze."""
    _sealed = False
    Unit = _SemanticUnit
    Quantity = _SemanticQuantity

    def __init__(self, *args, **kwargs):
        object.__setattr__(self, "_sealed", False)
        object.__setattr__(self, "_bootstrap_complete", False)
        object.__setattr__(self, "_in_currency_registration", False)
        object.__setattr__(self, "_in_internal_redefine", False)
        object.__setattr__(self, "_semantic_tag_cache", {})
        object.__setattr__(self, "_custom_registrations", ())
        object.__setattr__(self, "_definition_hash", "")
        object.__setattr__(self, "_active_ctx_var", ContextVar(f"commodex_active_ctx_{id(self)}", default=None),)
        super().__init__(*args, **kwargs)

    @property
    def _active_ctx(self):
        var = self.__dict__["_active_ctx_var"]
        chain = var.get()
        if chain is None:
            chain = pint.facets.context.objects.ContextChain()
            var.set(chain)
        return chain

    @_active_ctx.setter
    def _active_ctx(self, value):
        from core.units.exceptions import RegistryFrozenError
        if self.__dict__.get("_sealed", False):
            raise RegistryFrozenError("UnitRegistry is frozen; active contexts cannot be reassigned directly.")
        self.__dict__["_active_ctx_var"].set(value)

    def __setattr__(self, name, value):
        from core.units.exceptions import RegistryFrozenError
        if self.__dict__.get("_sealed", False):
            raise RegistryFrozenError("UnitRegistry is frozen; attributes and conversion policies cannot be reassigned.")
        super().__setattr__(name, value)

    def __delattr__(self, name):
        from core.units.exceptions import RegistryFrozenError
        if self.__dict__.get("_sealed", False):
            raise RegistryFrozenError("UnitRegistry is frozen; attributes cannot be deleted.")
        super().__delattr__(name)

    def _mutable(self):
        from core.units.exceptions import RegistryFrozenError
        if self.__dict__.get("_sealed", False):
            raise RegistryFrozenError("UnitRegistry is frozen.")

    def _redefine(self, definition):
        with _REGISTRY_LOCK:
            self._mutable()
            object.__setattr__(self, "_in_internal_redefine", True)
            try:
                return super()._redefine(definition)
            finally:
                object.__setattr__(self, "_in_internal_redefine", False)

    def define(self, definition):
        with _REGISTRY_LOCK:
            self._mutable()
            post_bootstrap = self.__dict__.get("_bootstrap_complete", False)
            in_redefine = self.__dict__.get("_in_internal_redefine", False)
            if post_bootstrap and not in_redefine:
                if not isinstance(definition, str):
                    raise TypeError("Post-bootstrap definitions must be validated unit definition strings.")
                allow_ccy = self.__dict__.get("_in_currency_registration", False)
                validate_unit_definition(definition, registry=self, allow_new_dimension=allow_ccy, allow_alias_directive=True, _allow_currency_dimension=allow_ccy,)
            result = super().define(definition)
            self.__dict__["_semantic_tag_cache"].clear()
            if post_bootstrap and not in_redefine and isinstance(definition, str):
                updated = (*self.__dict__["_custom_registrations"], definition.strip())
                object.__setattr__(self, "_custom_registrations", updated)
                object.__setattr__(self, "_definition_hash", _compute_definition_hash(COMMODITY_UNIT_DEFINITIONS, updated),)
            return result

    def load_definitions(self, file, is_resource=False):
        with _REGISTRY_LOCK:
            self._mutable()
            if not self.__dict__.get("_bootstrap_complete", False):
                return super().load_definitions(file, is_resource)
            if isinstance(file, (str, Path)) or is_resource:
                raise ValueError("Loading external Pint definition files or resources at runtime is prohibited; register validated definitions explicitly.")
            cleaned: list[str] = []
            for line in file:
                if not isinstance(line, str):
                    raise ValueError("All loaded definitions must be validated unit definition strings.")
                stripped = line.strip()
                if not stripped:
                    continue
                validate_unit_definition(stripped, registry=self, allow_new_dimension=False, allow_alias_directive=False,)
                cleaned.append(stripped)
            for stripped in cleaned:
                self.define(stripped)

    def add_context(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            self._mutable()
            return super().add_context(*args, **kwargs)

    def remove_context(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            self._mutable()
            return super().remove_context(*args, **kwargs)

    def _prepare_contexts(self, names_or_contexts, kwargs):
        from core.units.exceptions import RegistryFrozenError
        var = self.__dict__["_active_ctx_var"]
        current = var.get()
        if current is not None and current.defaults:
            merged_kwargs = dict(current.defaults, **kwargs)
        else:
            merged_kwargs = dict(kwargs)
        with _REGISTRY_LOCK:
            ctxs = [self._contexts[name] if isinstance(name, str) else name for name in names_or_contexts]
            for ctx in ctxs:
                if not isinstance(ctx, pint.facets.context.objects.Context):
                    raise TypeError(f"Expected pint.Context or registered context name, got {type(ctx).__name__}.")
                if ctx.redefinitions and self.__dict__.get("_sealed", False):
                    raise RegistryFrozenError("UnitRegistry is frozen; contexts with unit redefinitions are prohibited.")
                if not ctx.checked:
                    funcs_copy = dict(ctx.funcs)
                    for (src, dst), func in funcs_copy.items():
                        src_ = self._get_dimensionality(src)
                        dst_ = self._get_dimensionality(dst)
                        if src != src_ or dst != dst_:
                            ctx.remove_transformation(src, dst)
                            ctx.add_transformation(src_, dst_, func)
                    ctx.checked = True
            return tuple(pint.facets.context.objects.Context.from_context(ctx, **merged_kwargs) for ctx in ctxs)

    def enable_contexts(self, *names_or_contexts, **kwargs):
        with _REGISTRY_LOCK:
            self._mutable()
            prepared = self._prepare_contexts(names_or_contexts, kwargs)
            var = self.__dict__["_active_ctx_var"]
            previous = var.get()
            new_chain = (_clone_context_chain(previous) if previous is not None else pint.facets.context.objects.ContextChain())
            new_chain.insert_contexts(*prepared)
            var.set(new_chain)
            self._switch_context_cache_and_units()

    def disable_contexts(self, n=None):
        with _REGISTRY_LOCK:
            self._mutable()
            var = self.__dict__["_active_ctx_var"]
            previous = var.get()
            new_chain = (_clone_context_chain(previous) if previous is not None else pint.facets.context.objects.ContextChain())
            new_chain.remove_contexts(n)
            var.set(new_chain)
            self._switch_context_cache_and_units()

    @contextmanager
    def context(self, *names, **kwargs):
        prepared = self._prepare_contexts(names, kwargs)
        var = self.__dict__["_active_ctx_var"]
        previous = var.get()
        new_chain = (_clone_context_chain(previous) if previous is not None else pint.facets.context.objects.ContextChain())
        new_chain.insert_contexts(*prepared)
        has_redefs = any(ctx.redefinitions for ctx in new_chain.contexts) or (previous is not None and any(ctx.redefinitions for ctx in previous.contexts))
        token = var.set(new_chain)
        if has_redefs:
            with _REGISTRY_LOCK:
                self._switch_context_cache_and_units()
        try:
            yield self
        finally:
            var.reset(token)
            if has_redefs:
                with _REGISTRY_LOCK:
                    self._switch_context_cache_and_units()

    def get_name(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().get_name(*args, **kwargs)

    def get_dimensionality(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().get_dimensionality(*args, **kwargs)

    def parse_units(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().parse_units(*args, **kwargs)

    def parse_expression(self, *args, **kwargs):
        with _REGISTRY_LOCK:
            return super().parse_expression(*args, **kwargs)

    def _semantic_tags(self, container):
        cache = self.__dict__["_semantic_tag_cache"]
        if container in cache:
            return cache[container]
        tags = {}
        def walk(name, exponent, seen):
            if name in seen or name not in self._units:
                return
            definition = self._units[name]
            canonical = definition.name
            basis = {"wmt": "wet", "dmt": "dry", "dmtu": "contained", "contained_metal_tonne": "contained", "payable_metal_tonne": "payable", "admt": "air_dry", "bdmt": "bone_dry"}.get(canonical)
            pressure = canonical.lower()
            if basis:
                tags[basis] = tags.get(basis, 0) + exponent
                return
            if pressure in {"psig", "barg", "kpag", "mpag", "pag"}:
                tags["gauge"] = tags.get("gauge", 0) + exponent
                return
            for child, power in definition.reference.items():
                walk(child, exponent * power, seen | {name})
        for name, exponent in container.items():
            walk(name, exponent, set())
        result = {key: value for key, value in tags.items() if value}
        cache[container] = result
        return result

    def _project_mass_units(self, container):
        """Internal physical scaling for an already-authorized assay operation."""
        factor, roots = self._get_root_units(container)
        projected = dict(roots)
        for name in ("wmt", "dmt", "contained_metal_tonne", "payable_metal_tonne"):
            exponent = projected.pop(name, 0)
            if exponent:
                factor *= 1000 ** exponent
                projected["kilogram"] = projected.get("kilogram", 0) + exponent
        return factor, self.UnitsContainer({name: power for name, power in projected.items() if power})

    def _convert_mass_scale(self, value, src, dst):
        source_factor, source = self._project_mass_units(src)
        target_factor, target = self._project_mass_units(dst)
        number = Decimal(str(source_factor)) if isinstance(value, Decimal) else float(source_factor)
        result = super()._convert(value * number, source, target)
        divisor = Decimal(str(target_factor)) if isinstance(value, Decimal) else float(target_factor)
        return result / divisor

    def _convert(self, value, src, dst, *args, **kwargs):
        from core.units.validation import _safety_SEMANTIC_SCOPE
        with _REGISTRY_LOCK:
            if not _safety_SEMANTIC_SCOPE.get() and self._semantic_tags(src) != self._semantic_tags(dst):
                raise pint.DimensionalityError(src, dst, extra_msg="; pressure datum or mass-basis changes require a validated conversion/assay API")
            if _safety_SEMANTIC_SCOPE.get():
                semantic = {"wet", "dry", "contained", "payable"}
                source_tags = set(self._semantic_tags(src)) & semantic
                target_tags = set(self._semantic_tags(dst)) & semantic
                if source_tags and target_tags and source_tags != target_tags:
                    raise pint.DimensionalityError(src, dst, extra_msg="; use explicit moisture/grade/payability conversion")
                if source_tags or target_tags:
                    return self._convert_mass_scale(value, src, dst)
            return super()._convert(value, src, dst, *args, **kwargs)

def _new_registry() -> pint.UnitRegistry:
    registry = _GuardedUnitRegistry(autoconvert_offset_to_baseunit=True, on_redefinition="ignore")
    try:
        OffsetConverter(1, 273.15).to_reference(Decimal(0))
    except TypeError:
        # Replace only this registry's immutable definitions; preserve shared aliases.
        replacements = {}
        for name, definition in tuple(registry._units.items()):
            converter = definition.converter
            if isinstance(converter, OffsetConverter):
                if id(definition) not in replacements:
                    replacements[id(definition)] = replace(definition, converter=_DecimalOffsetConverter(converter.scale, converter.offset))
                registry._units[name] = replacements[id(definition)]
    for definition in COMMODITY_UNIT_DEFINITIONS:
        before = dict(registry._units)
        registry.define(definition)
        changed = {name for name, previous in before.items() if registry._units.get(name) is not previous}
        canonical = definition.split("=", 1)[0].strip()
        unexpected = changed if definition in CUSTOM_COMMODITY_UNIT_DEFINITIONS else changed - _BOOTSTRAP_COLLISIONS.get(canonical, frozenset())
        if unexpected:
            raise pint.errors.RedefinitionError(", ".join(sorted(unexpected)), "unit")
    registry._on_redefinition = "raise"
    object.__setattr__(registry, "_definition_hash", _compute_definition_hash(COMMODITY_UNIT_DEFINITIONS, ()))
    object.__setattr__(registry, "_bootstrap_complete", True)
    return registry

_MASTER_REGISTRY: pint.UnitRegistry | None = None

def get_registry() -> pint.UnitRegistry:
    """Atomically construct, bind, and publish the process-wide UnitRegistry singleton."""
    global _MASTER_REGISTRY
    published = _MASTER_REGISTRY
    if published is not None:
        return published
    with _REGISTRY_LOCK:
        if _MASTER_REGISTRY is None:
            candidate = _new_registry()
            pint.set_application_registry(candidate)
            _MASTER_REGISTRY = candidate
        return _MASTER_REGISTRY

def _locked_registry(function):
    @wraps(function)
    def locked(*args, **kwargs):
        with _REGISTRY_LOCK:
            return function(*args, **kwargs)
    return locked

# Module-level singleton instance and Quantity constructor export
ureg: Final[pint.UnitRegistry] = get_registry()
Q_ = ureg.Quantity

# ---------------------------------------------------------------------------
# Parsing and Normalisation Infrastructure
# ---------------------------------------------------------------------------

def _has_special_symbol_chars(symbol: str) -> bool:
    """Returns True if symbol contains characters other than letters, numbers, or underscores."""
    for char in symbol: # For every character in the symbol
        if not char.isalnum() and char != "_": # Check if it's non-acceptable characters
            return True # Return True if it is
    return False # Return False if all characters are valid identifier characters

_SAFE_SYMBOLS: Final[dict[str, str]] = {} # Initialise empty dictionary

for symbol, code in CURRENCY_SYMBOL_TO_CODE.items(): # Unpacks pairs from master currency dictionary (e.g. "$", "USD")
    if _has_special_symbol_chars(symbol) or symbol in {"kr", "Rs", "Fr", "R"}: # True for raw currency glyphs (e.g. $, £, C$), False for alphanumeric identifiers (USD, GBP)
        _SAFE_SYMBOLS[symbol] = code # Symbol code pair into _SAFE_SYMBOLS

_SORTED_SAFE_SYMBOLS: tuple[str, ...] = tuple(sorted(_SAFE_SYMBOLS.keys(), key=len, reverse=True)) # Sorts symbols by length descending to match multi-character prefixes before single glyphs

_CURRENCY_CODES_SET: Final[set[str]] = set(CURRENCY_REFERENCE_CODES)
_CURRENCY_CANONICAL_CODES: Final[dict[str, str]] = {code: code for code in CURRENCY_REFERENCE_CODES}
_CURRENCY_DEFINITION_LINES: Final[frozenset[str]] = (frozenset(CURRENCY_UNIT_DEFINITIONS) | frozenset(CURRENCY_MINOR_UNIT_DEFINITIONS)) #Composite set that merges all primary currency dimensional definitions with their corresponding fractional minor unit definitions

#translation table mapping colloquial trading screen shorthands, sub-unit glyphs, and vendor ticker conventions to their canonical Pint minor currency units
minor_mapping: Final[dict[str, str]] = _FreezableDict({
    "zac": "ZAR_cent",
    "zarc": "ZAR_cent",
    "c": "USD_cent",
    "¢": "USD_cent",
    "usc": "USD_cent",
    "us¢": "USD_cent",
    "p": "GBP_penny",
    "penny": "GBP_penny",
    "gbx": "GBP_penny",
    "€c": "EUR_cent",
    "eurc": "EUR_cent",
    "ac": "AUD_cent",
    "a¢": "AUD_cent",
    "cc": "CAD_cent",
    "c¢": "CAD_cent",
})

def _register_minor_unit_definition(definition: str) -> None:
    """Populates minor_mapping from a Pint currency minor-unit definition string."""
    parts = [part.strip() for part in definition.split("=")] #Splits input definition string e.g. "USD_cent = 0.01 * USD = US_cent = cent", parts becomes ["USD_cent", "0.01 * USD", "US_cent", "cent"]
    if not parts or not parts[0]: #evaluates whether the list is empty or if the leading element is empty string
        return #terminates execution if the definition contains no valid primary identifier
    canonical_identifier = parts[0] #Assigns first segment to the variable canonical
    minor_mapping[canonical_identifier] = canonical_identifier #Self-referencing key-value entry into minor_mapping
    minor_mapping[canonical_identifier.lower()] = canonical_identifier #Inserts lowercase representation of the canonical name
    for alias in parts[2:]: #Iterates over all segments from index 2 onward skipping 1 intentionally
        if not alias or not alias.isidentifier():  #Evaluates whethoer alias is empty or fails isidentifier
            continue
        minor_mapping[alias] = canonical_identifier #Maps exact-cased alias to canonical minor unit identifier
        if alias.upper() not in _CURRENCY_CODES_SET: #Checks whether the uppercase form of the alias colides with sovereign currency code
            minor_mapping[alias.lower()] = canonical_identifier #Maps lowercase form of the alias

# Establishing the two foundational lookup tables required to decompose and normalise <numerator>/<denominator> broker market quotes
for _minor_definition in CURRENCY_MINOR_UNIT_DEFINITIONS:
    _register_minor_unit_definition(_minor_definition)

irregular_denominators: Final[dict[str, str]] = _FreezableDict()

def _register_commodity_denominator_definition(definition: str) -> None:
    """This function extract valid unit identifiers and aliases from Pint definition strings
    (e.g. COMMODITY_UNIT_DEFINITIONS) and populates irregular_denominators so any commodity
    product denominator dynamically maps to its canonical unit"""
    parts = [part.strip() for part in definition.split("=")]  # Splits definition into components across '=' signs and trims whitespace
    if not parts or not parts[0]:  # Guard clause: returns early if definition string is empty or missing a leading token
        return  # Aborts processing for invalid or blank definition strings
    canonical_identifier = parts[0].split(":")[0].strip()  # Isolates base unit token, stripping away any dimensional metadata trailing a colon
    if not canonical_identifier.isidentifier():  # Ensures the extracted canonical token is a syntactically valid Python identifier
        return  # Skips registration if the token contains mathematical operators or invalid characters
    irregular_denominators[canonical_identifier] = canonical_identifier  # Maps exact-cased canonical identifier to itself for O(1) direct lookup
    irregular_denominators.setdefault(canonical_identifier.lower(), canonical_identifier)  # Registers lowercase variant without overwriting existing entries
    for alias in parts[2:]:  # Iterates through aliases starting at index 2, skipping index 1 (the mathematical conversion ratio)
        if not alias or not alias.isidentifier():  # Validates that alias token exists and is a clean alphanumeric identifier
            continue  # Discards blank entries and complex arithmetic expressions containing exponents or operators
        irregular_denominators[alias] = canonical_identifier  # Maps exact-cased alias directly back to canonical identifier
        irregular_denominators.setdefault(alias.lower(), canonical_identifier)  # Preserves case-insensitive lookup for broker screen shorthand inputs

for _commodity_definition in COMMODITY_UNIT_DEFINITIONS:  # Iterates sequentially through every unit definition string in the master commodity definitions tuple
    if _commodity_definition not in _CURRENCY_DEFINITION_LINES:  # Filters out currency and minor currency definitions using the composite set
        _register_commodity_denominator_definition(_commodity_definition)  # Parses and registers valid commodity identifiers and aliases into irregular_denominators

# Canonical market denominator overrides aligned with _QUOTE_REPLACEMENTS and _UNIT_TOKEN_REPLACEMENTS
for _denomination_aliases, _canonical_denomination in (  # Iterates through pairs of colloquial broker aliases and their canonical Pint units
    (("bbl", "BBL", "barrel_oil", "bbls", "barrel", "barrels"), "oil_barrel"),  # Standardises 42-gallon petroleum barrel variants, preventing Pint's dry-barrel collision
    (("gal", "Gal", "us_gal", "gallon_petroleum", "gallons"), "gallon"),  # Maps US liquid gallon screen abbreviations and plurals to standard volume
    (("MMBTU", "mmbtu", "MMBtu", "mmBtu", "MMBtu_IT"), "MMBtu_it"),  # Anchors natural gas benchmark tickers to the International Table energy basis
    (("MWh", "mwh"), "megawatt_hour"),  # Normalises wholesale electricity baseload contract blocks
    (("kWh", "kwh"), "kilowatt_hour"),  # Normalises retail and small-scale power generation units
    (("GWh", "gwh"), "gigawatt_hour"),  # Normalises transmission grid and utility-scale balance units
    (("th", "therm", "Therm", "THERM", "th_it"), "therm_it"),  # Maps retail/pipeline gas therm abbreviations directly to IT therms
    (("MT", "mt", "t", "metric_tonne", "tonne", "Tonne", "TONNE", "tonnes"), "tonne"),  # Normalises metric tonne bulk freight, grain, and metals mass units
    (("bu", "bushel", "bushels"), "bushel"),  # Maps standard agricultural grain contract volume/mass abbreviations
    (("lb", "lbs", "pound", "pounds"), "pound"),  # Maps imperial avoirdupois pound mass tickers and plurals
    (("toz", "oz", "ozt", "OZT", "TOZ", "troy_ounce", "troy_oz"), "toz"),  # Forces generic bullion 'oz' and troy ounce variants to LBMA/COMEX troy ounces
    (("cbm", "CBM", "cbm_oil", "m3_oil", "m3"), "meter ** 3"),  # Expands maritime freight and metric storage volume shorthand into valid algebraic syntax
    (("kg", "kgs"), "kilogram"),  # Maps SI base mass abbreviations and informal plural forms
    (("g", "gm"), "gram"),  # Normalises precious/minor metals gram abbreviations
    (("l", "ltr", "litre", "litres", "liters"), "liter"),  # Normalises metric liquid volume abbreviations across US and UK spellings
):
    for _alias in _denomination_aliases:  # Iterates across each individual alias string defined for the current commodity
        irregular_denominators[_alias] = _canonical_denomination  # Maps exact-cased alias token directly to canonical Pint unit for O(1) lookup
        if _alias.lower() != "mt" or _alias in ("mt", "MT"):  # Case collision guard: prevents mixed-case 'Mt' (megatonnes/million tonnes) from corrupting metric tonnes ('mt')
            irregular_denominators[_alias.lower()] = _canonical_denomination  # Registers lowercase version to enable case-insensitive screen quote matching

# ---------------------------------------------------------------------------
# Dynamic Quote Resolution Engine
# ---------------------------------------------------------------------------

def _resolve_quote_expression(quote: str) -> str | None:
    """Converts a quote suchas "USD/BBL" or "gbx per therm" into a parsable Pint unit expression,
      such as "USD / barrel". If it cannot identify or validate the quote it returns "None" """
    
    clean_quote = quote.strip() #Removes leading/trailing whitespaces from the string

    # Split quote into numerator and denominator across '/' or 'per'
    if "/" in clean_quote: 
        raw_numerator, raw_denomination = clean_quote.split("/", 1) #If "/" is in the quote, split it by "/" and set first section as number, second as denomination
    elif " per " in clean_quote.lower(): #If it uses "per", same thing as "/"
        parts = re.split(r"\s+per\s+", clean_quote, maxsplit=1, flags=re.I) #splits clean_quote wherever it finds the word per with whitespace on both sides, ignoring capitalisation.
        raw_numerator, raw_denomination = parts[0], parts[1] #Sets numerator and denominator
    else:
        return None

    num = raw_numerator.strip() #Strips incase of any whitespaces
    denomination = raw_denomination.strip() #Strips incase of any whitespaces
    # -----------------------------------------------------------------------
    # 1. Resolve Currency Numerator
    # -----------------------------------------------------------------------
    num_norm: str | None = None
    num_lower = num.lower()
    if num in minor_mapping:
        # Exact match in sub-unit mapping (e.g. "gbx" -> "GBP_penny")
        num_norm = minor_mapping[num]
    elif num_lower in minor_mapping:
        # Case-insensitive sub-unit match (e.g. "cents" -> "USD_cent")
        num_norm = minor_mapping[num_lower]
    elif num in CURRENCY_SYMBOL_TO_CODE:
        # Single-glyph or composite currency symbol (e.g. "$" -> "USD", "€" -> "EUR")
        num_norm = CURRENCY_SYMBOL_TO_CODE[num]
    elif num.upper() in _CURRENCY_CODES_SET:
        # Standard ISO-4217 3-letter currency code (e.g. "usd" -> "USD")
        num_norm = _CURRENCY_CANONICAL_CODES[num.upper()]
    elif num == "¥":
        # Disambiguate generic yen/yuan glyph to Japanese Yen by default
        num_norm = "JPY"
    elif num.endswith(("_cent", "_penny", "_fils", "_minor")):
        # Pre-canonicalised minor currency units pass through directly
        num_norm = num
    # Abort resolution if the numerator token cannot be identified
    if num_norm is None:
        return None
    # -----------------------------------------------------------------------
    # 2. Resolve Commodity Denominator
    # -----------------------------------------------------------------------
    # Lookup in irregular denominations table with case-insensitive fallback
    denomination_norm = irregular_denominators.get(denomination, irregular_denominators.get(denomination.lower(), denomination),)
    for pattern, replacement in _UNIT_TOKEN_REPLACEMENTS:
        denomination_norm = pattern.sub(replacement, denomination_norm)
    # Reject injection of invalid characters outside standard Pint identifiers,
    # spaces, and exponent/multiplication operators (e.g. "meter ** 3")
    sanitised_check = (denomination_norm.replace("_", "").replace(" ", "").replace("*", ""))

    if not sanitised_check.isalnum():
        return None
    expression = f"{num_norm} / {denomination_norm}"
    try:
        ureg.Unit(expression)
    except (pint.PintError, ValueError, TypeError):
        return None
    return expression

class QuoteReplacementsMap(dict[str, str]):
    """Dynamic, self-caching dictionary that can recognise quote keys it hasn’t stored yet. 
    It uses _resolve_quote_expression from the previous example to parse and normalise them

    __missing__: Resolves and caches a quote when its key is not already in the dictionary.

    __contains__: Checks whether a quote is stored or can be resolved.

    get: Returns a quote’s value, or a default if it cannot be resolved.
    """

    def _ensure_mutable(self) -> None:
        if _REGISTRY_FROZEN:
            raise RuntimeError("UnitRegistry is frozen; quote replacement mappings cannot be mutated.")

    def _restore(self, snapshot: Mapping[str, str]) -> None:
        dict.clear(self)
        dict.update(self, snapshot)

    def __setitem__(self, key: str, value: str) -> None:
        self._ensure_mutable()
        super().__setitem__(key, value)

    def __delitem__(self, key: str) -> None:
        self._ensure_mutable()
        super().__delitem__(key)

    def pop(self, *args: Any) -> Any:
        self._ensure_mutable()
        return super().pop(*args)

    def popitem(self) -> tuple[str, str]:
        self._ensure_mutable()
        return super().popitem()

    def clear(self) -> None:
        self._ensure_mutable()
        super().clear()

    def update(self, *args: Any, **kwargs: Any) -> None:
        self._ensure_mutable()
        super().update(*args, **kwargs)

    def setdefault(self, key: str, default: Any = None) -> Any:
        if super().__contains__(key):
            return dict.__getitem__(self, key)
        self._ensure_mutable()
        return super().setdefault(key, default)

    def __missing__(self, key: str) -> str:
        with _REGISTRY_LOCK: #Locks the registry while it resolves and possibly stores the quote, helping prevent simultaneous lookups from interfering.
            if dict.__contains__(self, key): #Checks whether the key is stored. Another lookup may have added it while this call was waiting for the lock.
                return dict.__getitem__(self, key) #If it is now stored, returns its value directly.
            resolved = _resolve_quote_expression(key) #Tries to parse and normalise the quote key.
            if resolved is not None: #Checks whether parsing succeeded.
                if not _REGISTRY_FROZEN:
                    dict.__setitem__(self, key, resolved) #Stores the resolved expression in the dictionary while unfrozen.
                return resolved #Returns the resolved expression.
            raise KeyError(key) 

    def __contains__(self, key: object) -> bool:
        if super().__contains__(key): #Checks whether the key is already stored in the dictionary.
            return True #Returns True if it is stored.
        return isinstance(key, str) and _resolve_quote_expression(key) is not None 
        #Otherwise, checks that the key is a string and tries to resolve it. 
        #Returns True if it can be resolved. This check does not store the result.

    def get(self, key: str, default: Any = None) -> Any:
        try: #Looks up the key. If it is missing, this invokes __missing__, which tries to resolve and cache it.
            return self[key]
        except KeyError:
            return default #Returns the supplied default value instead of raising KeyError.

# Master lookup map populated with baseline institutional benchmarks
_QUOTE_REPLACEMENTS: Final[QuoteReplacementsMap] = QuoteReplacementsMap({
    # Liquids & Refined Products
    "$/BBL": "USD / oil_barrel",
    "$/bbl": "USD / oil_barrel",
    "$/gal": "USD / gallon",
    "c/gal": "USD_cent / gallon",
    "cents/gal": "USD_cent / gallon",
    # Gas & Power
    "$/MMBTU": "USD / MMBtu_it",
    "$/mmbtu": "USD / MMBtu_it",
    "$/MWh": "USD / megawatt_hour",
    "p/th": "GBP_penny / therm_it",
    "p/therm": "GBP_penny / therm_it",
    "£/MWh": "GBP / megawatt_hour",
    "€/MWh": "EUR / megawatt_hour",
    # Agriculture & Freight
    "$/MT": "USD / tonne",
    "$/mt": "USD / tonne",
    "c/bu": "USD_cent / bushel",
    "cents/bu": "USD_cent / bushel",
    "c/lb": "USD_cent / pound",
    "cents/lb": "USD_cent / pound",
    "$/lb": "USD / pound",
    # Bullion & Precious Metals
    "$/toz": "USD / toz",
    "$/oz": "USD / toz",
    # Compliance Carbon & Emissions
    "€/t": "EUR / tonne_CO2",
    "€/EUA": "EUR / EU_allowance",
    "£/t": "GBP / tonne_CO2",
    "£/UKA": "GBP / UK_allowance",
    "$/CCA": "USD / California_carbon_allowance",
})

# ---------------------------------------------------------------------------
# Pre-Parsing Token Interventions & Regex Matchers
# ---------------------------------------------------------------------------
# Strict negative lookbehind/lookahead patterns to bypass Pint's prefix engine
# without corrupting compound variable names (e.g. daily_bbl_rate)
_UNIT_TOKEN_REPLACEMENTS: Final[list[tuple[re.Pattern[str], str]]] = [
    (re.compile(r"(?<![A-Za-z0-9_])psia(?![A-Za-z0-9_])", re.I), "psia"),
    (re.compile(r"(?<![A-Za-z0-9_])psig(?![A-Za-z0-9_])", re.I), "psig"),
    # Composite expressions must precede their individual token replacements.
    (re.compile(r"(?<![A-Za-z0-9_])MMBtu/MWh(?![A-Za-z0-9_])", re.I), "MMBtu_it / megawatt_hour"),
    (re.compile(r"(?<![A-Za-z0-9_])Btu/kWh(?![A-Za-z0-9_])", re.I), "btu_it / kilowatt_hour"),
    (re.compile(r"(?<![A-Za-z0-9_])bbl(?![A-Za-z0-9_])", re.I), "oil_barrel"),
    (re.compile(r"(?<![A-Za-z0-9_])MMbbl(?![A-Za-z0-9_])", re.I), "MMbbl"),
    (re.compile(r"(?<![A-Za-z0-9_])Mbbl(?![A-Za-z0-9_])", re.I), "kbbl"),
    (re.compile(r"(?<![A-Za-z0-9_])MMBtu(?![A-Za-z0-9_])", re.I), "MMBtu_it"),
    (re.compile(r"(?<![A-Za-z0-9_])MMcf(?![A-Za-z0-9_])", re.I), "MMscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])Mcf(?![A-Za-z0-9_])", re.I), "Mscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])Bcf(?![A-Za-z0-9_])", re.I), "Bscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])Tcf(?![A-Za-z0-9_])", re.I), "Tscf_60F_14_73psia"),
    (re.compile(r"(?<![A-Za-z0-9_])MMboe(?![A-Za-z0-9_])", re.I), "MMboe"),
    (re.compile(r"(?<![A-Za-z0-9_])Mboe(?![A-Za-z0-9_])", re.I), "kboe"),
    (re.compile(r"(?<![A-Za-z0-9_])cwt(?![A-Za-z0-9_])", re.I), "commodity_short_hundredweight"),
    (re.compile(r"(?<![A-Za-z0-9_])cbm(?![A-Za-z0-9_])", re.I), "meter ** 3"),
]

# Require digits and valid three-digit comma groups, retaining decimal/exponent support.
#This builds a regular expression for matching a numeric magnitude, such as 1234, 1,234.50 or 2.5e-3.
_MAGNITUDE_PATTERN = r"(?:(?:\d{1,3}(?:,\d{3})+|\d+)(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?"

#Builds and stores a regular expression for matching quotes that start with a 
#currency symbol, such as -$1,234.50/bbl or ($1,234.50/bbl). 
#It prepares the matcher; it does not match a quote itself
# Matches "$1,250.50/bbl", "($45.20/bbl)", and "$.50/gal".
# Longer currency symbols must precede their single-character suffixes.
def _recompile_prefix_currency_pattern() -> None:
    """Refresh prefix matching after symbol additions; caller holds the registry lock."""
    global _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN
    symbols = tuple(sorted(_SAFE_SYMBOLS, key=len, reverse=True)) #Takes the symbols in _SAFE_SYMBOLS sorting them from longest to shortest and stores the rorted results in tuple called symbols
    #Compiles a regex that splits a currency quote into its optional brackets, sign, symbol amount and denomination
    pattern = re.compile(
        rf"^\s*(?P<bracket>\()?\s*(?P<sign>[+-]?)\s*"
        rf"(?P<symbol>{'|'.join(re.escape(s) for s in symbols)})\s*"
        rf"(?P<magnitude>{_MAGNITUDE_PATTERN})\s*"
        r"(?:/\s*(?P<denomination>.+?))?\s*(?(bracket)\))\s*$"
    )
    _SORTED_SAFE_SYMBOLS = symbols
    _PREFIX_CURRENCY_PATTERN = pattern

_recompile_prefix_currency_pattern()

# Matches "1,500 MWh", "(2,500) bbl", "2.5%", and "15 bps".
_NUMERIC_MAGNITUDE_PATTERN = re.compile(
    r"^\s*(?P<bracket>\()?\s*"
    rf"(?P<magnitude>[+-]?{_MAGNITUDE_PATTERN})"
    r"(?(bracket)\s*\))"
    r"(?:\s+|(?=[^\d\s,.)+\-]))(?P<unit>.+?)\s*$"
)

# Separate batch magnitudes must not silently discard an embedded unit.
_BARE_MAGNITUDE_PATTERN = re.compile(
    r"^\s*(?P<bracket>\()?\s*"
    rf"(?P<magnitude>[+-]?{_MAGNITUDE_PATTERN})"
    r"\s*(?(bracket)\))\s*$"
)

from core.units.currencies import CURRENCY_CONVENTIONS

# Symbols that can represent several currencies.
# Add every shared glyph or market alias your registry accepts.
_SHARED_SYMBOL_CURRENCIES: Final[Mapping[str, frozenset[str]]] = MappingProxyType({symbol: frozenset(codes) for symbol, codes in CURRENCY_CONVENTIONS["shared_symbols"].items()})

def _resolve_currency_symbol(symbol: str, default_currency: str | None,) -> str:
    """Resolve any registered currency symbol, using context for shared symbols."""
    try:
        code = CURRENCY_SYMBOL_TO_CODE[symbol]
    except KeyError:
        raise ValueError(f"Unknown currency symbol: {symbol!r}") from None
    
    if default_currency is None:
        return code
    if default_currency.upper() not in _CURRENCY_CANONICAL_CODES:
        raise ValueError(f"Unknown ISO 4217 currency code: {default_currency}")
    
    default_currency = _CURRENCY_CANONICAL_CODES[default_currency.upper()]
    valid_codes = _SHARED_SYMBOL_CURRENCIES.get(symbol, frozenset({code}))
    
    if default_currency not in valid_codes:
        raise ValueError(f"{symbol!r} is not a valid currency symbol for {default_currency}")
    return default_currency

# ---------------------------------------------------------------------------
# Foreign-Registry Validation, Snapshots, Compatibility & Definition Guards
# ---------------------------------------------------------------------------

def validate_registry_quantity(value: Any, *, allow_unit: bool = True, registry: pint.UnitRegistry | None = None,) -> Any:
    """Verify that a Pint Quantity or Unit belongs to the active Commodex UnitRegistry.

    Raises ForeignRegistryError when an object was constructed by a foreign
    ``pint.UnitRegistry()`` or an isolated ``sandbox_registry()``.
    """
    target = ureg if registry is None else registry
    if isinstance(value, pint.Quantity):
        if getattr(value, "_REGISTRY", None) is not target:
            raise ForeignRegistryError("Quantity belongs to a foreign UnitRegistry; all quantities must originate from the active Commodex registry.")
    elif allow_unit and isinstance(value, pint.Unit):
        if getattr(value, "_REGISTRY", None) is not target:
            raise ForeignRegistryError("Unit belongs to a foreign UnitRegistry; all units must originate from the active Commodex registry.")
    return value

def _validate_identifier_token(token: str, *, label: str) -> None:
    from core.units.exceptions import AmbiguousUnitError
    if not token or not token.isidentifier():
        raise ValueError(f"{label} {token!r} must be a valid identifier.")
    if keyword.iskeyword(token) or token in {"True", "False", "None"}:
        raise ValueError(f"{label} {token!r} cannot be a Python keyword.")
    if token.startswith("_"):
        raise ValueError(f"{label} {token!r} cannot start with an underscore.")
    if token in _FORBIDDEN_DEFINITION_TOKENS:
        raise AmbiguousUnitError(f"{label} {token!r} is a reserved ambiguous market token.")

def _check_token_collision(token: str, target_registry: pint.UnitRegistry, *, is_currency_registration: bool = False,) -> None:
    if token in target_registry._units:
        raise pint.errors.RedefinitionError(token, "unit")
    if target_registry is ureg and not is_currency_registration:
        if token.upper() in _CURRENCY_CODES_SET or token in CURRENCY_SYMBOL_TO_CODE or token in minor_mapping:
            raise pint.errors.RedefinitionError(token, "unit")

def _validate_definition_rhs_ast(rhs: str, target_registry: pint.UnitRegistry,) -> None:
    if any(ch in rhs for ch in ("^", "[", "]", ";", "#", "@", "=", "'", '"', "`", "\\")):
        raise ValueError(f"Unit definition RHS {rhs!r} contains prohibited syntax.")
    try:
        tree = ast.parse(rhs, mode="eval")
    except (SyntaxError, ValueError) as error:
        raise ValueError(f"Invalid algebraic unit definition RHS {rhs!r}.") from error

    unit_refs: list[str] = []

    def eval_scale(node: ast.AST) -> float:
        if isinstance(node, ast.Constant):
            if isinstance(node.value, bool) or not isinstance(node.value, (int, float)):
                raise ValueError("Unit definition constants must be real numbers.")
            val = float(node.value)
            if not math.isfinite(val) or val <= 0.0 or val > 1e30:
                raise ValueError(f"Numeric factor in unit definition must be finite and in (0, 1e30], got {node.value!r}.")
            return val
        if isinstance(node, ast.Name):
            name = node.id
            if keyword.iskeyword(name) or name.startswith("_"):
                raise ValueError(f"Invalid unit identifier {name!r} in definition RHS.")
            if name not in target_registry._units and not target_registry.parse_unit_name(name):
                raise ValueError(f"Unknown reference unit {name!r} in unit definition.")
            unit_refs.append(name)
            return 1.0
        if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.UAdd):
            return eval_scale(node.operand)
        if isinstance(node, ast.BinOp) and isinstance(node.op, (ast.Mult, ast.Div)):
            left = eval_scale(node.left)
            right = eval_scale(node.right)
            if right == 0.0:
                raise ValueError("Division by zero in unit definition.")
            return left * right if isinstance(node.op, ast.Mult) else left / right
        if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Pow):
            base_scale = eval_scale(node.left)
            exp_node = node.right
            sign = 1
            if isinstance(exp_node, ast.UnaryOp) and isinstance(exp_node.op, (ast.UAdd, ast.USub)):
                sign = -1 if isinstance(exp_node.op, ast.USub) else 1
                exp_node = exp_node.operand
            if (not isinstance(exp_node, ast.Constant) or isinstance(exp_node.value, bool) or not isinstance(exp_node.value, int)):
                raise ValueError("Exponents in custom unit definitions must be integer literals.")
            exponent = sign * exp_node.value
            if exponent == 0 or abs(exponent) > _MAX_UNIT_EXPONENT:
                raise ValueError(f"Unit exponent must be a non-zero integer within [-{_MAX_UNIT_EXPONENT}, {_MAX_UNIT_EXPONENT}].")
            return base_scale ** exponent
        raise ValueError("Custom unit definitions allow only positive numbers, existing units, *, /, and bounded integer ** powers.")

    net_scale = eval_scale(tree.body)
    if not unit_refs:
        raise ValueError("Custom unit definitions must reference at least one existing physical or financial unit.")
    if not math.isfinite(net_scale) or not (1e-36 <= net_scale <= 1e36):
        raise ValueError(f"Net scaling factor {net_scale} in unit definition is out of safe bounds [1e-36, 1e36].")

def validate_unit_definition(definition: str, *, registry: pint.UnitRegistry | None = None, allow_new_dimension: bool = False, allow_alias_directive: bool = False, _allow_currency_dimension: bool = False,) -> tuple[str, str, tuple[str, ...]]:
    """Validate a user- or config-supplied Pint definition string before registration.

    Blocks multi-line payloads, ``@import``/``@context``/``@system``/``@group`` directives,
    prefix definitions (``foo- = ...``), affine offset units (``; offset:``), comments,
    unapproved root dimensions, non-positive/extreme scales, and collisions with existing
    units or aliases. Returns ``(canonical_name, rhs, aliases)``.
    """
    if not isinstance(definition, str):
        raise ValueError("Unit definition must be a string.")
    if len(definition) > _MAX_UNIT_STRING_LENGTH:
        raise ValueError(f"Unit definition exceeds maximum length of {_MAX_UNIT_STRING_LENGTH} characters.")
    if any(ch in ("\n", "\r", "\x00") or (ord(ch) < 32 and ch != "\t") or ord(ch) == 127 for ch in definition):
        raise ValueError("Unit definition must be a single line without control characters.")
    stripped = definition.strip()
    if not stripped:
        raise ValueError("Unit definition must be a non-empty string.")
    if "#" in stripped:
        raise ValueError("Inline '#' comments are prohibited in runtime unit definitions.")
    if ";" in stripped or "offset:" in stripped.lower():
        raise ValueError("Affine/offset unit definitions ('; offset:') are prohibited at runtime.")

    target_registry = ureg if registry is None else registry

    if stripped.startswith("@"):
        if not allow_alias_directive or not stripped.startswith("@alias "):
            raise ValueError("Pint '@' directives (@import, @context, @alias, @system, @group) are prohibited in custom unit definitions.")
        body = stripped[len("@alias "):].strip()
        parts = [part.strip() for part in body.split("=")]
        if len(parts) < 2 or any(not p for p in parts):
            raise ValueError(f"Invalid @alias definition: {definition!r}")
        existing = parts[0]
        _validate_identifier_token(existing, label="Alias target")
        if existing not in target_registry._units:
            raise pint.errors.UndefinedUnitError(existing)
        aliases = tuple(parts[1:])
        seen_aliases: set[str] = {existing}
        for alias in aliases:
            _validate_identifier_token(alias, label="Unit alias")
            if alias in seen_aliases:
                raise ValueError(f"Duplicate identifier {alias!r} in @alias definition.")
            seen_aliases.add(alias)
            _check_token_collision(alias, target_registry)
        return existing, existing, aliases

    if "@" in stripped:
        raise ValueError("Pint '@' directives are prohibited in custom unit definitions.")

    parts = [part.strip() for part in stripped.split("=")]
    if len(parts) < 2 or len(parts) > 12 or any(not part for part in parts):
        raise ValueError(f"Invalid unit definition format: {definition!r}")

    canonical = parts[0]
    if "-" in canonical or ":" in canonical:
        raise ValueError("Prefix definitions and colon metadata are prohibited in custom unit definitions.")
    _validate_identifier_token(canonical, label="Unit name")
    _check_token_collision(canonical, target_registry, is_currency_registration=_allow_currency_dimension,)

    rhs = parts[1]
    raw_aliases = parts[2:]
    aliases_out: list[str] = []
    seen: set[str] = {canonical}
    for idx, alias in enumerate(raw_aliases):
        if idx == 0 and alias == "_":
            continue
        _validate_identifier_token(alias, label="Unit alias")
        if alias in seen:
            raise ValueError(f"Duplicate identifier {alias!r} in unit definition.")
        seen.add(alias)
        _check_token_collision(alias, target_registry, is_currency_registration=_allow_currency_dimension,)
        aliases_out.append(alias)

    if "[" in rhs or "]" in rhs:
        if not allow_new_dimension:
            raise ValueError("Defining new root dimensions ('[...]') is prohibited in custom unit definitions.")
        dim_match = re.fullmatch(r"\[([A-Za-z][A-Za-z0-9_]{0,63})\]", rhs)
        if dim_match is None:
            raise ValueError(f"Invalid root dimension syntax {rhs!r}.")
        dim_name = dim_match.group(1)
        if dim_name.startswith("currency_") and not _allow_currency_dimension:
            raise ValueError("Currency dimensions can only be created via register_currency().")
    else:
        _validate_definition_rhs_ast(rhs, target_registry)

    return canonical, rhs, tuple(aliases_out)

def _validate_external_unit_text(raw: str) -> None:
    """Reject directive, definition, multi-line, and overlong external unit/quote strings."""
    if not isinstance(raw, str):
        raise TypeError(f"Expected string unit expression, got {type(raw).__name__}.")
    if len(raw) > _MAX_UNIT_STRING_LENGTH:
        raise ValueError(f"Unit expression exceeds maximum length of {_MAX_UNIT_STRING_LENGTH} characters.")
    if any(ch in ("\n", "\r", "\x00") or (ord(ch) < 32 and ch != "\t") or ord(ch) == 127 for ch in raw):
        raise ValueError("Unit expression must not contain newlines or control characters.")
    if any(ch in raw for ch in ("@", "=", ";", "#", "[", "]", "'", '"', "`", "\\", "{", "}")):
        raise ValueError(f"Unit expression {raw!r} contains prohibited directive or definition syntax.")
    if "__" in raw:
        raise ValueError(f"Unit expression {raw!r} contains prohibited dunder identifier syntax.")

def _validate_sanitized_unit_expression(expr: str) -> None:
    """Reject function calls, attribute access, and unbounded exponents in unit expressions."""
    _validate_external_unit_text(expr)
    if re.search(r"[A-Za-z0-9_]\s*\(", expr):
        raise ValueError(f"Function call syntax is prohibited in unit expression {expr!r}.")
    if re.search(r"[A-Za-z_]\s*\.|\.\s*[A-Za-z_]", expr):
        raise ValueError(f"Attribute access syntax is prohibited in unit expression {expr!r}.")
    for match in re.finditer(r"(?:\*\*|\^)\s*([+-]?\d+)?", expr):
        digits = match.group(1)
        if digits is None:
            raise ValueError(f"Unit exponent in {expr!r} must be an integer literal.")
        if abs(int(digits)) > _MAX_UNIT_EXPONENT:
            raise ValueError(f"Unit exponent {digits} in {expr!r} exceeds maximum magnitude {_MAX_UNIT_EXPONENT}.")
    if re.search(r"(?:\*\*|\^)\s*[+-]?\d+\s*\.", expr):
        raise ValueError(f"Non-integer unit exponent is prohibited in {expr!r}.")

@dataclass(frozen=True, slots=True)
class RegistrySnapshot:
    """Immutable point-in-time audit snapshot of a UnitRegistry and its reference data."""
    schema_version: int
    registry_version: str
    currency_data_version: str
    baseline_definition_hash: str
    definition_hash: str
    is_frozen: bool
    unit_count: int
    currency_codes: tuple[str, ...]
    custom_registrations: tuple[str, ...]
    collision_allowlist: Mapping[str, frozenset[str]]
    reference_data_hashes: Mapping[str, str]

def get_registry_definition_hash(registry: pint.UnitRegistry | None = None) -> str:
    """Return the SHA-256 definition hash of the active or supplied UnitRegistry."""
    target = ureg if registry is None else registry
    with _REGISTRY_LOCK:
        stored = getattr(target, "_definition_hash", "")
        if stored:
            return stored
        custom = getattr(target, "_custom_registrations", ())
        return _compute_definition_hash(COMMODITY_UNIT_DEFINITIONS, custom)

def get_registry_snapshot(registry: pint.UnitRegistry | None = None) -> RegistrySnapshot:
    """Capture an immutable audit snapshot of registry version, hashes, and state."""
    from core.units.currencies import CURRENCY_CATALOG
    target = ureg if registry is None else registry
    with _REGISTRY_LOCK:
        is_master = target is ureg
        frozen = _REGISTRY_FROZEN if is_master else bool(getattr(target, "_sealed", False))
        custom = tuple(getattr(target, "_custom_registrations", ()))
        codes = tuple(sorted(_CURRENCY_CODES_SET)) if is_master else CURRENCY_REFERENCE_CODES
        unit_count = len(target._units)
        def_hash = get_registry_definition_hash(target)
    return RegistrySnapshot(schema_version=REGISTRY_SCHEMA_VERSION, registry_version=REGISTRY_VERSION, currency_data_version=CURRENCY_CATALOG.version, baseline_definition_hash=BASELINE_DEFINITION_HASH, definition_hash=def_hash, is_frozen=frozen, unit_count=unit_count, currency_codes=codes, custom_registrations=custom, collision_allowlist=PINT_COLLISION_ALLOWLIST, reference_data_hashes=REFERENCE_DATA_HASHES,)

def assert_registry_compatible(*, expected_schema_version: int | None = None, expected_registry_version: str | None = None, expected_definition_hash: str | None = None, allow_superset_extensions: bool = False, required_units: Mapping[str, str] | Iterable[str] = (), registry: pint.UnitRegistry | None = None,) -> None:
    """Assert that a registry satisfies version, definition-hash, and unit-contract checks. Raises RegistryCompatibilityError on any mismatch."""
    target = ureg if registry is None else registry
    with _REGISTRY_LOCK:
        if expected_schema_version is not None and expected_schema_version != REGISTRY_SCHEMA_VERSION:
            raise RegistryCompatibilityError(f"Registry schema_version mismatch: expected {expected_schema_version}, active is {REGISTRY_SCHEMA_VERSION}.")
        if expected_registry_version is not None and expected_registry_version != REGISTRY_VERSION:
            raise RegistryCompatibilityError(f"Registry version mismatch: expected {expected_registry_version!r}, active is {REGISTRY_VERSION!r}.")
        actual_hash = get_registry_definition_hash(target)
        if expected_definition_hash is not None and actual_hash != expected_definition_hash:
            if not (allow_superset_extensions and expected_definition_hash == BASELINE_DEFINITION_HASH):
                raise RegistryCompatibilityError(f"Registry definition_hash mismatch: expected {expected_definition_hash!r}, active is {actual_hash!r}.")
        if isinstance(required_units, Mapping):
            for unit_expr, expected_dim_expr in required_units.items():
                try:
                    resolved = target.Unit(unit_expr)
                    expected_dim = target.get_dimensionality(expected_dim_expr)
                except (pint.PintError, ValueError, TypeError) as error:
                    raise RegistryCompatibilityError(f"Required unit contract failed for {unit_expr!r}: {error}") from error
                if resolved.dimensionality != expected_dim:
                    raise RegistryCompatibilityError(f"Required unit {unit_expr!r} has dimensionality {resolved.dimensionality}, expected {expected_dim}.")
        else:
            for unit_expr in required_units:
                try:
                    target.Unit(unit_expr)
                except (pint.PintError, ValueError, TypeError) as error:
                    raise RegistryCompatibilityError(f"Required unit {unit_expr!r} is missing or invalid: {error}") from error

def check_registry_compatibility(*, expected_schema_version: int | None = None, expected_registry_version: str | None = None, expected_definition_hash: str | None = None, allow_superset_extensions: bool = False, required_units: Mapping[str, str] | Iterable[str] = (), registry: pint.UnitRegistry | None = None,) -> bool:
    """Return True when the registry satisfies all requested compatibility constraints."""
    try:
        assert_registry_compatible(expected_schema_version=expected_schema_version, expected_registry_version=expected_registry_version, expected_definition_hash=expected_definition_hash, allow_superset_extensions=allow_superset_extensions, required_units=required_units, registry=registry,)
        return True
    except RegistryCompatibilityError:
        return False

# ---------------------------------------------------------------------------
# Public Sanitisation & Ingestion APIs
# ---------------------------------------------------------------------------

@contextmanager
def sandbox_registry() -> Iterator[pint.UnitRegistry]:
    """Yield an independent registry containing the baseline master definitions.

    Runtime desk additions and master freeze state are not copied. Define test
    units directly on the yielded registry; module APIs still use the master.
    The application registry is never rebound, even when the context raises.
    """
    sandbox = _new_registry()
    yield sandbox

@contextmanager
def _registry_transaction():
    global _REGISTRY_FROZEN, _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN
    with _REGISTRY_LOCK:
        _check_frozen()
        ctx_var = ureg.__dict__["_active_ctx_var"]
        memo = {
            id(ureg): ureg,
            id(ureg.Quantity): ureg.Quantity,
            id(ureg.Unit): ureg.Unit,
            id(ctx_var): ctx_var,
        }
        state = copy.deepcopy(ureg.__dict__, memo)
        names = ("CURRENCY_MINOR_DIGITS", "CURRENCY_SYMBOL_TO_CODE", "_SAFE_SYMBOLS", "_CURRENCY_CODES_SET", "_CURRENCY_CANONICAL_CODES", "minor_mapping", "irregular_denominators", "_QUOTE_REPLACEMENTS")
        snapshots = {name: dict(globals()[name]) if isinstance(globals()[name], dict) else globals()[name].copy() for name in names}
        tokens = _UNIT_TOKEN_REPLACEMENTS[:]
        frozen, symbols, pattern = _REGISTRY_FROZEN, _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN
        try:
            yield
        except BaseException:
            _REGISTRY_FROZEN = frozen
            object.__setattr__(ureg, "_sealed", False)
            ureg.__dict__.clear()
            ureg.__dict__.update(state)
            for name, snapshot in snapshots.items():
                target = globals()[name]
                if hasattr(target, "_restore"):
                    target._restore(snapshot)
                else:
                    target.clear()
                    target.update(snapshot)
            _UNIT_TOKEN_REPLACEMENTS[:] = tokens
            _SORTED_SAFE_SYMBOLS, _PREFIX_CURRENCY_PATTERN = symbols, pattern
            raise

def freeze_registry() -> None:
    """Seal public registration after service startup configuration.

    Register approved currencies and commodity units during process
    initialisation, then call ``freeze_registry()`` before serving requests.
    Freezing is permanent for this process and blocks the public registration
    APIs, registry mutators, and module-level metadata mappings. Reading units
    and activating an existing context via ``with ureg.context(...):`` remain available.
    """
    global _REGISTRY_FROZEN
    with _REGISTRY_LOCK:
        _REGISTRY_FROZEN = True
        object.__setattr__(ureg, "_sealed", True)

def is_registry_frozen() -> bool:
    """Return whether public runtime registration has been disabled."""
    with _REGISTRY_LOCK:
        return _REGISTRY_FROZEN

def _check_frozen() -> None:
    """Check the registration guard while the caller holds _REGISTRY_LOCK."""
    from core.units.exceptions import RegistryFrozenError
    if _REGISTRY_FROZEN:
        raise RegistryFrozenError("UnitRegistry is frozen. Dynamic registration of custom units or currencies is prohibited in locked production environments.")

def register_currency(code: str, *, decimals: int, symbols: tuple[str, ...] = (), aliases: tuple[str, ...] = (), sub_units: dict[str, str] | None = None, quote_aliases: dict[str, str] | None = None,) -> None:
    """
    Registers a fiat, crypto, or central bank digital currency (CBDC) into the registry.

        --- Parameters ---

        code : str
            Canonical currency symbol (e.g. "BTC", "ETH", "USDT").

        decimals : int
            Settlement precision exponent (e.g. 8 for BTC, 18 for ETH, 6 for USDC).
            Used downstream by Decimal Banker's Rounding engines.

        symbols : tuple[str, ...]
            Graphical or compound screen symbols (e.g. a cryptocurrency glyph).

        aliases : tuple[str, ...]
            Alternative screen symbols or ticker tags (e.g. ("XBT", "bitcoin")).

        sub_units : dict[str, str] | None
            Named fractional denominations mapped to their scaling definition relative to the canonical unit.
            e.g. {"satoshi": "1e-8 * BTC = sat = sats"}

        quote_aliases : dict[str, str] | None
            Broker/exchange screen quote normalisations (e.g. {"BTC/USD": "BTC / USD"}).
    """

    cleaned_currency_code_attempt = code.strip() # Remove leading/trailing whitespace from the currency code
    if not cleaned_currency_code_attempt or not cleaned_currency_code_attempt.isalnum(): # If the cleaned code is empty or contains non-alphanumeric characters, raise an error
        raise ValueError(f"Currency code '{code}' must be a non-empty alphanumeric string.")
    if isinstance(decimals, bool) or not isinstance(decimals, int) or decimals < 0:
        raise ValueError(f"Decimals for currency '{code}' must be a non-negative integer, got {decimals}.")

    # Construct dimensional definition isolating this currency
    if aliases: # If there are aliases provided, join them into a string for the definition
        joined_aliases = " = ".join(aliases) # Join aliases with " = " to create a definition string
        alias_attempt = f" = {joined_aliases}" # Create the alias attempt string
    else:
        alias_attempt = "" # If no aliases are provided, set it to an empty string

    base_definition = f"{cleaned_currency_code_attempt} = [currency_{cleaned_currency_code_attempt}]{alias_attempt}" # Construct the base definition string for the currency, including its dimension and any aliases

    with _registry_transaction(): # Roll back partial definitions and metadata on failure
        _check_frozen()
        for symbol in symbols:
            if not isinstance(symbol, str) or not symbol or symbol != symbol.strip() or any(char.isspace() for char in symbol):
                raise ValueError("Currency symbols must be non-empty strings without whitespace.")
            if not _has_special_symbol_chars(symbol) or any(char in symbol for char in "/=()*0123456789"):
                raise ValueError(f"Invalid graphical currency symbol: {symbol!r}")
            if symbol in CURRENCY_SYMBOL_TO_CODE:
                raise ValueError(f"Currency symbol {symbol!r} is already registered.")
        if cleaned_currency_code_attempt in ureg: # If the currency code is already registered in the unit registry, raise an error
            raise ValueError(f"Currency '{cleaned_currency_code_attempt}' is already registered in the unit registry.")

        # 1. Register canonical base currency dimension
        object.__setattr__(ureg, "_in_currency_registration", True)
        try:
            ureg.define(base_definition)
        finally:
            object.__setattr__(ureg, "_in_currency_registration", False)

        # 2. Register minor / sub-denominations (e.g., satoshis, gwei, wei)
        if sub_units:
            for sub_name, sub_expressions in sub_units.items():
                if not isinstance(sub_name, str) or not isinstance(sub_expressions, str):
                    raise ValueError("Currency sub_units keys and expressions must be strings.")
                sub_definition = f"{sub_name} = {sub_expressions}"
                validate_unit_definition(sub_definition, registry=ureg, allow_new_dimension=False, allow_alias_directive=False,)
                ureg.define(sub_definition)
                _register_minor_unit_definition(sub_definition)

        # 3. Synchronise downstream settlement engine metadata
        _CURRENCY_CODES_SET.add(cleaned_currency_code_attempt.upper())
        _CURRENCY_CODES_SET.update(alias.upper() for alias in aliases)
        _CURRENCY_CANONICAL_CODES[cleaned_currency_code_attempt.upper()] = cleaned_currency_code_attempt
        for alias in aliases:
            _CURRENCY_CANONICAL_CODES[alias.upper()] = cleaned_currency_code_attempt
        CURRENCY_MINOR_DIGITS[cleaned_currency_code_attempt] = decimals
        for alias in aliases:
            CURRENCY_MINOR_DIGITS[alias] = decimals

        # Publish symbols before the new regex so every regex match is resolvable.
        if symbols:
            for symbol in symbols:
                CURRENCY_SYMBOL_TO_CODE[symbol] = cleaned_currency_code_attempt
                _SAFE_SYMBOLS[symbol] = cleaned_currency_code_attempt
            _recompile_prefix_currency_pattern()

        # 4. Register vendor market data quote substitutions
        if quote_aliases:
            for raw_quote, normalised_expressions in quote_aliases.items():
                _validate_external_unit_text(raw_quote)
                _validate_sanitized_unit_expression(normalised_expressions)
                ureg.Unit(normalised_expressions)
                _QUOTE_REPLACEMENTS[raw_quote] = normalised_expressions


def register_commodity_unit(definition: str, *, quote_aliases: dict[str, str] | None = None,) -> None:
    """
    Registers a custom or desk-specific commodity unit into the active registry.

    --- Parameters ---

    definition : str
        A valid Pint unit definition string (e.g. "drum_200l = 200 * liter = steel_drum").

    quote_aliases : dict[str, str] | None
        Optional raw broker screen quote substitutions to map to this unit
        (e.g. {"$/drum": "USD / drum_200l"}).
    """

    if not isinstance(definition, str) or not definition.strip(): # If the definition is not a string or is empty after stripping whitespace
        raise ValueError("Unit definition must be a non-empty string.") # Raise an error

    clean_definition = definition.strip()

    # Check if the unit is already defined and strictly conforms to the safe definition grammar
    with _registry_transaction(): # Roll back partial definitions and metadata on failure
        _check_frozen()
        canonical_name, _, validated_aliases = validate_unit_definition(clean_definition, registry=ureg, allow_new_dimension=False, allow_alias_directive=False,)

        ureg.define(clean_definition) # Register into active Pint registry
        _register_commodity_denominator_definition(clean_definition)

        # Resolve custom symbols and aliases before Pint can interpret them as
        # prefixed units. Keep matching case-sensitive for desk-defined names.
        for token in validated_aliases:
            _UNIT_TOKEN_REPLACEMENTS.insert(0, (re.compile(r"(?<![A-Za-z0-9_])" + re.escape(token) + r"(?![A-Za-z0-9_])"), canonical_name,))

        if quote_aliases: # Register vendor screen quote replacements if provided
            for screen_quote, normalised_expressions in quote_aliases.items():
                _validate_external_unit_text(screen_quote)
                _validate_sanitized_unit_expression(normalised_expressions)
                ureg.Unit(normalised_expressions)
                _QUOTE_REPLACEMENTS[screen_quote] = normalised_expressions

@_locked_registry
def sanitize_unit_string(raw_unit: str, default_currency: str | None = None, *, context: QuoteParsingContext | None = None) -> str:
    """Normalize unit text; context selects source-aware commodity quote parsing."""
    from core.units.exceptions import AmbiguousUnitError
    _validate_external_unit_text(raw_unit)
    if context is not None:
        from core.units.quote_parser import _quote_parsing_normalize_quote_expression as _normalize_quote_expression
        result = _normalize_quote_expression(raw_unit, context, default_currency)
        _validate_sanitized_unit_expression(result)
        return result
    text = re.sub(r"\s+per\s+", " / ", raw_unit.strip(), flags=re.I)
    ambiguous = re.search(r"(?<![A-Za-z0-9_])(ppt|mu|ton_mile)(?![A-Za-z0-9_])", text)
    if ambiguous:
        raise AmbiguousUnitError(f"Ambiguous market unit {ambiguous[0]!r}; use an explicit ratio, area or freight convention.")

    # 1. Fast-path lookup against cached quote mappings
    if default_currency is None and text in _QUOTE_REPLACEMENTS:
        return _QUOTE_REPLACEMENTS[text]

    # 2. Extract leading currency glyphs and convert to ISO codes (e.g. "$/barrel" -> "USD/barrel")
    for symbol in _SORTED_SAFE_SYMBOLS:
        # Ensure symbol is standalone or followed immediately by a slash or space
        if text.startswith(symbol) and (len(text) == len(symbol) or text[len(symbol)] in {"/", " "}):
            resolved_code = _resolve_currency_symbol(symbol, default_currency)
            text = resolved_code + text[len(symbol):]
            break

    if text in _QUOTE_REPLACEMENTS:
        return _QUOTE_REPLACEMENTS[text]

    # 3. Normalise high-frequency exchange minor-unit prefixes
    if text.startswith("p/"):
        text = "GBP_penny/" + text[2:]
    elif text.startswith("c/"):
        text = "USD_cent/" + text[2:]
    elif text.startswith("cents/"):
        text = "USD_cent/" + text[6:]

    # 4. Apply strict-boundary regex replacements for problematic tokens (e.g. bbl, MMBtu)
    for pattern, replacement in _UNIT_TOKEN_REPLACEMENTS:
        text = pattern.sub(replacement, text)

    _validate_sanitized_unit_expression(text)
    return text

@_locked_registry
def parse_unit(unit: str | pint.Unit, *, default_currency: str | None = None, context: QuoteParsingContext | None = None) -> pint.Unit:
    """Parse and normalise a raw unit expression into a validated Pint Unit."""
    if isinstance(unit, pint.Unit):
        validate_registry_quantity(unit)
        if context is not None:
            raise TypeError("Source-aware quote parsing requires original text, not a pre-parsed Unit.")
        return unit
    if not isinstance(unit, str) or not unit.strip():
        raise ValueError(f"Unit expression must be a non-empty string, got {unit!r}")

    sanitized = sanitize_unit_string(unit, default_currency, context=context)
    return ureg.Unit(sanitized)

def _parse_magnitude(value, as_decimal):
    if isinstance(value, str):
        match = _BARE_MAGNITUDE_PATTERN.fullmatch(value)
        if match is None:
            raise ValueError(f"Invalid numeric magnitude: {value!r}")
        value = match.group("magnitude").replace(",", "")
        if match.group("bracket"):
            value = "-" + value.lstrip("+-")
    if isinstance(value, bool):
        raise TypeError("Boolean values are not quantity magnitudes.")
    magnitude = Decimal(str(value)) if as_decimal else float(value)
    if not Decimal(str(magnitude)).is_finite():
        raise ValueError("Quantity magnitude must be finite and representable.")
    return magnitude

@_locked_registry
def parse_quantity(value: int | float | Decimal | str, unit: str | None = None, *, default_currency: str | None = None, as_decimal: bool = True, context: QuoteParsingContext | None = None,) -> pint.Quantity:
    """Parses numeric expressions or quotation strings into validated Pint Quantities.

    Handles currency-prefixed strings ("$75.50/bbl"), spaced numbers ("75.50 USD/bbl"),
    accounting negatives ("(1,500) MWh"), percentages ("2.5%"), and scientific
    notation ("1.5e6 scfd"). Preserves high-precision arithmetic by
    defaulting to `Decimal` magnitude storage.

    With context, require a commodity quote and resolve its original tokens
    using explicit source conventions before any legacy ticker normalization.
    """
    if isinstance(value, (pint.Quantity, pint.Unit)):
        validate_registry_quantity(value)
    if isinstance(unit, (pint.Quantity, pint.Unit)):
        validate_registry_quantity(unit)
    if isinstance(value, str) and unit is None:
        _validate_external_unit_text(value)
        # Pattern A: Leading currency symbol expressions (e.g. "$75.50/bbl", "-£45.20 / MWh")
        prefix_match = _PREFIX_CURRENCY_PATTERN.match(value)
        if prefix_match:
            sign = "-" if prefix_match.group("bracket") else (prefix_match.group("sign") or "")
            symbol = prefix_match.group("symbol")
            mag_str = sign + prefix_match.group("magnitude").replace(",", "")
            denomination = prefix_match.group("denomination")

            currency_code = _resolve_currency_symbol(symbol, default_currency) if context is None else symbol
            # Keep the screen ticker intact so explicit mappings such as EUR
            # carbon quotes retain their commodity identity on round trips.
            unit_str = f"{symbol}/{denomination}" if denomination else currency_code
            sanitized = sanitize_unit_string(unit_str, default_currency, context=context)
            magnitude = _parse_magnitude(mag_str, as_decimal)
            return ureg.Quantity(magnitude, sanitized)

        # Pattern B: Standard numeric value with trailing unit string (e.g. "75.50 USD/bbl")
        num_match = _NUMERIC_MAGNITUDE_PATTERN.match(value)
        if num_match:
            mag_str = num_match.group("magnitude").replace(",", "")
            if num_match.group("bracket"):
                mag_str = "-" + mag_str.lstrip("+-")
            unit_str = num_match.group("unit")
            sanitized = sanitize_unit_string(unit_str, default_currency, context=context)
            magnitude = _parse_magnitude(mag_str, as_decimal)
            return ureg.Quantity(magnitude, sanitized)

        if _BARE_MAGNITUDE_PATTERN.fullmatch(value):
            if context is not None:
                raise ValueError("A commodity quote requires currency/commodity units.")
            return ureg.Quantity(_parse_magnitude(value, as_decimal), ureg.dimensionless)

        # Pattern C: Standalone dimensional unit strings with implicit unit magnitude (e.g. "USD/bbl")
        parsed = ureg.Quantity(sanitize_unit_string(value, default_currency, context=context))
        magnitude = _parse_magnitude(parsed.magnitude, as_decimal)
        return ureg.Quantity(magnitude, parsed.units)

    # Pattern D: Numerical value supplied without an explicit unit (dimensionless quantity)
    if unit is None:
        if context is not None:
            raise ValueError("A commodity quote requires currency/commodity units.")
        return ureg.Quantity(_parse_magnitude(value, as_decimal), ureg.dimensionless)

    # Pattern E: Explicit separate value and unit parameters passed
    sanitized = sanitize_unit_string(unit, default_currency, context=context)
    return ureg.Quantity(_parse_magnitude(value, as_decimal), sanitized)

# ---------------------------------------------------------------------------
# Institutional Validation & Unit Enforcement Helpers
# ---------------------------------------------------------------------------

def validate_dimensionality(quantity: pint.Quantity, expected_dimension: str | pint.Unit | pint.Quantity,) -> bool:
    """Verifies that a quantity strictly conforms to an expected dimensional context.

    Prevents execution errors such as pricing volumetric barrels against mass-denominated
    hedging curves.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, received {type(quantity).__name__}")
    validate_registry_quantity(quantity)

    if isinstance(expected_dimension, (pint.Unit, pint.Quantity)):
        validate_registry_quantity(expected_dimension)
        return quantity.check(expected_dimension.dimensionality)

    # Resolve raw string through sanitiser prior to registry dimensionality verification
    parsed_expected = ureg(sanitize_unit_string(expected_dimension))
    return quantity.check(parsed_expected.dimensionality)

def ensure_unit(quantity: pint.Quantity, target_unit: str | pint.Unit,) -> pint.Quantity:
    """Converts a quantity to a specified unit, raising Pint's DimensionalityError on mismatch."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, received {type(quantity).__name__}")
    validate_registry_quantity(quantity)
    if isinstance(target_unit, (pint.Unit, pint.Quantity)):
        validate_registry_quantity(target_unit)
    sanitized_target = (sanitize_unit_string(target_unit) if isinstance(target_unit, str) else target_unit)
    return quantity.to(sanitized_target)

def are_convertible(source_unit: str | pint.Unit | pint.Quantity, target_unit: str | pint.Unit | pint.Quantity, *, default_currency: str | None = None,) -> bool:
    """Return whether ordinary conversion is dimensionally and semantically safe."""
    validate_registry_quantity(source_unit)
    validate_registry_quantity(target_unit)
    def resolve_unit(value: str | pint.Unit | pint.Quantity) -> pint.Unit:
        if isinstance(value, pint.Quantity):
            return value.units
        return parse_unit(value, default_currency=default_currency)
    try:
        source, target = resolve_unit(source_unit), resolve_unit(target_unit)
        if source.dimensionality != target.dimensionality:
            return False
        return ureg._semantic_tags(source._units) == ureg._semantic_tags(target._units)
    except ForeignRegistryError:
        raise
    except (pint.PintError, TypeError, ValueError):
        return False

# ---------------------------------------------------------------------------
# Reverse Market Formatting
# ---------------------------------------------------------------------------

_CANONICAL_TO_MARKET_TICKER: Final[Mapping[str, str]] = MappingProxyType({
    "USD / oil_barrel": "$/bbl",
    "USD / gallon": "$/gal",
    "USD_cent / gallon": "c/gal",
    "USD / MMBtu_it": "$/MMBtu",
    "USD / megawatt_hour": "$/MWh",
    "GBP_penny / therm_it": "p/th",
    "GBP / megawatt_hour": "\u00a3/MWh",
    "EUR / megawatt_hour": "\u20ac/MWh",
    "USD / tonne": "$/mt",
    "USD_cent / bushel": "c/bu",
    "USD_cent / pound": "c/lb",
    "USD / pound": "$/lb",
    "USD / toz": "$/toz",
    "EUR / tonne_CO2": "\u20ac/t",
    "EUR / EU_allowance": "\u20ac/EUA",
    "GBP / tonne_CO2": "\u00a3/t",
    "GBP / UK_allowance": "\u00a3/UKA",
    "USD / California_carbon_allowance": "$/CCA",
    "USD / gigawatt_hour": "$/GWh",
    "EUR / gigawatt_hour": "\u20ac/GWh",
    "GBP / gigawatt_hour": "\u00a3/GWh",
    "USD / kilowatt_hour": "$/kWh",
    "USD_cent / kilowatt_hour": "c/kWh",
    "GBP_penny / kilowatt_hour": "p/kWh",
    "USD / Mscf_60F_14_73psia": "$/Mcf",
    "USD / dth_it": "$/Dth",
    "USD / boe_58": "$/boe",
    "USD / commodity_short_hundredweight": "$/cwt",
    "USD_cent / commodity_short_hundredweight": "c/cwt",
})

# Resolve aliases through Pint once, so its preferred spelling cannot prevent a match.
_MARKET_TICKERS: Final[Mapping[pint.Unit, str]] = MappingProxyType({ureg.Unit(canonical): ticker for canonical, ticker in _CANONICAL_TO_MARKET_TICKER.items()})
_CANONICAL_UNIT_DISPLAY: Final[Mapping[str, str]] = MappingProxyType({
    "gCO2e / megajoule": "CI",
    "CI_unit": "CI",
    # Render complete ratios before replacing individual energy identifiers.
    "btu_it / kilowatt_hour": "Btu/kWh",
    "MMBtu_it / megawatt_hour": "MMBtu/MWh",
    "heat_rate_us": "Btu/kWh",
    "heat_rate_macro": "MMBtu/MWh",
    "heat_rate_metric": "MJ/kWh",
    "heat_rate_gas_turbine": "GJ/MWh",
    "oil_barrel": "bbl",
    "MMBtu_it": "MMBtu",
    "meter ** 3": "cbm",
    "therm_it": "th",
    "metric_tonne": "mt",
    "tonne": "mt",
    "megawatt_hour": "MWh",
    "kilowatt_hour": "kWh",
    "gigawatt_hour": "GWh",
    "dth_it": "Dth",
    "commodity_short_hundredweight": "cwt",
    "Mscf_60F_14_73psia": "Mcf",
    "MMscf_60F_14_73psia": "MMcf",
    "Bscf_60F_14_73psia": "Bcf",
    "Tscf_60F_14_73psia": "Tcf",
    "scf_60F_14_73psia": "scf",
    "boe_58": "boe",
    "basis_point": "bps",
    "percent": "%",
})

def _format_quantity_latex(quantity: pint.Quantity, magnitude: str) -> str:
    """Typeset clean unit names using Pint's structural factors and exponents."""
    escapes = {"\\": r"\backslash{}", "{": r"\{", "}": r"\}", "_": r"\_", "%": r"\%", "$": r"\$", "&": r"\&", "#": r"\#", "^": r"\hat{}", "~": r"\sim{}",}

    def display_latex(display: str) -> str:
        if "/" in display:
            numerator, denominator = display.split("/", 1)
            return rf"\frac{{{display_latex(numerator.strip())}}}{{{display_latex(denominator.strip())}}}"
        escaped = "".join(escapes.get(char, char) for char in display)
        return rf"\mathrm{{{escaped}}}"

    # A complete display override includes composites such as carbon intensity.
    whole_display = _CANONICAL_UNIT_DISPLAY.get(str(quantity.units))
    if whole_display is not None:
        unit_latex = display_latex(whole_display)
    else:
        numerator_terms: list[str] = []
        denominator_terms: list[str] = []
        for name, exponent in quantity.to_tuple()[1]:
            if exponent == 0:
                continue
            display = _CANONICAL_UNIT_DISPLAY.get(name, name)
            term = display_latex(display)
            power = abs(exponent)
            if power != 1:
                if "/" in display:
                    term = rf"\left({term}\right)"
                term = rf"{term}^{{{power}}}"
            (numerator_terms if exponent > 0 else denominator_terms).append(term)
        numerator = r"\,".join(numerator_terms) or "1"
        denominator = r"\,".join(denominator_terms)
        unit_latex = rf"\frac{{{numerator}}}{{{denominator}}}" if denominator else numerator
        if not numerator_terms and not denominator_terms:
            unit_latex = ""

    scientific = re.fullmatch(r"([+-]?(?:\d+(?:\.\d*)?|\.\d+))[eE]([+-]?\d+)", magnitude)
    if scientific:
        magnitude = rf"{scientific[1]}\times 10^{{{int(scientific[2])}}}"
    return rf"{magnitude}\,{unit_latex}" if unit_latex else magnitude

@_locked_registry
def format_quantity(quantity: pint.Quantity, *, style: str = "market", precision: int | None = None,) -> str:
    """Format a quantity as market, standard, canonical, or LaTeX notation.

    Precision is a non-negative number of decimal places; None preserves the
    magnitude's existing representation without currency-based rounding.
    Unknown market tickers fall back to standard unit display without conversion.
    LaTeX returns a math fragment without dollar delimiters.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)
    if style not in {"market", "standard", "canonical", "latex"}:
        raise ValueError(f"Unknown quantity format style: {style!r}")
    if precision is not None and (isinstance(precision, bool) or not isinstance(precision, int) or precision < 0):
        raise ValueError("precision must be a non-negative integer or None.")

    magnitude = str(quantity.magnitude) if precision is None else f"{quantity.magnitude:.{precision}f}"
    unit_string = str(quantity.units)
    if style == "canonical":
        return f"{magnitude} {unit_string}"
    if style == "latex":
        return _format_quantity_latex(quantity, magnitude)
    if style == "market":
        ticker = _MARKET_TICKERS.get(quantity.units)
        if ticker is not None:
            if ticker.startswith(("$", "\u00a3", "\u20ac", "\u00a5")):
                return f"{ticker[0]}{magnitude}{ticker[1:]}"
            return f"{magnitude} {ticker}"

    for canonical, display in _CANONICAL_UNIT_DISPLAY.items():
        unit_string = re.sub(r"(?<![A-Za-z0-9_])" + re.escape(canonical) + r"(?![A-Za-z0-9_])", display, unit_string,)
    return f"{magnitude} {unit_string}"

# ---------------------------------------------------------------------------
# Semantic Classification & Batch Ingestion
# ---------------------------------------------------------------------------

@cache
def _classify_dimensions(dimensions) -> tuple[bool, bool, bool, bool]:
    """Cache classification by immutable Pint dimensionality, including desk units."""
    currencies = {dimension: exponent for dimension, exponent in dimensions.items() if dimension.startswith("[currency_")}
    isolated_currency = len(dimensions) == 1 and list(currencies.values()) == [1]
    physical = {dimension: exponent for dimension, exponent in dimensions.items() if dimension not in currencies}
    commodity_price = (len(currencies) == 1 and list(currencies.values()) == [1] and any(exponent < 0 for exponent in physical.values()))
    flow_rate = not currencies and dimensions.get("[time]") == -1 and any(dimension != "[time]" and exponent > 0 for dimension, exponent in physical.items())
    fx_rate = (not physical and len(currencies) == 2 and sorted(currencies.values()) == [-1, 1])
    return isolated_currency, commodity_price, flow_rate, fx_rate

def is_currency(quantity: pint.Quantity) -> bool:
    """Return True for an isolated currency or its fractional denominations."""
    if not isinstance(quantity, pint.Quantity):
        return False
    validate_registry_quantity(quantity)
    return _classify_dimensions(quantity.dimensionality)[0]

def is_commodity_price(quantity: pint.Quantity) -> bool:
    """Return True for one currency over physical dimensions, excluding FX rates."""
    if not isinstance(quantity, pint.Quantity):
        return False
    validate_registry_quantity(quantity)
    return _classify_dimensions(quantity.dimensionality)[1]

def is_flow_rate(quantity: pint.Quantity) -> bool:
    """Return True for non-financial quantities per time, such as bpd or tph."""
    if not isinstance(quantity, pint.Quantity):
        return False
    validate_registry_quantity(quantity)
    return _classify_dimensions(quantity.dimensionality)[2]

def is_fx_rate(quantity_or_unit: pint.Quantity | pint.Unit) -> bool:
    """Return True for pure currency ratios such as EUR/USD or GBP/USD."""
    if not isinstance(quantity_or_unit, (pint.Quantity, pint.Unit)):
        return False
    validate_registry_quantity(quantity_or_unit)
    return _classify_dimensions(quantity_or_unit.dimensionality)[3]

@_locked_registry
def parse_quantities(values: list[str | int | float | Decimal], units: str | list[str] | None = None, *, default_currency: str | None = None, as_decimal: bool = True, context: QuoteParsingContext | None = None,) -> list[pint.Quantity]:
    """Batch parse expressions or magnitudes with shared or paired units.

    Repeated units are sanitised once per call. Paired lists must have equal
    lengths. Separate magnitudes support grouped numbers and accounting signs.
    """
    if units is None:
        return [parse_quantity(value, default_currency=default_currency, as_decimal=as_decimal, context=context) for value in values]
    if not isinstance(units, (str, list)):
        raise TypeError("units must be a string, a list of strings, or None.")
    if isinstance(units, list) and len(values) != len(units):
        raise ValueError("values and units must have equal lengths.")

    sanitised_units: dict[str, pint.Unit] = {}
    def resolve_unit(raw_unit: str) -> pint.Unit:
        if raw_unit not in sanitised_units:
            sanitised_units[raw_unit] = ureg.Unit(sanitize_unit_string(raw_unit, default_currency, context=context))
        return sanitised_units[raw_unit]

    shared_unit = resolve_unit(units) if isinstance(units, str) else None
    results: list[pint.Quantity] = []
    for index, value in enumerate(values):
        magnitude = _parse_magnitude(value, as_decimal)
        parsed_unit = shared_unit if shared_unit is not None else resolve_unit(units[index])
        results.append(ureg.Quantity(magnitude, parsed_unit))
    return results

# ---------------------------------------------------------------------------
# Quote Structure, Scalar Serialisation & Startup Diagnostics
# ---------------------------------------------------------------------------

@dataclass(frozen=True, slots=True)
class QuoteComponents:
    """Immutable financial and physical metadata for a two-sided market quote."""
    currency: str
    commodity_unit: str
    is_minor_currency: bool
    settlement_decimals: int | None
    canonical_expression: str

@_locked_registry
def decompose_quote(quote: str, default_currency: str | None = None, *, context: QuoteParsingContext | None = None) -> QuoteComponents:
    """Validate and decompose a currency/commodity quote without price conversion.

    Currency dimensions identify registered proprietary sub-units as well as
    sovereign minor units. Denominators must be physical, not currencies or FX.
    """
    if not isinstance(quote, str):
        raise TypeError("quote must be a string.")
    expression = sanitize_unit_string(quote, default_currency, context=context)
    if "/" not in expression:
        raise ValueError(f"Expression {quote!r} is not a two-sided quote.")
    numerator_string, denomination_string = (part.strip() for part in expression.split("/", 1))
    try:
        numerator = ureg.Unit(numerator_string)
        # Derive the denominator from the complete expression so chained
        # divisions preserve Pint's associativity rather than reversing a rate.
        complete_unit = ureg.Unit(expression)
        denomination = numerator / complete_unit
        dimensions = numerator.dimensionality
        if len(dimensions) != 1:
            raise ValueError("Quote numerator must be an isolated currency.")
        dimension, exponent = next(iter(dimensions.items()))
        if exponent != 1 or not dimension.startswith("[currency_"):
            raise ValueError("Quote numerator must be an isolated currency.")
        currency = dimension[len("[currency_"):-1]
        if currency.upper() not in _CURRENCY_CANONICAL_CODES:
            raise ValueError("Quote numerator currency is not registered.")
        currency = _CURRENCY_CANONICAL_CODES[currency.upper()]
        if not denomination.dimensionality or any(name.startswith("[currency_") for name in denomination.dimensionality):
            raise ValueError("Quote denominator must represent a physical commodity unit.")
        minor = ureg.Quantity(1, numerator).to(currency).magnitude < 1
    except (pint.PintError, ValueError, TypeError) as error:
        raise ValueError(f"Invalid commodity quote {quote!r}: {error}") from error
    return QuoteComponents(currency=currency, commodity_unit=str(denomination), is_minor_currency=minor, settlement_decimals=CURRENCY_MINOR_DIGITS.get(currency), canonical_expression=f"{numerator} / ({denomination})" if "/" in str(denomination) else f"{numerator} / {denomination}",)

@_locked_registry
def quantity_to_dict(quantity: pint.Quantity, *, include_metadata: bool = False,) -> dict[str, Any]:
    """Serialise a finite scalar quantity to JSON-safe magnitude/unit strings.

    When ``include_metadata=True``, embeds ``schema_version``, ``registry_version``,
    and ``definition_hash`` for cross-process compatibility verification.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)
    magnitude = quantity.magnitude
    if isinstance(magnitude, bool) or not isinstance(magnitude, (int, float, Decimal)):
        raise TypeError("Only scalar int, float, and Decimal magnitudes are supported.")
    if not Decimal(str(magnitude)).is_finite():
        raise ValueError("Quantity magnitude must be finite.")
    payload: dict[str, Any] = {"m": str(magnitude), "u": str(quantity.units)}
    if include_metadata:
        payload["schema_version"] = REGISTRY_SCHEMA_VERSION
        payload["registry_version"] = REGISTRY_VERSION
        payload["definition_hash"] = get_registry_definition_hash()
    return payload

@_locked_registry
def quantity_from_dict(payload: dict[str, Any], *, as_decimal: bool = True, expected_schema_version: int | None = None, expected_registry_version: str | None = None, expected_definition_hash: str | None = None, allow_superset_extensions: bool = False,) -> pint.Quantity:
    """Restore exact stored unit identifiers, verifying registry compatibility metadata when present."""
    if not isinstance(payload, dict) or not {"m", "u"} <= payload.keys():
        raise ValueError("Invalid quantity payload; must contain 'm' and 'u'.")

    schema_ver = payload.get("schema_version")
    if schema_ver is not None:
        if isinstance(schema_ver, bool) or not isinstance(schema_ver, int):
            raise RegistryCompatibilityError("Payload 'schema_version' must be an integer.")
        if expected_schema_version is not None and schema_ver != expected_schema_version:
            raise RegistryCompatibilityError(f"Payload schema_version {schema_ver} conflicts with expected_schema_version {expected_schema_version}.")
        expected_schema_version = schema_ver

    reg_ver = payload.get("registry_version")
    if reg_ver is not None:
        if not isinstance(reg_ver, str):
            raise RegistryCompatibilityError("Payload 'registry_version' must be a string.")
        if expected_registry_version is not None and reg_ver != expected_registry_version:
            raise RegistryCompatibilityError(f"Payload registry_version {reg_ver!r} conflicts with expected_registry_version {expected_registry_version!r}.")
        expected_registry_version = reg_ver

    def_hash = payload.get("definition_hash")
    if def_hash is not None:
        if not isinstance(def_hash, str):
            raise RegistryCompatibilityError("Payload 'definition_hash' must be a string.")
        if expected_definition_hash is not None and def_hash != expected_definition_hash:
            raise RegistryCompatibilityError(f"Payload definition_hash {def_hash!r} conflicts with expected_definition_hash {expected_definition_hash!r}.")
        expected_definition_hash = def_hash

    assert_registry_compatible(expected_schema_version=expected_schema_version, expected_registry_version=expected_registry_version, expected_definition_hash=expected_definition_hash, allow_superset_extensions=allow_superset_extensions,)

    raw_magnitude, raw_unit = payload["m"], payload["u"]
    if isinstance(raw_magnitude, bool) or not isinstance(raw_magnitude, (str, int, float, Decimal)):
        raise ValueError("Payload magnitude must be a numeric scalar or string.")
    if not isinstance(raw_unit, str) or not raw_unit.strip():
        raise ValueError("Payload unit must be a non-empty string.")
    _validate_sanitized_unit_expression(raw_unit)
    try:
        decimal_magnitude = Decimal(str(raw_magnitude))
        if not decimal_magnitude.is_finite():
            raise ValueError("Payload magnitude must be finite.")
        magnitude = decimal_magnitude if as_decimal else float(decimal_magnitude)
        if not Decimal(str(magnitude)).is_finite():
            raise ValueError("Payload magnitude overflows float precision.")
        return ureg.Quantity(magnitude, ureg.Unit(raw_unit))
    except (InvalidOperation, pint.PintError, TypeError, OverflowError) as error:
        raise ValueError(f"Invalid quantity payload: {error}") from error


from core.units.validation import authorized_conversion

@authorized_conversion
def verify_registry_integrity() -> None:
    """Audit static definitions, dimensional isolation, shields and cached quotes.

    Raises AssertionError with the failing invariant, including under python -O.
    This diagnostic validates existing baselines; it registers no conversions.
    """
    def check_scale(name: str, quantity: pint.Quantity, target: str, expected: float) -> None:
        actual = float(quantity.to(target).magnitude)
        if not abs(actual - expected) <= max(1e-9, abs(expected) * 1e-12):
            raise AssertionError(f"{name}: expected {expected} {target}, got {actual}.")

    check_scale("Petroleum barrel", ureg.Quantity(1, "bbl"), "gallon", 42)
    check_scale("IT Btu", ureg.Quantity(1, "btu_it"), "joule", 1055.05585262)
    try:
        ureg.Quantity(1, "USD").to("EUR")
    except pint.DimensionalityError:
        pass
    else:
        raise AssertionError("USD converted to EUR without an explicit FX context.")
    for token, expected in (("Mcf", 1e3), ("MMcf", 1e6), ("Bcf", 1e9), ("Tcf", 1e12)):
        check_scale(token, ureg.Quantity(1, sanitize_unit_string(token)), "scf_60F_14_73psia", expected)
    check_scale("Basis point", ureg.Quantity(100, "basis_point"), "percent", 1)
    weather_index = ureg.Quantity(10, "delta_degF") * ureg.Quantity(1, "day")
    check_scale("Weather degree day", weather_index, "degree_day", 10)
    check_scale("Weather degree hour", ureg.Quantity(1, "degree_day_C"), "degree_hour_C", 24)
    heat_rate = ureg.Quantity(7.5, "MMBtu_th / MWh_e")
    expected_dimensions = (ureg.joule_thermal / ureg.joule_electric).dimensionality

    if heat_rate.dimensionality != expected_dimensions:
        raise AssertionError("Heat rate must retain energy-per-electrical-output dimensions.")

    check_scale("Heat rate quoting scale", ureg.Quantity(1, "heat_rate_macro"), "heat_rate_us", 1000)
    check_scale("Metric heat rate quoting scale", ureg.Quantity(1, "heat_rate_gas_turbine"), "heat_rate_metric", 1)
    check_scale("Hydrogen mass", ureg.Quantity(1, "tonne_H2"), "kg_H2", 1000)
    check_scale("Hydrogen normal volume", ureg.Quantity(1, "Nm3_H2"), "Nm3_0C_101325kPa_dry", 1)
    check_scale("SAF gallon", ureg.Quantity(1, "gallon_SAF"), "gallon", 1)

    ci_dimensions = (ureg.gram / ureg.megajoule).dimensionality
    if ureg.Quantity(35.5, "CI").dimensionality != ci_dimensions:
        raise AssertionError("Carbon intensity must have mass-per-energy dimensions.")

    check_scale("Carbon intensity mass basis", ureg.Quantity(1, "gCO2e"), "gram", 1)
    check_scale("Absolute pipeline pressure", ureg.Quantity(14.73, "psia"), "psi", 14.73)
    check_scale("Gauge pipeline pressure scale", ureg.Quantity(600, "psig"), "psi", 600)
    check_scale("Absolute bar pressure", ureg.Quantity(1, "bara"), "bar", 1)
    check_scale("Gauge bar pressure scale", ureg.Quantity(1, "barg"), "bar", 1)
    check_scale("Gold kilobar", ureg.Quantity(1, "kilobar_gold"), "toz", 1000 / 31.1034768)
    check_scale("Ten tola bar", ureg.Quantity(1, "ten_tola_bar"), "tola", 10)
    check_scale("Nominal gold delivery bar", ureg.Quantity(1, "good_delivery_gold"), "toz", 400)
    check_scale("Nominal silver delivery bar", ureg.Quantity(1, "good_delivery_silver"), "toz", 1000)
    check_scale("Power ramp scale", ureg.Quantity(1, "MW_per_minute"), "kW_per_minute", 1000)
    check_scale("Battery C-rate", ureg.Quantity(1, "C_rate"), "1 / hour", 1)
    check_scale("Air-dry pulp label", ureg.Quantity(1, "admt"), "metric_ton", 1)
    check_scale("Bone-dry pulp label", ureg.Quantity(1, "bdmt"), "metric_ton", 1)
    check_scale("Corn crop bushel", ureg.Quantity(1, "bushel_corn"), "pound", 56)
    check_scale("Canola crop bushel", ureg.Quantity(1, "bushel_canola"), "pound", 50)
    check_scale("Rice crop bushel", ureg.Quantity(1, "bushel_rice"), "pound", 45)
    check_scale("Metric absolute pressure", ureg.Quantity(101.325, "kpaa"), "kilopascal", 101.325)
    check_scale("Metric gauge pressure", ureg.Quantity(500, "kpag"), "kilopascal", 500)
    check_scale("Uranium tonne", ureg.Quantity(1, "mtU"), "kgU", 1000)
    check_scale("Metric tonne unit isolation", ureg.Quantity(1, "mtu"), "kilogram", 10)

    for unit, dimension in (("m3_20C_101325kPa_dry", "[standard_gas_volume_20C_101325kPa_dry]"), ("scf_60F_14696psi_dry", "[standard_gas_volume_60F_1atm_dry]"),):
        if ureg.Unit(unit).dimensionality != {dimension: 1}:
            raise AssertionError(f"Gas reference dimensional isolation failed for {unit}.")
    with _REGISTRY_LOCK:
        quotes = tuple(_QUOTE_REPLACEMENTS.items())
    for screen_quote, expression in quotes:
        try:
            ureg.Unit(expression)
        except (pint.PintError, ValueError, TypeError) as error:
            raise AssertionError(f"Invalid cached quote {screen_quote!r}: {expression!r}") from error

# ---------------------------------------------------------------------------
# Module Type Aliases & Public Exports
# ---------------------------------------------------------------------------

Unit: Final = ureg.Unit
Quantity: Final = ureg.Quantity

__all__ = ('BASELINE_DEFINITION_HASH', 'COMMODITY_UNIT_DEFINITIONS', 'CURRENCY_CODES', 'CURRENCY_MINOR_DIGITS', 'CURRENCY_MINOR_UNIT_EXPONENT', 'CURRENCY_REFERENCE_CODES', 'CURRENCY_SYMBOL_TO_CODE', 'ForeignRegistryError', 'PINT_COLLISION_ALLOWLIST', 'Q_', 'Quantity', 'QuoteComponents', 'QuoteParsingContext', 'QuoteTokenRule', 'REFERENCE_DATA_HASHES', 'REGISTRY_SCHEMA_VERSION', 'REGISTRY_VERSION', 'RegistryCompatibilityError', 'RegistrySnapshot', 'Unit', 'are_convertible', 'assert_registry_compatible', 'check_registry_compatibility', 'decompose_quote', 'ensure_unit', 'format_quantity', 'freeze_registry', 'get_registry', 'get_registry_definition_hash', 'get_registry_snapshot', 'is_commodity_price', 'is_currency', 'is_flow_rate', 'is_fx_rate', 'is_registry_frozen', 'parse_quantities', 'parse_quantity', 'parse_unit', 'quantity_from_dict', 'quantity_to_dict', 'register_commodity_unit', 'register_currency', 'sandbox_registry', 'sanitize_unit_string', 'ureg', 'validate_dimensionality', 'validate_registry_quantity', 'validate_unit_definition', 'verify_registry_integrity')
_PUBLIC_ALIASES = {'RegistryFrozenError': ('exceptions', 'RegistryFrozenError'), 'AmbiguousUnitError': ('exceptions', 'AmbiguousUnitError'), 'CURRENCY_CATALOG': ('currencies', 'CURRENCY_CATALOG'), '_SEMANTIC_SCOPE': ('validation', '_safety_SEMANTIC_SCOPE'), 'RegistryVersionMismatchError': ('exceptions', 'RegistryVersionMismatchError'), 'QuoteTokenRule': ('quote_parser', 'QuoteTokenRule'), 'CURRENCY_CONVENTIONS': ('currencies', 'CURRENCY_CONVENTIONS'), '_CURRENCY_DATA_DIRECTORY': ('currencies', 'DATA_DIRECTORY'), '_reference_document': ('standards', '_reference_reference_document'), 'CommodityUnitError': ('exceptions', 'CommodityUnitError'), 'QuoteParsingContext': ('quote_parser', 'QuoteParsingContext'), 'authorized_conversion': ('validation', 'authorized_conversion')}

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
    from core.units.quote_parser import QuoteParsingContext
    from core.units.quote_parser import QuoteTokenRule