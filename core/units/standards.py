"""Commodex standards implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from datetime import date
from decimal import Decimal
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from types import MappingProxyType
import json
import re



REFERENCE_VERSION = "2026-11-27"



REFERENCE_DIRECTORY = (Path(__file__).resolve().parents[0] / "reference_data")



_reference_TOKEN = re.compile(r"[A-Za-z0-9_.-]+\Z")





def _reference_token(value):
    if not isinstance(value, str) or not _reference_TOKEN.fullmatch(value) or value in {".", ".."}:
        raise ValueError("Reference names and versions must be simple identifiers.")
    return value





def _reference_freeze(value):
    if isinstance(value, dict):
        return MappingProxyType({key: _reference_freeze(item) for key, item in value.items()})
    if isinstance(value, list):
        return tuple(_reference_freeze(item) for item in value)
    return value





def _reference_decode(value):
    if isinstance(value, Mapping):
        if set(value) == {"decimal"}:
            return Decimal(value["decimal"])
        if set(value) == {"quantity", "unit"}:
            from core.units.registry import ureg
            return ureg.Quantity(_reference_decode(value["quantity"]), value["unit"])
        if set(value) == {"tuple"}:
            return tuple(_reference_decode(item) for item in value["tuple"])
        return {key: _reference_decode(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_reference_decode(item) for item in value]
    return value





@dataclass(frozen=True, slots=True)
class ReferenceRecord:
    id: str
    version: str
    category: str
    source: Mapping
    valid_from: date | None
    valid_to: date | None
    units: Mapping
    uncertainty: Mapping
    value_kind: str
    data: Mapping

    def constructor_arguments(self):
        """Return fresh arguments so callers cannot mutate the cached catalogue."""
        payload = _reference_decode(self.data)
        return payload["args"], payload["kwargs"]





@lru_cache(maxsize=None)
def _reference_reference_document(section, category, version=REFERENCE_VERSION):
    if section not in {"assays", "units", "standards", "currencies", "aliases"}:
        raise ValueError("Unknown reference section.")
    _reference_token(category)
    _reference_token(version)
    directory = REFERENCE_DIRECTORY
    manifest = json.loads((directory / "standards" / version / "manifest.json").read_text(encoding="utf-8"))
    if manifest["schema_version"] != 1 or manifest["version"] != version:
        raise ValueError("Unsupported reference manifest schema or version.")
    relative = f"{section}/{category}.json"
    payload = (directory / section / version / (category + ".json")).read_bytes()
    if sha256(payload).hexdigest() != manifest["files"][relative]:
        raise ValueError(f"Reference data integrity failure: {relative}")
    document = json.loads(payload)
    if document["schema_version"] != 1 or document["version"] != version:
        raise ValueError("Reference catalogue schema/version mismatch.")
    return _reference_freeze(document)





@lru_cache(maxsize=None)
def _reference_catalogue(category, version):
    document = _reference_reference_document("assays", category, version)
    records = {}
    for raw in document["records"]:
        required = {"id", "source", "valid_from", "valid_to", "units", "uncertainty", "value_kind", "data"}
        if not required <= raw.keys() or raw["id"] in records:
            raise ValueError("Incomplete or duplicate reference record.")
        if raw["value_kind"] not in {"measured_value", "commercial_convention", "indicative_benchmark"}:
            raise ValueError("Unknown reference value kind.")
        if not raw["source"].get("citation") or not raw["uncertainty"].get("status"):
            raise ValueError("Reference source and uncertainty status are required.")
        start = date.fromisoformat(raw["valid_from"]) if raw["valid_from"] else None
        end = date.fromisoformat(raw["valid_to"]) if raw["valid_to"] else None
        if start and end and start > end:
            raise ValueError("Invalid reference validity interval.")
        records[raw["id"]] = ReferenceRecord(raw["id"], version, category,
            _reference_freeze(raw["source"]), start, end, _reference_freeze(raw["units"]),
            _reference_freeze(raw["uncertainty"]), raw["value_kind"], _reference_freeze(raw["data"]))
    return MappingProxyType(records)





def reference_record(category, identifier, *, version=REFERENCE_VERSION, as_of=None):
    """Find a pinned record. Validity bounds are inclusive; null means unknown.

    An as-of lookup rejects unknown validity rather than inventing market history.
    """
    record = _reference_catalogue(category, version)[identifier]
    if as_of is not None:
        when = as_of if type(as_of) is date else date.fromisoformat(as_of)
        if record.valid_from is None or record.valid_to is None:
            raise ValueError("Reference validity is unknown; cannot guarantee an as-of value.")
        if not record.valid_from <= when <= record.valid_to:
            raise ValueError("Reference record is outside its validity period.")
    return record





def reference_records(category, *, version=REFERENCE_VERSION):
    return tuple(_reference_catalogue(category, version).values())





def construct_reference(cls, category, identifier, *, version=REFERENCE_VERSION, overrides=None):
    record = reference_record(category, identifier, version=version)
    args, kwargs = record.constructor_arguments()
    if overrides:
        kwargs.update(overrides)
    return cls(*args, **kwargs)




STANDARD_CATEGORIES = ("metrology", "petroleum", "natural_gas", "agriculture", "metals")





def recommended_standards(category=None, *, version=REFERENCE_VERSION):
    categories = STANDARD_CATEGORIES if category is None else (category,)
    if any(item not in STANDARD_CATEGORIES for item in categories):
        raise ValueError("Unknown standards category.")
    return tuple(record for item in categories
                 for record in _reference_reference_document("standards", item, version)["records"])


__all__ = ('REFERENCE_VERSION', 'ReferenceRecord', 'STANDARD_CATEGORIES', 'recommended_standards', 'reference_record', 'reference_records')
_PUBLIC_ALIASES = {'_reference_document': ('standards', '_reference_reference_document')}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
