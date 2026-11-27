"""Commodex logistics_units implementation and public interfaces."""
from __future__ import annotations





LOGISTICS_UNITS = ("item", "dozen", "contract", "lot", "delivery_unit", "crate", "sack", "gross_tonnage", "net_tonnage")





def logistics_unit(name):
    from core.units.registry import parse_unit
    if name not in LOGISTICS_UNITS:
        raise ValueError(f"Unknown logistics unit: {name!r}")
    return parse_unit(name)


__all__ = ('LOGISTICS_UNITS', 'logistics_unit')
_PUBLIC_ALIASES = {'parse_unit': ('registry', 'parse_unit')}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
