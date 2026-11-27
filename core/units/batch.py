"""Commodex batch implementation and public interfaces."""
from __future__ import annotations






def convert_many(quantities, target_unit, **options):
    """Convert in input order; failures propagate without returning partial data."""
    from core.units.conversions import convert
    return tuple(convert(quantity, target_unit, **options) for quantity in quantities)


__all__ = ('convert_many',)
_PUBLIC_ALIASES = {'convert': ('conversions', 'convert')}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
