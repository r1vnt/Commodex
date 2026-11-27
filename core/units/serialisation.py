"""Commodex serialisation implementation and public interfaces."""
from __future__ import annotations



__all__ = ('quantity_from_dict', 'quantity_to_dict')
_PUBLIC_ALIASES = {'quantity_from_dict': ('registry', 'quantity_from_dict'), 'quantity_to_dict': ('registry', 'quantity_to_dict')}

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
    from core.units.registry import quantity_from_dict
    from core.units.registry import quantity_to_dict
