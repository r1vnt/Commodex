"""Commodex types implementation and public interfaces."""
from __future__ import annotations

__all__ = ('Quantity', 'Unit')
_PUBLIC_ALIASES = {'Quantity': ('registry', 'Quantity'), 'Unit': ('registry', 'Unit')}

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
    from core.units.registry import Quantity
    from core.units.registry import Unit
