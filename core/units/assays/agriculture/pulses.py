"""Commodex assays agriculture pulses implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from typing import ClassVar
import pint

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class PulseAssay(_agricultural_families_AgriculturalFamilyAssay):
    """Dry pulses with optional whole-to-split processing yield."""
    split_yield_fraction: Decimal | None = None
    _commodities = frozenset({"peas", "lentils", "chickpeas", "dry-beans"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "split_yield_fraction", _agricultural_families_yield(self.split_yield_fraction, "split yield"))

    def split_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.split_yield_fraction, "split yield")

    def whole_pulse_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.split_yield_fraction, "split yield")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (PulseAssay,))
__all__ = ('PulseAssay', 'convert_pulse', 'create_pulse_context')
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
    from core.units.assays.agriculture.pulses import convert_pulse
    from core.units.assays.agriculture.pulses import create_pulse_context
