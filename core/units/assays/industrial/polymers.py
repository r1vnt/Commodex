"""Commodex assays industrial polymers implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class PolymerAssay(MaterialSpecification):
    polymer: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        MaterialSpecification.__post_init__(self)
        if self.polymer not in {"PE", "HDPE", "LDPE", "LLDPE", "PP", "PVC", "PET", "PS", "ABS"}:
            raise ValueError("Unknown polymer family.")

    def pellet_mass(self, quantity, *, process_yield_fraction):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        return scale_mass(quantity, _converters_petroleum_fraction(process_yield_fraction, "polymer process yield"))





@dataclass(frozen=True, slots=True, kw_only=True)
class RecycledPolymerAssay(PolymerAssay):
    contamination_fraction: Decimal
    wash_recovery_fraction: Decimal
    melt_recovery_fraction: Decimal

    def __post_init__(self):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        PolymerAssay.__post_init__(self)
        for name in ("contamination_fraction", "wash_recovery_fraction", "melt_recovery_fraction"):
            object.__setattr__(self, name, _converters_petroleum_fraction(getattr(self, name), name))

    def clean_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, (1 - self.contamination_fraction) * self.wash_recovery_fraction)

    def pellet_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(self.clean_mass(quantity), self.melt_recovery_fraction)

    def recycled_feed_required(self, quantity):
        from core.units.assays.base import scale_mass
        fraction = (1 - self.contamination_fraction) * self.wash_recovery_fraction * self.melt_recovery_fraction
        if fraction == 0:
            raise ValueError("Cannot invert zero recycled-polymer recovery.")
        return scale_mass(quantity, 1 / fraction)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (PolymerAssay, RecycledPolymerAssay,))
__all__ = ('PolymerAssay', 'RecycledPolymerAssay', 'convert_polymer', 'convert_recycled_polymer', 'create_polymer_context', 'create_recycled_polymer_context')
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
    from core.units.assays.industrial.polymers import convert_polymer
    from core.units.assays.industrial.polymers import convert_recycled_polymer
    from core.units.assays.industrial.polymers import create_polymer_context
    from core.units.assays.industrial.polymers import create_recycled_polymer_context
