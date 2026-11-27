"""Commodex assays industrial industrial_minerals implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class IndustrialMineralAssay(MaterialSpecification):
    mineral: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        MaterialSpecification.__post_init__(self)
        object.__setattr__(self, "mineral", _converters_petroleum_text(self.mineral, "mineral"))

    def commercial_component_mass(self, quantity, *, component, mass_basis="dry"):
        return self.contained_mass(quantity, component=component, mass_basis=mass_basis)


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class GemstoneAssay(MaterialSpecification):
    def polished_mass(self, quantity, *, recovery_fraction):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        return scale_mass(quantity, _converters_petroleum_fraction(recovery_fraction, "polishing recovery"))

    def rough_mass(self, quantity, *, recovery_fraction):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        fraction = _converters_petroleum_fraction(recovery_fraction, "polishing recovery")
        if fraction == 0:
            raise ValueError("Cannot invert zero polishing recovery.")
        return scale_mass(quantity, 1 / fraction)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (IndustrialMineralAssay, GemstoneAssay,))
__all__ = ('GemstoneAssay', 'IndustrialMineralAssay', 'convert_gemstone', 'convert_industrial_mineral', 'create_gemstone_context', 'create_industrial_mineral_context')
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
    from core.units.assays.industrial.industrial_minerals import convert_gemstone
    from core.units.assays.industrial.industrial_minerals import convert_industrial_mineral
    from core.units.assays.industrial.industrial_minerals import create_gemstone_context
    from core.units.assays.industrial.industrial_minerals import create_industrial_mineral_context
