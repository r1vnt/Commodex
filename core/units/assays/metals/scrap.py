"""Commodex assays metals scrap implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class MetalScrapAssay(MaterialSpecification):
    scrap_grade: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        MaterialSpecification.__post_init__(self)
        object.__setattr__(self, "scrap_grade", _converters_petroleum_text(self.scrap_grade, "scrap grade"))
        if self.composition is None:
            raise ValueError("Scrap requires its multi-element composition.")

    def unrecovered_mass(self, quantity, *, mass_basis="dry"):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        dry = self._dry_feed(quantity, mass_basis)
        outputs = sum((mass(self.recovered_mass(dry, component=k)) for k in self.composition.mass_fractions), Decimal(0))
        return _converters_petroleum_output(mass(dry) - outputs, "kg", quantity)


from core.units.assays.metals.concentrates import MultiMetalConcentrateAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class BatteryRecyclingAssay(MultiMetalConcentrateAssay):
    def black_mass_from_pack(self, quantity, *, black_mass_yield_fraction):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        return scale_mass(quantity, _converters_petroleum_fraction(black_mass_yield_fraction, "black mass yield"))


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (MetalScrapAssay, BatteryRecyclingAssay,))
__all__ = ('BatteryRecyclingAssay', 'MetalScrapAssay', 'convert_battery_recycling', 'convert_metal_scrap', 'create_battery_recycling_context', 'create_metal_scrap_context')
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
    from core.units.assays.metals.scrap import convert_battery_recycling
    from core.units.assays.metals.scrap import convert_metal_scrap
    from core.units.assays.metals.scrap import create_battery_recycling_context
    from core.units.assays.metals.scrap import create_metal_scrap_context
