"""Commodex assays environmental water implementation and public interfaces."""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class WaterAssay(MaterialSpecification):
    salinity_mass_fraction: Decimal | None = None
    dissolved_solids_fraction: Decimal | None = None

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        for name in ("salinity_mass_fraction", "dissolved_solids_fraction"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, _converters_petroleum_fraction(getattr(self, name), name))
        if self.salinity_mass_fraction is not None and self.dissolved_solids_fraction is not None and self.salinity_mass_fraction > self.dissolved_solids_fraction:
            raise ValueError("Salt mass cannot exceed total dissolved solids.")

    def salt_mass(self, quantity):
        from core.units.assays.base import scale_mass
        if self.salinity_mass_fraction is None:
            raise ValueError("Explicit salinity is required.")
        return scale_mass(quantity, self.salinity_mass_fraction)

    def dissolved_solids_mass(self, quantity):
        from core.units.assays.base import scale_mass
        if self.dissolved_solids_fraction is None:
            raise ValueError("Explicit dissolved-solids fraction is required.")
        return scale_mass(quantity, self.dissolved_solids_fraction)

    def solvent_water_mass(self, quantity):
        from core.units.assays.base import scale_mass
        if self.dissolved_solids_fraction is None:
            raise ValueError("Explicit dissolved-solids fraction is required.")
        return scale_mass(quantity, 1 - self.dissolved_solids_fraction)

    def dilution_water_required(self, quantity, *, target_salinity, diluent_salinity=Decimal(0)):
        from core.units.assays.industrial.chemicals import dilution_addition
        if self.salinity_mass_fraction is None:
            raise ValueError("Explicit salinity is required.")
        return dilution_addition(quantity, self.salinity_mass_fraction, target_salinity, diluent_salinity)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (WaterAssay,))
__all__ = ('WaterAssay', 'convert_water', 'create_water_context')
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
    from core.units.assays.environmental.water import convert_water
    from core.units.assays.environmental.water import create_water_context
