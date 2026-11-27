"""Commodex assays energy synthetic_fuels implementation and public interfaces."""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class SyntheticFuelAssay(MaterialSpecification):
    fuel: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        MaterialSpecification.__post_init__(self)
        if self.fuel not in {"e-methanol", "e-diesel", "e-kerosene", "synthetic-methane"}:
            raise ValueError("Declare e-methanol, e-diesel, e-kerosene or synthetic-methane.")

    def fuel_energy(self, quantity, *, basis="lhv", target_unit="MWh"):
        from core.units.heat_content import convert_heat_content
        return convert_heat_content(quantity, target_unit, assay=self.as_heat_content(), basis=basis)


from core.units.assays.industrial.chemicals import ChemicalSolutionAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class MethanolFuelAssay(ChemicalSolutionAssay):
    solute_formula: str = "CH3OH"

    def __post_init__(self):
        from core.units.assays.industrial.chemicals import ChemicalSolutionAssay
        ChemicalSolutionAssay.__post_init__(self)
        if self.solute_formula not in {"CH3OH", "CH4O"}:
            raise ValueError("Methanol formula must be CH3OH or CH4O.")

    def fuel_energy(self, quantity, *, basis="lhv", target_unit="MWh"):
        from core.units.heat_content import convert_heat_content
        return convert_heat_content(quantity, target_unit, assay=self.as_heat_content(), basis=basis)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (SyntheticFuelAssay, MethanolFuelAssay,))
__all__ = ('MethanolFuelAssay', 'SyntheticFuelAssay', 'convert_methanol_fuel', 'convert_synthetic_fuel', 'create_methanol_fuel_context', 'create_synthetic_fuel_context')
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
    from core.units.assays.energy.synthetic_fuels import convert_methanol_fuel
    from core.units.assays.energy.synthetic_fuels import convert_synthetic_fuel
    from core.units.assays.energy.synthetic_fuels import create_methanol_fuel_context
    from core.units.assays.energy.synthetic_fuels import create_synthetic_fuel_context
