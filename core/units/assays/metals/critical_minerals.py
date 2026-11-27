"""Commodex assays metals critical_minerals implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, field
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class ExplorationBrineAssay(MaterialSpecification):
    concentrations: Mapping[str, pint.Quantity]
    extraction_recoveries: Mapping[str, Decimal] = field(default_factory=dict)

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.base import _spec_engine_QuantityMap
        from core.units.assays.base import fractions
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        MaterialSpecification.__post_init__(self)
        if not self.concentrations:
            raise ValueError("Supply measured brine concentrations with mass/volume units.")
        concentrations = {_converters_petroleum_text(k, "brine component"): _converters_petroleum_quantity_number(v, "kg/m^3", "brine concentration", minimum=0) for k, v in self.concentrations.items()}
        recoveries = fractions(self.extraction_recoveries, "extraction recoveries")
        if recoveries.keys() - concentrations.keys():
            raise ValueError("Extraction recovery refers to an unassayed brine component.")
        if self.density is not None and sum(concentrations.values(), Decimal(0)) > self.density.magnitude:
            raise ValueError("Dissolved mass exceeds total brine density.")
        object.__setattr__(self, "concentrations", _spec_engine_QuantityMap(concentrations, "kg/m^3"))
        object.__setattr__(self, "extraction_recoveries", recoveries)

    def dissolved_mass(self, quantity, *, component):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if component not in self.concentrations:
            raise ValueError("Unknown brine component.")
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "brine volume", minimum=0)
        return _converters_petroleum_output(volume * self.concentrations[component].magnitude, "kg", quantity)

    def recoverable_mass(self, quantity, *, component):
        from core.units.assays.base import scale_mass
        if component not in self.extraction_recoveries:
            raise ValueError("Explicit extraction recovery is required.")
        return scale_mass(self.dissolved_mass(quantity, component=component), self.extraction_recoveries[component])

    def brine_required(self, quantity, *, component):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        if component not in self.extraction_recoveries or component not in self.concentrations:
            raise ValueError("Explicit concentration and extraction recovery are required.")
        factor = self.concentrations[component].magnitude * self.extraction_recoveries[component]
        if factor == 0:
            raise ValueError("Cannot invert zero extraction.")
        return _converters_petroleum_output(mass(quantity) / factor, "m^3", quantity)


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class CriticalMineralAssay(MaterialSpecification):
    commodity: str
    reported_compound: str | None = None

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.composition import StoichiometricCompound
        MaterialSpecification.__post_init__(self)
        object.__setattr__(self, "commodity", _converters_petroleum_text(self.commodity, "critical mineral"))
        if self.reported_compound is not None:
            StoichiometricCompound(self.reported_compound)

    def reported_to_element(self, quantity, *, element):
        from core.units.composition import StoichiometricCompound
        if self.reported_compound is None:
            raise ValueError("Declare the reported chemical compound.")
        return StoichiometricCompound(self.reported_compound).contained_element(quantity, element=element)

    def purified_mass(self, quantity, *, component, target_fraction, recovery_fraction, mass_basis="dry"):
        from core.units.assays.base import beneficiated_mass
        return beneficiated_mass(self, quantity, component=component, target_fraction=target_fraction,
                                  recovery_fraction=recovery_fraction, mass_basis=mass_basis)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (ExplorationBrineAssay, CriticalMineralAssay,))
__all__ = ('CriticalMineralAssay', 'ExplorationBrineAssay', 'convert_critical_mineral', 'convert_exploration_brine', 'create_critical_mineral_context', 'create_exploration_brine_context')
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
    from core.units.assays.metals.critical_minerals import convert_critical_mineral
    from core.units.assays.metals.critical_minerals import convert_exploration_brine
    from core.units.assays.metals.critical_minerals import create_critical_mineral_context
    from core.units.assays.metals.critical_minerals import create_exploration_brine_context
