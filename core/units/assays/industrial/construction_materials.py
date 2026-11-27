"""Commodex assays industrial construction_materials implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class RefractoryMaterialAssay(MaterialSpecification):
    loss_on_ignition_fraction: Decimal | None = None

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        if self.loss_on_ignition_fraction is not None:
            object.__setattr__(self, "loss_on_ignition_fraction", _converters_petroleum_fraction(self.loss_on_ignition_fraction, "loss on ignition", allow_one=False))

    def calcined_mass(self, quantity, *, mass_basis="dry"):
        from core.units.assays.base import scale_mass
        if self.loss_on_ignition_fraction is None:
            raise ValueError("Declare a loss on ignition for this material/process.")
        return scale_mass(self._dry_feed(quantity, mass_basis), 1 - self.loss_on_ignition_fraction)

    def uncalcined_mass(self, quantity):
        from core.units.assays.base import scale_mass
        if self.loss_on_ignition_fraction is None:
            raise ValueError("Declare loss on ignition.")
        return scale_mass(quantity, 1 / (1 - self.loss_on_ignition_fraction))


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class CementClinkerAssay(MaterialSpecification):
    clinker_fraction: Decimal | None = None

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        if self.clinker_fraction is not None:
            fraction = _converters_petroleum_fraction(self.clinker_fraction, "clinker ratio")
            if fraction == 0:
                raise ValueError("Clinker ratio must be positive.")
            object.__setattr__(self, "clinker_fraction", fraction)

    def cement_from_clinker(self, quantity):
        from core.units.assays.base import scale_mass
        if self.clinker_fraction is None:
            raise ValueError("Declare the cement clinker ratio.")
        return scale_mass(quantity, 1 / self.clinker_fraction)

    def clinker_from_cement(self, quantity):
        from core.units.assays.base import scale_mass
        if self.clinker_fraction is None:
            raise ValueError("Declare the cement clinker ratio.")
        return scale_mass(quantity, self.clinker_fraction)

    def calcination_co2_mass(self, quantity, *, limestone_purity_fraction, conversion_fraction, mass_basis="dry"):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.composition import StoichiometricCompound
        reacted = scale_mass(self._dry_feed(quantity, mass_basis), _converters_petroleum_fraction(limestone_purity_fraction, "limestone purity") * _converters_petroleum_fraction(conversion_fraction, "calcination conversion"))
        # CaCO3 -> CaO + CO2. Remaining feed is unchanged inert/unreacted solid.
        ratio = StoichiometricCompound("CO2").molar_mass.magnitude / StoichiometricCompound("CaCO3").molar_mass.magnitude
        return scale_mass(reacted, ratio)

    def calcined_solid_mass(self, quantity, **conditions):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        mass_basis = conditions.get("mass_basis", "dry")
        dry = self._dry_feed(quantity, mass_basis)
        released = self.calcination_co2_mass(quantity, **conditions)
        return _converters_petroleum_output(mass(dry) - mass(released), "kg", quantity).to(quantity.units)


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class AggregateAssay(MaterialSpecification):
    def compacted_volume(self, quantity, *, volume_reduction_fraction):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "loose volume", minimum=0)
        fraction = _converters_petroleum_fraction(volume_reduction_fraction, "volume reduction", allow_one=False)
        return _converters_petroleum_output(volume * (1 - fraction), "m^3", quantity)

    def loose_volume(self, quantity, *, volume_reduction_fraction):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "compacted volume", minimum=0)
        fraction = _converters_petroleum_fraction(volume_reduction_fraction, "volume reduction", allow_one=False)
        return _converters_petroleum_output(volume / (1 - fraction), "m^3", quantity)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (RefractoryMaterialAssay, CementClinkerAssay, AggregateAssay,))
__all__ = ('AggregateAssay', 'CementClinkerAssay', 'RefractoryMaterialAssay', 'convert_aggregate', 'convert_cement_clinker', 'convert_refractory_material', 'create_aggregate_context', 'create_cement_clinker_context', 'create_refractory_material_context')
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
    from core.units.assays.industrial.construction_materials import convert_aggregate
    from core.units.assays.industrial.construction_materials import convert_cement_clinker
    from core.units.assays.industrial.construction_materials import convert_refractory_material
    from core.units.assays.industrial.construction_materials import create_aggregate_context
    from core.units.assays.industrial.construction_materials import create_cement_clinker_context
    from core.units.assays.industrial.construction_materials import create_refractory_material_context
