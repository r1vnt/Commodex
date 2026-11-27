"""Commodex assays industrial industrial_gases implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class IndustrialGasAssay(MaterialSpecification):
    gas_composition: GasComposition
    reference: GasReferenceConditions
    backend: GasPropertyBackend

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.composition import GasComposition
        from core.units.reference_conditions import GasReferenceConditions
        from core.units.thermodynamics import GasPropertyBackend
        from core.units.thermodynamics import gas_density
        if not isinstance(self.gas_composition, GasComposition) or not isinstance(self.reference, GasReferenceConditions) or not isinstance(self.backend, GasPropertyBackend):
            raise TypeError("Supply gas composition, reference conditions and an explicit backend.")
        evaluated = gas_density(self.gas_composition, self.reference, backend=self.backend)
        if self.density is not None:
            measured = _converters_petroleum_positive_quantity(self.density, "kg/m^3", "gas density")
            if abs(measured / evaluated.magnitude - 1) > Decimal("1e-8"):
                raise ValueError("Declared gas density conflicts with the selected EOS/state.")
        object.__setattr__(self, "density", evaluated)
        object.__setattr__(self, "reference_temperature", self.reference.temperature)
        object.__setattr__(self, "reference_pressure", self.reference.absolute_pressure)
        MaterialSpecification.__post_init__(self)

    def gas_volume_to_mass(self, quantity):
        from core.units.assays.energy.lng import _converters_gas_geometric_volume
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        return _converters_petroleum_output(_converters_gas_geometric_volume(quantity, self.reference) * self.density.magnitude, "kg", quantity)

    def mass_to_gas_volume(self, quantity, *, target_unit="m^3"):
        from core.units.assays.base import mass
        from core.units.assays.energy.lng import _converters_gas_gas_volume_output
        return _converters_gas_gas_volume_output(mass(quantity) / self.density.magnitude, target_unit, self.reference, quantity)

    def as_context(self, *, name=None):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.base import add_molar_context
        from core.units.registry import ureg
        from core.units.conversions import _GAS_REFERENCE_PRESETS
        context = MaterialSpecification.as_context(self, name=name)
        for tagged, _geometric, _temperature, _pressure in _GAS_REFERENCE_PRESETS.values():
            unit = ureg.Unit(tagged)
            context.add_transformation(unit.dimensionality, "[mass]", lambda _r, q: self.gas_volume_to_mass(q))

            def inverse(destination):
                return lambda _r, q: self.mass_to_gas_volume(q, target_unit=destination)

            context.add_transformation("[mass]", unit.dimensionality, inverse(unit))
        return add_molar_context(context, self.gas_composition.molar_mass)

    def as_heat_content(self, **overrides):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.lng import _converters_gas_geometric_volume
        from core.units.registry import ureg
        if self.heat_content is not None or {"gcv", "ncv"} & overrides.keys():
            return MaterialSpecification.as_heat_content(self, **overrides)
        from core.units.heat_content import HeatContentAssay
        values = self.gas_composition.calorific_values()
        if values.hhv.magnitude == 0:
            raise ValueError("This industrial gas has no combustion calorific value.")
        standard_density = None
        from core.units.conversions import _GAS_REFERENCE_PRESETS
        for tagged, geometric, _temperature, _pressure in _GAS_REFERENCE_PRESETS.values():
            try:
                reference_volume = _converters_gas_geometric_volume(ureg.Quantity(Decimal(1), tagged), self.reference)
            except ValueError:
                continue
            standard_density = ureg.Quantity(reference_volume * self.density.magnitude, ureg.kg / ureg.Unit(tagged))
            break
        specification = HeatContentAssay(name=self.name, gcv=values.hhv, ncv=values.lhv, density=self.density, standard_density=standard_density,
                                molar_mass=self.gas_composition.molar_mass,
                                reference_temperature=self.reference.temperature, reference_pressure=self.reference.absolute_pressure)
        from dataclasses import replace
        return replace(specification, **overrides)


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class RefrigerantAssay(MaterialSpecification):
    blend_mass_fractions: Mapping[str, Decimal]
    component_molar_masses: Mapping[str, pint.Quantity]

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.base import _spec_engine_QuantityMap
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.composition import ElementComposition
        MaterialSpecification.__post_init__(self)
        blend = ElementComposition(self.blend_mass_fractions)
        if blend.unassayed_fraction != 0:
            raise ValueError("A complete refrigerant blend must sum to one.")
        if set(self.component_molar_masses) != set(blend.mass_fractions):
            raise ValueError("Provide molar mass for every refrigerant component.")
        masses = {k: _converters_petroleum_positive_quantity(v, "g/mol", "refrigerant molar mass") for k, v in self.component_molar_masses.items()}
        object.__setattr__(self, "blend_mass_fractions", blend.mass_fractions)
        # Canonical scalar values prevent mutation of nested quantity metadata.
        object.__setattr__(self, "component_molar_masses", _spec_engine_QuantityMap(masses, "g/mol"))

    @property
    def molar_mass(self):
        from core.units.registry import ureg
        reciprocal = sum((w / self.component_molar_masses[k].magnitude for k, w in self.blend_mass_fractions.items()), Decimal(0))
        return ureg.Quantity(1 / reciprocal, "g/mol")

    @property
    def mole_fractions(self):
        amounts = {k: w / self.component_molar_masses[k].magnitude for k, w in self.blend_mass_fractions.items()}
        total = sum(amounts.values(), Decimal(0))
        return MappingProxyType({k: v / total for k, v in amounts.items()})

    def component_mass(self, quantity, *, component):
        from core.units.assays.base import scale_mass
        if component not in self.blend_mass_fractions:
            raise ValueError("Unknown refrigerant blend component.")
        return scale_mass(quantity, self.blend_mass_fractions[component])

    def as_context(self, *, name=None):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.base import add_molar_context
        return add_molar_context(MaterialSpecification.as_context(self, name=name), self.molar_mass)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (IndustrialGasAssay, RefrigerantAssay,))
__all__ = ('IndustrialGasAssay', 'RefrigerantAssay', 'convert_industrial_gas', 'convert_refrigerant', 'create_industrial_gas_context', 'create_refrigerant_context')
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
    from core.units.composition import GasComposition
    from core.units.thermodynamics import GasPropertyBackend
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.assays.industrial.industrial_gases import convert_industrial_gas
    from core.units.assays.industrial.industrial_gases import convert_refrigerant
    from core.units.assays.industrial.industrial_gases import create_industrial_gas_context
    from core.units.assays.industrial.industrial_gases import create_refrigerant_context
