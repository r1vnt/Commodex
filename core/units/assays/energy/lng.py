"""Commodex assays energy lng implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass, field, replace
from decimal import Decimal
from types import MappingProxyType
import pint




_converters_gas_AIR_MOLAR_MASS = Decimal("28.96546")  # same conventional dry-air mass as conversions.py



_converters_gas_CP_NAMES = {
    "CH4": "Methane", "C2H6": "Ethane", "C3H8": "Propane", "iC4H10": "IsoButane", "nC4H10": "n-Butane",
    "iC5H12": "Isopentane", "nC5H12": "n-Pentane", "C6H14": "n-Hexane", "C7H16": "n-Heptane",
    "C8H18": "n-Octane", "C9H20": "n-Nonane", "C10H22": "n-Decane", "H2": "Hydrogen",
    "N2": "Nitrogen", "CO2": "CarbonDioxide", "CO": "CarbonMonoxide", "O2": "Oxygen",
    "H2O": "Water", "H2S": "HydrogenSulfide", "He": "Helium", "Ar": "Argon",
}





def _converters_gas_geometric_volume(volume: pint.Quantity, reference: GasReferenceConditions) -> Decimal:
    from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
    from core.units.heat_content import _converters_energy_resolve_gas_reference_label
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.conversions import _GAS_REFERENCE_PRESETS
    if not isinstance(volume, pint.Quantity):
        raise TypeError("Gas volume requires an explicit quantity.")
    validate_registry_quantity(volume)
    if volume.dimensionality == (ureg.meter**3).dimensionality:
        return _converters_petroleum_quantity_number(volume, "meter**3", "gas volume", minimum=0)
    key = _converters_energy_resolve_gas_reference_label(volume.units, reference)
    tagged, geometric, _, _ = _GAS_REFERENCE_PRESETS[key]
    magnitude = _converters_petroleum_quantity_number(volume, tagged, "tagged gas volume", minimum=0)
    return _converters_petroleum_quantity_number(ureg.Quantity(magnitude, geometric), "meter**3", "gas volume", minimum=0)





def _converters_gas_gas_volume_output(volume_m3: Decimal, unit: str, reference: GasReferenceConditions, source: pint.Quantity) -> pint.Quantity:
    from core.units.assays.energy.petroleum import _converters_petroleum_output
    from core.units.heat_content import _converters_energy_resolve_gas_reference_label
    from core.units.registry import ureg
    from core.units.conversions import _GAS_REFERENCE_PRESETS
    unit = ureg.Unit(unit)
    if unit.dimensionality == (ureg.meter**3).dimensionality:
        result = ureg.Quantity(volume_m3, "meter**3").to(unit)
    else:
        key = _converters_energy_resolve_gas_reference_label(unit, reference)
        tagged, geometric, _, _ = _GAS_REFERENCE_PRESETS[key]
        number = ureg.Quantity(volume_m3, "meter**3").to(geometric).magnitude
        result = ureg.Quantity(number, tagged).to(unit)
    return _converters_petroleum_output(result.magnitude, result.units, source)





def composition_gas_energy(volume: pint.Quantity, composition: GasComposition, reference: GasReferenceConditions, *,
                            backend: GasPropertyBackend, basis: str,
                            calorific_values: GasCalorificValues | None = None, target_unit: str = "MJ") -> pint.Quantity:
    from core.units.assays.energy.petroleum import _converters_petroleum_output
    from core.units.heat_content import _converters_energy_normalize_basis
    from core.units.registry import ureg
    from core.units.thermodynamics import gas_density
    basis = _converters_energy_normalize_basis(basis)
    values = calorific_values if calorific_values is not None else composition.calorific_values()
    if reference.combustion_reference_temperature is not None and abs(reference.combustion_reference_temperature.to("K").magnitude - values.combustion_temperature.magnitude) > Decimal("1e-8"):
        raise ValueError("Gas energy and combustion reference temperatures differ.")
    mass = _converters_gas_geometric_volume(volume, reference) * gas_density(composition, reference, backend=backend).magnitude
    specific = values.hhv if basis == "HHV" else values.lhv
    energy = ureg.Quantity(mass * specific.magnitude, "MJ").to(target_unit)
    return _converters_petroleum_output(energy.magnitude, energy.units, volume)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LNGCargo(_safety_QuantitySnapshot):
    name: str
    composition: GasComposition
    liquid_state: LiquidState
    measured_density: pint.Quantity | None = None
    density_source: str | None = None
    density_backend: GasPropertyBackend | None = None
    calorific_values: GasCalorificValues | None = None

    def __post_init__(self) -> None:
        from core.units.assays.energy.petroleum import LiquidState
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.composition import GasCalorificValues
        from core.units.composition import GasComposition
        from core.units.registry import ureg
        from core.units.thermodynamics import GasPropertyBackend
        object.__setattr__(self, "name", _converters_petroleum_text(self.name, "cargo name"))
        if not isinstance(self.composition, GasComposition) or not isinstance(self.liquid_state, LiquidState):
            raise TypeError("LNG cargo requires composition and liquid state.")
        if (self.measured_density is None) == (self.density_backend is None):
            raise ValueError("Supply measured cargo density or a composition/temperature density backend, exactly one.")
        if self.measured_density is not None:
            density = _converters_petroleum_positive_quantity(self.measured_density, "kg/meter**3", "measured LNG density")
            object.__setattr__(self, "measured_density", ureg.Quantity(density, "kg/meter**3"))
            object.__setattr__(self, "density_source", _converters_petroleum_text(self.density_source, "density measurement source"))
        elif not isinstance(self.density_backend, GasPropertyBackend):
            raise TypeError("density_backend must be a GasPropertyBackend.")
        else:
            if self.density_source is not None and self.density_source != self.density_backend.source:
                raise ValueError("Cargo density source conflicts with its selected backend.")
            object.__setattr__(self, "density_source", self.density_backend.source)
        if self.composition.water_mole_fraction:
            raise ValueError("LNG recipe must be water-free; ice/water requires a dedicated multiphase model.")
        if self.calorific_values is None:
            object.__setattr__(self, "calorific_values", self.composition.calorific_values())
        elif not isinstance(self.calorific_values, GasCalorificValues):
            raise TypeError("Supply measured HHV/LHV as GasCalorificValues.")
        if self.calorific_values.hhv.magnitude <= 0:
            raise ValueError("LNG cargo requires positive calorific value; noncombustible cryogenic gases are separate materials.")

    @property
    def liquid_density(self) -> pint.Quantity:
        if self.measured_density is not None:
            return self.measured_density
        return self.density_backend.liquid_density(self.composition, self.liquid_state)

    def mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        volume_m3 = _converters_petroleum_quantity_number(volume, "meter**3", "LNG liquid volume", minimum=0)
        return _converters_petroleum_output(volume_m3 * self.liquid_density.magnitude, "kg", volume)

    def energy(self, volume: pint.Quantity, *, basis: str, target_unit: str = "MMBtu_it") -> pint.Quantity:
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.heat_content import _converters_energy_normalize_basis
        from core.units.registry import ureg
        basis = _converters_energy_normalize_basis(basis)
        specific = self.calorific_values.hhv if basis == "HHV" else self.calorific_values.lhv
        mass = _converters_petroleum_quantity_number(self.mass(volume), "kg", "cargo mass", minimum=0)
        energy = ureg.Quantity(mass * specific.magnitude, "MJ").to(target_unit)
        return _converters_petroleum_output(energy.magnitude, energy.units, volume)

    def gas_equivalent(self, volume: pint.Quantity, reference: GasReferenceConditions, *, backend: GasPropertyBackend,
                       target_unit: str = "meter**3") -> pint.Quantity:
        """Complete vaporization without fractionation, water addition or component loss."""
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.thermodynamics import gas_density
        mass = _converters_petroleum_quantity_number(self.mass(volume), "kg", "cargo mass", minimum=0)
        density = gas_density(self.composition, reference, backend=backend).magnitude
        return _converters_gas_gas_volume_output(mass / density, target_unit, reference, volume)

    def liquid_equivalent(self, gas_volume: pint.Quantity, reference: GasReferenceConditions, *, backend: GasPropertyBackend) -> pint.Quantity:
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.thermodynamics import gas_density
        mass = _converters_gas_geometric_volume(gas_volume, reference) * gas_density(self.composition, reference, backend=backend).magnitude
        return _converters_petroleum_output(mass / self.liquid_density.magnitude, "meter**3", gas_volume)

    def expansion_ratio(self, reference: GasReferenceConditions, *, backend: GasPropertyBackend) -> Decimal:
        from core.units.thermodynamics import gas_density
        return self.liquid_density.magnitude / gas_density(self.composition, reference, backend=backend).magnitude


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LNGInventory(_safety_QuantitySnapshot):
    cargo: LNGCargo
    liquid_volume: pint.Quantity

    def __post_init__(self) -> None:
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.registry import ureg
        if not isinstance(self.cargo, LNGCargo):
            raise TypeError("Inventory requires an LNGCargo specification.")
        volume = _converters_petroleum_quantity_number(self.liquid_volume, "meter**3", "liquid inventory", minimum=0)
        object.__setattr__(self, "liquid_volume", ureg.Quantity(volume, "meter**3"))


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LNGTransferAdjustment(_safety_QuantitySnapshot):
    name: str
    energy: pint.Quantity
    basis: str
    direction: str
    source: str
    combustion_temperature: pint.Quantity

    def __post_init__(self) -> None:
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.heat_content import _converters_energy_normalize_basis
        from core.units.registry import ureg
        object.__setattr__(self, "name", _converters_petroleum_text(self.name, "adjustment name"))
        object.__setattr__(self, "source", _converters_petroleum_text(self.source, "adjustment source"))
        object.__setattr__(self, "basis", _converters_energy_normalize_basis(self.basis))
        object.__setattr__(self, "combustion_temperature", _converters_petroleum_temperature(self.combustion_temperature))
        if self.direction not in {"credit", "debit"}:
            raise ValueError("An adjustment requires explicit credit/debit direction.")
        energy = _converters_petroleum_quantity_number(self.energy, "MJ", "adjustment energy", minimum=0)
        object.__setattr__(self, "energy", ureg.Quantity(energy, "MJ"))


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LNGTransferResult(_safety_QuantitySnapshot):
    opening_energy: pint.Quantity
    closing_energy: pint.Quantity
    gross_transferred_energy: pint.Quantity
    net_transferred_energy: pint.Quantity
    direction: str
    basis: str
    adjustments: tuple[LNGTransferAdjustment, ...]
    opening_inventory: LNGInventory
    closing_inventory: LNGInventory
    opening_density: pint.Quantity
    closing_density: pint.Quantity





def calculate_lng_transfer(opening: LNGInventory, closing: LNGInventory, *, direction: str, basis: str,
                            adjustments: tuple[LNGTransferAdjustment, ...], target_unit: str = "MMBtu_it") -> LNGTransferResult:
    """Survey energy difference with explicitly signed contractual adjustments.

    Supply adjustments=() to declare none. Vapor return, BOG, fuel, heel,
    flare and other adjustments are never assumed, double-counted or netted
    automatically. Each survey carries its own density/composition/temperature.
    """
    from core.units.heat_content import _converters_energy_normalize_basis
    from core.units.registry import ureg
    basis = _converters_energy_normalize_basis(basis)
    if direction not in {"loading", "unloading"}:
        raise ValueError("Transfer direction must be loading or unloading.")
    if not isinstance(opening, LNGInventory) or not isinstance(closing, LNGInventory):
        raise TypeError("Supply opening and closing LNG inventories.")
    if opening.cargo.calorific_values.combustion_temperature != closing.cargo.calorific_values.combustion_temperature:
        raise ValueError("Survey combustion reference temperatures differ.")
    adjustments = tuple(adjustments)
    if not all(isinstance(item, LNGTransferAdjustment) for item in adjustments):
        raise TypeError("Adjustments must be LNGTransferAdjustment records.")
    if len({item.name for item in adjustments}) != len(adjustments):
        raise ValueError("Duplicate LNG transfer adjustment names.")
    first_density, last_density = opening.cargo.liquid_density, closing.cargo.liquid_density

    def energy(inventory, density):
        from core.units.registry import ureg
        values = inventory.cargo.calorific_values
        specific = values.hhv if basis == "HHV" else values.lhv
        return ureg.Quantity(inventory.liquid_volume.magnitude * density.magnitude * specific.magnitude, "MJ")

    first, last = energy(opening, first_density), energy(closing, last_density)
    gross = last.magnitude - first.magnitude if direction == "loading" else first.magnitude - last.magnitude
    if gross < 0:
        raise ValueError("Survey energy difference conflicts with loading/unloading direction.")
    net = gross
    for adjustment in adjustments:
        if adjustment.basis != basis:
            raise ValueError("Adjustment HHV/LHV basis conflicts with cargo energy basis.")
        if adjustment.combustion_temperature != opening.cargo.calorific_values.combustion_temperature:
            raise ValueError("Adjustment combustion reference temperature conflicts with cargo data.")
        net += adjustment.energy.magnitude if adjustment.direction == "credit" else -adjustment.energy.magnitude
    if net < 0:
        raise ValueError("LNG adjustments produce negative delivered energy.")
    return LNGTransferResult(first.to(target_unit), last.to(target_unit), ureg.Quantity(gross, "MJ").to(target_unit),
                             ureg.Quantity(net, "MJ").to(target_unit), direction, basis, adjustments,
                             opening, closing, first_density, last_density)


__all__ = ('GasAssay', 'GasReferenceConditions', 'LNGCargo', 'LNGInventory', 'LNGTransferAdjustment', 'LNGTransferResult', 'calculate_lng_transfer', 'composition_gas_energy')
_PUBLIC_ALIASES = {'GasAssay': ('assays.energy.natural_gas', 'GasAssay'), 'GasReferenceConditions': ('reference_conditions', 'GasReferenceConditions')}

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
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.composition import GasCalorificValues
    from core.units.composition import GasComposition
    from core.units.thermodynamics import GasPropertyBackend
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.assays.energy.petroleum import LiquidState
