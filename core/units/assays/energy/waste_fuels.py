"""Commodex assays energy waste_fuels implementation and public interfaces."""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import pint

from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GasUpgradeResult(_safety_QuantitySnapshot):
    product_composition: GasComposition
    product_volume: pint.Quantity
    feed_mass: pint.Quantity
    product_mass: pint.Quantity
    offgas_mass: pint.Quantity
    removed_water_mass: pint.Quantity


from core.units.assays.industrial.industrial_gases import IndustrialGasAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class BiogasAssay(IndustrialGasAssay):
    def __post_init__(self):
        from core.units.assays.industrial.industrial_gases import IndustrialGasAssay
        IndustrialGasAssay.__post_init__(self)
        if self.gas_composition.mole_fractions.get("CH4", 0) <= 0:
            raise ValueError("Biogas requires positive methane content.")

    def upgrade(self, quantity, *, methane_recovery_fraction, target_methane_fraction,
                target_reference, target_backend, target_unit="m^3"):
        """Dry product with retained impurity proportions, conserving components.

        Recovery is methane moles retained. Nonmethane components retain their
        relative dry proportions; this is an explicit separation assumption,
        not a prediction of an upgrading plant or chemical reaction.
        """
        from core.units.assays.energy.lng import _converters_gas_gas_volume_output
        from core.units.assays.energy.lng import _converters_gas_geometric_volume
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.composition import GasComposition
        from core.units.thermodynamics import gas_density
        recovery = _converters_petroleum_fraction(methane_recovery_fraction, "methane recovery")
        purity = _converters_petroleum_fraction(target_methane_fraction, "product methane purity")
        dry = self.gas_composition.dry()
        source_purity = dry.mole_fractions.get("CH4", Decimal(0))
        if purity < source_purity or purity == 0:
            raise ValueError("Upgrading target cannot reduce dry methane purity.")
        values = {"CH4": purity}
        if purity < 1:
            if source_purity == 1:
                raise ValueError("No impurities exist to supply the requested blend.")
            values.update({k: v * (1 - purity) / (1 - source_purity) for k, v in dry.mole_fractions.items() if k != "CH4"})
        product = GasComposition(values, normalize=True)
        volume_m3 = _converters_gas_geometric_volume(quantity, self.reference)
        feed_mass = volume_m3 * self.density.magnitude
        feed_moles = feed_mass / self.gas_composition.molar_mass.to("kg/mol").magnitude
        water_moles = feed_moles * self.gas_composition.water_mole_fraction
        dry_moles = feed_moles * (1 - self.gas_composition.water_mole_fraction)
        product_moles = dry_moles * source_purity * recovery / purity
        product_mass = product_moles * product.molar_mass.to("kg/mol").magnitude
        water_mass = water_moles * GasComposition({"H2O": 1}).molar_mass.to("kg/mol").magnitude
        offgas = feed_mass - product_mass - water_mass
        if offgas < -max(feed_mass, Decimal(1)) * Decimal("1e-20"):
            raise ArithmeticError("Upgrading mass balance produced negative offgas.")
        density = gas_density(product, target_reference, backend=target_backend)
        product_volume = _converters_gas_gas_volume_output(product_mass / density.magnitude, target_unit, target_reference, quantity)
        return GasUpgradeResult(product, product_volume, _converters_petroleum_output(feed_mass, "kg", quantity),
                                _converters_petroleum_output(product_mass, "kg", quantity), _converters_petroleum_output(max(offgas, Decimal(0)), "kg", quantity),
                                _converters_petroleum_output(water_mass, "kg", quantity))

    def upgraded_volume(self, quantity, **conditions):
        return self.upgrade(quantity, **conditions).product_volume





@dataclass(frozen=True, slots=True, kw_only=True)
class BiomethaneAssay(BiogasAssay):
    def injection_quality(self, *, minimum_methane, maximum_co2, maximum_h2s, maximum_water):
        return self.gas_composition.quality(minimum_methane=minimum_methane, maximum_co2=maximum_co2,
                                            maximum_h2s=maximum_h2s, maximum_water=maximum_water)


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class WasteDerivedFuelAssay(MaterialSpecification):
    inert_dry_fraction: Decimal
    fuel: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        if self.fuel not in {"RDF", "SRF", "municipal-waste", "other"}:
            raise ValueError("Declare the waste-derived fuel category.")
        object.__setattr__(self, "inert_dry_fraction", _converters_petroleum_fraction(self.inert_dry_fraction, "dry inert fraction"))

    def combustible_dry_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(self.wet_to_dry_mass(quantity), 1 - self.inert_dry_fraction)

    def fuel_energy(self, quantity, *, basis="lhv", target_unit="MWh"):
        from core.units.heat_content import convert_heat_content
        return convert_heat_content(quantity, target_unit, assay=self.as_heat_content(), basis=basis)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (BiogasAssay, BiomethaneAssay, WasteDerivedFuelAssay,))
__all__ = ('BiogasAssay', 'BiomethaneAssay', 'GasUpgradeResult', 'WasteDerivedFuelAssay', 'convert_biogas', 'convert_biomethane', 'convert_waste_derived_fuel', 'create_biogas_context', 'create_biomethane_context', 'create_waste_derived_fuel_context')
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
    from core.units.assays.energy.waste_fuels import convert_biogas
    from core.units.assays.energy.waste_fuels import convert_biomethane
    from core.units.assays.energy.waste_fuels import convert_waste_derived_fuel
    from core.units.assays.energy.waste_fuels import create_biogas_context
    from core.units.assays.energy.waste_fuels import create_biomethane_context
    from core.units.assays.energy.waste_fuels import create_waste_derived_fuel_context
