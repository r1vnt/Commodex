"""Commodex processing implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, field
from decimal import Decimal
from enum import Enum
from itertools import permutations
from types import MappingProxyType
from typing import Any, Final, Literal
import math
import pint




def _energy_environment_specs_efficiency(value, label):
    from core.units.assays.energy.petroleum import _converters_petroleum_fraction
    from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
    if isinstance(value, pint.Quantity):
        value = _converters_petroleum_quantity_number(value, "dimensionless", label, minimum=0)
    fraction = _converters_petroleum_fraction(value, label)
    if fraction == 0:
        raise ValueError(f"{label} must be positive.")
    return fraction


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True, kw_only=True)
class EnergyStorageSpecification(_safety_QuantitySnapshot):
    name: str
    charge_efficiency: Decimal
    discharge_efficiency: Decimal
    standing_loss_fraction_per_hour: Decimal | None = None
    capacity: pint.Quantity | None = None

    def __post_init__(self):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.registry import ureg
        object.__setattr__(self, "name", _converters_petroleum_text(self.name, "storage name"))
        for name in ("charge_efficiency", "discharge_efficiency"):
            object.__setattr__(self, name, _energy_environment_specs_efficiency(getattr(self, name), name))
        if self.standing_loss_fraction_per_hour is not None:
            object.__setattr__(self, "standing_loss_fraction_per_hour", _converters_petroleum_fraction(self.standing_loss_fraction_per_hour, "hourly standing loss", allow_one=False))
        if self.capacity is not None:
            object.__setattr__(self, "capacity", ureg.Quantity(_converters_petroleum_positive_quantity(self.capacity, "MWh", "stored energy capacity"), "MWh"))

    def _retention(self, duration):
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if duration is None:
            if self.standing_loss_fraction_per_hour not in {None, 0}:
                raise ValueError("Declare storage duration when standing losses are specified.")
            return Decimal(1)
        hours = _converters_petroleum_quantity_number(duration, "hour", "storage duration", minimum=0)
        if self.standing_loss_fraction_per_hour is None and hours != 0:
            raise ValueError("Declare standing loss, including explicit zero, for a holding period.")
        return (1 - (self.standing_loss_fraction_per_hour or 0)) ** hours

    def stored_energy(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        input_energy = _converters_petroleum_quantity_number(quantity, "MWh", "charge energy", minimum=0)
        stored = input_energy * self.charge_efficiency
        if self.capacity is not None and stored > self.capacity.magnitude:
            raise ValueError("Charge exceeds declared stored-energy capacity.")
        return _converters_petroleum_output(stored, "MWh", quantity).to(quantity.units)

    def output_energy(self, quantity, *, duration=None):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        stored = _converters_petroleum_quantity_number(self.stored_energy(quantity), "MWh", "stored energy", minimum=0)
        return _converters_petroleum_output(stored * self._retention(duration) * self.discharge_efficiency, "MWh", quantity).to(quantity.units)

    def input_required(self, quantity, *, duration=None):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        output = _converters_petroleum_quantity_number(quantity, "MWh", "required output", minimum=0)
        factor = self.charge_efficiency * self.discharge_efficiency * self._retention(duration)
        result = _converters_petroleum_output(output / factor, "MWh", quantity).to(quantity.units)
        self.stored_energy(result)  # Capacity validation for the inverse route.
        return result

    def as_context(self, *, name=None):
        # Energy -> energy losses cannot be a Pint dimension-changing edge.
        return pint.Context(name or self.name)

    def convert_quantity(self, quantity, target_unit, *, basis=None):
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.registry import parse_unit
        _converters_petroleum_quantity_number(quantity, "MWh", "storage energy", minimum=0)
        return quantity.to(parse_unit(target_unit))





def _energy_systems_fraction(value, label, *, positive=False):
    from core.units.heat_content import _calorific_models_dec
    from core.units.heat_content import _calorific_models_q
    if isinstance(value, pint.Quantity):
        value = _calorific_models_q(value, "dimensionless", label).magnitude
    result = _calorific_models_dec(value, label)
    if not 0 <= result <= 1 or positive and result == 0:
        raise ValueError(f"{label} must be {'positive and ' if positive else ''}within [0, 1].")
    return result





def _energy_systems_energy(value, *, role=None):
    from core.units.dimensions import physical_quantity
    from core.units.dimensions import roles
    from core.units.heat_content import _calorific_models_dec
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(value, pint.Quantity):
        raise TypeError("Energy or power requires a Pint quantity.")
    validate_registry_quantity(value)
    signature = roles(value.units)
    role = role or (next(iter(signature)) if len(signature) == 1 else "thermal")
    value = physical_quantity(value, role=role)
    number = _calorific_models_dec(value.magnitude, "energy/power")
    if number < 0 or value.dimensionality not in {ureg.J.dimensionality, ureg.W.dimensionality}:
        raise ValueError("Supply finite nonnegative energy or power.")
    return ureg.Quantity(number, value.units)





def _energy_systems_output(value, source, target=None, *, role=None):
    """Keep caller float/Decimal arithmetic after Decimal-safe unit conversion."""
    from core.units.dimensions import restore
    from core.units.dimensions import roles
    from core.units.dimensions import typed_unit
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    signature = roles(source.units)
    if role is None and signature and value.dimensionality in {ureg.J.dimensionality, ureg.W.dimensionality}:
        role = next(iter(signature))
        target = target or source.units
    if role is not None:
        value = restore(value, target or typed_unit(value, role), role=role)
    elif target is not None:
        value = value.to(parse_unit(target))
    return ureg.Quantity(value.magnitude if isinstance(source.magnitude, Decimal) else float(value.magnitude), value.units)





def _energy_systems_kind(value):
    from core.units.heat_content import _kind
    if value.lower().strip() in {"gar", "nar"}:
        raise ValueError("Use HHV/LHV for system energy basis; declare analytical basis on the assay.")
    return _kind(value)





def convert_thermal_efficiency(efficiency, *, assay, source_basis, target_basis):
    """Same output and fuel amount: eta_target = eta_source * H_source/H_target.

    LHV efficiencies can exceed 1 for condensing equipment. This conversion
    changes the fuel-energy denominator; it does not change gross/net output.
    """
    from core.units.heat_content import _calorific_models_dec
    from core.units.heat_content import _calorific_models_q
    from core.units.heat_content import heat_content_from_assay
    source_basis, target_basis = _energy_systems_kind(source_basis), _energy_systems_kind(target_basis)
    if isinstance(efficiency, pint.Quantity):
        efficiency = _calorific_models_q(efficiency, "dimensionless", "efficiency").magnitude
    eta = _calorific_models_dec(efficiency, "efficiency")
    if eta <= 0:
        raise ValueError("Efficiency must be positive.")
    specification = heat_content_from_assay(assay)
    if specification.energy_role != "combustion":
        raise ValueError("Thermal efficiency requires fuel combustion heat, not recovered-carrier hydrogen energy.")
    gross_heat = specification.calorific_value("hhv")
    net_heat = specification.energy_density(gross_heat.units, kind="lhv")
    heat = {"gcv": gross_heat.magnitude, "ncv": net_heat.magnitude}
    source, target = heat[source_basis], heat[target_basis]
    # Pure thermal fuel efficiency must be physically admissible on HHV basis.
    gross = gross_heat.magnitude
    if eta * source > gross + Decimal("1e-12") * gross:
        raise ValueError("Efficiency exceeds the fuel's HHV energy input.")
    return eta * source / target





def _energy_systems_parasitic(gross, *, parasitic_load=None, parasitic_fraction=None):
    if (parasitic_load is None) == (parasitic_fraction is None):
        raise ValueError("Select exactly one parasitic_load or parasitic_fraction, including explicit zero.")
    if parasitic_load is not None:
        load = _energy_systems_energy(parasitic_load, role="electrical").to(gross.units)
    else:
        load = gross * _energy_systems_fraction(parasitic_fraction, "parasitic fraction of gross")
    if load > gross:
        raise ValueError("Parasitic load exceeds gross generation.")
    return load





def net_generation(gross, *, parasitic_load=None, parasitic_fraction=None, target_unit=None):
    """Net output = gross output - explicit plant auxiliary consumption."""
    value = _energy_systems_energy(gross, role="electrical")
    result = value - _energy_systems_parasitic(value, parasitic_load=parasitic_load, parasitic_fraction=parasitic_fraction)
    return _energy_systems_output(result, gross, target_unit, role="electrical")





def gross_generation(net, *, parasitic_load=None, parasitic_fraction=None, target_unit=None):
    value = _energy_systems_energy(net, role="electrical")
    if (parasitic_load is None) == (parasitic_fraction is None):
        raise ValueError("Select exactly one parasitic_load or parasitic_fraction.")
    if parasitic_load is not None:
        result = value + _energy_systems_energy(parasitic_load, role="electrical").to(value.units)
    else:
        fraction = _energy_systems_fraction(parasitic_fraction, "parasitic fraction of gross")
        if fraction == 1:
            raise ValueError("Cannot invert 100% parasitic consumption.")
        result = value / (1 - fraction)
    return _energy_systems_output(result, net, target_unit, role="electrical")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class CHPBalance(_safety_QuantitySnapshot):
    fuel_input: pint.Quantity
    gross_electricity: pint.Quantity
    parasitic_load: pint.Quantity
    net_electricity: pint.Quantity
    useful_heat: pint.Quantity
    rejected_heat: pint.Quantity
    fuel_basis: str

    @property
    def gross_electrical_efficiency(self):
        from core.units.dimensions import physical_quantity
        return (physical_quantity(self.gross_electricity, role="electrical") / physical_quantity(self.fuel_input, role="thermal")).to("").magnitude

    @property
    def net_electrical_efficiency(self):
        from core.units.dimensions import physical_quantity
        return (physical_quantity(self.net_electricity, role="electrical") / physical_quantity(self.fuel_input, role="thermal")).to("").magnitude

    @property
    def total_useful_efficiency(self):
        from core.units.dimensions import physical_quantity
        useful = physical_quantity(self.net_electricity, role="electrical") + physical_quantity(self.useful_heat, role="thermal")
        return (useful / physical_quantity(self.fuel_input, role="thermal")).to("").magnitude

    @property
    def energy_balance_residual(self):
        """First-law accounting inside this explicitly declared CHP process."""
        from core.units.dimensions import physical_quantity
        from core.units.dimensions import restore
        residual = physical_quantity(self.fuel_input, role="thermal")
        for quantity, role in ((self.gross_electricity, "electrical"), (self.useful_heat, "thermal"), (self.rejected_heat, "thermal")):
            residual -= physical_quantity(quantity, role=role)
        return restore(residual, self.fuel_input.units, role="thermal")





@dataclass(frozen=True, slots=True)
class CHPSpecification:
    """Fractions of input thermal energy, with parasitics a fraction of gross electricity.

    Efficiencies sum to at most 1 on the selected basis. Condensing systems whose
    combined LHV efficiency exceeds 1 should declare their HHV efficiencies.
    EPA: https://www.epa.gov/chp/methods-calculating-chp-efficiency
    """
    electrical_efficiency: Decimal | float
    useful_heat_efficiency: Decimal | float
    fuel_basis: str
    parasitic_fraction: Decimal | float = 0

    def __post_init__(self):
        for key in ("electrical_efficiency", "useful_heat_efficiency", "parasitic_fraction"):
            object.__setattr__(self, key, _energy_systems_fraction(getattr(self, key), key))
        if self.electrical_efficiency + self.useful_heat_efficiency > 1:
            raise ValueError("Gross electrical and useful-heat efficiencies cannot exceed input energy.")
        object.__setattr__(self, "fuel_basis", _energy_systems_kind(self.fuel_basis))

    def outputs(self, fuel_input, *, fuel_basis):
        if _energy_systems_kind(fuel_basis) != self.fuel_basis:
            raise ValueError("Fuel input HHV/LHV basis must match the CHP specification.")
        energy = _energy_systems_energy(fuel_input, role="thermal")
        if energy.magnitude == 0:
            raise ValueError("CHP efficiency balance requires positive fuel input.")
        gross = energy * self.electrical_efficiency
        load = gross * self.parasitic_fraction
        heat = energy * self.useful_heat_efficiency
        rejected = energy - gross - heat
        values = (energy, gross, load, gross-load, heat, rejected)
        roles = ("thermal", "electrical", "electrical", "electrical", "thermal", "thermal")
        return CHPBalance(*(_energy_systems_output(value, fuel_input, role=role) for value, role in zip(values, roles)), self.fuel_basis)





@dataclass(frozen=True, slots=True)
class HeatPumpSpecification:
    coefficient_of_performance: Decimal | float

    def __post_init__(self):
        from core.units.heat_content import _calorific_models_dec
        value = _calorific_models_dec(self.coefficient_of_performance, "COP")
        if value <= 0:
            raise ValueError("COP must be positive; it is not a combustion efficiency.")
        object.__setattr__(self, "coefficient_of_performance", value)

    def useful_heat(self, electrical_input, *, target_unit=None):
        return _energy_systems_output(_energy_systems_energy(electrical_input, role="electrical") * self.coefficient_of_performance, electrical_input, target_unit, role="thermal")

    def electrical_input(self, useful_heat, *, target_unit=None):
        return _energy_systems_output(_energy_systems_energy(useful_heat, role="thermal") / self.coefficient_of_performance, useful_heat, target_unit, role="electrical")


from core.units.validation import _safety_QuantitySnapshot
from core.units.validation import authorized_conversion


@dataclass(frozen=True, slots=True)
class ElectrolysisSpecification(_safety_QuantitySnapshot):
    """Declared system electricity consumption, including selected balance-of-plant scope.

    Specific electricity can also be built from explicitly chosen H2 HHV/LHV
    efficiency. Additional thermal input is separate for heat-assisted systems.
    https://www.energy.gov/cmei/fuels/technical-targets-proton-exchange-membrane-electrolysis
    """
    electricity_per_kg_h2: pint.Quantity
    source: str
    additional_heat_per_kg_h2: pint.Quantity | None = None

    def __post_init__(self):
        from core.units.dimensions import physical_quantity
        from core.units.dimensions import restore
        from core.units.heat_content import _calorific_models_q
        from core.units.heat_content import _calorific_models_text
        object.__setattr__(self, "source", _calorific_models_text(self.source, "electrolysis source/system boundary"))
        electricity = _calorific_models_q(physical_quantity(self.electricity_per_kg_h2, role="electrical"), "kWh/kg", "electrolysis electricity")
        if electricity.magnitude <= 0:
            raise ValueError("Electrolysis electricity must be positive.")
        object.__setattr__(self, "electricity_per_kg_h2", restore(electricity, "kWh_e/kg", role="electrical"))
        if self.additional_heat_per_kg_h2 is not None:
            heat = _calorific_models_q(physical_quantity(self.additional_heat_per_kg_h2, role="thermal"), "kWh/kg", "additional heat", nonnegative=True)
            object.__setattr__(self, "additional_heat_per_kg_h2", restore(heat, "kWh_th/kg", role="thermal"))

    @classmethod
    def from_efficiency(cls, *, efficiency, hydrogen_heat_content, basis, source):
        cls._validate_hydrogen_heat(hydrogen_heat_content)
        eta = _energy_systems_fraction(efficiency, "electrolysis efficiency", positive=True)
        heat = hydrogen_heat_content.energy_density("kWh/kg", kind=_energy_systems_kind(basis))
        return cls(heat / eta, source)

    @authorized_conversion
    def hydrogen_output(self, electrical_input, *, target_unit="kg"):
        from core.units.dimensions import physical_quantity
        result = _energy_systems_energy(electrical_input, role="electrical") / physical_quantity(self.electricity_per_kg_h2, role="electrical")
        return _energy_systems_output(result, electrical_input, target_unit)

    @authorized_conversion
    def electricity_required(self, hydrogen_output, *, target_unit="MWh_e"):
        from core.units.dimensions import physical_quantity
        from core.units.heat_content import _calorific_models_dec
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        target = parse_unit(target_unit)
        if not isinstance(hydrogen_output, pint.Quantity):
            raise TypeError("Hydrogen output requires a mass or mass flow quantity.")
        validate_registry_quantity(hydrogen_output)
        if _calorific_models_dec(hydrogen_output.magnitude, "hydrogen output") < 0 or not any(hydrogen_output.is_compatible_with(unit) for unit in ("kg", "kg/hour")):
            raise ValueError("Hydrogen output must be nonnegative mass or mass flow.")
        value = ureg.Quantity(_calorific_models_dec(hydrogen_output.magnitude, "hydrogen output"), hydrogen_output.units)
        return _energy_systems_output(value * physical_quantity(self.electricity_per_kg_h2, role="electrical"), hydrogen_output, target, role="electrical")

    def efficiency(self, *, hydrogen_heat_content, basis, include_heat_input=False):
        from core.units.dimensions import physical_quantity
        self._validate_hydrogen_heat(hydrogen_heat_content)
        heat = hydrogen_heat_content.energy_density("kWh/kg", kind=_energy_systems_kind(basis))
        input_energy = physical_quantity(self.electricity_per_kg_h2, role="electrical")
        if include_heat_input:
            if self.additional_heat_per_kg_h2 is None:
                raise ValueError("Declare additional heat, including zero, for total-input efficiency.")
            input_energy += physical_quantity(self.additional_heat_per_kg_h2, role="thermal")
        return (heat / input_energy).to("").magnitude

    @staticmethod
    def _validate_hydrogen_heat(specification):
        from core.units.heat_content import HeatContentAssay
        if not isinstance(specification, HeatContentAssay) or specification.energy_role != "combustion":
            raise TypeError("Supply the pure-H2 combustion HeatContentAssay.")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class HydrogenCrackingBalance(_safety_QuantitySnapshot):
    carrier_mass: pint.Quantity
    delivered_hydrogen_mass: pint.Quantity
    delivered_hydrogen_energy: pint.Quantity
    carrier_combustion_energy: pint.Quantity | None
    cracking_energy_input: pint.Quantity
    heating_basis: str


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class HydrogenCrackingSpecification(_safety_QuantitySnapshot):
    hydrogen_mass_fraction: Decimal | float
    recovery_fraction: Decimal | float
    hydrogen_heat_content: object
    cracking_energy_per_kg_carrier: pint.Quantity
    carrier_heat_content: object | None = None

    def __post_init__(self):
        from core.units.dimensions import physical_quantity
        from core.units.dimensions import restore
        from core.units.heat_content import _calorific_models_q
        from core.units.heat_content import HeatContentAssay
        for key in ("hydrogen_mass_fraction", "recovery_fraction"):
            object.__setattr__(self, key, _energy_systems_fraction(getattr(self, key), key))
        if not isinstance(self.hydrogen_heat_content, HeatContentAssay):
            raise TypeError("Supply the pure-H2 heat-content specification.")
        if self.carrier_heat_content is not None and not isinstance(self.carrier_heat_content, HeatContentAssay):
            raise TypeError("Carrier combustion heat must be a separate HeatContentAssay.")
        if self.hydrogen_heat_content.energy_role != "combustion" or self.carrier_heat_content is not None and self.carrier_heat_content.energy_role != "combustion":
            raise ValueError("Supply actual H2/carrier combustion values, not already recovered-hydrogen equivalents.")
        energy = _calorific_models_q(physical_quantity(self.cracking_energy_per_kg_carrier, role="thermal"), "MJ/kg", "cracking energy", nonnegative=True)
        object.__setattr__(self, "cracking_energy_per_kg_carrier", restore(energy, "MJ_th/kg", role="thermal"))

    def balance(self, carrier_mass, *, basis):
        from core.units.dimensions import physical_quantity
        from core.units.heat_content import _calorific_models_q
        mass = _calorific_models_q(carrier_mass, "kg", "carrier mass", nonnegative=True)
        delivered = mass * self.hydrogen_mass_fraction * self.recovery_fraction
        hydrogen_energy = delivered * self.hydrogen_heat_content.energy_density("MJ/kg", kind=_energy_systems_kind(basis))
        combustion = None if self.carrier_heat_content is None else mass * self.carrier_heat_content.energy_density("MJ/kg", kind=_energy_systems_kind(basis))
        demand = mass * physical_quantity(self.cracking_energy_per_kg_carrier, role="thermal")
        return HydrogenCrackingBalance(_energy_systems_output(mass, carrier_mass), _energy_systems_output(delivered, carrier_mass), _energy_systems_output(hydrogen_energy, carrier_mass, "MJ_th", role="thermal"),
            None if combustion is None else _energy_systems_output(combustion, carrier_mass, "MJ_th", role="thermal"), _energy_systems_output(demand, carrier_mass, "MJ_th", role="thermal"), _energy_systems_kind(basis))





@dataclass(frozen=True, slots=True)
class EnergyLoss:
    name: str
    fraction: Decimal | float

    def __post_init__(self):
        from core.units.heat_content import _calorific_models_text
        object.__setattr__(self, "name", _calorific_models_text(self.name, "loss stage"))
        object.__setattr__(self, "fraction", _energy_systems_fraction(self.fraction, "energy loss"))


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class EnergyLossBalance(_safety_QuantitySnapshot):
    input_energy: pint.Quantity
    delivered_energy: pint.Quantity
    stage_losses: tuple[tuple[str, pint.Quantity], ...]





def _energy_systems_stages(stages):
    result = tuple(stages)
    if not all(isinstance(stage, EnergyLoss) for stage in result):
        raise TypeError("Loss stages must be explicit EnergyLoss records.")
    if len({stage.name for stage in result}) != len(result):
        raise ValueError("Loss stages must have distinct names.")
    return result





def apply_energy_losses(energy, stages):
    current = _energy_systems_energy(energy)
    losses = []
    for stage in _energy_systems_stages(stages):
        loss = current * stage.fraction
        losses.append((stage.name, _energy_systems_output(loss, energy)))
        current -= loss
    return EnergyLossBalance(_energy_systems_output(_energy_systems_energy(energy), energy), _energy_systems_output(current, energy), tuple(losses))





def energy_before_losses(delivered, stages):
    retention = Decimal(1)
    for stage in _energy_systems_stages(stages):
        retention *= 1 - stage.fraction
    if retention == 0:
        raise ValueError("Cannot invert complete energy loss.")
    return _energy_systems_output(_energy_systems_energy(delivered) / retention, delivered)





def blend_heat_content(assays, weights, *, weighting, name="fuel blend", density=None, gas_state=None):
    """Mass- or mole-weighted heat at compatible analytical/combustion bases.

    Weights are either fractions totalling one or positive amount quantities
    all expressed in the selected mass/mole dimension. Gas mixtures require an
    explicitly calculated blend GasMoistureState; Z is never linearly blended.
    Mixing/reaction enthalpy is excluded. Use GasComposition.blend to recompute
    a composition and its ISO/component-derived calorific values instead.
    """
    from core.units.heat_content import _calorific_models_q
    from core.units.heat_content import _calorific_models_text
    from core.units.registry import ureg
    from core.units.heat_content import HeatContentAssay
    from core.units.heat_content import GasMoistureState
    from core.units.heat_content import heat_content_from_assay
    from core.units.heat_content import CalorificMethod
    specs, weights = tuple(heat_content_from_assay(item) for item in assays), tuple(weights)
    if not specs or len(specs) != len(weights) or weighting not in {"mass", "molar"}:
        raise ValueError("Supply matching nonempty assays/weights and mass or molar weighting.")
    first = specs[0]
    from core.units.heat_content import ThermodynamicBasis
    if first.thermodynamic_basis == ThermodynamicBasis.UNSPECIFIED or first.combustion_temperature is None:
        raise ValueError("Fuel blending requires explicit thermodynamic basis and combustion temperature.")
    if first.basis in {"DAF", "DMMF"}:
        raise ValueError("Normalize excluded-constituent analyses to a whole-fuel basis before blending.")
    if any(spec.energy_role != "combustion" for spec in specs):
        raise ValueError("Fuel blending requires combustion values, not recovered-hydrogen equivalents.")
    for spec in specs[1:]:
        if spec.basis != first.basis or spec.thermodynamic_basis != first.thermodynamic_basis or spec.combustion_temperature != first.combustion_temperature:
            raise ValueError("Blend components must share analytical and combustion reference bases.")
    if any(isinstance(weight, pint.Quantity) for weight in weights):
        if not all(isinstance(weight, pint.Quantity) for weight in weights):
            raise ValueError("Do not mix fractional weights with amount quantities.")
        numbers = [_calorific_models_q(weight, "kg" if weighting == "mass" else "mol", "blend weight", nonnegative=True).magnitude for weight in weights]
        total = sum(numbers, Decimal(0))
        if total <= 0:
            raise ValueError("Blend weights must have positive total.")
        fractions = [weight / total for weight in numbers]
    else:
        fractions = [_energy_systems_fraction(weight, "blend fraction") for weight in weights]
        if abs(sum(fractions, Decimal(0)) - 1) > Decimal("1e-12"):
            raise ValueError("Fractional blend weights must total one; normalization is explicit.")
    values = {}
    unit = "MJ/kg" if weighting == "mass" else "MJ/mol"
    for kind in ("gcv", "ncv"):
        # Never derive a missing counterpart merely because blending is requested.
        if all(getattr(spec, kind) is not None for spec in specs):
            values[kind] = sum((spec.energy_density(unit, kind=kind) * weight for spec, weight in zip(specs, fractions)), ureg.Quantity(Decimal(0), unit))
    if not values:
        raise ValueError("Blend components have no common supplied calorific value.")
    metadata = dict(basis=first.basis, thermodynamic_basis=first.thermodynamic_basis,
                    combustion_temperature=first.combustion_temperature, density=density,
                    method=CalorificMethod.COMPONENT_SUM, method_source=f"{weighting}-weighted blend: " + "; ".join(spec.name for spec in specs))
    molar_mass = None
    if all(spec.molar_mass is not None for spec in specs):
        if weighting == "molar":
            molar_mass = sum((spec.molar_mass * weight for spec, weight in zip(specs, fractions)), ureg.Quantity(Decimal(0), "kg/mol"))
        else:
            molar_mass = ureg.Quantity(1 / sum((weight / spec.molar_mass.magnitude for spec, weight in zip(specs, fractions)), Decimal(0)), "kg/mol")
        metadata["molar_mass"] = molar_mass
    if first.basis in {"DRY_GAS", "WET_GAS"}:
        if not isinstance(gas_state, GasMoistureState) or weighting == "mass" and molar_mass is None:
            raise ValueError("Gas blending requires molar masses and an explicit blend GasMoistureState.")
        if weighting == "molar":
            mole_weights = fractions
        else:
            moles = [weight / spec.molar_mass.magnitude for spec, weight in zip(specs, fractions)]
            total = sum(moles, Decimal(0))
            mole_weights = [value / total for value in moles]
        water = sum((weight * spec._gas_state(spec.basis).water_mole_fraction for spec, weight in zip(specs, mole_weights)), Decimal(0))
        if abs(water - gas_state.water_mole_fraction) > Decimal("1e-12"):
            raise ValueError("Declared blend-state water fraction disagrees with weighted gas composition.")
        metadata["dry_gas_state" if first.basis == "DRY_GAS" else "wet_gas_state"] = gas_state
    else:
        if weighting == "mass":
            mass_weights = fractions
        else:
            if molar_mass is None:
                if first.basis == "DB":
                    return HeatContentAssay(name=_calorific_models_text(name, "blend name"), **values, **metadata)
                raise ValueError("Molar blending of non-gas fuel requires each molar mass to preserve analytical moisture.")
            mass_weights = [weight * spec.molar_mass.magnitude / molar_mass.magnitude for spec, weight in zip(specs, fractions)]
        # Measured values can be blended on a declared common AR/AD basis
        # without knowing sample moisture. Leave those metadata unknown.
        if first.basis == "AR" and any(spec.moisture_pct is None and spec.as_received_moisture_pct is None for spec in specs):
            return HeatContentAssay(name=_calorific_models_text(name, "blend name"), **values, **metadata)
        if first.basis == "AD" and any(spec.moisture_pct is None and spec.air_dried_moisture_pct is None for spec in specs):
            return HeatContentAssay(name=_calorific_models_text(name, "blend name"), **values, **metadata)
        moisture = sum((weight * spec._moisture(spec.basis) for spec, weight in zip(specs, mass_weights)), Decimal(0))
        metadata["moisture_pct"] = moisture * 100
        for field in ("ash_db_pct", "mineral_matter_db_pct"):
            if all(getattr(spec, field) is not None for spec in specs):
                metadata[field] = sum((weight * (1 - spec._moisture(spec.basis)) * getattr(spec, field) for spec, weight in zip(specs, mass_weights)), Decimal(0)) / (1 - moisture)
        # Named fuel conditions cannot be reconstructed from moisture alone.
        if first.basis in {"AS_FIRED", "AS_DELIVERED", "REFERENCE_MOISTURE"}:
            from core.units.heat_content import FuelCondition
            from core.units.heat_content import _CONDITION_FIELDS
            metadata[_CONDITION_FIELDS[first.basis]] = FuelCondition(moisture * 100, density=density)
    return HeatContentAssay(name=_calorific_models_text(name, "blend name"), **values, **metadata)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GenerationProcess(_safety_QuantitySnapshot):
    """Explicit thermal-to-electrical process for amounts, rates and quotes.

    Declare exactly one efficiency or heat rate and the gross/net output scope.
    This process represents generation parity; it supplies no fuel assay, FX,
    carbon allowance conversion, or tolling-contract assumptions.
    """
    efficiency: object = None
    heat_rate: pint.Quantity | None = None
    output_basis: str = "net"
    fuel_basis: str = "lhv"
    name: str = "declared generation process"

    def __post_init__(self):
        from core.units.heat_content import _calorific_models_text
        from core.units.registry import ureg
        from core.units.heat_content import _efficiency
        from core.units.heat_content import efficiency_from_heat_rate
        if (self.efficiency is None) == (self.heat_rate is None):
            raise ValueError("Generation requires exactly one efficiency or heat_rate.")
        eta = _efficiency(self.efficiency) if self.heat_rate is None else efficiency_from_heat_rate(self.heat_rate)
        if self.efficiency is not None:
            object.__setattr__(self, "efficiency", eta)
        else:
            from core.units.dimensions import physical_heat_rate
            rate = physical_heat_rate(self.heat_rate).magnitude
            object.__setattr__(self, "heat_rate", ureg.Quantity(rate, ureg.thermal_energy_kind / ureg.electrical_energy_kind))
        if self.output_basis not in {"gross", "net"}:
            raise ValueError("Declare gross or net generation output basis.")
        object.__setattr__(self, "fuel_basis", _energy_systems_kind(self.fuel_basis))
        object.__setattr__(self, "name", _calorific_models_text(self.name, "generation process name"))

    @property
    def thermal_efficiency(self):
        from core.units.heat_content import _efficiency
        from core.units.heat_content import efficiency_from_heat_rate
        return _efficiency(self.efficiency) if self.heat_rate is None else efficiency_from_heat_rate(self.heat_rate)

    def convert(self, quantity, target_unit):
        from core.units.dimensions import physical_quantity
        from core.units.dimensions import physical_unit
        from core.units.dimensions import restore
        from core.units.dimensions import roles
        from core.units.heat_content import _calorific_models_dec
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        if not isinstance(quantity, pint.Quantity):
            raise TypeError("Generation conversion requires a Pint quantity.")
        validate_registry_quantity(quantity)
        target = parse_unit(target_unit)
        source_signature, target_signature = roles(quantity.units), roles(target)
        if not source_signature and not target_signature:
            raise ValueError("Generation endpoints need an explicit thermal or electrical role.")
        source_role = next(iter(source_signature)) if len(source_signature) == 1 else None
        target_role = next(iter(target_signature)) if len(target_signature) == 1 else None
        if source_role is None:
            source_role = "thermal" if target_role == "electrical" else "electrical"
        if target_role is None:
            target_role = "electrical" if source_role == "thermal" else "thermal"
        physical = physical_quantity(quantity, role=source_role)
        physical_target = physical_unit(target, role=target_role)
        value = _calorific_models_dec(quantity.magnitude, "generation amount/quote")
        eta = self.thermal_efficiency
        from core.units.conversions import _numerator_currency_iso
        source_currency = _numerator_currency_iso(quantity) if any(key.startswith("[currency_") for key in quantity.dimensionality) else None
        target_currency = _numerator_currency_iso(ureg.Quantity(1, target)) if any(key.startswith("[currency_") for key in target.dimensionality) else None
        quote = source_currency is not None
        if quote:
            if source_currency != target_currency:
                raise ValueError("Apply FX explicitly before generation-price conversion.")
            denomination = ureg.Unit(source_currency) / physical.units
            if denomination.dimensionality not in {ureg.J.dimensionality, ureg.W.dimensionality}:
                raise ValueError("Generation quotes require currency per energy or power.")
        elif physical.dimensionality not in {ureg.J.dimensionality, ureg.W.dimensionality}:
            raise ValueError("Generation converts thermal/electrical energy, power or quotes; convert fuel quantities to thermal energy first.")
        factor = Decimal(1)
        if source_role != target_role:
            factor = eta if source_role == "thermal" else 1 / eta
            if quote:
                factor = 1 / factor
        result = ureg.Quantity(value * factor, physical.units).to(physical_target)
        result = restore(result, target, role=target_role)
        return ureg.Quantity(result.magnitude if isinstance(quantity.magnitude, Decimal) else float(result.magnitude), target)


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# 3. ProcessingSpecification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ProcessingSpecification:
    """Immutable feedstock-to-product yield, metallurgical recovery, loss, and margin model."""

    input_material: str
    output_yields: Mapping[str, Decimal] = field(default_factory=dict)
    recovery_rates: Mapping[str, Decimal] = field(default_factory=dict)
    process_losses: Mapping[str, Decimal] = field(default_factory=dict)

    def __post_init__(self) -> None:
        from core.units.composition import _materials_clean_name
        from core.units.composition import _materials_validate_fraction
        in_mat = _materials_clean_name(self.input_material, "input_material")
        if (
            not isinstance(self.output_yields, Mapping)
            or not isinstance(self.recovery_rates, Mapping)
            or not isinstance(self.process_losses, Mapping)
        ):
            raise TypeError("output_yields, recovery_rates, and process_losses must be mappings.")

        clean_yields: dict[str, Decimal] = {}
        for raw_k, raw_v in self.output_yields.items():
            k = _materials_clean_name(raw_k, "output product name")
            if k in clean_yields:
                raise ValueError(f"Duplicate output product key: {k!r}.")
            clean_yields[k] = _materials_validate_fraction(raw_v, f"output_yields[{k!r}]")

        clean_rec: dict[str, Decimal] = {}
        for raw_k, raw_v in self.recovery_rates.items():
            k = _materials_clean_name(raw_k, "recovery_rates key")
            if k in clean_rec:
                raise ValueError(f"Duplicate recovery_rates key: {k!r}.")
            clean_rec[k] = _materials_validate_fraction(raw_v, f"recovery_rates[{k!r}]")

        clean_loss: dict[str, Decimal] = {}
        for raw_k, raw_v in self.process_losses.items():
            k = _materials_clean_name(raw_k, "process_losses key")
            if k in clean_loss:
                raise ValueError(f"Duplicate process_losses key: {k!r}.")
            clean_loss[k] = _materials_validate_fraction(raw_v, f"process_losses[{k!r}]")

        total_accounted = sum(clean_yields.values(), Decimal(0)) + sum(clean_loss.values(), Decimal(0))
        if total_accounted > Decimal("1.000000000001"):
            raise ValueError(
                f"Feedstock mass balance violated for {in_mat!r}: output_yields + process_losses = {total_accounted} > 1.0."
            )

        object.__setattr__(self, "input_material", in_mat)
        object.__setattr__(self, "output_yields", MappingProxyType(clean_yields))
        object.__setattr__(self, "recovery_rates", MappingProxyType(clean_rec))
        object.__setattr__(self, "process_losses", MappingProxyType(clean_loss))

    def __hash__(self) -> int:
        return hash((
            self.input_material,
            tuple(sorted(self.output_yields.items())),
            tuple(sorted(self.recovery_rates.items())),
            tuple(sorted(self.process_losses.items())),
        ))

    @property
    def total_output_yield(self) -> Decimal:
        return sum(self.output_yields.values(), Decimal(0))

    @property
    def total_process_loss(self) -> Decimal:
        return sum(self.process_losses.values(), Decimal(0))

    @property
    def unaccounted_mass_fraction(self) -> Decimal:
        return max(Decimal(0), Decimal(1) - self.total_output_yield - self.total_process_loss)

    def effective_yield(self, output: str, *, apply_recovery: bool = True) -> Decimal:
        """Net mass fraction of ``output`` per unit of ``input_material``."""
        from core.units.composition import _materials_clean_name
        key = _materials_clean_name(output, "output")
        if key not in self.output_yields:
            raise KeyError(f"Output {key!r} is not defined in output_yields for {self.input_material!r}.")
        yld = self.output_yields[key]
        if apply_recovery and key in self.recovery_rates:
            yld = yld * self.recovery_rates[key]
        return yld

    def input_to_output(
        self,
        input_mass: pint.Quantity,
        output: str,
        *,
        apply_recovery: bool = True,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Convert feedstock mass to produced mass of ``output``."""
        from core.units.composition import _materials_kg_to_target
        from core.units.composition import _materials_mass_in_kg
        from core.units.registry import parse_unit
        kg_in, _ = _materials_mass_in_kg(input_mass, "input_mass")
        yld = self.effective_yield(output, apply_recovery=apply_recovery)
        dst = parse_unit(target_unit) if target_unit is not None else input_mass.units
        return _materials_kg_to_target(kg_in * yld, dst, input_mass)

    def output_to_input(
        self,
        output_mass: pint.Quantity,
        output: str,
        *,
        apply_recovery: bool = True,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Convert target output mass back to required feedstock mass."""
        from core.units.composition import _materials_kg_to_target
        from core.units.composition import _materials_mass_in_kg
        from core.units.registry import parse_unit
        kg_out, _ = _materials_mass_in_kg(output_mass, "output_mass")
        yld = self.effective_yield(output, apply_recovery=apply_recovery)
        if yld <= Decimal(0):
            raise ValueError(f"Effective yield for {output!r} is zero; required feedstock mass is undefined.")
        dst = parse_unit(target_unit) if target_unit is not None else output_mass.units
        return _materials_kg_to_target(kg_out / yld, dst, output_mass)

    def all_outputs(
        self,
        input_mass: pint.Quantity,
        *,
        apply_recovery: bool = True,
        target_unit: str | Unit | None = None,
    ) -> dict[str, pint.Quantity]:
        """Compute produced masses for all co-products from ``input_mass``."""
        return {
            name: self.input_to_output(
                input_mass, name, apply_recovery=apply_recovery, target_unit=target_unit
            )
            for name in self.output_yields
        }

    def recovered_component_mass(
        self,
        input_mass: pint.Quantity,
        composition: MaterialComposition,
        component: str,
        *,
        input_basis: CompositionBasis | str | None = None,
        target_unit: str | Unit = "payable_metal_tonne",
    ) -> pint.Quantity:
        """Compute recovered/payable component mass combining ``composition`` grade and ``recovery_rates``."""
        from core.units.composition import _materials_clean_name
        from core.units.composition import _materials_kg_to_target
        from core.units.composition import _materials_to_decimal
        from core.units.registry import parse_unit
        key = _materials_clean_name(component, "component")
        if key not in self.recovery_rates:
            raise KeyError(f"No recovery_rate defined for component {key!r} on {self.input_material!r}.")
        contained = composition.contained_mass(
            input_mass, key, input_basis=input_basis, target_unit="kilogram"
        )
        rec = self.recovery_rates[key]
        dst = parse_unit(target_unit)
        return _materials_kg_to_target(_materials_to_decimal(contained.magnitude, "contained") * rec, dst, input_mass)

    @authorized_conversion
    def gross_processing_worth(
        self,
        output_prices: Mapping[str, pint.Quantity],
        *,
        target_unit: str | Unit = "USD / metric_tonne",
        apply_recovery: bool = True,
    ) -> pint.Quantity:
        """Compute co-product Gross Processing Worth (GPW) per unit of feedstock mass."""
        from core.units.composition import _materials_MASS_FAMILY_DIMS
        from core.units.composition import _materials_preserve_numeric_type
        from core.units.composition import _materials_to_decimal
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        from core.units.conversions import _split_ratio_units

        if not self.output_yields:
            raise ValueError("ProcessingSpecification defines no output_yields.")
        dst = parse_unit(target_unit)
        dst_num, dst_den = _split_ratio_units(dst)
        if dst_den.dimensionality not in _materials_MASS_FAMILY_DIMS:
            raise ValueError("gross_processing_worth target_unit must be currency per mass.")

        total_per_kg_dst_ccy = Decimal(0)
        sample_sources: list[Any] = []
        for prod in self.output_yields:
            if prod not in output_prices:
                raise KeyError(f"Missing market price for output product {prod!r}.")
            px = output_prices[prod]
            validate_registry_quantity(px)
            sample_sources.append(px)
            px_per_kg = px.to(dst_num / ureg.kilogram)
            yld = self.effective_yield(prod, apply_recovery=apply_recovery)
            total_per_kg_dst_ccy += _materials_to_decimal(px_per_kg.magnitude, f"price[{prod!r}]") * yld

        gpw = ureg.Quantity(total_per_kg_dst_ccy, dst_num / ureg.kilogram).to(dst)
        return ureg.Quantity(_materials_preserve_numeric_type(gpw.magnitude, *sample_sources), dst)

    @authorized_conversion
    def processing_margin(
        self,
        input_price: pint.Quantity,
        output_prices: Mapping[str, pint.Quantity],
        *,
        processing_cost: pint.Quantity | None = None,
        apply_recovery: bool = True,
    ) -> pint.Quantity:
        """Compute net crush / refining / smelting margin in ``input_price.units``."""
        from core.units.composition import _materials_preserve_numeric_type
        from core.units.composition import _materials_to_decimal
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        validate_registry_quantity(input_price)
        gpw = self.gross_processing_worth(
            output_prices, target_unit=input_price.units, apply_recovery=apply_recovery
        )
        in_dec = _materials_to_decimal(input_price.magnitude, "input_price")
        gpw_dec = _materials_to_decimal(gpw.magnitude, "gpw")
        cost_dec = Decimal(0)
        if processing_cost is not None:
            validate_registry_quantity(processing_cost)
            cost_dec = _materials_to_decimal(processing_cost.to(input_price.units).magnitude, "processing_cost")
        net = gpw_dec - in_dec - cost_dec
        return ureg.Quantity(
            _materials_preserve_numeric_type(net, input_price, processing_cost, *output_prices.values()),
            input_price.units,
        )

    def blend_with(
        self,
        other: ProcessingSpecification,
        *,
        self_input_fraction: Decimal | float,
    ) -> ProcessingSpecification:
        """Blend two ProcessingSpecifications by feedstock mass fraction, conserving recovered outputs."""
        from core.units.composition import _materials_validate_fraction
        if not isinstance(other, ProcessingSpecification):
            raise TypeError(f"Expected ProcessingSpecification, got {type(other).__name__}.")
        w1 = _materials_validate_fraction(self_input_fraction, "self_input_fraction")
        w2 = Decimal(1) - w1

        all_outs = sorted(set(self.output_yields) | set(other.output_yields))
        b_yields = {
            k: (w1 * self.output_yields.get(k, Decimal(0))) + (w2 * other.output_yields.get(k, Decimal(0)))
            for k in all_outs
        }

        all_losses = sorted(set(self.process_losses) | set(other.process_losses))
        b_losses = {
            k: (w1 * self.process_losses.get(k, Decimal(0))) + (w2 * other.process_losses.get(k, Decimal(0)))
            for k in all_losses
        }

        all_recs = sorted(set(self.recovery_rates) | set(other.recovery_rates))
        b_recs: dict[str, Decimal] = {}
        for k in all_recs:
            r1 = self.recovery_rates.get(k, Decimal(1) if k in self.output_yields else Decimal(0))
            r2 = other.recovery_rates.get(k, Decimal(1) if k in other.output_yields else Decimal(0))
            if k in b_yields and b_yields[k] > Decimal(0):
                u1 = w1 * self.output_yields.get(k, Decimal(0))
                u2 = w2 * other.output_yields.get(k, Decimal(0))
                b_recs[k] = ((u1 * r1) + (u2 * r2)) / (u1 + u2)
            else:
                b_recs[k] = (w1 * r1) + (w2 * r2)

        label = self.input_material if self.input_material == other.input_material else f"{self.input_material}/{other.input_material}"
        return ProcessingSpecification(
            input_material=label,
            output_yields=b_yields,
            recovery_rates=b_recs,
            process_losses=b_losses,
        )


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (EnergyStorageSpecification,))
__all__ = ('CHPBalance', 'CHPSpecification', 'ElectrolysisSpecification', 'EnergyLoss', 'EnergyLossBalance', 'EnergyStorageSpecification', 'GenerationProcess', 'HeatPumpSpecification', 'HydrogenCrackingBalance', 'HydrogenCrackingSpecification', 'ProcessingSpecification', 'apply_energy_losses', 'blend_heat_content', 'convert_energy_storage', 'convert_thermal_efficiency', 'create_energy_storage_context', 'energy_before_losses', 'gross_generation', 'net_generation')
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
    from core.units.composition import CompositionBasis
    from core.units.composition import MaterialComposition
    from core.units.registry import Unit
    from core.units.processing import convert_energy_storage
    from core.units.processing import create_energy_storage_context
