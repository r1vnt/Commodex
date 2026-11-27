"""Commodex heat_content implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Iterator, Mapping
from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, field
from decimal import Decimal
from enum import StrEnum
from itertools import permutations
from pathlib import Path
from types import MappingProxyType
from typing import Any, Final
from typing import Any, Final, Literal
import json
import pint

class CalorificMethod(StrEnum):
    MEASURED = "measured"
    COMPONENT_SUM = "component_sum"
    LATENT_WATER_APPROXIMATION = "latent_water"
    EMPIRICAL_CORRELATION = "empirical"
    EXPLICIT_RATIO = "explicit_ratio"
    STANDARD_LABORATORY = "standard_laboratory"

class ThermodynamicBasis(StrEnum):
    UNSPECIFIED = "unspecified"
    CONSTANT_PRESSURE = "constant_pressure"
    CONSTANT_VOLUME = "constant_volume"

def _calorific_models_dec(value, label):
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise TypeError(f"{label} must be a numeric scalar.")
    result = Decimal(str(value))
    if not result.is_finite():
        raise ValueError(f"{label} must be finite.")
    return result

def _calorific_models_text(value, label):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a nonempty string.")
    return value.strip()

def _calorific_models_q(value, unit, label, *, nonnegative=False):
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"{label} requires a quantity with units.")
    validate_registry_quantity(value)
    result = ureg.Quantity(_calorific_models_dec(value.magnitude, label), value.units).to(unit)
    if nonnegative and result.magnitude < 0:
        raise ValueError(f"{label} must be nonnegative.")
    return result

def _calorific_models_temp(value):
    result = _calorific_models_q(value, "K", "combustion temperature")
    if result.magnitude <= 0:
        raise ValueError("Combustion temperature must exceed absolute zero.")
    return result





def _calorific_models_analysis_basis(value):
    from core.units.heat_content import _basis
    result = _basis(value)
    if result not in {"AR", "AD", "DB", "DAF", "DMMF"}:
        raise ValueError("Fuel analysis requires AR, AD, DB, DAF or DMMF basis.")
    return result





def _calorific_models_percentages(instance, names):
    total = Decimal(0)
    for name in names:
        value = _calorific_models_dec(getattr(instance, name), name)
        if not 0 <= value <= 100:
            raise ValueError(f"{name} must be in [0, 100].")
        object.__setattr__(instance, name, value)
        total += value
    if abs(total - 100) > Decimal("0.000001"):
        raise ValueError("Declared fuel analysis must total 100 percent; normalization is explicit.")
    if instance.basis in {"DB", "DAF", "DMMF"} and instance.moisture_pct != 0:
        raise ValueError("Dry analytical bases require zero sample moisture.")
    if instance.basis in {"DAF", "DMMF"} and instance.ash_pct != 0:
        raise ValueError("Excluded mineral/ash-free analysis must declare zero residual ash.")
    if instance.moisture_pct == 100:
        raise ValueError("Fuel analysis must contain dry matter.")





@dataclass(frozen=True, slots=True)
class UltimateAnalysis:
    """Mass percentages on one basis; H and O exclude the sample's water.

    All CHONS, ash, sample water and explicitly unreported material total 100%.
    Moisture-bound hydrogen is reported separately, never used as fuel hydrogen.
    """
    carbon_pct: Decimal | float
    hydrogen_fuel_bound_pct: Decimal | float
    oxygen_pct: Decimal | float
    nitrogen_pct: Decimal | float
    sulphur_pct: Decimal | float
    ash_pct: Decimal | float
    moisture_pct: Decimal | float
    other_pct: Decimal | float = 0
    basis: str = "AR"

    def __post_init__(self):
        object.__setattr__(self, "basis", _calorific_models_analysis_basis(self.basis))
        _calorific_models_percentages(self, ("carbon_pct", "hydrogen_fuel_bound_pct", "oxygen_pct", "nitrogen_pct",
                            "sulphur_pct", "ash_pct", "moisture_pct", "other_pct"))

    @property
    def moisture_bound_hydrogen_pct(self):
        return self.moisture_pct / Decimal("8.936")

    @property
    def total_hydrogen_pct(self):
        return self.hydrogen_fuel_bound_pct + self.moisture_bound_hydrogen_pct





@dataclass(frozen=True, slots=True)
class ProximateAnalysis:
    """Moisture, ash, volatile matter and fixed carbon, all on the same basis."""
    moisture_pct: Decimal | float
    ash_pct: Decimal | float
    volatile_matter_pct: Decimal | float
    fixed_carbon_pct: Decimal | float
    basis: str = "AR"

    def __post_init__(self):
        object.__setattr__(self, "basis", _calorific_models_analysis_basis(self.basis))
        _calorific_models_percentages(self, ("moisture_pct", "ash_pct", "volatile_matter_pct", "fixed_carbon_pct"))


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class EmpiricalHeatCorrelation(_safety_QuantitySnapshot):
    """An explicitly supplied linear ultimate-analysis HHV correlation.

    Coefficients are MJ/kg per unit mass fraction (not percentage point).
    No built-in coefficient set or applicability claim is assumed. The source
    and accepted element ranges are part of the user's chosen correlation.
    """
    name: str
    source: str
    coefficients: Mapping[str, Decimal]
    basis: str
    combustion_temperature: pint.Quantity
    thermodynamic_basis: ThermodynamicBasis
    valid_ranges_pct: Mapping[str, tuple[Decimal, Decimal]]
    intercept_mj_per_kg: Decimal = Decimal(0)

    def __post_init__(self):
        for name in ("name", "source"):
            object.__setattr__(self, name, _calorific_models_text(getattr(self, name), name))
        object.__setattr__(self, "basis", _calorific_models_analysis_basis(self.basis))
        object.__setattr__(self, "combustion_temperature", _calorific_models_temp(self.combustion_temperature))
        thermo = ThermodynamicBasis(self.thermodynamic_basis)
        if thermo == ThermodynamicBasis.UNSPECIFIED:
            raise ValueError("Correlation requires its thermodynamic basis.")
        object.__setattr__(self, "thermodynamic_basis", thermo)
        allowed = set(UltimateAnalysis.__dataclass_fields__) - {"basis"}
        if not self.coefficients or not set(self.coefficients) <= allowed or not set(self.valid_ranges_pct) <= allowed:
            raise ValueError("Correlation coefficients and validity ranges must name ultimate-analysis fields.")
        if not set(self.coefficients) <= set(self.valid_ranges_pct):
            raise ValueError("Declare validity ranges for every correlation coefficient.")
        coefficients = {key: _calorific_models_dec(value, key) for key, value in self.coefficients.items()}
        ranges = {}
        for key, (low, high) in self.valid_ranges_pct.items():
            low, high = _calorific_models_dec(low, key), _calorific_models_dec(high, key)
            if not 0 <= low <= high <= 100:
                raise ValueError("Correlation percentage validity ranges must lie within [0, 100].")
            ranges[key] = (low, high)
        object.__setattr__(self, "coefficients", MappingProxyType(coefficients))
        object.__setattr__(self, "valid_ranges_pct", MappingProxyType(ranges))
        object.__setattr__(self, "intercept_mj_per_kg", _calorific_models_dec(self.intercept_mj_per_kg, "intercept"))

    def evaluate(self, analysis: UltimateAnalysis):
        from core.units.registry import ureg
        if not isinstance(analysis, UltimateAnalysis) or analysis.basis != self.basis:
            raise ValueError("Correlation and ultimate analysis must share their declared analytical basis.")
        for key, (low, high) in self.valid_ranges_pct.items():
            if not low <= getattr(analysis, key) <= high:
                raise ValueError(f"{key} is outside correlation applicability.")
        value = self.intercept_mj_per_kg + sum((coefficient * getattr(analysis, key) / 100 for key, coefficient in self.coefficients.items()), Decimal(0))
        if value <= 0:
            raise ValueError("Correlation produced a nonpositive calorific value.")
        return ureg.Quantity(value, "MJ/kg")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LaboratoryCalorificResult(_safety_QuantitySnapshot):
    """External laboratory-standard result, with report/method traceability.

    This record imports a calculation or report. It does not execute ISO 1928,
    ISO 18125, acid corrections or bomb-calorimeter calibration itself.
    """
    standard: str
    report_id: str
    source: str
    basis: str
    thermodynamic_basis: ThermodynamicBasis
    combustion_temperature: pint.Quantity
    gcv: pint.Quantity | None = None
    ncv: pint.Quantity | None = None

    def __post_init__(self):
        for name in ("standard", "report_id", "source"):
            object.__setattr__(self, name, _calorific_models_text(getattr(self, name), name))
        object.__setattr__(self, "basis", _calorific_models_analysis_basis(self.basis))
        thermo = ThermodynamicBasis(self.thermodynamic_basis)
        if thermo == ThermodynamicBasis.UNSPECIFIED:
            raise ValueError("Laboratory report must declare constant pressure or volume.")
        object.__setattr__(self, "thermodynamic_basis", thermo)
        object.__setattr__(self, "combustion_temperature", _calorific_models_temp(self.combustion_temperature))
        if self.gcv is None and self.ncv is None:
            raise ValueError("Laboratory result requires at least one measured/calculated value.")
        for name in ("gcv", "ncv"):
            value = getattr(self, name)
            if value is not None:
                value = _calorific_models_q(value, "MJ/kg", name)
                if value.magnitude <= 0:
                    raise ValueError("Laboratory calorific values must be positive.")
                object.__setattr__(self, name, value)
        if self.gcv is not None and self.ncv is not None and self.ncv > self.gcv:
            raise ValueError("NCV cannot exceed GCV.")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class ISO6976Table(_safety_QuantitySnapshot):
    """User-supplied ISO 6976 component data at one combustion temperature.

    Each component maps to (gross, net) kJ/mol. Table edition/source are required;
    no unprovided temperatures, pseudo-components or uncertainties are inferred.
    Metering density/Z remains an independent explicit gas-state calculation.
    """
    edition: str
    source: str
    combustion_temperature: pint.Quantity
    component_values_kj_per_mol: Mapping[str, tuple[Decimal, Decimal]]

    def __post_init__(self):
        for name in ("edition", "source"):
            object.__setattr__(self, name, _calorific_models_text(getattr(self, name), name))
        object.__setattr__(self, "combustion_temperature", _calorific_models_temp(self.combustion_temperature))
        values = {}
        from core.units.composition import GAS_MOLAR_MASSES
        for name, (gross, net) in self.component_values_kj_per_mol.items():
            if name not in GAS_MOLAR_MASSES:
                raise ValueError(f"Unknown canonical gas component {name!r}.")
            if name in values:
                raise ValueError("Duplicate component aliases in ISO 6976 table.")
            gross, net = _calorific_models_dec(gross, "gross component value"), _calorific_models_dec(net, "net component value")
            if not 0 <= net <= gross:
                raise ValueError("Component data must satisfy 0 <= net <= gross.")
            if name == "H2O" and (gross != 0 or net != 0):
                raise ValueError("Spectator water requires zero component heat.")
            values[name] = (gross, net)
        if not values:
            raise ValueError("ISO 6976 table cannot be empty.")
        object.__setattr__(self, "component_values_kj_per_mol", MappingProxyType(values))

    def calculate(self, composition):
        from core.units.registry import ureg
        from core.units.composition import GasComposition
        from core.units.composition import GasCalorificValues
        if not isinstance(composition, GasComposition):
            raise TypeError("ISO 6976 calculation requires GasComposition.")
        missing = set(composition.mole_fractions) - set(self.component_values_kj_per_mol) - {"H2O"}
        if missing:
            raise ValueError(f"Missing ISO 6976 component values: {sorted(missing)}.")
        totals = [sum((fraction * self.component_values_kj_per_mol[name][index]
                       for name, fraction in composition.mole_fractions.items() if name != "H2O"), Decimal(0)) for index in (0, 1)]
        return GasCalorificValues(*(ureg.Quantity(total / composition.molar_mass.magnitude, "MJ/kg") for total in totals),
            self.combustion_temperature, f"ISO 6976:{self.edition} component summation; {self.source}",
            method=CalorificMethod.COMPONENT_SUM, thermodynamic_basis=ThermodynamicBasis.CONSTANT_PRESSURE)





def water_latent_heat(temperature, *, method):
    """Explicit saturation latent heat from CoolProp HEOS (IAPWS-95) or IF97.

    This is a saturation-property approximation for a fuel latent-water model,
    not the standard-state ideal vapour enthalpy or a laboratory calculation.
    https://coolprop.org/coolprop/HighLevelAPI.html#vapor-liquid-and-saturation-states
    """
    from core.units.registry import ureg
    temperature = _calorific_models_temp(temperature)
    if method not in {"iapws95", "if97"}:
        raise ValueError("Latent-heat method must be 'iapws95' or 'if97'.")
    if not Decimal("273.16") <= temperature.magnitude < Decimal("647.096"):
        raise ValueError("Liquid/vapour saturation latent heat requires 273.16 <= T < 647.096 K.")
    try:
        from CoolProp.CoolProp import PropsSI
    except ImportError as exc:
        raise ImportError("Temperature-dependent water latent heat requires optional CoolProp.") from exc
    fluid = "HEOS::Water" if method == "iapws95" else "IF97::Water"
    vapor = PropsSI("H", "T", float(temperature.magnitude), "Q", 1, fluid)
    liquid = PropsSI("H", "T", float(temperature.magnitude), "Q", 0, fluid)
    return ureg.Quantity(Decimal(str(vapor - liquid)), "J/kg").to("MJ/kg")





def ideal_gas_combustion_correction(delta_gas_moles, temperature):
    """Signed delta(PV) per fuel amount = (n_products - n_reactants) R T.

    Stoichiometry and product-water phase must match the specific HHV/LHV
    reaction. No condensed-phase PV or real-gas correction is inferred.
    """
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(delta_gas_moles, pint.Quantity):
        raise TypeError("Supply signed stoichiometric delta gas moles per fuel amount.")
    validate_registry_quantity(delta_gas_moles)
    amount = ureg.Quantity(_calorific_models_dec(delta_gas_moles.magnitude, "delta gas moles"), delta_gas_moles.units)
    if not any(amount.is_compatible_with(unit) for unit in ("mol/kg", "mol/mol", "mol/m^3")):
        raise ValueError("Delta gas moles must be per kg, mole or physical volume of fuel.")
    return amount * ureg.Quantity(Decimal("8.314462618"), "J/mol/K") * _calorific_models_temp(temperature)





def convert_combustion_heat(heat, *, source, target, pressure_volume_correction):
    """Convert positive released heat: Qp = Qv - delta(PV).

    The signed reaction correction must be explicitly supplied in matching
    specific-energy units. This thermodynamic identity does not perform bomb
    calorimeter acid, fuse-wire or laboratory-standard corrections.
    """
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    source, target = ThermodynamicBasis(source), ThermodynamicBasis(target)
    if ThermodynamicBasis.UNSPECIFIED in {source, target}:
        raise ValueError("Constant-pressure/volume heat conversion requires explicit bases.")
    if source == target:
        raise ValueError("Source and target thermodynamic bases must differ.")
    if not isinstance(heat, pint.Quantity):
        raise TypeError("Combustion heat requires a quantity.")
    validate_registry_quantity(heat)
    value = _calorific_models_dec(heat.magnitude, "combustion heat")
    if value <= 0 or not any(heat.is_compatible_with(unit) for unit in ("MJ/kg", "MJ/mol", "MJ/m^3")):
        raise ValueError("Combustion heat must be positive specific energy.")
    correction = _calorific_models_q(pressure_volume_correction, heat.units, "reaction PV correction").magnitude
    result = value - correction if target == ThermodynamicBasis.CONSTANT_PRESSURE else value + correction
    if result <= 0:
        raise ValueError("Reaction correction produced nonpositive heat.")
    return ureg.Quantity(result if isinstance(heat.magnitude, Decimal) else float(result), heat.units)




ANALYTICAL_BASES: Final[tuple[str, ...]] = (
    "AR", "AD", "DB", "DAF", "DMMF", "AS_FIRED", "AS_DELIVERED",
    "REFERENCE_MOISTURE", "DRY_GAS", "WET_GAS",
)



_GAS_BASES = frozenset({"DRY_GAS", "WET_GAS"})



_DRY_SOLID_BASES = frozenset({"DB", "DAF", "DMMF"})



_CONDITION_FIELDS = {"AS_FIRED": "as_fired_state", "AS_DELIVERED": "as_delivered_state",
                     "REFERENCE_MOISTURE": "reference_moisture_state"}





def _number(value: Any, label: str, *, positive: bool = False) -> Decimal:
    from core.units.conversions import _as_decimal
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise TypeError(f"{label} must be a numeric scalar.")
    result = _as_decimal(value)
    if not result.is_finite() or (result <= 0 if positive else result < 0):
        raise ValueError(f"{label} must be finite and {'positive' if positive else 'non-negative'}.")
    return result





def _fraction(value: Any, label: str) -> Decimal:
    result = _number(value, label) / 100
    if result >= 1:
        raise ValueError(f"{label} must be in [0, 100).")
    return result





def _basis(value: str) -> str:
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    if not isinstance(value, str):
        raise TypeError("Analytical basis must be a string.")
    result = value.upper().strip().replace("-", "_").replace(" ", "_")
    result = {"DRY": "DB", "D": "DB", "ADB": "AD", "MAF": "DAF", "AF": "AS_FIRED",
              "FIRED": "AS_FIRED", "DELIVERED": "AS_DELIVERED", "REFERENCE": "REFERENCE_MOISTURE",
              "RM": "REFERENCE_MOISTURE"}.get(result, result)
    if result not in ANALYTICAL_BASES:
        raise IncompatibleAnalyticalBasisError(f"Analytical basis must be one of {ANALYTICAL_BASES} (ADB and MAF are supported aliases).")
    return result





def _kind(value: str) -> str:
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    result = value.lower().strip()
    if result in {"hhv", "gcv", "gar", "gross"}:
        return "gcv"
    if result in {"lhv", "ncv", "nar", "net"}:
        return "ncv"
    raise IncompatibleAnalyticalBasisError("Heating-value basis must be GCV/HHV or NCV/LHV.")





def _report_basis(kind, basis):
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    if kind.lower().strip() in {"gar", "nar"}:
        if basis is not None and _basis(basis) != "AR":
            raise IncompatibleAnalyticalBasisError("GAR/NAR mean gross/net AS RECEIVED; analytical basis must be AR.")
        return "AR"
    return basis





def _quantity(value: pint.Quantity, label: str) -> pint.Quantity:
    from core.units.conversions import _as_decimal
    from core.units.dimensions import physical_quantity
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"{label} must be a Pint quantity with units.")
    validate_registry_quantity(value)
    value = physical_quantity(value, role="thermal")
    _number(value.magnitude, label, positive=True)
    # Use Decimal internally without converting a caller's float arithmetic
    # into an assertion of additional measurement precision.
    return ureg.Quantity(_as_decimal(value.magnitude), value.units)





def _is_volume_dimension(dimensionality: Any) -> bool:
    from core.units.registry import ureg
    return dimensionality == (ureg.meter ** 3).dimensionality or (
        len(dimensionality) == 1 and "gas_volume" in next(iter(dimensionality))
        and next(iter(dimensionality.values())) == 1
    )





def _input_type(result: pint.Quantity, source: pint.Quantity) -> pint.Quantity:
    from core.units.precision import pad_ieee754_decimal
    from core.units.registry import ureg
    if isinstance(source.magnitude, Decimal):
        return ureg.Quantity(pad_ieee754_decimal(result.magnitude), result.units)
    return ureg.Quantity(float(result.magnitude), result.units)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class FuelCondition(_safety_QuantitySnapshot):
    """Declared firing, delivery or contractual moisture condition.

    Moisture is percent of total condition mass. Optional density and T/P apply
    to this condition only. Temperature is metadata, not a sensible-heat or
    variable-latent-heat correction to the combustion reference calculation.
    """
    moisture_pct: Decimal | float
    temperature: pint.Quantity | None = None
    absolute_pressure: pint.Quantity | None = None
    density: pint.Quantity | None = None
    label: str | None = None

    def __post_init__(self):
        object.__setattr__(self, "moisture_pct", _fraction(self.moisture_pct, "condition moisture") * 100)
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature as _temperature
        from core.units.assays.energy.petroleum import _converters_petroleum_pressure as _pressure
        if self.temperature is not None:
            object.__setattr__(self, "temperature", _temperature(self.temperature))
        if self.absolute_pressure is not None:
            object.__setattr__(self, "absolute_pressure", _pressure(self.absolute_pressure))
        if self.density is not None:
            object.__setattr__(self, "density", _quantity(self.density, "condition density").to("kg/m^3"))
        if self.label is not None and (not isinstance(self.label, str) or not self.label.strip()):
            raise ValueError("Condition label must be a nonempty string when supplied.")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GasMoistureState(_safety_QuantitySnapshot):
    """Gas metering state and its water-vapour fraction, independently of fuel moisture.

    A complete water mole fraction is required in GasReferenceConditions.
    Dry-gas molar mass is required for mass-basis conversions, while positive
    declared Z is required for volume/mole conversions. No EOS is inferred.
    """
    reference: Any
    dry_gas_molar_mass: pint.Quantity | None = None
    gas_constant: pint.Quantity = field(default_factory=lambda: ureg.Quantity(Decimal("8.314462618"), "J/mol/K"))
    source: str = "declared gas state"

    def __post_init__(self):
        from core.units.exceptions import MissingReferenceConditionsError
        from core.units.reference_conditions import GasReferenceConditions
        if not isinstance(self.reference, GasReferenceConditions):
            raise TypeError("Gas state reference must be GasReferenceConditions.")
        if self.reference.water_mole_fraction is None:
            raise MissingReferenceConditionsError("Gas state requires explicit water_mole_fraction, including zero for dry gas.")
        if self.dry_gas_molar_mass is not None:
            object.__setattr__(self, "dry_gas_molar_mass", _quantity(self.dry_gas_molar_mass, "dry-gas molar mass").to("kg/mol"))
        object.__setattr__(self, "gas_constant", _quantity(self.gas_constant, "gas constant").to("J/mol/K"))
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("Gas-state source must be a nonempty string.")

    @property
    def water_mole_fraction(self):
        return self.reference.water_mole_fraction

    @property
    def molar_mass(self):
        from core.units.exceptions import MissingAssayError
        from core.units.registry import ureg
        if self.dry_gas_molar_mass is None:
            raise MissingAssayError("Mass-basis gas conversion requires dry_gas_molar_mass.")
        from core.units.composition import GAS_MOLAR_MASSES
        y = self.water_mole_fraction
        water_mass = GAS_MOLAR_MASSES["H2O"] / 1000
        return ureg.Quantity((1 - y) * self.dry_gas_molar_mass.magnitude + y * water_mass, "kg/mol")

    @property
    def molar_volume(self):
        from core.units.exceptions import MissingReferenceConditionsError
        from core.units.registry import ureg
        z = self.reference.compressibility_factor
        if z is None:
            raise MissingReferenceConditionsError("Gas volume conversion requires an explicit compressibility factor Z.")
        volume = z * self.gas_constant.magnitude * self.reference.temperature.magnitude / self.reference.absolute_pressure.magnitude
        return ureg.Quantity(volume, "m^3/mol")

    @property
    def density(self):
        return (self.molar_mass / self.molar_volume).to("kg/m^3")

    def volume_per_mole(self, unit):
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
        parsed = parse_unit(unit)
        geometric = _geometric_volume(ureg.Quantity(Decimal(1), parsed), self.reference)
        return ureg.Quantity(self.molar_volume.magnitude / geometric, parsed / ureg.mol)

    @classmethod
    def from_composition(cls, composition, reference, **metadata):
        from core.units.composition import GasComposition
        from core.units.reference_conditions import GasReferenceConditions
        from dataclasses import replace
        if not isinstance(composition, GasComposition) or not isinstance(reference, GasReferenceConditions):
            raise TypeError("Supply GasComposition and GasReferenceConditions.")
        if reference.water_mole_fraction is not None and reference.water_mole_fraction != composition.water_mole_fraction:
            raise ValueError("Gas-state water fraction conflicts with composition.")
        reference = replace(reference, water_mole_fraction=composition.water_mole_fraction)
        return cls(reference, composition.dry().molar_mass, **metadata)





@dataclass(frozen=True, slots=True)
class HeatContentAssay:
    """Calorific data usable alongside any commodity assay.

    ``gcv``/``ncv`` are energy per mass, volume or mole quantities. Supply
    density and/or molar_mass to bridge those bases. ``standard_density`` can
    additionally bridge a registry-tagged standard gas volume. All densities
    describe this declared analytical basis and reference state.

    Hydrogen-carrier adapters report delivered H2 energy after cracking,
    matching HydrogenAssay, rather than the carrier's combustion heat.

    DMMF uses its own measured/declared mineral matter. Positive mineral
    calorific correction terms represent energy included in dry whole fuel
    but excluded from the mineral-free value; signed endothermic corrections
    are also supported. Zero/omitted corrections give mass normalization only.
    FuelCondition temperatures describe physical fuel states, independently of
    the combustion reference used by the calorific data and latent-heat model.
    """

    name: str = "custom_heat_content"
    gcv: pint.Quantity | None = None
    ncv: pint.Quantity | None = None
    basis: str = "AR"
    moisture_pct: Decimal | float | None = None
    as_received_moisture_pct: Decimal | float | None = None
    air_dried_moisture_pct: Decimal | float | None = None
    hydrogen_pct: Decimal | float | None = None
    ash_db_pct: Decimal | float | None = None
    density: pint.Quantity | None = None
    standard_density: pint.Quantity | None = None
    molar_mass: pint.Quantity | None = None
    reference_temperature: pint.Quantity | None = None
    reference_pressure: pint.Quantity | None = None
    ncv_to_gcv_ratio: Decimal | float | None = None
    latent_heat_mj_per_kg_water: Decimal | float = Decimal("2.442")
    water_per_kg_hydrogen: Decimal | float = Decimal("8.936")
    # Standard-state parity inherited only from an explicitly supplied GasAssay.
    standard_volume_factors: tuple[tuple[str, Decimal], ...] = ()
    gcv_reference_temperature: pint.Quantity | None = None
    ncv_reference_temperature: pint.Quantity | None = None
    gcv_reference_pressure: pint.Quantity | None = None
    ncv_reference_pressure: pint.Quantity | None = None
    mineral_matter_db_pct: Decimal | float | None = None
    mineral_matter_gcv_correction: pint.Quantity | None = None
    mineral_matter_ncv_correction: pint.Quantity | None = None
    as_fired_state: FuelCondition | None = None
    as_delivered_state: FuelCondition | None = None
    reference_moisture_state: FuelCondition | None = None
    dry_gas_state: GasMoistureState | None = None
    wet_gas_state: GasMoistureState | None = None
    method: CalorificMethod = CalorificMethod.MEASURED
    method_source: str | None = None
    gcv_method: CalorificMethod | None = None
    ncv_method: CalorificMethod | None = None
    thermodynamic_basis: ThermodynamicBasis = ThermodynamicBasis.UNSPECIFIED
    combustion_temperature: pint.Quantity | None = None
    latent_heat_method: str = "fixed"
    moisture_correction_method: CalorificMethod | None = None
    ultimate_analysis: UltimateAnalysis | None = None
    proximate_analysis: ProximateAnalysis | None = None
    laboratory_standard: str | None = None
    laboratory_report_id: str | None = None
    derivations: tuple[str, ...] = ()
    energy_role: str = "combustion"

    def __post_init__(self) -> None:
        from core.units.conversions import _as_decimal
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        from core.units.exceptions import InvalidPhysicalStateError
        from core.units.exceptions import MissingAssayError
        from core.units.exceptions import MissingReferenceConditionsError
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        object.__setattr__(self, "basis", _basis(self.basis))
        object.__setattr__(self, "method", CalorificMethod(self.method))
        if self.energy_role not in {"combustion", "recovered_hydrogen"}:
            raise ValueError("energy_role must distinguish combustion from recovered_hydrogen.")
        object.__setattr__(self, "thermodynamic_basis", ThermodynamicBasis(self.thermodynamic_basis))
        if self.moisture_correction_method is not None:
            correction_method = CalorificMethod(self.moisture_correction_method)
            if correction_method != CalorificMethod.LATENT_WATER_APPROXIMATION:
                raise ValueError("Sample-water correction requires explicit latent_water approximation.")
            object.__setattr__(self, "moisture_correction_method", correction_method)
        for key in ("gcv_method", "ncv_method"):
            if getattr(self, key) is not None:
                object.__setattr__(self, key, CalorificMethod(getattr(self, key)))
            elif getattr(self, key[:3]) is not None:
                selected = self.method if self.method in {CalorificMethod.COMPONENT_SUM, CalorificMethod.EMPIRICAL_CORRELATION, CalorificMethod.STANDARD_LABORATORY} else CalorificMethod.MEASURED
                object.__setattr__(self, key, selected)
        for key in ("method_source", "laboratory_standard", "laboratory_report_id"):
            value = getattr(self, key)
            if value is not None and (not isinstance(value, str) or not value.strip()):
                raise ValueError(f"{key} must be a nonempty string.")
        if self.method in {CalorificMethod.COMPONENT_SUM, CalorificMethod.EMPIRICAL_CORRELATION, CalorificMethod.STANDARD_LABORATORY} and not self.method_source:
            raise ValueError("Calculated calorific methods require method_source provenance.")
        if self.method == CalorificMethod.STANDARD_LABORATORY and (not self.laboratory_standard or not self.laboratory_report_id or self.thermodynamic_basis == ThermodynamicBasis.UNSPECIFIED or self.combustion_temperature is None):
            raise ValueError("Standard laboratory data require standard, report, thermodynamic basis and combustion temperature.")
        if self.method == CalorificMethod.EXPLICIT_RATIO and self.ncv_to_gcv_ratio is None:
            raise ValueError("explicit_ratio method requires ncv_to_gcv_ratio.")
        if self.combustion_temperature is not None:
            from core.units.heat_content import _calorific_models_temp as _temp
            object.__setattr__(self, "combustion_temperature", _temp(self.combustion_temperature))
        if self.latent_heat_method not in {"fixed", "iapws95", "if97"}:
            raise ValueError("latent_heat_method must be fixed, iapws95 or if97.")
        if self.latent_heat_method != "fixed":
            if (self.method != CalorificMethod.LATENT_WATER_APPROXIMATION and self.moisture_correction_method != CalorificMethod.LATENT_WATER_APPROXIMATION) or self.combustion_temperature is None:
                raise ValueError("Temperature-dependent latent heat requires latent_water method and combustion_temperature.")
            latent = water_latent_heat(self.combustion_temperature, method=self.latent_heat_method)
            object.__setattr__(self, "latent_heat_mj_per_kg_water", latent.magnitude)
        elif self.method == CalorificMethod.LATENT_WATER_APPROXIMATION and self.combustion_temperature is None:
            object.__setattr__(self, "combustion_temperature", ureg.Quantity(Decimal("298.15"), "K"))
        if not isinstance(self.derivations, tuple) or not all(isinstance(item, str) and item.strip() for item in self.derivations):
            raise ValueError("derivations must be a tuple of nonempty descriptions.")
        for key, model in (("ultimate_analysis", UltimateAnalysis), ("proximate_analysis", ProximateAnalysis)):
            analysis = getattr(self, key)
            if analysis is None:
                continue
            if not isinstance(analysis, model) or analysis.basis != self.basis:
                raise IncompatibleAnalyticalBasisError("Fuel analysis must match the declared calorific analytical basis.")
            if self.moisture_pct is not None and _as_decimal(self.moisture_pct) != analysis.moisture_pct:
                raise ValueError("Fuel analysis moisture conflicts with calorific moisture.")
            object.__setattr__(self, "moisture_pct", analysis.moisture_pct)
            ash_db = analysis.ash_pct / (1 - analysis.moisture_pct / 100)
            if self.ash_db_pct is not None and abs(_as_decimal(self.ash_db_pct) - ash_db) > Decimal("0.000001"):
                raise ValueError("Fuel analysis ash conflicts with dry-basis ash metadata.")
            object.__setattr__(self, "ash_db_pct", ash_db)
            if key == "ultimate_analysis":
                if self.hydrogen_pct is not None and _as_decimal(self.hydrogen_pct) != analysis.hydrogen_fuel_bound_pct:
                    raise ValueError("Fuel-bound hydrogen conflicts with ultimate analysis.")
                object.__setattr__(self, "hydrogen_pct", analysis.hydrogen_fuel_bound_pct)
        if self.gcv is None and self.ncv is None:
            raise MissingAssayError("Supply at least one measured GCV/HHV or NCV/LHV quantity.")
        for field_name in ("moisture_pct", "as_received_moisture_pct", "air_dried_moisture_pct", "ash_db_pct", "mineral_matter_db_pct"):
            value = getattr(self, field_name)
            if value is not None:
                _fraction(value, field_name)
        for field_name in _CONDITION_FIELDS.values():
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, FuelCondition):
                raise TypeError(f"{field_name} must be FuelCondition.")
        for field_name in ("dry_gas_state", "wet_gas_state"):
            value = getattr(self, field_name)
            if value is not None and not isinstance(value, GasMoistureState):
                raise TypeError(f"{field_name} must be GasMoistureState.")
        if self.dry_gas_state is not None and self.dry_gas_state.water_mole_fraction != 0:
            raise ValueError("dry_gas_state must declare zero water mole fraction.")
        if self.dry_gas_state is not None and self.wet_gas_state is not None:
            a, b = self.dry_gas_state.dry_gas_molar_mass, self.wet_gas_state.dry_gas_molar_mass
            if a is not None and b is not None and abs(a.magnitude - b.magnitude) > max(a.magnitude, b.magnitude) * Decimal("1e-12"):
                raise ValueError("Wet and dry states must describe the same dry-gas composition/molar mass.")
        for field_name in ("mineral_matter_gcv_correction", "mineral_matter_ncv_correction"):
            value = getattr(self, field_name)
            if value is not None:
                if not isinstance(value, pint.Quantity):
                    raise TypeError("Mineral-matter calorific corrections require energy-per-dry-mass quantities.")
                validate_registry_quantity(value)
                number = _as_decimal(value.magnitude)
                if not number.is_finite():
                    raise ValueError("Mineral-matter calorific correction must be finite.")
                object.__setattr__(self, field_name, ureg.Quantity(number, value.units).to("MJ/kg"))
        if self.basis in _CONDITION_FIELDS:
            condition = self._condition(self.basis)
            if condition is None:
                raise MissingReferenceConditionsError(f"{self.basis} requires an explicit {_CONDITION_FIELDS[self.basis]}.")
            if self.moisture_pct is not None and _as_decimal(self.moisture_pct) != condition.moisture_pct:
                raise IncompatibleAnalyticalBasisError("Sample moisture conflicts with the declared fuel condition.")
            object.__setattr__(self, "moisture_pct", condition.moisture_pct)
            self._inherit_state_metadata(condition.temperature, condition.absolute_pressure, condition.density)
        if self.basis in _GAS_BASES:
            state = self._gas_state(self.basis)
            if self.moisture_pct not in (None, 0):
                raise IncompatibleAnalyticalBasisError("Gas water vapour uses water mole fractions, not solid-fuel moisture_pct.")
            self._inherit_state_metadata(state.reference.temperature, state.reference.absolute_pressure,
                                        state.density if state.dry_gas_molar_mass is not None and state.reference.compressibility_factor is not None else None)
            if state.dry_gas_molar_mass is not None:
                if self.molar_mass is not None:
                    declared = _quantity(self.molar_mass, "gas molar mass").to("kg/mol").magnitude
                    if abs(declared / state.molar_mass.magnitude - 1) > Decimal("1e-12"):
                        raise ValueError("Declared molar_mass conflicts with gas moisture state.")
                else:
                    object.__setattr__(self, "molar_mass", state.molar_mass)
            # Named standard gas tags declare water-free volumes. Wet gas must
            # use a geometric volume at the explicit wet metering state.
            from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
            for value in (self.gcv, self.ncv):
                if isinstance(value, pint.Quantity):
                    denominator = ureg.MJ / value.units
                    if _is_volume_dimension(denominator.dimensionality):
                        _geometric_volume(ureg.Quantity(Decimal(1), denominator), state.reference)
        if self.hydrogen_pct is not None and _number(self.hydrogen_pct, "hydrogen_pct") > 100:
            raise ValueError("hydrogen_pct must be in [0, 100].")
        for field in ("latent_heat_mj_per_kg_water", "water_per_kg_hydrogen"):
            _number(getattr(self, field), field, positive=True)
        if self.ncv_to_gcv_ratio is not None:
            if _number(self.ncv_to_gcv_ratio, "ncv_to_gcv_ratio", positive=True) > 1:
                raise ValueError("ncv_to_gcv_ratio must be in (0, 1].")
        if self.basis in _DRY_SOLID_BASES and self.moisture_pct not in (None, 0):
            raise IncompatibleAnalyticalBasisError("Dry-basis assays cannot declare nonzero sample moisture.")
        declared_moisture = self.as_received_moisture_pct if self.basis == "AR" else self.air_dried_moisture_pct if self.basis == "AD" else None
        if self.moisture_pct is not None and declared_moisture is not None:
            if _as_decimal(self.moisture_pct) != _as_decimal(declared_moisture):
                raise IncompatibleAnalyticalBasisError("Sample moisture must agree with the moisture declared for its analytical basis.")
        dimensions = set()
        for field in ("density", "standard_density"):
            value = getattr(self, field)
            if value is None:
                continue
            value = _quantity(value, field)
            volume_dim = ureg.kg.dimensionality / value.dimensionality
            if not _is_volume_dimension(volume_dim):
                raise ValueError(f"{field} must have mass per volume units.")
            key = str(volume_dim)
            if key in dimensions:
                raise ValueError("Supply only one density for each volume dimensionality.")
            dimensions.add(key)
            if self.basis in _GAS_BASES and field == "standard_density":
                from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
                state = self._gas_state(self.basis)
                volume_unit = ureg.kg / value.units
                volume_m3 = _geometric_volume(ureg.Quantity(Decimal(1), volume_unit), state.reference)
                if state.dry_gas_molar_mass is not None and state.reference.compressibility_factor is not None:
                    expected = state.density.magnitude * volume_m3
                    if abs(value.magnitude / expected - 1) > Decimal("1e-12"):
                        raise ValueError("standard_density conflicts with the declared gas state.")
            object.__setattr__(self, field, value)
        if self.molar_mass is not None:
            object.__setattr__(self, "molar_mass", _quantity(self.molar_mass, "molar_mass").to("kg/mol"))
        if self.reference_temperature is not None:
            _quantity(self.reference_temperature.to("kelvin"), "reference_temperature")
        if self.reference_pressure is not None:
            pressure = self.reference_pressure
            from core.units.conversions import _classify_pressure_datum
            if _classify_pressure_datum(pressure.units) == "gauge":
                raise InvalidPhysicalStateError("Reference pressure must be absolute, not gauge.")
            _quantity(pressure.to("pascal"), "reference_pressure")
        # Measured values must already describe one reference state. No hidden
        # ideal-gas correction or density extrapolation is performed here.
        from core.units.conversions import _classify_pressure_datum
        for kind, unit in (("temperature", "kelvin"), ("pressure", "pascal")):
            references = []
            for prefix in ("", "gcv_", "ncv_"):
                reference = getattr(self, prefix + "reference_" + kind)
                if reference is None:
                    continue
                if prefix and getattr(self, prefix[:-1]) is None:
                    raise MissingAssayError("A calorific reference state requires its measured value.")
                if kind == "pressure" and _classify_pressure_datum(reference.units) == "gauge":
                    raise InvalidPhysicalStateError("Calorific reference pressure must be absolute.")
                references.append(_quantity(reference.to(unit), "calorific reference " + kind).magnitude)
            if references and any(abs(value - references[0]) > Decimal("1e-8") for value in references[1:]):
                raise MissingReferenceConditionsError("Measured calorific reference states differ; normalize explicitly first.")
        if self.gcv is not None and self.ncv is not None:
            gross_denominator = ureg.MJ.dimensionality / self.gcv.dimensionality
            net_denominator = ureg.MJ.dimensionality / self.ncv.dimensionality
            if _is_volume_dimension(gross_denominator) and _is_volume_dimension(net_denominator) and gross_denominator != net_denominator:
                raise MissingReferenceConditionsError("Measured GCV and NCV use different volume reference states; normalize them explicitly first.")
        seen = set()
        for unit, factor in self.standard_volume_factors:
            if unit in seen or not _is_volume_dimension(ureg.Unit(unit).dimensionality):
                raise ValueError("Invalid/duplicate standard gas volume reference.")
            seen.add(unit)
            _number(factor, "standard volume factor", positive=True)
        for field in ("gcv", "ncv"):
            value = getattr(self, field)
            if value is not None:
                object.__setattr__(self, field, self._to_mass_specific(_quantity(value, field)))
        if self.gcv is not None and self.ncv is not None and self.gcv.dimensionality != self.ncv.dimensionality:
            raise MissingAssayError("Measured calorific values require matching density/molar-mass metadata and reference states.")
        if self.gcv is not None and self.ncv is not None and self.ncv > self.gcv:
            raise InvalidPhysicalStateError("NCV/LHV cannot exceed GCV/HHV.")

    def _densities(self) -> tuple[pint.Quantity, ...]:
        return tuple(q for q in (self.density, self.standard_density) if q is not None)

    def _condition(self, basis):
        field_name = _CONDITION_FIELDS.get(basis)
        return getattr(self, field_name) if field_name is not None else None

    def _gas_state(self, basis):
        from core.units.exceptions import MissingReferenceConditionsError
        state = self.dry_gas_state if basis == "DRY_GAS" else self.wet_gas_state
        if state is None:
            raise MissingReferenceConditionsError(f"{basis} requires an explicit {'dry' if basis == 'DRY_GAS' else 'wet'}_gas_state.")
        return state

    def _inherit_state_metadata(self, temperature, pressure, density):
        from core.units.conversions import _as_decimal
        for field_name, value, unit in (("reference_temperature", temperature, "K"),
                                        ("reference_pressure", pressure, "Pa"), ("density", density, "kg/m^3")):
            if value is None:
                continue
            declared = getattr(self, field_name)
            if declared is not None:
                a = _as_decimal(declared.to(unit).magnitude)
                b = _as_decimal(value.to(unit).magnitude)
                if abs(a - b) > max(abs(a), abs(b), Decimal(1)) * Decimal("1e-12"):
                    raise ValueError(f"{field_name} conflicts with the declared fuel/gas condition.")
            else:
                object.__setattr__(self, field_name, value)

    def _to_mass_specific(self, value: pint.Quantity) -> pint.Quantity:
        from core.units.exceptions import MissingAssayError
        from core.units.registry import ureg
        volumetric = [value, *self._standard_specific_values(value)]
        candidates = list(volumetric)
        candidates.extend(candidate / density for candidate in volumetric for density in self._densities())
        if self.molar_mass is not None:
            candidates.append(value / self.molar_mass)
        for candidate in candidates:
            if candidate.is_compatible_with("MJ/kg"):
                return candidate.to("MJ/kg")
        # Volumetric-only gas and molar-only specifications do not imply a
        # mass density. They still support energy parity on their native basis.
        denominator_dim = ureg.MJ.dimensionality / value.dimensionality
        if _is_volume_dimension(denominator_dim) or denominator_dim == ureg.mol.dimensionality:
            return value
        raise MissingAssayError("Calorific units need energy per mass, or matching density/molar_mass metadata.")

    def _standard_specific_values(self, value: pint.Quantity) -> tuple[pint.Quantity, ...]:
        """Energy densities using inherited scf-per-unit gas reference factors.

        This maps a supplied calorific density, rather than substituting the
        assay benchmark. Groningen retains its contractual energy parity.
        Without an inherited mapping, unrelated gas states stay unrelated.
        """
        from core.units.conversions import _as_decimal
        from core.units.registry import ureg
        if not self.standard_volume_factors or self.basis == "WET_GAS":
            return ()
        reference_heat = None
        if value.is_compatible_with(ureg.MJ / ureg.scf):
            reference_heat = value.to(ureg.MJ / ureg.scf).magnitude
        else:
            for unit, factor in self.standard_volume_factors:
                target = ureg.MJ / ureg.Unit(unit)
                if value.is_compatible_with(target):
                    reference_heat = value.to(target).magnitude / _as_decimal(factor)
                    break
        if reference_heat is None:
            return ()
        return tuple(ureg.Quantity(reference_heat * _as_decimal(factor), ureg.MJ / ureg.Unit(unit))
                     for unit, factor in self.standard_volume_factors)

    def _moisture(self, basis: str) -> Decimal:
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        from core.units.exceptions import MissingAssayError
        from core.units.exceptions import MissingReferenceConditionsError
        if basis in _DRY_SOLID_BASES:
            return Decimal(0)
        if basis in _CONDITION_FIELDS:
            state = self._condition(basis)
            if state is None:
                raise MissingReferenceConditionsError(f"{basis} requires an explicit {_CONDITION_FIELDS[basis]}.")
            return _fraction(state.moisture_pct, "condition moisture")
        if basis in _GAS_BASES:
            raise IncompatibleAnalyticalBasisError("Gas basis uses water mole fractions, not solid-fuel moisture.")
        value = self.moisture_pct if basis == self.basis else None
        if basis == "AR":
            value = self.as_received_moisture_pct if self.as_received_moisture_pct is not None else value
        if basis == "AD":
            value = self.air_dried_moisture_pct if self.air_dried_moisture_pct is not None else value
        if value is None:
            raise MissingAssayError(f"Moisture is required for {basis} analytical basis.")
        return _fraction(value, "moisture_pct")

    def _mass_per_dry_kg(self, basis):
        from core.units.analytical_basis import dry_matter_fraction
        from core.units.exceptions import MissingAssayError
        if basis == "DAF":
            if self.ash_db_pct is None:
                raise MissingAssayError("ash_db_pct is required for DAF conversion.")
            return 1 - _fraction(self.ash_db_pct, "ash_db_pct")
        if basis == "DMMF":
            if self.mineral_matter_db_pct is None:
                raise MissingAssayError("mineral_matter_db_pct is required for DMMF; ash is not mineral matter.")
            return 1 - _fraction(self.mineral_matter_db_pct, "mineral_matter_db_pct")
        return 1 / dry_matter_fraction(self._moisture(basis))

    def _mineral_correction(self, kind):
        correction = getattr(self, "mineral_matter_" + kind + "_correction")
        return Decimal(0) if correction is None else correction.magnitude

    def _latent_loss(self) -> Decimal:
        from core.units.conversions import _as_decimal
        from core.units.exceptions import MissingAssayError
        if self.hydrogen_pct is None:
            raise MissingAssayError("Missing opposite heating value: supply hydrogen_pct or an explicit NCV/GCV ratio.")
        # Gas spectator water remains vapour on both HHV and LHV conventions;
        # only chemically formed water contributes to their difference.
        moisture = Decimal(0) if self.basis in _GAS_BASES else self._moisture(self.basis)
        water = moisture + _as_decimal(self.water_per_kg_hydrogen) * _as_decimal(self.hydrogen_pct) / 100
        return _as_decimal(self.latent_heat_mj_per_kg_water) * water

    def _gas_molar_value(self, value, state):
        """MJ per total (wet or dry) mole at this metering state."""
        from core.units.registry import ureg
        if value.dimensionality == (ureg.MJ / ureg.kg).dimensionality:
            return (value.to("MJ/kg") * state.molar_mass).to("MJ/mol")
        if value.dimensionality == (ureg.MJ / ureg.mol).dimensionality:
            return value.to("MJ/mol")
        denominator = ureg.MJ / value.units
        if _is_volume_dimension(denominator.dimensionality):
            return (value * state.volume_per_mole(denominator)).to("MJ/mol")
        raise ValueError("Gas calorific value needs energy per mass, mole or metered volume.")

    def _gas_value_from_molar(self, value, state, unit):
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        target = parse_unit(unit)
        if target.dimensionality == (ureg.MJ / ureg.kg).dimensionality:
            return (value / state.molar_mass).to(target)
        if target.dimensionality == (ureg.MJ / ureg.mol).dimensionality:
            return value.to(target)
        denominator = ureg.MJ / target
        if _is_volume_dimension(denominator.dimensionality):
            return (value / state.volume_per_mole(denominator)).to(target)
        raise ValueError("Target gas calorific units must be energy per mass, mole or metered volume.")

    def _gas_basis_value(self, value, source_basis, target_basis, unit):
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        if source_basis not in _GAS_BASES or target_basis not in _GAS_BASES:
            raise IncompatibleAnalyticalBasisError("Wet/dry gas conversions must use gas bases on both sides; solid-fuel moisture is a separate basis system.")
        source = self._gas_state(source_basis)
        target = self._gas_state(target_basis)
        source_denominator = ureg.MJ / value.units
        target_unit = parse_unit(unit)
        target_denominator = ureg.MJ / target_unit
        if source_basis == target_basis and _is_volume_dimension(source_denominator.dimensionality) and _is_volume_dimension(target_denominator.dimensionality):
            # Same-state tagged/geometric volumes need no EOS evaluation: the
            # heating value was already measured per unit metered volume.
            from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
            source_volume = _geometric_volume(ureg.Quantity(Decimal(1), source_denominator), source.reference)
            target_volume = _geometric_volume(ureg.Quantity(Decimal(1), target_denominator), target.reference)
            return ureg.Quantity(value.magnitude / source_volume * target_volume, target_unit)
        per_dry_mole = self._gas_molar_value(value, source) / (1 - source.water_mole_fraction)
        per_target_mole = per_dry_mole * (1 - target.water_mole_fraction)
        return self._gas_value_from_molar(per_target_mole, target, unit)

    def calorific_value(self, kind: str = "ncv", *, basis: str | None = None) -> pint.Quantity:
        """Returns MJ/kg (or native volume/molar units); supplied pairs remain exact."""
        from core.units.conversions import _as_decimal
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        from core.units.exceptions import MissingAssayError
        from core.units.registry import ureg
        basis = _report_basis(kind, basis)
        kind = _kind(kind)
        value = getattr(self, kind)
        if value is None:
            other = self.gcv if kind == "ncv" else self.ncv
            if self.method == CalorificMethod.EXPLICIT_RATIO:
                ratio = _as_decimal(self.ncv_to_gcv_ratio)
                value = other * (ratio if kind == "ncv" else 1 / ratio)
            elif self.method == CalorificMethod.LATENT_WATER_APPROXIMATION:
                if self.thermodynamic_basis == ThermodynamicBasis.CONSTANT_VOLUME:
                    raise IncompatibleAnalyticalBasisError("Latent-water approximation does not perform constant-volume/pressure corrections; supply matching values or explicit reaction corrections.")
                if not other.is_compatible_with("MJ/kg"):
                    raise ValueError("Hydrogen/moisture estimation requires mass-specific calorific data.")
                loss = ureg.Quantity(self._latent_loss(), "MJ/kg")
                value = other - loss if kind == "ncv" else other + loss
            else:
                raise MissingAssayError("Missing calorific value: select method='latent_water' or method='explicit_ratio', or supply both measured values. No approximation is selected implicitly.")
            _number(value.magnitude, "Derived calorific value", positive=True)
        target = self.basis if basis is None else _basis(basis)
        if target == self.basis:
            return value
        if self.basis in _GAS_BASES or target in _GAS_BASES:
            if value.dimensionality == (ureg.MJ / ureg.kg).dimensionality:
                unit = "MJ/kg"
            elif value.dimensionality == (ureg.MJ / ureg.mol).dimensionality:
                unit = "MJ/mol"
            else:
                # Wet gas cannot be reported against a dry-only volume tag.
                unit = "MJ/m^3"
            return self._gas_basis_value(value, self.basis, target, unit)
        if not value.is_compatible_with("MJ/kg"):
            raise ValueError("Analytical basis changes require mass-specific calorific data.")
        # Remove sample-water evaporation loss before scaling dry matter;
        # add the target moisture loss afterwards. This also preserves a
        # supplied measured pair without replacing it with an approximation.
        source_m = self._moisture(self.basis)
        target_m = self._moisture(target)
        if kind == "ncv" and source_m != target_m:
            self._require_moisture_correction()
        latent = _as_decimal(self.latent_heat_mj_per_kg_water)
        magnitude = value.magnitude
        if kind == "ncv":
            magnitude += latent * source_m
        # Resolve through one kg of the original dry whole material. DAF and
        # DMMF remove distinct measured masses; neither substitutes for the other.
        if self.basis in {"DAF", "DMMF"}:
            magnitude *= self._mass_per_dry_kg(self.basis)
        else:
            magnitude /= 1 - source_m
        if self.basis == "DMMF":
            magnitude += self._mineral_correction(kind)
        if target == "DMMF":
            magnitude -= self._mineral_correction(kind)
        if target in {"DAF", "DMMF"}:
            magnitude /= self._mass_per_dry_kg(target)
        else:
            magnitude *= 1 - target_m
        if kind == "ncv":
            magnitude -= latent * target_m
        _number(magnitude, "Converted calorific value", positive=True)
        return ureg.Quantity(magnitude, "MJ/kg")

    def energy_density(self, unit: str | pint.Unit = "MJ/kg", *, kind: str = "ncv", basis: str | None = None) -> pint.Quantity:
        """Report fuel heat, preserving an explicitly thermal energy denomination."""
        from core.units.dimensions import physical_unit
        from core.units.dimensions import restore
        from core.units.registry import parse_unit
        target = parse_unit(unit)
        physical = physical_unit(target, role="thermal")
        result = self._energy_density(physical, kind=kind, basis=basis)
        return restore(result, target, role="thermal")

    def _energy_density(self, unit: str | pint.Unit = "MJ/kg", *, kind: str = "ncv", basis: str | None = None) -> pint.Quantity:
        """Converts mass, molar or reference-volume energy density explicitly.

        Density metadata applies only to the declared basis; basis-normalized
        mass values need a new density before volumetric conversion.
        """
        from core.units.conversions import _as_decimal
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        from core.units.exceptions import MissingAssayError
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        basis = _report_basis(kind, basis)
        target = parse_unit(unit)
        target_basis = self.basis if basis is None else _basis(basis)
        value = self.calorific_value(kind, basis=basis)
        if value.is_compatible_with(target):
            return value.to(target)
        if self.basis in _GAS_BASES or target_basis in _GAS_BASES:
            if target_basis == "DRY_GAS" and self.standard_volume_factors:
                # Explicit legacy dry-gas reference parity survives moisture
                # projections; it never labels wet gas as a dry standard volume.
                from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
                state = self._gas_state("DRY_GAS")
                for tagged, _ in self.standard_volume_factors:
                    try:
                        _geometric_volume(ureg.Quantity(Decimal(1), tagged), state.reference)
                    except ValueError:
                        continue
                    reference_value = self._gas_basis_value(self.calorific_value(kind), self.basis, "DRY_GAS", ureg.MJ / ureg.Unit(tagged))
                    # The map is defined on the dry profile even when this
                    # object currently represents its wet counterpart.
                    for unit_name, factor in self.standard_volume_factors:
                        if reference_value.is_compatible_with(ureg.MJ / ureg.Unit(unit_name)):
                            per_scf = reference_value.to(ureg.MJ / ureg.Unit(unit_name)).magnitude / _as_decimal(factor)
                            for destination, destination_factor in self.standard_volume_factors:
                                candidate = ureg.Quantity(per_scf * _as_decimal(destination_factor), ureg.MJ / ureg.Unit(destination))
                                if candidate.is_compatible_with(target):
                                    return candidate.to(target)
                            break
                    break
            return self._gas_basis_value(self.calorific_value(kind), self.basis, target_basis, target)
        if target_basis != self.basis:
            condition = self._condition(target_basis)
            if condition is not None and condition.density is not None:
                candidate = value * condition.density
                if candidate.is_compatible_with(target):
                    return candidate.to(target)
            raise IncompatibleAnalyticalBasisError("Volumetric/molar conversion after a basis change requires new density/molar metadata at the target condition.")
        candidates = [value * density for density in self._densities()]
        for candidate in (value, *candidates):
            mapped = self._standard_specific_values(candidate)
            if mapped:
                candidates.extend(mapped)
                break
        if self.molar_mass is not None:
            candidates.append(value * self.molar_mass)
        for candidate in candidates:
            if candidate.is_compatible_with(target):
                return candidate.to(target)
        raise MissingAssayError("Target energy-density units require matching density or molar_mass metadata.")

    def as_context(self, *, basis: str = "lhv", name: str | None = None) -> pint.Context:
        return create_heat_content_context(self, basis=basis, name=name)

    def for_basis(self, basis: str, *, density=None, molar_mass=None, kind: str | None = None):
        """Return an independent spec for quantities measured on the chosen basis.

        A source density is never reused after a basis change. Named target
        conditions can supply their measured density; gas densities require
        the target P/T/Z and mixture molar mass. Known HHV/LHV values are
        rebased separately; a missing value is resolved only if its original
        hydrogen/ratio data can resolve it before the change.
        Select kind to project only HHV or LHV, without requesting an unrelated
        counterpart or its moisture correction.
        """
        from core.units.registry import ureg
        from dataclasses import replace
        requested = None if kind is None else _kind(kind)
        target = _basis(basis if kind is None else _report_basis(kind, basis))
        if target == self.basis:
            updates = {}
            if density is not None:
                updates["density"] = density
            if molar_mass is not None:
                updates["molar_mass"] = molar_mass
            if requested is not None:
                updates[requested] = self.calorific_value(requested)
                updates[requested + "_method"] = self.calorific_method(requested)
                updates["ncv" if requested == "gcv" else "gcv"] = None
            return replace(self, **updates)
        values = {}
        can_resolve = (self.method == CalorificMethod.LATENT_WATER_APPROXIMATION and self.hydrogen_pct is not None
                       or self.method == CalorificMethod.EXPLICIT_RATIO)
        for kind in ("gcv", "ncv"):
            selected = requested is None or requested == kind
            values[kind] = self.calorific_value(kind, basis=target) if selected and (getattr(self, kind) is not None or can_resolve) else None
        condition = self._condition(target)
        temperature = condition.temperature if condition is not None else None
        pressure = condition.absolute_pressure if condition is not None else None
        if density is None and condition is not None:
            density = condition.density
        standard_density = None
        if target in _GAS_BASES:
            state = self._gas_state(target)
            temperature, pressure = state.reference.temperature, state.reference.absolute_pressure
            if molar_mass is None and state.dry_gas_molar_mass is not None:
                molar_mass = state.molar_mass
            if density is None and state.dry_gas_molar_mass is not None and state.reference.compressibility_factor is not None:
                density = state.density
            if state.water_mole_fraction == 0 and state.dry_gas_molar_mass is not None and state.reference.compressibility_factor is not None:
                from core.units.conversions import _GAS_REFERENCE_PRESETS
                from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
                for tagged, _, _, _ in _GAS_REFERENCE_PRESETS.values():
                    try:
                        volume = _geometric_volume(ureg.Quantity(Decimal(1), tagged), state.reference)
                    except ValueError:
                        continue
                    standard_density = ureg.Quantity(state.density.magnitude * volume, ureg.kg / ureg.Unit(tagged))
                    break
            moisture = None
        else:
            moisture = self._moisture(target) * 100
        retained_moistures = {}
        if self.basis == "AR" and self.as_received_moisture_pct is None:
            retained_moistures["as_received_moisture_pct"] = self._moisture("AR") * 100
        if self.basis == "AD" and self.air_dried_moisture_pct is None:
            retained_moistures["air_dried_moisture_pct"] = self._moisture("AD") * 100
        return replace(self, **values, **retained_moistures, basis=target, moisture_pct=moisture, density=density,
                       standard_density=standard_density, molar_mass=molar_mass,
                       reference_temperature=temperature, reference_pressure=pressure,
                       hydrogen_pct=None, standard_volume_factors=self.standard_volume_factors if self.basis in _GAS_BASES and target in _GAS_BASES else (),
                       gcv_reference_temperature=None, ncv_reference_temperature=None,
                       gcv_reference_pressure=None, ncv_reference_pressure=None,
                       ultimate_analysis=None, proximate_analysis=None,
                       gcv_method=self.calorific_method("gcv") if values["gcv"] is not None else None,
                       ncv_method=self.calorific_method("ncv") if values["ncv"] is not None else None,
                       derivations=self.derivations + (f"Analytical basis {self.basis} -> {target}",))

    def _require_moisture_correction(self):
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        from core.units.exceptions import MissingAssayError
        if self.thermodynamic_basis == ThermodynamicBasis.CONSTANT_VOLUME:
            raise IncompatibleAnalyticalBasisError("NCV latent-water moisture correction requires constant-pressure data or a prior explicit reaction PV correction.")
        if self.method != CalorificMethod.LATENT_WATER_APPROXIMATION and self.moisture_correction_method != CalorificMethod.LATENT_WATER_APPROXIMATION:
            raise MissingAssayError("NCV moisture rebasing requires explicit moisture_correction_method='latent_water' (or method='latent_water').")

    def calorific_method(self, kind="ncv"):
        """Provenance for a supplied or explicitly derived calorific value."""
        from core.units.exceptions import MissingAssayError
        key = _kind(kind)
        declared = getattr(self, key + "_method")
        if declared is not None and getattr(self, key) is not None:
            return declared
        if getattr(self, key) is not None:
            return self.method if self.method in {CalorificMethod.COMPONENT_SUM, CalorificMethod.EMPIRICAL_CORRELATION, CalorificMethod.STANDARD_LABORATORY} else CalorificMethod.MEASURED
        if self.method not in {CalorificMethod.LATENT_WATER_APPROXIMATION, CalorificMethod.EXPLICIT_RATIO}:
            raise MissingAssayError("Missing value has no selected calculation method.")
        return self.method

    def for_thermodynamic_basis(self, target, *, gcv_pv_correction=None, ncv_pv_correction=None):
        """Project Qp/Qv with explicit signed reaction delta(PV) corrections."""
        from dataclasses import replace
        target = ThermodynamicBasis(target)
        values = {}
        for kind, correction in (("gcv", gcv_pv_correction), ("ncv", ncv_pv_correction)):
            value = getattr(self, kind)
            if value is None and correction is not None:
                raise ValueError(f"A {kind} PV correction requires its calorific value.")
            if value is not None:
                if correction is None:
                    raise ValueError(f"Supply {kind}_pv_correction for its product-water phase.")
                values[kind] = convert_combustion_heat(value, source=self.thermodynamic_basis, target=target,
                                                     pressure_volume_correction=correction)
        return replace(self, **values, thermodynamic_basis=target,
                       method=CalorificMethod.MEASURED, latent_heat_method="fixed",
                       derivations=self.derivations + (f"Reaction PV correction: {self.thermodynamic_basis} -> {target}",))

    @classmethod
    def from_laboratory(cls, result: LaboratoryCalorificResult, **metadata):
        if not isinstance(result, LaboratoryCalorificResult):
            raise TypeError("Supply LaboratoryCalorificResult with report provenance.")
        data = dict(gcv=result.gcv, ncv=result.ncv, basis=result.basis,
                    thermodynamic_basis=result.thermodynamic_basis, combustion_temperature=result.combustion_temperature,
                    method=CalorificMethod.STANDARD_LABORATORY, method_source=result.source,
                    laboratory_standard=result.standard, laboratory_report_id=result.report_id)
        if set(data) & metadata.keys():
            raise ValueError("Laboratory result identity cannot be overridden; transform the resulting assay explicitly.")
        return cls(**data, **metadata)

    @classmethod
    def from_ultimate_analysis(cls, analysis, *, correlation, **metadata):
        if not isinstance(analysis, UltimateAnalysis) or not isinstance(correlation, EmpiricalHeatCorrelation):
            raise TypeError("Supply UltimateAnalysis and an explicit EmpiricalHeatCorrelation.")
        return cls(gcv=correlation.evaluate(analysis), basis=analysis.basis, ultimate_analysis=analysis,
                   method=CalorificMethod.EMPIRICAL_CORRELATION, method_source=f"{correlation.name}: {correlation.source}",
                   thermodynamic_basis=correlation.thermodynamic_basis, combustion_temperature=correlation.combustion_temperature, **metadata)

    @classmethod
    def from_gas_composition(cls, composition, *, dry_reference, wet_reference, name="composition gas", basis="DRY_GAS", gas_constant=None, component_table=None):
        """Derive both heating values from dry composition, then add spectator water.

        Reference states carry explicitly selected Z. The optional gas constant
        records the chosen EOS convention; the default is CODATA 2019.
        """
        from core.units.composition import GasComposition
        from core.units.reference_conditions import GasReferenceConditions
        if not isinstance(composition, GasComposition) or not isinstance(dry_reference, GasReferenceConditions) or not isinstance(wet_reference, GasReferenceConditions):
            raise TypeError("Supply GasComposition and explicit dry/wet GasReferenceConditions.")
        if composition.water_mole_fraction and composition.water_mole_fraction != wet_reference.water_mole_fraction:
            raise ValueError("Wet composition conflicts with wet_reference water fraction.")
        if wet_reference.water_mole_fraction is None:
            raise ValueError("wet_reference requires an explicit water mole fraction.")
        dry = composition.dry()
        wet = dry.with_water(wet_reference.water_mole_fraction)
        options = {} if gas_constant is None else {"gas_constant": gas_constant}
        dry_state = GasMoistureState.from_composition(dry, dry_reference, **options)
        wet_state = GasMoistureState.from_composition(wet, wet_reference, **options)
        if component_table is not None and not isinstance(component_table, ISO6976Table):
            raise TypeError("component_table must be ISO6976Table.")
        heating = dry.calorific_values() if component_table is None else component_table.calculate(dry)
        specification = cls(name=name, gcv=heating.hhv, ncv=heating.lhv, basis="DRY_GAS",
                            dry_gas_state=dry_state, wet_gas_state=wet_state,
                            method=heating.method, method_source=heating.source,
                            thermodynamic_basis=heating.thermodynamic_basis, combustion_temperature=heating.combustion_temperature)
        target = _basis(basis)
        if target not in _GAS_BASES:
            raise ValueError("Composition-derived gas needs DRY_GAS or WET_GAS basis.")
        return specification if target == "DRY_GAS" else specification.for_basis(target)

    @classmethod
    def from_assay(cls, assay: Any, **overrides: Any) -> HeatContentAssay:
        return heat_content_from_assay(assay, **overrides)





def heat_content_from_assay(assay: Any, **overrides: Any) -> HeatContentAssay:
    """Adapts fuel assays, or adds explicit calorific data to any other assay.

    No heat is inferred for crude/products, metals, ores, crops or nuclear fuel. Their callers
    supply gcv/ncv quantities. Existing carrier recovery and measured/estimated
    pairs are preserved. Overrides are validated by HeatContentAssay.
    """
    from core.units.conversions import _as_decimal
    from core.units.registry import ureg
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.assays.energy.coal import CoalAssay
    from core.units.assays.energy.biofuels import CleanFuelAssay
    from core.units.assays.energy.biofuels import BiofuelAssay
    from core.units.assays.energy.ngl import NGLAssay
    from core.units.assays.energy.biomass import BiomassAssay
    from core.units.assays.energy.petcoke import PetcokeAssay
    from core.units.assays.energy.hydrogen import HydrogenAssay
    from core.units.reference_conditions import GasReferenceConditions
    if isinstance(assay, HeatContentAssay):
        from dataclasses import replace
        return replace(assay, **overrides)
    adapter = getattr(assay, "as_heat_content", None)
    if callable(adapter):
        from dataclasses import replace
        import inspect
        accepts_overrides = any(p.kind == p.VAR_KEYWORD for p in inspect.signature(adapter).parameters.values())
        specification = adapter(**overrides) if accepts_overrides else adapter()
        if not isinstance(specification, HeatContentAssay):
            raise TypeError("assay.as_heat_content() must return HeatContentAssay.")
        return specification if accepts_overrides else replace(specification, **overrides)
    data: dict[str, Any] = {"name": getattr(assay, "name", type(assay).__name__)}
    if isinstance(assay, GasAssay):
        data.update(gcv=ureg.Quantity(assay.gcv_btu_per_scf, "btu_it/scf"),
                    ncv=ureg.Quantity(assay.ncv_btu_per_scf, "btu_it/scf"))
        # GasAssay defines scf at 60 F / 14.73 psia. Its conversion methods
        # carry this reference correction; derive a matching mass density.
        if assay.density_kg_per_scf is not None:
            data["standard_density"] = ureg.Quantity(assay.density_kg_per_scf, "kg/scf")
        # Use the same canonical reference constants as the original gas
        # context, rather than round-tripping its rounded psi label.
        reference = GasReferenceConditions.from_standard("ferc")
        data.update(reference_temperature=reference.temperature, reference_pressure=reference.absolute_pressure)
        factors = []
        with ureg.context(assay.as_context(basis="gcv")):
            for unit in ("scf", "scf_tx", "scf_60F_14696psi_dry", "Sm3", "Nm3", "Sm3_20C", "m3_Groningen"):
                factors.append((unit, _as_decimal(ureg.Quantity(Decimal(1), unit).to("scf").magnitude)))
        data["standard_volume_factors"] = tuple(factors)
    elif isinstance(assay, CoalAssay):
        data.update(gcv=ureg.Quantity(assay.gj_per_mt_gar, "MJ/kg"), ncv=ureg.Quantity(assay.gj_per_mt_nar, "MJ/kg"),
                    moisture_pct=assay.total_moisture_pct, air_dried_moisture_pct=assay.inherent_moisture_pct,
                    hydrogen_pct=assay.hydrogen_ar_pct)
        if assay.ash_pct is not None:
            data["ash_db_pct"] = _as_decimal(assay.ash_pct) / (1 - _as_decimal(assay.total_moisture_pct) / 100)
    elif isinstance(assay, (CleanFuelAssay, BiofuelAssay, NGLAssay, BiomassAssay, PetcokeAssay, HydrogenAssay)):
        mass = ureg.Quantity(Decimal(1), "kg")
        data.update(gcv=assay.mass_to_energy(mass, basis="hhv") / mass,
                    ncv=assay.mass_to_energy(mass, basis="lhv") / mass)
        if isinstance(assay, (BiomassAssay, PetcokeAssay)):
            data.update(moisture_pct=assay.moisture_pct, ash_db_pct=assay.ash_db_pct,
                        density=ureg.Quantity(assay.bulk_density_kg_per_m3, "kg/m^3"))
        if isinstance(assay, BiomassAssay):
            data["hydrogen_pct"] = _as_decimal(assay.hydrogen_db_pct) * (1 - _as_decimal(assay.moisture_pct) / 100)
        if isinstance(assay, CleanFuelAssay):
            if assay.liquid_density_kg_per_m3 is not None:
                data["density"] = ureg.Quantity(assay.liquid_density_kg_per_m3, "kg/m^3")
            if assay.density_kg_per_nm3 is not None:
                data["standard_density"] = ureg.Quantity(assay.density_kg_per_nm3, "kg/Nm3")
        if isinstance(assay, (BiofuelAssay, NGLAssay, HydrogenAssay)):
            data["density"] = ureg.Quantity(assay.density_kg_per_m3, "kg/m^3")
        if isinstance(assay, NGLAssay) and assay.molecular_weight_g_per_mol is not None:
            data["molar_mass"] = ureg.Quantity(assay.molecular_weight_g_per_mol, "g/mol")
        if isinstance(assay, HydrogenAssay):
            data.update(reference_temperature=ureg.Quantity(assay.reference_temp_c, "degC"),
                        reference_pressure=ureg.Quantity(assay.reference_pressure_bar, "bar"))
            if assay.state == "carrier":
                data["energy_role"] = "recovered_hydrogen"
    else:
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialProperties
        from core.units.assays.base import MaterialSpec
        if isinstance(assay, (MaterialProperties, MaterialSpec)):
            props = assay.properties if isinstance(assay, MaterialSpec) else assay
            comp = assay.composition if isinstance(assay, MaterialSpec) else None
            dens = props.effective_true_density or props.effective_bulk_density
            if dens is not None:
                data["density"] = dens
            if props.molar_mass is not None:
                data["molar_mass"] = props.molar_mass
            moist = props.moisture if props.moisture is not None else (comp.moisture_fraction if comp is not None else None)
            if moist is not None:
                data["moisture_pct"] = moist * Decimal(100)
            if comp is not None and comp.ash_fraction is not None:
                if comp.basis in {CompositionBasis.DRY, CompositionBasis.DRY_ASH_FREE}:
                    data["ash_db_pct"] = comp.ash_fraction * Decimal(100)
                elif comp.basis in {CompositionBasis.AS_RECEIVED, CompositionBasis.AIR_DRIED} and moist is not None:
                    data["ash_db_pct"] = (comp.ash_fraction / (Decimal(1) - moist)) * Decimal(100)
    data.update(overrides)
    return HeatContentAssay(**data)





def gcv_to_ncv(gcv: pint.Quantity, **metadata: Any) -> pint.Quantity:
    """Returns supplied NCV or an explicit hydrogen/moisture/ratio estimate."""
    return _input_type(HeatContentAssay(gcv=gcv, **metadata).energy_density(gcv.units, kind="ncv"), gcv)





def ncv_to_gcv(ncv: pint.Quantity, **metadata: Any) -> pint.Quantity:
    return _input_type(HeatContentAssay(ncv=ncv, **metadata).energy_density(ncv.units, kind="gcv"), ncv)





def as_received_to_dry(
    value: pint.Quantity,
    moisture_pct: Any,
    *,
    kind: str = "gcv",
    latent_heat_mj_per_kg_water: Decimal | float = Decimal("2.442"),
    method: CalorificMethod = CalorificMethod.MEASURED,
) -> pint.Quantity:
    data = {_kind(kind): value}
    spec = HeatContentAssay(**data, moisture_pct=moisture_pct, latent_heat_mj_per_kg_water=latent_heat_mj_per_kg_water, method=method)
    return _input_type(spec.energy_density(value.units, kind=kind, basis="DB"), value)





def dry_to_as_received(
    value: pint.Quantity,
    moisture_pct: Any,
    *,
    kind: str = "gcv",
    latent_heat_mj_per_kg_water: Decimal | float = Decimal("2.442"),
    method: CalorificMethod = CalorificMethod.MEASURED,
) -> pint.Quantity:
    """Reverse dry-basis heating values using an explicit sample-water latent heat."""
    from core.units.exceptions import MissingAssayError
    from core.units.registry import ureg
    latent_heat = _number(latent_heat_mj_per_kg_water, "latent_heat_mj_per_kg_water", positive=True)
    m = _fraction(moisture_pct, "moisture_pct")
    dry = _quantity(value, "calorific value").to("MJ/kg")
    result = dry * (1 - m)
    if _kind(kind) == "ncv":
        if CalorificMethod(method) != CalorificMethod.LATENT_WATER_APPROXIMATION:
            raise MissingAssayError("NCV moisture rebasing requires method='latent_water'.")
        result -= ureg.Quantity(latent_heat * m, "MJ/kg")
    _number(result.magnitude, "Converted calorific value", positive=True)
    return _input_type(result.to(value.units), value)





def create_heat_content_context(assay: HeatContentAssay, *, basis: str = "lhv", name: str | None = None) -> pint.Context:
    """Bridges energy, mass, molar and explicitly declared volume states."""
    from core.units.conversions import _as_decimal
    from core.units.conversions import _scale_factor
    from core.units.registry import ureg
    if _report_basis(basis, None) == "AR" and assay.basis != "AR":
        assay = assay.for_basis("AR", kind=basis)
    value = assay.calorific_value(basis)
    if not value.is_compatible_with("MJ/kg"):
        denominator = ureg.MJ / value.units
        nodes = [(denominator, value.magnitude), (ureg.MJ, Decimal(1))]
        heat = None
    else:
        heat = value.to("MJ/kg").magnitude
        nodes = [(ureg.kg, heat), (ureg.MJ, Decimal(1))]
    for density in assay._densities() if heat is not None else ():
        # Retain the density's declared volume unit (including standard tags).
        volume = ureg.kg / density.units
        rho = density.to(ureg.kg / volume).magnitude
        nodes.append((volume, heat * rho))
    if assay.molar_mass is not None and heat is not None:
        nodes.append((ureg.mol, heat * assay.molar_mass.to("kg/mol").magnitude))
    if assay.standard_volume_factors and assay.basis != "WET_GAS":
        reference_heat = assay.energy_density("MJ/scf", kind=basis).magnitude
        present = {unit.dimensionality for unit, _ in nodes}
        for unit, factor in assay.standard_volume_factors:
            parsed = ureg.Unit(unit)
            if parsed.dimensionality not in present:
                nodes.append((parsed, reference_heat * _as_decimal(factor)))
                present.add(parsed.dimensionality)
    if assay.basis in _GAS_BASES:
        state = assay._gas_state(assay.basis)
        native = ureg.MJ / value.units
        known_mass = native.dimensionality == ureg.kg.dimensionality
        known_molar = native.dimensionality == ureg.mol.dimensionality
        known_volume = _is_volume_dimension(native.dimensionality)
        has_mass = state.dry_gas_molar_mass is not None
        has_z = state.reference.compressibility_factor is not None
        can_molar = known_molar or known_mass and has_mass or known_volume and has_z
        can_mass = known_mass or can_molar and has_mass
        can_volume = known_volume or can_molar and has_z
        # Gas-state transformations use that state's mixture molar mass and Z,
        # even when only a subset of mass/molar/volume bridges is identifiable.
        nodes = [(native, value.magnitude), (ureg.MJ, Decimal(1))]
        present = {unit.dimensionality for unit, _ in nodes}

        def add(unit):
            from core.units.registry import ureg
            if unit.dimensionality not in present:
                factor = assay.energy_density(ureg.MJ / unit, kind=basis).magnitude
                nodes.append((unit, factor))
                present.add(unit.dimensionality)

        if can_mass:
            add(ureg.kg)
        if can_molar:
            add(ureg.mol)
        if can_volume:
            add(ureg.meter ** 3)
            from core.units.conversions import _GAS_REFERENCE_PRESETS
            from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
            for tagged, _, _, _ in _GAS_REFERENCE_PRESETS.values():
                try:
                    _geometric_volume(ureg.Quantity(Decimal(1), tagged), state.reference)
                except ValueError:
                    continue
                add(ureg.Unit(tagged))
        if assay.basis == "DRY_GAS" and assay.standard_volume_factors:
            reference_heat = assay.energy_density("MJ/scf", kind=basis).magnitude
            for tagged, factor in assay.standard_volume_factors:
                unit = ureg.Unit(tagged)
                if unit.dimensionality not in present:
                    nodes.append((unit, reference_heat * _as_decimal(factor)))
                    present.add(unit.dimensionality)
    # A declared fuel context can express its heat in thermal units, never
    # electrical output. Semantic markers still isolate electricity.
    nodes.append((ureg.megajoule_thermal, Decimal(1)))
    context = pint.Context(name or f"heat_content_{assay.name}_{_kind(basis)}")

    def transform(source: pint.Unit, destination: pint.Unit, factor: Decimal):
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        def apply(_registry: pint.UnitRegistry, quantity: pint.Quantity) -> pint.Quantity:
            from core.units.conversions import _scale_factor
            from core.units.registry import ureg
            q = quantity.to(source)
            return ureg.Quantity(q.magnitude * _scale_factor(q.magnitude, factor), destination)
        return apply

    for (source, source_heat), (destination, destination_heat) in permutations(nodes, 2):
        context.add_transformation(source.dimensionality, destination.dimensionality,
                                   transform(source, destination, source_heat / destination_heat))
    return context


from core.units.validation import authorized_conversion


@authorized_conversion
def convert_heat_content(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: Any,
    basis: str = "lhv",
    analytical_basis: str | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
    **metadata: Any,
) -> Any:
    """Converts quantities, rates, quotes or a specific calorific value.

    Commodity assays are adapted automatically; nonfuel assays need explicit
    gcv/ncv metadata. ``basis`` means HHV/LHV; ``analytical_basis`` means AR/DB.
    Basis normalization of calorific densities is explicit, never a hidden
    conversion of two identically dimensioned mass units.
    """
    from core.units.conversions import _as_decimal
    from core.units.dimensions import physical_quantity
    from core.units.dimensions import physical_unit
    from core.units.dimensions import restore
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError("Expected a Pint quantity.")
    validate_registry_quantity(quantity)
    original_quantity = quantity
    original_target = parse_unit(target_unit)
    quantity = physical_quantity(quantity, role="thermal")
    target_unit = physical_unit(original_target, role="thermal")
    if not _as_decimal(quantity.magnitude).is_finite():
        raise ValueError("Heat-content conversion input must be finite.")
    spec = heat_content_from_assay(assay, **metadata)
    if strict:
        sel_method = spec.calorific_method(_kind(basis))
        if sel_method in {CalorificMethod.LATENT_WATER_APPROXIMATION, CalorificMethod.EXPLICIT_RATIO, CalorificMethod.EMPIRICAL_CORRELATION}:
            raise StandardsApplicabilityError(
                f"strict=True prohibits estimated calorific method {sel_method.value!r}; supply measured/tabulated calorific data."
            )
    target = parse_unit(target_unit)
    # A calorific-density input represents the selected heating value on the
    # declared analytical basis. Preserve the supplied value rather than
    # silently substituting the assay's benchmark.
    density_dims = [ureg.MJ.dimensionality / unit.dimensionality for unit, _ in _heat_nodes(spec)]
    if quantity.dimensionality in density_dims:
        from dataclasses import replace
        supplied = replace(spec, gcv=quantity if _kind(basis) == "gcv" else None,
                           ncv=quantity if _kind(basis) == "ncv" else None)
        result = supplied.energy_density(target, kind=basis, basis=analytical_basis)
        converted = _input_type(restore(result, original_target, role="thermal"), original_quantity)
    else:
        if analytical_basis is not None and _basis(analytical_basis) != spec.basis:
            raise ValueError("Use calorific-density quantities for analytical basis changes, or assay.for_basis(...) before converting amounts.")
        # Execute the thermal context here rather than calling the master router:
        # conversions.convert delegates energy crossings to this function.
        from core.units.conversions import _convert_contextual_quantity
        context = create_heat_content_context(spec, basis=basis)
        with ureg.context(context):
            result = _convert_contextual_quantity(quantity, target, [context])
        converted = _input_type(restore(result, original_target, role="thermal"), original_quantity)

    if not return_metadata:
        return converted
    from core.units.provenance import build_conversion_result
    return build_conversion_result(
        converted,
        original_quantity,
        assay=spec,
        heating_basis=_kind(basis),
        analytical_basis=analytical_basis or spec.basis,
        conversion_path=(str(original_quantity.units), f"heat_content:{_kind(basis)}", str(converted.units)),
        method="convert_heat_content",
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
    )





def _heat_nodes(assay: HeatContentAssay) -> list[tuple[pint.Unit, None]]:
    from core.units.registry import ureg
    nodes = [(ureg.kg, None)]
    nodes.extend((ureg.kg / density.units, None) for density in assay._densities())
    if assay.molar_mass is not None:
        nodes.append((ureg.mol, None))
    if assay.basis != "WET_GAS":
        nodes.extend((ureg.Unit(unit), None) for unit, _ in assay.standard_volume_factors)
    for value in (assay.gcv, assay.ncv):
        if value is not None:
            nodes.append((ureg.MJ / value.units, None))
    if assay.basis in _GAS_BASES:
        nodes.append((ureg.meter ** 3, None))
        from core.units.conversions import _GAS_REFERENCE_PRESETS
        from core.units.assays.energy.lng import _converters_gas_geometric_volume as _geometric_volume
        for tagged, _, _, _ in _GAS_REFERENCE_PRESETS.values():
            try:
                _geometric_volume(ureg.Quantity(Decimal(1), tagged), assay._gas_state(assay.basis).reference)
            except ValueError:
                continue
            nodes.append((ureg.Unit(tagged), None))
    return nodes





def _efficiency(value: Any) -> Decimal:
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if isinstance(value, pint.Quantity):
        validate_registry_quantity(value)
        value = value.to(ureg.dimensionless).magnitude
    result = _number(value, "Efficiency", positive=True)
    if result > 1:
        raise ValueError("Efficiency must be in (0, 1]; use a percent quantity for percentages.")
    return result





def efficiency_from_heat_rate(heat_rate: pint.Quantity) -> Decimal:
    """Heat rate is fuel energy / output energy, e.g. Btu_it/kWh."""
    from core.units.dimensions import physical_heat_rate
    rate = physical_heat_rate(heat_rate).magnitude
    return _efficiency(1 / rate)





def heat_rate_from_efficiency(efficiency: Any, unit: str | pint.Unit = "heat_rate") -> pint.Quantity:
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.dimensions import roles
    target = parse_unit(unit)
    signature = roles(target)
    if signature and signature != {"thermal": 1, "electrical": -1}:
        raise ValueError("Heat rate must be thermal input per electrical output.")
    unit = ureg.thermal_energy_kind / ureg.electrical_energy_kind if signature else ureg.dimensionless
    return ureg.Quantity(1 / _efficiency(efficiency), unit).to(target)





def thermal_to_output(energy: pint.Quantity, *, efficiency: Any = None,
                      heat_rate: pint.Quantity | None = None, target_unit: str | pint.Unit = "MWh_e") -> pint.Quantity:
    """Fuel thermal energy to output; use exactly one efficiency or heat rate."""
    return _thermal_parity(energy, efficiency, heat_rate, target_unit, inverse=False)





def output_to_thermal(energy: pint.Quantity, *, efficiency: Any = None,
                      heat_rate: pint.Quantity | None = None, target_unit: str | pint.Unit = "MWh_th") -> pint.Quantity:
    return _thermal_parity(energy, efficiency, heat_rate, target_unit, inverse=True)





def _thermal_parity(energy: pint.Quantity, efficiency: Any, heat_rate: pint.Quantity | None,
                    target_unit: str | pint.Unit, *, inverse: bool) -> pint.Quantity:
    from core.units.conversions import _scale_factor
    from core.units.dimensions import physical_quantity
    from core.units.dimensions import restore
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(energy, pint.Quantity):
        raise TypeError("Expected a Pint quantity for energy.")
    validate_registry_quantity(energy)
    source_role, target_role = ("electrical", "thermal") if inverse else ("thermal", "electrical")
    original_energy = energy
    energy = physical_quantity(energy, role=source_role)
    if (efficiency is None) == (heat_rate is None):
        raise ValueError("Supply exactly one efficiency or heat_rate.")
    eta = _efficiency(efficiency) if heat_rate is None else efficiency_from_heat_rate(heat_rate)
    if energy.dimensionality not in (ureg.joule.dimensionality, ureg.watt.dimensionality):
        raise ValueError("Thermal parity requires energy or power units.")
    _number(energy.magnitude, "Energy")
    result = ureg.Quantity(energy.magnitude * _scale_factor(energy.magnitude, 1 / eta if inverse else eta), energy.units)
    return _input_type(restore(result, parse_unit(target_unit), role=target_role), original_energy)




_converters_energy_DATA_FILE: Final[Path] = (Path(__file__).resolve().parents[0] / "reference_data") / "units" / "energy_specs.json"



_converters_energy_BTU_IT_JOULES: Final[Decimal] = Decimal("1055.05585262")



_converters_energy_P_GPA_14696_KPA: Final[Decimal] = Decimal("101.325")




HeatingBasis = Literal["LHV", "HHV"]



GallonEquivalentKind = Literal["GGE", "DGE"]





def _converters_energy_normalize_basis(basis: str) -> HeatingBasis:
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    if not isinstance(basis, str):
        raise TypeError(f"basis must be 'LHV' or 'HHV', got {type(basis).__name__}.")
    key = basis.strip().upper()
    if key in {"LHV", "NCV", "NET"}:
        return "LHV"
    if key in {"HHV", "GCV", "GROSS"}:
        return "HHV"
    raise IncompatibleAnalyticalBasisError(f"Heating-value basis must be 'LHV' or 'HHV', got {basis!r}.")





@dataclass(frozen=True, slots=True)
class EnergySpec:
    """Named or custom energy-equivalent gallon specification (GGE or DGE)."""

    name: str
    energy_per_gallon_btu: float
    basis: HeatingBasis
    source: str
    fuel_type: GallonEquivalentKind = "GGE"

    def __post_init__(self) -> None:
        from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
        if not isinstance(self.name, str) or not self.name.strip():
            raise ValueError("EnergySpec.name must be a non-empty string.")
        if not isinstance(self.source, str) or not self.source.strip():
            raise ValueError("EnergySpec.source must be a non-empty string.")
        btu_dec = _converters__interpolation_validate_real_scalar(self.energy_per_gallon_btu, "energy_per_gallon_btu")
        if btu_dec <= 0:
            raise ValueError(f"energy_per_gallon_btu must be strictly positive, got {self.energy_per_gallon_btu}.")
        norm_basis = _converters_energy_normalize_basis(self.basis)
        if not isinstance(self.fuel_type, str) or self.fuel_type.strip().upper() not in {"GGE", "DGE"}:
            raise ValueError(f"fuel_type must be 'GGE' or 'DGE', got {self.fuel_type!r}.")
        object.__setattr__(self, "name", self.name.strip())
        object.__setattr__(self, "energy_per_gallon_btu", float(btu_dec))
        object.__setattr__(self, "basis", norm_basis)
        object.__setattr__(self, "source", self.source.strip())
        object.__setattr__(self, "fuel_type", self.fuel_type.strip().upper())

    @property
    def energy_per_gallon_joule(self) -> float:
        """Exact energy per equivalent gallon in Joules (IT Btu basis)."""
        return float(Decimal(str(self.energy_per_gallon_btu)) * _converters_energy_BTU_IT_JOULES)

    @property
    def energy_per_gallon_mj(self) -> float:
        """Energy per equivalent gallon in megajoules (MJ)."""
        return float((Decimal(str(self.energy_per_gallon_btu)) * _converters_energy_BTU_IT_JOULES) / Decimal("1e6"))

    @property
    def energy_per_gallon_kwh(self) -> float:
        """Energy per equivalent gallon in kilowatt-hours (kWh)."""
        return float((Decimal(str(self.energy_per_gallon_btu)) * _converters_energy_BTU_IT_JOULES) / Decimal("3.6e6"))

    def as_quantity(self) -> pint.Quantity:
        """Return energy density as a Pint Quantity in ``btu_it / gallon``."""
        from core.units.registry import ureg
        return ureg.Quantity(Decimal(str(self.energy_per_gallon_btu)), "btu_it / gallon")





@dataclass(frozen=True, slots=True)
class GallonEquivalent:
    """Specification-tagged gallon-equivalent result (never a physical tank volume)."""

    value: float | Decimal
    fuel_type: GallonEquivalentKind
    spec: str
    basis: HeatingBasis
    energy_per_gallon_btu: float
    source: str

    def __post_init__(self) -> None:
        from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
        val_dec = _converters__interpolation_validate_real_scalar(self.value, "GallonEquivalent.value")
        if val_dec < 0:
            raise ValueError("GallonEquivalent.value must be non-negative.")
        if self.fuel_type not in {"GGE", "DGE"}:
            raise ValueError(f"fuel_type must be 'GGE' or 'DGE', got {self.fuel_type!r}.")
        object.__setattr__(self, "basis", _converters_energy_normalize_basis(self.basis))





@dataclass(frozen=True, slots=True)
class GasEnergyResult:
    """Natural-gas energy result retaining explicit heating-value basis and metering state."""

    energy: pint.Quantity
    basis: HeatingBasis
    reference_conditions: str

    def __post_init__(self) -> None:
        from core.units.registry import ureg
        if not isinstance(self.energy, pint.Quantity):
            raise TypeError("GasEnergyResult.energy must be a pint.Quantity.")
        if self.energy.dimensionality != ureg.joule.dimensionality:
            raise pint.DimensionalityError(
                self.energy.units, ureg.joule, self.energy.dimensionality, ureg.joule.dimensionality
            )
        object.__setattr__(self, "basis", _converters_energy_normalize_basis(self.basis))
        if not isinstance(self.reference_conditions, str) or not self.reference_conditions.strip():
            raise ValueError("reference_conditions must be a non-empty identifier.")

    def __iter__(self) -> Iterator[Any]:
        """Support ``energy_qty, basis = gas_volume_to_energy(...)`` unpacking."""
        yield self.energy
        yield self.basis





def _converters_energy_load_energy_specs(
    path: Path = _converters_energy_DATA_FILE,
) -> tuple[Mapping[str, EnergySpec], Mapping[str, EnergySpec], Mapping[str, Mapping[str, float]]]:
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported energy_specs.json schema_version: {payload.get('schema_version')!r}")

    gge_map: dict[str, EnergySpec] = {}
    for key, item in payload["gge_specs"].items():
        gge_map[key.upper()] = EnergySpec(
            name=item["name"],
            energy_per_gallon_btu=item["energy_per_gallon_btu"],
            basis=item["basis"],
            source=item["source"],
            fuel_type=item.get("fuel_type", "GGE"),
        )

    dge_map: dict[str, EnergySpec] = {}
    for key, item in payload["dge_specs"].items():
        dge_map[key.upper()] = EnergySpec(
            name=item["name"],
            energy_per_gallon_btu=item["energy_per_gallon_btu"],
            basis=item["basis"],
            source=item["source"],
            fuel_type=item.get("fuel_type", "DGE"),
        )

    comp_map: dict[str, Mapping[str, float]] = {
        comp: MappingProxyType({k: float(v) for k, v in vals.items()})
        for comp, vals in payload["gas_components_iso6976_gpa2172"]["components"].items()
    }
    return MappingProxyType(gge_map), MappingProxyType(dge_map), MappingProxyType(comp_map)





GGE_SPECS, DGE_SPECS, ISO6976_COMPONENTS = _converters_energy_load_energy_specs()





def resolve_energy_spec(
    spec: str | EnergySpec,
    *,
    fuel_type: GallonEquivalentKind = "GGE",
) -> EnergySpec:
    """Resolve a named GGE/DGE standard from ``energy_specs.json`` or validate an ``EnergySpec``."""
    from core.units.exceptions import IncompatibleCommodityError
    expected_kind = fuel_type.strip().upper()
    if expected_kind not in {"GGE", "DGE"}:
        raise ValueError(f"fuel_type must be 'GGE' or 'DGE', got {fuel_type!r}.")
    registry_map = GGE_SPECS if expected_kind == "GGE" else DGE_SPECS

    if isinstance(spec, EnergySpec):
        if spec.fuel_type != expected_kind:
            raise IncompatibleCommodityError(
                f"EnergySpec {spec.name!r} has fuel_type={spec.fuel_type!r}, expected {expected_kind!r}."
            )
        return spec
    if not isinstance(spec, str) or not spec.strip():
        raise TypeError("spec must be a named standard string or an EnergySpec instance.")

    key = spec.strip().upper()
    # Allow e.g. "NIST_HB44_DGE" -> "NIST_HB44" lookup on DGE_SPECS
    if key not in registry_map and expected_kind == "DGE" and key.endswith("_DGE"):
        key = key[:-4]
    if key not in registry_map:
        raise ValueError(
            f"Unknown {expected_kind} specification {spec!r}. "
            f"Available {expected_kind} standards: {tuple(registry_map.keys())}."
        )
    return registry_map[key]





def _converters_energy_resolve_basis_multiplier(
    source_basis: HeatingBasis,
    target_basis: HeatingBasis,
    basis_conversion_ratio: float | Decimal | None,
) -> Decimal:
    """Validate and return exact Decimal multiplier E_target / E_source across LHV <-> HHV."""
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if source_basis == target_basis:
        if basis_conversion_ratio is not None:
            ratio = _converters__interpolation_validate_real_scalar(basis_conversion_ratio, "basis_conversion_ratio")
            if ratio != Decimal(1):
                raise ValueError(
                    f"basis_conversion_ratio must be omitted or 1 when source and target share {source_basis} basis."
                )
        return Decimal(1)

    if basis_conversion_ratio is None:
        raise IncompatibleAnalyticalBasisError(
            f"Energy basis mismatch ({source_basis} -> {target_basis}): "
            "refusing to convert without an explicit basis_conversion_ratio."
        )
    ratio = _converters__interpolation_validate_real_scalar(basis_conversion_ratio, "basis_conversion_ratio")
    if ratio <= 0:
        raise ValueError("basis_conversion_ratio must be strictly positive.")
    if source_basis == "HHV" and target_basis == "LHV" and ratio > Decimal(1):
        raise ValueError(
            f"Converting HHV to LHV requires basis_conversion_ratio (LHV/HHV) in (0, 1], got {basis_conversion_ratio}."
        )
    if source_basis == "LHV" and target_basis == "HHV" and ratio < Decimal(1):
        raise ValueError(
            f"Converting LHV to HHV requires basis_conversion_ratio (HHV/LHV) >= 1, got {basis_conversion_ratio}."
        )
    return ratio





def _converters_energy_infer_or_validate_energy_basis(
    energy: pint.Quantity,
    basis: str | None,
) -> HeatingBasis:
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    from core.units.exceptions import MissingAssayError
    unit_name = str(energy.units)
    unit_basis: HeatingBasis | None = None
    if unit_name == "MWh_lhv":
        unit_basis = "LHV"
    elif unit_name == "MWh_hhv":
        unit_basis = "HHV"

    if basis is None:
        if unit_basis is not None:
            return unit_basis
        raise MissingAssayError(
            "Explicit fuel energy basis ('LHV' or 'HHV') is required when converting to GGE/DGE."
        )
    resolved = _converters_energy_normalize_basis(basis)
    if unit_basis is not None and unit_basis != resolved:
        raise IncompatibleAnalyticalBasisError(
            f"Quantity unit {unit_name!r} ({unit_basis}) contradicts basis={basis!r}."
        )
    return resolved





def _converters_energy_energy_to_equivalent_gallon(
    energy: pint.Quantity,
    spec: str | EnergySpec,
    *,
    fuel_type: GallonEquivalentKind,
    basis: str | None = None,
    basis_conversion_ratio: float | Decimal | None = None,
) -> GallonEquivalent:
    from core.units.conversions import _as_decimal
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if not isinstance(energy, pint.Quantity):
        raise TypeError(f"energy must be a pint.Quantity, got {type(energy).__name__}.")
    from core.units.registry import validate_registry_quantity
    validate_registry_quantity(energy)
    from core.units.dimensions import physical_quantity
    physical = physical_quantity(energy, role="thermal")
    if physical.dimensionality != ureg.joule.dimensionality:
        raise pint.DimensionalityError(
            energy.units, ureg.btu_it, energy.dimensionality, ureg.btu_it.dimensionality
        )
    mag_dec = _converters__interpolation_validate_real_scalar(energy.magnitude, "energy")
    if mag_dec < 0:
        raise ValueError("Fuel energy must be non-negative.")

    source_basis = _converters_energy_infer_or_validate_energy_basis(energy, basis)
    resolved_spec = resolve_energy_spec(spec, fuel_type=fuel_type)
    basis_mult = _converters_energy_resolve_basis_multiplier(source_basis, resolved_spec.basis, basis_conversion_ratio)

    energy_btu = _as_decimal(physical.to(ureg.btu_it).magnitude) * basis_mult
    spec_btu = Decimal(str(resolved_spec.energy_per_gallon_btu))
    eq_val = _converters__interpolation_match_scalar_type(energy_btu / spec_btu, energy)
    return GallonEquivalent(
        value=eq_val,
        fuel_type=fuel_type,
        spec=resolved_spec.name,
        basis=resolved_spec.basis,
        energy_per_gallon_btu=resolved_spec.energy_per_gallon_btu,
        source=resolved_spec.source,
    )





def _converters_energy_equivalent_gallon_to_energy(
    equivalent: GallonEquivalent | float | int | Decimal,
    spec: str | EnergySpec | None,
    *,
    fuel_type: GallonEquivalentKind,
    target_unit: str | pint.Unit = "btu_it",
    target_basis: str | None = None,
    basis_conversion_ratio: float | Decimal | None = None,
) -> pint.Quantity:
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    from core.units.exceptions import IncompatibleCommodityError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if isinstance(equivalent, GallonEquivalent):
        if equivalent.fuel_type != fuel_type:
            raise IncompatibleCommodityError(
                f"Expected {fuel_type} GallonEquivalent, got {equivalent.fuel_type!r}."
            )
        if spec is not None:
            checked = resolve_energy_spec(spec, fuel_type=fuel_type)
            if (
                checked.name != equivalent.spec
                or checked.basis != equivalent.basis
                or abs(checked.energy_per_gallon_btu - equivalent.energy_per_gallon_btu) > 1e-9
            ):
                raise ValueError(
                    f"Supplied spec {checked.name!r} conflicts with GallonEquivalent metadata {equivalent.spec!r}."
                )
            resolved_spec = checked
        else:
            resolved_spec = EnergySpec(
                name=equivalent.spec,
                energy_per_gallon_btu=equivalent.energy_per_gallon_btu,
                basis=equivalent.basis,
                source=equivalent.source,
                fuel_type=equivalent.fuel_type,
            )
        raw_val = equivalent.value
    else:
        if spec is None:
            raise ValueError(f"spec is required when converting a numeric {fuel_type} value to energy.")
        resolved_spec = resolve_energy_spec(spec, fuel_type=fuel_type)
        raw_val = equivalent

    val_dec = _converters__interpolation_validate_real_scalar(raw_val, f"{fuel_type} amount")
    if val_dec < 0:
        raise ValueError(f"{fuel_type} amount must be non-negative.")

    dst_unit = parse_unit(target_unit)
    from core.units.dimensions import physical_unit
    if physical_unit(dst_unit, role="thermal").dimensionality != ureg.joule.dimensionality:
        raise pint.DimensionalityError(
            ureg.btu_it, dst_unit, ureg.btu_it.dimensionality, dst_unit.dimensionality
        )

    dst_unit_name = str(dst_unit)
    implied_dst_basis: HeatingBasis | None = (
        "LHV" if dst_unit_name == "MWh_lhv" else "HHV" if dst_unit_name == "MWh_hhv" else None
    )
    if target_basis is None:
        effective_target_basis = implied_dst_basis or resolved_spec.basis
    else:
        effective_target_basis = _converters_energy_normalize_basis(target_basis)
        if implied_dst_basis is not None and implied_dst_basis != effective_target_basis:
            raise IncompatibleAnalyticalBasisError(
                f"target_unit {dst_unit_name!r} ({implied_dst_basis}) contradicts target_basis={target_basis!r}."
            )

    basis_mult = _converters_energy_resolve_basis_multiplier(resolved_spec.basis, effective_target_basis, basis_conversion_ratio)
    total_btu = val_dec * Decimal(str(resolved_spec.energy_per_gallon_btu)) * basis_mult
    from core.units.dimensions import restore
    converted = restore(ureg.Quantity(total_btu, ureg.btu_it), dst_unit, role="thermal")
    out_mag = _converters__interpolation_match_scalar_type(converted.magnitude, raw_val)
    return ureg.Quantity(out_mag, dst_unit)





def energy_to_gge(
    energy: pint.Quantity,
    spec: str | EnergySpec,
    *,
    basis: str | None = None,
    basis_conversion_ratio: float | Decimal | None = None,
) -> GallonEquivalent:
    """Convert fuel energy to Gasoline Gallon Equivalent (GGE) under an explicit spec."""
    return _converters_energy_energy_to_equivalent_gallon(
        energy,
        spec,
        fuel_type="GGE",
        basis=basis,
        basis_conversion_ratio=basis_conversion_ratio,
    )





def gge_to_energy(
    gge: GallonEquivalent | float | int | Decimal,
    spec: str | EnergySpec | None = None,
    *,
    target_unit: str | pint.Unit = "btu_it",
    target_basis: str | None = None,
    basis_conversion_ratio: float | Decimal | None = None,
) -> pint.Quantity:
    """Convert Gasoline Gallon Equivalent (GGE) back to physical energy."""
    return _converters_energy_equivalent_gallon_to_energy(
        gge,
        spec,
        fuel_type="GGE",
        target_unit=target_unit,
        target_basis=target_basis,
        basis_conversion_ratio=basis_conversion_ratio,
    )





def energy_to_dge(
    energy: pint.Quantity,
    spec: str | EnergySpec,
    *,
    basis: str | None = None,
    basis_conversion_ratio: float | Decimal | None = None,
) -> GallonEquivalent:
    """Convert fuel energy to Diesel Gallon Equivalent (DGE) under an explicit spec."""
    return _converters_energy_energy_to_equivalent_gallon(
        energy,
        spec,
        fuel_type="DGE",
        basis=basis,
        basis_conversion_ratio=basis_conversion_ratio,
    )





def dge_to_energy(
    dge: GallonEquivalent | float | int | Decimal,
    spec: str | EnergySpec | None = None,
    *,
    target_unit: str | pint.Unit = "btu_it",
    target_basis: str | None = None,
    basis_conversion_ratio: float | Decimal | None = None,
) -> pint.Quantity:
    """Convert Diesel Gallon Equivalent (DGE) back to physical energy."""
    return _converters_energy_equivalent_gallon_to_energy(
        dge,
        spec,
        fuel_type="DGE",
        target_unit=target_unit,
        target_basis=target_basis,
        basis_conversion_ratio=basis_conversion_ratio,
    )





# ---------------------------------------------------------------------------
# Natural Gas Composition Heating Value (ISO 6976 / GPA 2172) & Volume <-> Energy
# ---------------------------------------------------------------------------


def calculate_gas_calorific_value(
    molar_composition: Mapping[str, float | int | Decimal],
    *,
    standard: Literal["ISO_6976", "GPA_2172"] = "ISO_6976",
    reference_conditions: str | GasReferenceConditions | None = None,
    basis: str = "HHV",
    real_gas: bool = True,
    normalize_composition: bool = False,
) -> pint.Quantity:
    """Compute reference-component sums from the bundled rounded ISO/GPA data.

    This is not a complete standards implementation. Custom real-gas states
    require explicit Z; dry metering and combustion reference states are checked.

    - ``standard="ISO_6976"``: Returns ``MJ / Sm3`` (default ``reference_conditions="iso15"``),
      ``MJ / Nm3`` (``"din0"``), ``MJ / Sm3_20C`` (``"gost20"``), or ``MJ / m**3`` for custom
      dry ``GasReferenceConditions`` (15 °C combustion reference).
    - ``standard="GPA_2172"``: Returns ``btu_it / scf_60F_14696psi_dry`` (default
      ``reference_conditions="atm60f"`` at 60 °F / 14.696 psia), ``btu_it / scf`` (``"ferc"``
      at 14.73 psia), or ``btu_it / scf_tx`` (``"tx"`` at 14.65 psia).
    """
    from core.units.conversions import UNIVERSAL_GAS_CONSTANT_J_MOL_K
    from core.units.conversions import _GAS_REFERENCE_PRESETS
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if not isinstance(molar_composition, Mapping) or not molar_composition:
        raise ValueError("molar_composition must be a non-empty mapping of component -> mole fraction.")
    if not isinstance(real_gas, bool):
        raise TypeError("real_gas must be a boolean.")
    norm_basis = _converters_energy_normalize_basis(basis)
    std_key = standard.strip().upper()
    if std_key not in {"ISO_6976", "GPA_2172"}:
        raise ValueError(f"standard must be 'ISO_6976' or 'GPA_2172', got {standard!r}.")

    cleaned: dict[str, Decimal] = {}
    use_decimal = False
    for raw_comp, raw_frac in molar_composition.items():
        if not isinstance(raw_comp, str) or raw_comp.strip() not in ISO6976_COMPONENTS:
            raise ValueError(
                f"Unsupported gas component {raw_comp!r}. "
                f"Supported components: {tuple(ISO6976_COMPONENTS.keys())}."
            )
        if isinstance(raw_frac, Decimal):
            use_decimal = True
        frac = _converters__interpolation_validate_real_scalar(raw_frac, f"mole fraction for {raw_comp}")
        if frac < 0:
            raise ValueError(f"Mole fraction for {raw_comp!r} cannot be negative.")
        comp_name = raw_comp.strip()
        cleaned[comp_name] = cleaned.get(comp_name, Decimal(0)) + frac

    total_mol = sum(cleaned.values(), Decimal(0))
    if not isinstance(normalize_composition, bool):
        raise TypeError("normalize_composition must be a boolean.")
    if total_mol <= 0 or total_mol != 1 and not normalize_composition:
        raise ValueError(f"Molar composition fractions must sum to one (got {total_mol}); request normalize_composition=True explicitly.")
    norm_comp = {k: v / total_mol if normalize_composition else v for k, v in cleaned.items() if v > 0}
    if norm_comp.get("H2O", 0) > 0:
        raise ValueError("This reference-data function requires dry gas; use gas_properties for explicit wet composition.")

    if std_key == "ISO_6976":
        molar_key = "molar_hhv_kj_mol_15c" if norm_basis == "HHV" else "molar_lhv_kj_mol_15c"
        h_molar_kj = sum(
            x * Decimal(str(ISO6976_COMPONENTS[c][molar_key])) for c, x in norm_comp.items()
        )
        if h_molar_kj <= 0:
            raise ValueError("Gas composition has zero calorific value.")

        sum_sqrt_b = sum(
            x * Decimal(str(ISO6976_COMPONENTS[c]["sqrt_b_15c"])) for c, x in norm_comp.items()
        )
        z_15c = Decimal(1) - (sum_sqrt_b * sum_sqrt_b) if real_gas else Decimal(1)

        if reference_conditions is None:
            ref_key = "iso15"
        elif isinstance(reference_conditions, str):
            ref_key = reference_conditions.strip().lower()
        else:
            ref_key = ""

        if isinstance(reference_conditions, GasReferenceConditions):
            if reference_conditions.water_mole_fraction != Decimal(0):
                raise StandardsApplicabilityError("ISO 6976 composition heating value here requires dry metering conditions (water_mole_fraction=0).")
            if reference_conditions.combustion_reference_temperature is not None and abs(reference_conditions.combustion_reference_temperature.to("K").magnitude - Decimal("288.15")) > Decimal("1e-8"):
                raise StandardsApplicabilityError("Bundled component combustion data are referenced to 15 C; normalize explicitly first.")
            if real_gas and reference_conditions.compressibility_factor is None:
                raise MissingReferenceConditionsError("Custom real-gas metering conditions require explicit compressibility_factor; no 15 C Z extrapolation.")
            if not real_gas and reference_conditions.compressibility_factor not in (None, Decimal(1)):
                raise IncompatibleAnalyticalBasisError("real_gas=False conflicts with the supplied non-ideal metering Z.")
            t_k = reference_conditions.temperature.magnitude
            p_pa = reference_conditions.absolute_pressure.magnitude
            z_state = reference_conditions.compressibility_factor if real_gas else Decimal(1)
            vol_unit_str = "meter ** 3"
        else:
            if ref_key not in {"iso15", "din0", "gost20"}:
                raise ValueError(
                    f"ISO_6976 reference_conditions must be 'iso15', 'din0', 'gost20', or GasReferenceConditions, got {reference_conditions!r}."
                )
            if real_gas and ref_key != "iso15":
                raise StandardsApplicabilityError("Summation-factor data are at 15 C; supply explicit GasReferenceConditions with Z or choose real_gas=False for other metering temperatures.")
            tagged_unit, _, temp_k_flt, press_kpa_flt = _GAS_REFERENCE_PRESETS[ref_key]
            t_k = Decimal(str(temp_k_flt))
            p_pa = Decimal(str(press_kpa_flt)) * Decimal(1000)
            z_state = z_15c
            vol_unit_str = tagged_unit

        r_const = Decimal(str(UNIVERSAL_GAS_CONSTANT_J_MOL_K))
        # V_m (m^3/mol) = Z * R * T / P; MJ/m^3 = (kJ/mol * 1e-3) / V_m
        v_m_m3_per_mol = (z_state * r_const * t_k) / p_pa
        mj_per_vol = (h_molar_kj * Decimal("1e-3")) / v_m_m3_per_mol
        mag = mj_per_vol if use_decimal else float(mj_per_vol)
        return ureg.Quantity(mag, f"megajoule / {vol_unit_str}")

    # GPA_2172
    gpa_key = "gpa_hhv_btu_scf_60f_1atm" if norm_basis == "HHV" else "gpa_lhv_btu_scf_60f_1atm"
    ideal_btu_scf_1atm = sum(
        x * Decimal(str(ISO6976_COMPONENTS[c][gpa_key])) for c, x in norm_comp.items()
    )
    if ideal_btu_scf_1atm <= 0:
        raise ValueError("Gas composition has zero calorific value.")

    sum_sqrt_b_60f = sum(
        x * Decimal(str(ISO6976_COMPONENTS[c]["sqrt_b_60f"])) for c, x in norm_comp.items()
    )
    z_60f = Decimal(1) - (sum_sqrt_b_60f * sum_sqrt_b_60f) if real_gas else Decimal(1)

    ref_key = "atm60f" if reference_conditions is None else (
        reference_conditions.strip().lower() if isinstance(reference_conditions, str) else ""
    )
    if ref_key not in {"atm60f", "ferc", "tx"}:
        raise ValueError(
            f"GPA_2172 reference_conditions must be 'atm60f' (14.696 psia), 'ferc' (14.73 psia), or 'tx' (14.65 psia), got {reference_conditions!r}."
        )
    tagged_unit, _, _, press_kpa_flt = _GAS_REFERENCE_PRESETS[ref_key]
    p_ratio = Decimal(str(press_kpa_flt)) / _converters_energy_P_GPA_14696_KPA
    btu_per_scf = (ideal_btu_scf_1atm * p_ratio) / z_60f
    mag = btu_per_scf if use_decimal else float(btu_per_scf)
    return ureg.Quantity(mag, f"btu_it / {tagged_unit}")





def _converters_energy_resolve_gas_reference_label(
    volume_unit: pint.Unit,
    reference_conditions: str | GasReferenceConditions | None,
) -> str:
    """Verify consistency between a gas volume unit and declared reference_conditions."""
    from core.units.conversions import _GAS_REFERENCE_PRESETS
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.registry import ureg
    vol_dim = volume_unit.dimensionality
    matched_preset: str | None = None
    for key, (tagged, _, temp_k, press_kpa) in _GAS_REFERENCE_PRESETS.items():
        if vol_dim == ureg.Unit(tagged).dimensionality:
            matched_preset = key
            if isinstance(reference_conditions, str) and reference_conditions.strip().lower() != key:
                raise IncompatibleAnalyticalBasisError(
                    f"Tagged gas volume unit {volume_unit} ({key}) contradicts reference_conditions={reference_conditions!r}."
                )
            if isinstance(reference_conditions, GasReferenceConditions):
                if (
                    abs(reference_conditions.temperature.magnitude - Decimal(str(temp_k))) > Decimal("1e-8")
                    or abs(reference_conditions.absolute_pressure.magnitude - Decimal(str(press_kpa)) * 1000) > Decimal("1e-6")
                    or reference_conditions.water_mole_fraction not in (None, Decimal(0))
                ):
                    raise IncompatibleAnalyticalBasisError(
                        f"Tagged gas volume unit {volume_unit} contradicts the supplied GasReferenceConditions."
                    )
            return key

    if vol_dim == (ureg.meter ** 3).dimensionality:
        if reference_conditions is None:
            raise MissingReferenceConditionsError(
                "Geometric gas volumes (e.g. m^3, ft^3) require explicit reference_conditions "
                "(e.g. 'iso15', 'din0', 'ferc', 'atm60f', or GasReferenceConditions)."
            )
        if isinstance(reference_conditions, GasReferenceConditions):
            return reference_conditions.standard_id or (
                f"{reference_conditions.temperature.magnitude}K_{reference_conditions.absolute_pressure.magnitude}Pa"
            )
        if isinstance(reference_conditions, str) and reference_conditions.strip():
            return reference_conditions.strip().lower()
        raise ValueError("reference_conditions must be a non-empty string or GasReferenceConditions.")

    raise pint.DimensionalityError(
        volume_unit,
        ureg.Sm3_15C_101325kPa_dry,
        vol_dim,
        ureg.Sm3_15C_101325kPa_dry.dimensionality,
        extra_msg="; expected a gas reference volume or geometric volume with explicit reference_conditions",
    )


from core.units.validation import authorized_conversion


@authorized_conversion
def gas_volume_to_energy(
    volume: pint.Quantity | float | int | Decimal,
    calorific_value: pint.Quantity | float | int | Decimal,
    *,
    basis: str = "HHV",
    reference_conditions: str | GasReferenceConditions | None = None,
    target_unit: str | pint.Unit | None = None,
) -> GasEnergyResult:
    """Convert natural-gas volume at reference conditions to thermal energy using measured CV.

    Requires strictly positive ``calorific_value``, explicit heating-value ``basis``
    (``"HHV"`` or ``"LHV"``), and consistent metering reference conditions.
    """
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    norm_basis = _converters_energy_normalize_basis(basis)

    if isinstance(volume, pint.Quantity) and isinstance(calorific_value, pint.Quantity):
        from core.units.registry import validate_registry_quantity
        validate_registry_quantity(volume)
        validate_registry_quantity(calorific_value)
        vol_dec = _converters__interpolation_validate_real_scalar(volume.magnitude, "volume")
        cv_dec = _converters__interpolation_validate_real_scalar(calorific_value.magnitude, "calorific_value")
        if vol_dec < 0:
            raise ValueError("Gas volume must be non-negative.")
        if cv_dec <= 0:
            raise ValueError("Calorific value must be strictly positive.")

        ref_label = _converters_energy_resolve_gas_reference_label(volume.units, reference_conditions)
        expected_cv_dim = ureg.joule.dimensionality / volume.dimensionality
        if calorific_value.dimensionality != expected_cv_dim:
            raise pint.DimensionalityError(
                calorific_value.units,
                ureg.megajoule / volume.units,
                calorific_value.dimensionality,
                expected_cv_dim,
                extra_msg="; volume and calorific_value must share the exact same gas reference state",
            )

        product = ureg.Quantity(vol_dec, volume.units) * ureg.Quantity(cv_dec, calorific_value.units)
        if target_unit is not None:
            dst = parse_unit(target_unit)
            if dst.dimensionality != ureg.joule.dimensionality:
                raise pint.DimensionalityError(
                    product.units, dst, product.dimensionality, dst.dimensionality
                )
        else:
            one_cv = (
                ureg.Quantity(Decimal(1), calorific_value.units)
                * ureg.Quantity(Decimal(1), volume.units)
            ).to_reduced_units()
            dst = one_cv.units
            if dst.dimensionality != ureg.joule.dimensionality:
                dst = ureg.megajoule

        converted = product.to(dst)
        out_mag = _converters__interpolation_match_scalar_type(converted.magnitude, volume, calorific_value)
        return GasEnergyResult(
            energy=ureg.Quantity(out_mag, dst),
            basis=norm_basis,
            reference_conditions=ref_label,
        )

    if isinstance(volume, pint.Quantity) != isinstance(calorific_value, pint.Quantity):
        raise TypeError("Pass both volume and calorific_value as pint.Quantity, or both as numeric scalars with reference_conditions.")

    # Scalar path (m^3 and MJ/m^3 at explicitly declared reference_conditions)
    if reference_conditions is None:
        raise MissingReferenceConditionsError(
            "Numeric gas_volume_to_energy requires explicit reference_conditions "
            "(e.g. 'iso15' for 15 °C / 101.325 kPa or 'atm60f' for 60 °F / 14.696 psia)."
        )
    vol_dec = _converters__interpolation_validate_real_scalar(volume, "volume")
    cv_dec = _converters__interpolation_validate_real_scalar(calorific_value, "calorific_value")
    if vol_dec < 0:
        raise ValueError("Gas volume must be non-negative.")
    if cv_dec <= 0:
        raise ValueError("Calorific value must be strictly positive.")
    ref_label = _converters_energy_resolve_gas_reference_label(ureg.meter ** 3, reference_conditions)
    dst = parse_unit(target_unit) if target_unit is not None else ureg.megajoule
    energy_mj = ureg.Quantity(vol_dec * cv_dec, ureg.megajoule).to(dst)
    out_mag = _converters__interpolation_match_scalar_type(energy_mj.magnitude, volume, calorific_value)
    return GasEnergyResult(
        energy=ureg.Quantity(out_mag, dst),
        basis=norm_basis,
        reference_conditions=ref_label,
    )


from core.units.validation import authorized_conversion


@authorized_conversion
def energy_to_gas_volume(
    energy: pint.Quantity,
    calorific_value: pint.Quantity,
    *,
    basis: str = "HHV",
    reference_conditions: str | GasReferenceConditions | None = None,
    target_unit: str | pint.Unit | None = None,
) -> pint.Quantity:
    """Convert thermal energy to natural-gas volume using an explicit measured calorific value."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    _converters_energy_normalize_basis(basis)
    if not isinstance(energy, pint.Quantity) or not isinstance(calorific_value, pint.Quantity):
        raise TypeError("energy and calorific_value must be pint.Quantity instances.")
    from core.units.registry import validate_registry_quantity
    validate_registry_quantity(energy)
    validate_registry_quantity(calorific_value)
    if energy.dimensionality != ureg.joule.dimensionality:
        raise pint.DimensionalityError(
            energy.units, ureg.megajoule, energy.dimensionality, ureg.megajoule.dimensionality
        )
    e_dec = _converters__interpolation_validate_real_scalar(energy.magnitude, "energy")
    cv_dec = _converters__interpolation_validate_real_scalar(calorific_value.magnitude, "calorific_value")
    if e_dec < 0:
        raise ValueError("Energy must be non-negative.")
    if cv_dec <= 0:
        raise ValueError("Calorific value must be strictly positive.")

    vol_unit_from_cv = (ureg.Quantity(Decimal(1), energy.units) / ureg.Quantity(Decimal(1), calorific_value.units)).to_reduced_units().units
    _converters_energy_resolve_gas_reference_label(vol_unit_from_cv, reference_conditions)

    dst = parse_unit(target_unit) if target_unit is not None else vol_unit_from_cv
    if dst.dimensionality != vol_unit_from_cv.dimensionality:
        raise pint.DimensionalityError(
            vol_unit_from_cv, dst, vol_unit_from_cv.dimensionality, dst.dimensionality
        )
    raw_vol = (ureg.Quantity(e_dec, energy.units) / ureg.Quantity(cv_dec, calorific_value.units)).to(dst)
    out_mag = _converters__interpolation_match_scalar_type(raw_vol.magnitude, energy, calorific_value)
    return ureg.Quantity(out_mag, dst)


__all__ = ('ANALYTICAL_BASES', 'CalorificMethod', 'DGE_SPECS', 'EmpiricalHeatCorrelation', 'EnergySpec', 'FuelCondition', 'GGE_SPECS', 'GallonEquivalent', 'GasEnergyResult', 'GasMoistureState', 'HeatContentAssay', 'ISO6976Table', 'ISO6976_COMPONENTS', 'LaboratoryCalorificResult', 'ProximateAnalysis', 'ThermodynamicBasis', 'UltimateAnalysis', 'as_received_to_dry', 'calculate_gas_calorific_value', 'convert_combustion_heat', 'convert_heat_content', 'create_heat_content_context', 'dge_to_energy', 'dry_to_as_received', 'efficiency_from_heat_rate', 'energy_to_dge', 'energy_to_gas_volume', 'energy_to_gge', 'gas_volume_to_energy', 'gcv_to_ncv', 'gge_to_energy', 'heat_content_from_assay', 'heat_rate_from_efficiency', 'ideal_gas_combustion_correction', 'ncv_to_gcv', 'output_to_thermal', 'resolve_energy_spec', 'thermal_to_output', 'water_latent_heat')
_PUBLIC_ALIASES = {'IncompatibleAnalyticalBasisError': ('exceptions', 'IncompatibleAnalyticalBasisError'), 'authorized_conversion': ('validation', 'authorized_conversion'), 'pad_ieee754_decimal': ('precision', 'pad_ieee754_decimal'), '_scale_factor': ('conversions', '_scale_factor'), 'MissingReferenceConditionsError': ('exceptions', 'MissingReferenceConditionsError'), 'ureg': ('registry', 'ureg'), 'validate_registry_quantity': ('registry', 'validate_registry_quantity'), 'parse_unit': ('registry', 'parse_unit'), 'MissingAssayError': ('exceptions', 'MissingAssayError'), 'dry_matter_fraction': ('analytical_basis', 'dry_matter_fraction'), '_QuantitySnapshot': ('validation', '_safety_QuantitySnapshot'), '_energy_restore': ('dimensions', 'restore'), 'StandardsApplicabilityError': ('exceptions', 'StandardsApplicabilityError'), '_as_decimal': ('conversions', '_as_decimal'), '_energy_unit': ('dimensions', 'physical_unit'), 'InvalidPhysicalStateError': ('exceptions', 'InvalidPhysicalStateError'), '_energy_physical': ('dimensions', 'physical_quantity')}

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
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.registry import ureg
