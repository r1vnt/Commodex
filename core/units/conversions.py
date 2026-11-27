"""Commodex conversions implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from pathlib import Path
from types import MappingProxyType
from typing import Final
from typing import TYPE_CHECKING, Any, Final
import inspect
import json
import math
import pint
import warnings



# ---------------------------------------------------------------------------
# 1. Physical, Thermodynamic & Custody-Transfer Reference Constants
# ---------------------------------------------------------------------------

# --- Petroleum & Water Reference Densities (API MPMS Ch. 11.1 / ASTM D1250) ---
# Standard reference density of air-free water at 60°F (15.5556°C) in kg/m³ (US/API basis)
WATER_DENSITY_60F_KG_M3: Final[float] = 999.016



# Standard reference density of air-free water at 15°C (59°F) in kg/m³ (ISO/Metric basis)
WATER_DENSITY_15C_KG_M3: Final[float] = 999.102



# Standard reference density of air-free water at 60°F in lb/gal (US liquid)
WATER_DENSITY_60F_LB_GAL: Final[float] = 8.337183637



# Exact volume of one 42-US-gallon petroleum barrel in cubic metres
BARREL_TO_M3: Final[float] = 0.158987294928



# Mass of 1 bbl of water at 60°F in metric tonnes (~0.158830764 mt)
MT_PER_BBL_WATER: Final[float] = (BARREL_TO_M3 * WATER_DENSITY_60F_KG_M3) / 1000.0



# Exact cubic feet per cubic metre (1 / 0.3048³)
FT3_PER_M3: Final[float] = 35.31466672148859



# Exact International Table MJ per Btu_it (1055.05585262 J / 1e6)
MJ_PER_BTU_IT: Final[float] = 0.00105505585262




# Legacy temperature-only empirical constants; not a complete API MPMS 11.1 implementation.
# For bounded/product-specific CTL/CPL/CTPL use converters.petroleum methods.
API_VCF_CRUDE_K0: Final[float] = 341.0957



API_VCF_CRUDE_K1: Final[float] = 0.0



API_VCF_CRUDE_K2: Final[float] = 0.0



API_VCF_GASOLINE_K0: Final[float] = 192.4571



API_VCF_GASOLINE_K1: Final[float] = 0.2438



API_VCF_DISTILLATE_K0: Final[float] = 103.8720



API_VCF_DISTILLATE_K1: Final[float] = 0.2701




# --- Fundamental Gas & Ideal Gas Law Constants (CODATA 2019 / AGA Report No. 8) ---
# Exact CODATA molar gas constant R in J / (mol * K)
UNIVERSAL_GAS_CONSTANT_J_MOL_K: Final[float] = 8.314462618



# Ideal gas standard molar volume at 0°C (273.15 K) and 1 atm (101.325 kPa) in m³/kmol
MOLAR_VOLUME_0C_1ATM_M3_KMOL: Final[float] = 22.413969545014137



# Ideal gas standard molar volume at 15°C (288.15 K) and 1 atm (101.325 kPa) in m³/kmol
MOLAR_VOLUME_15C_1ATM_M3_KMOL: Final[float] = 23.644800214596617



# Standard molar mass of dry air in kg/kmol (AGA Report No. 8 / ISO 6976)
DRY_AIR_MOLAR_MASS_KG_KMOL: Final[float] = 28.96546



# Standard density of dry air at 0°C, 101.325 kPa in kg/m³
DRY_AIR_DENSITY_0C_KG_M3: Final[float] = 1.292923



# Standard density of dry air at 15°C, 101.325 kPa in kg/m³
DRY_AIR_DENSITY_15C_KG_M3: Final[float] = 1.225525




# --- Cryogenic LNG & Clean Energy Carrier Benchmarks (GIIGNL / ISO 6976) ---
# Typical liquefied natural gas density at atmospheric boiling point (-162°C) in kg/m³
LNG_LIQUID_DENSITY_KG_M3: Final[float] = 450.0



# Standard GIIGNL nominal LNG volume-to-gas expansion factor (Sm³ gas per m³ liquid LNG)
LNG_GAS_EXPANSION_SM3_PER_M3: Final[float] = 600.0



# Pure Hydrogen (H2) gravimetric calorific values (ISO 6976 / NIST) in MJ/kg
H2_LHV_MJ_PER_KG: Final[float] = 119.96



H2_HHV_MJ_PER_KG: Final[float] = 141.88



# Pure Hydrogen standard volumetric density at 0°C, 101.325 kPa in kg/Nm³
H2_DENSITY_0C_KG_NM3: Final[float] = 0.08988



# Anhydrous Ammonia (NH3) gravimetric calorific values in MJ/kg
NH3_LHV_MJ_PER_KG: Final[float] = 18.60



NH3_HHV_MJ_PER_KG: Final[float] = 22.50




# --- Atmospheric & Barometric Datum Constants (ISO 2533 Standard Atmosphere) ---
# Standard sea-level atmospheric pressure in Pascals (101.325 kPa = 14.6959488 psia = 1.01325 bara)
STANDARD_ATMOSPHERE_PA: Final[float] = 101325.0



# Standard sea-level temperature in Kelvin (15°C)
STANDARD_SEA_LEVEL_TEMP_K: Final[float] = 288.15



# Standard tropospheric temperature lapse rate in K/m
BAROMETRIC_LAPSE_RATE_K_M: Final[float] = 0.0065



# Standard gravitational acceleration in m/s²
STANDARD_GRAVITY_M_S2: Final[float] = 9.80665





# ---------------------------------------------------------------------------
# Refutas Viscosity Blending Index (Standard Refinery LP Approximation)
# ---------------------------------------------------------------------------

def _viscosity_blending_number(viscosity_cst: float) -> float:
    """Computes Refutas Viscosity Blending Number (VBN) from kinematic viscosity (cSt)."""
    if not math.isfinite(viscosity_cst) or viscosity_cst <= 0.2:
        raise ValueError(f"Refutas kinematic viscosity must exceed 0.2 cSt, got {viscosity_cst}")
    return 14.534 * math.log(math.log(viscosity_cst + 0.8)) + 10.975





def _viscosity_from_vbn(vbn: float) -> float:
    """Inverts Refutas VBN back to kinematic viscosity (cSt)."""
    return math.exp(math.exp((vbn - 10.975) / 14.534)) - 0.8





# ---------------------------------------------------------------------------
# Thermodynamic Gas Reference Constants (Ideal Gas Law: PV / T = const)
# ---------------------------------------------------------------------------

_T_60F_K: Final[float] = 288.70555555555555



_T_20C_K: Final[float] = 293.15  # GOST 2939-63 Eurasian/CIS reference temperature (20°C)



_T_15C_K: Final[float] = 288.15



_T_0C_K: Final[float] = 273.15




_P_14_73_KPA: Final[float] = 101.559775234



_P_14_65_KPA: Final[float] = 101.008194216



_P_1ATM_KPA: Final[float] = 101.325




# Ideal-gas molar density ratios relative to FERC 60°F / 14.73 psia (assumes Z_src / Z_dst = 1.0)
_DENSITY_RATIO_FERC: Final[float] = 1.0



_DENSITY_RATIO_TX: Final[float] = (_P_14_65_KPA / _T_60F_K) / (_P_14_73_KPA / _T_60F_K)



_DENSITY_RATIO_ISO15: Final[float] = (_P_1ATM_KPA / _T_15C_K) / (_P_14_73_KPA / _T_60F_K)



_DENSITY_RATIO_DIN0: Final[float] = (_P_1ATM_KPA / _T_0C_K) / (_P_14_73_KPA / _T_60F_K)



_DENSITY_RATIO_GOST20: Final[float] = (_P_1ATM_KPA / _T_20C_K) / (_P_14_73_KPA / _T_60F_K)




_GRONINGEN_MJ_PER_M3: Final[float] = 35.17





def _ideal_gas_density_per_scf(molar_mass_g_per_mol: float) -> float:
    """Mass per FERC scf for the compatibility gas model's ideal reference states."""
    return (_P_14_73_KPA * 1000 * (molar_mass_g_per_mol / 1000)
            / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * _T_60F_K) * (0.3048 ** 3))





def _scale_factor(magnitude: Any, factor: float | Decimal | int) -> Decimal | float:
    """Returns a Decimal or float multiplier matching the incoming Quantity magnitude."""
    if isinstance(magnitude, Decimal):
        return Decimal(str(factor))
    return float(factor)





_DAK_A: Final[tuple[float, ...]] = (
    0.3265,    # A1
    -1.0700,   # A2
    -0.5339,   # A3
    0.01569,   # A4
    -0.05165,  # A5
    0.5475,    # A6
    -0.7361,   # A7
    0.1844,    # A8
    0.1056,    # A9
    0.6134,    # A10
    0.7210,    # A11
)




# Canonical 21-component AGA-8 Detail / GERG-2008 molar ordering
_AGA8_COMPONENTS: Final[tuple[str, ...]] = (
    "CH4", "N2", "CO2", "C2H6", "C3H8", "iC4H10", "nC4H10",
    "iC5H12", "nC5H12", "C6H14", "C7H16", "C8H18", "C9H20", "C10H22",
    "H2", "O2", "CO", "H2O", "H2S", "He", "Ar",
)



_DAK_UNSUPPORTED_NON_HYDROCARBONS: Final[frozenset[str]] = frozenset({"H2", "O2", "CO", "H2O", "He", "Ar"})





def _validated_composition(composition, *, normalize=False):
    from core.units.thermodynamics import SolverDomainError
    if not composition:
        raise SolverDomainError("A non-empty molar composition is required.")
    unknown = set(composition) - set(_AGA8_COMPONENTS)
    if unknown:
        raise SolverDomainError(f"Unsupported gas components: {sorted(unknown)}")
    if not isinstance(normalize, bool) or any(isinstance(value, bool) for value in composition.values()):
        raise TypeError("Composition fractions and normalization policy must have their declared numeric/boolean types.")
    values = {name: float(value) for name, value in composition.items()}
    if any(not math.isfinite(value) or value < 0 for value in values.values()):
        raise SolverDomainError("Composition fractions must be finite and non-negative.")
    total = sum(values.values())
    if total <= 0:
        raise SolverDomainError("Composition must have positive total moles.")
    if not normalize and not math.isclose(total, 1., rel_tol=0., abs_tol=1e-12):
        raise SolverDomainError("Composition fractions must sum to one; request normalization explicitly.")
    return {name: value / total if normalize else value for name, value in values.items() if value > 0}





def _checked_backend_z(value: Any, backend_label: str) -> float:
    from core.units.thermodynamics import SolverBackendError
    try:
        z_val = float(value)
    except (TypeError, ValueError) as error:
        raise SolverBackendError(f"{backend_label} returned a non-numeric Z-factor: {value!r}") from error
    if not math.isfinite(z_val) or z_val <= 0.0:
        raise SolverBackendError(f"{backend_label} returned an invalid Z-factor: {z_val!r}")
    return z_val





def _detect_package_version(module_obj: Any, dist_name: str) -> str | None:
    ver = getattr(module_obj, "__version__", None)
    if isinstance(ver, str) and ver.strip():
        return ver.strip()
    import sys as _sys
    mod_file = getattr(module_obj, "__file__", None)
    if mod_file is not None and dist_name in _sys.modules:
        try:
            from importlib.metadata import version as _pkg_version
            resolved = _pkg_version(dist_name)
            if resolved and resolved.strip():
                return resolved.strip()
        except Exception:
            pass
    return None


from core.units.standards import _reference_reference_document
_GAS_REFERENCE_PRESETS = {
    record["id"]: (record["tagged_unit"], record["geometric_unit"], Decimal(record["temperature"]), Decimal(record["pressure"]))
    for record in _reference_reference_document("units", "gas_standards")["records"]
}





def _absolute_reference_temperature(value: pint.Quantity, label: str) -> pint.Quantity:
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"{label} must be a temperature quantity.")
    validate_registry_quantity(value)
    if any(str(name).startswith("delta_") for name in value.units._units):
        raise ValueError(f"{label} must be an absolute temperature, not an interval.")
    converted = ureg.Quantity(_as_decimal(value.magnitude), ureg.Unit(value.units._units)).to("kelvin")
    if not converted.magnitude.is_finite() or converted.magnitude <= 0:
        raise ValueError(f"{label} must be finite and greater than zero kelvin.")
    return converted





def _reference_geometric_unit(unit: pint.Unit, state: GasReferenceConditions) -> tuple[pint.Unit, pint.Unit]:
    from core.units.registry import ureg
    if unit.dimensionality == (ureg.meter**3).dimensionality:
        return unit, unit
    for key, (tagged, geometric, temperature, pressure) in _GAS_REFERENCE_PRESETS.items():
        tag = ureg.Unit(tagged)
        if unit.dimensionality != tag.dimensionality:
            continue
        if (abs(state.temperature.magnitude - Decimal(str(temperature))) > Decimal("1e-8")
                or abs(state.absolute_pressure.magnitude - Decimal(str(pressure)) * 1000) > Decimal("1e-6")
                or state.water_mole_fraction != 0):
            raise ValueError("Tagged dry gas-volume unit contradicts the supplied reference conditions.")
        return tag, ureg.Unit(geometric)
    raise ValueError("Expected geometric gas volume or a matching named reference volume; calorific-equivalent volumes need an energy model.")





def convert_gas_reference_volume(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    source: GasReferenceConditions,
    target: GasReferenceConditions,
    conserve: str = "dry_gas",
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Convert single-phase gas volumes between arbitrary contractual states.

    dry_gas conserves non-water moles and multiplies P1/P2*T2/T1*Z2/Z1 by
    (1-y1)/(1-y2). It requires both water fractions (zero explicitly means dry).
    total_gas uses the fixed-total-mole equation and never performs drying;
    conflicting known water fractions are rejected. It may use unknown fractions
    only when the caller explicitly chooses this mode.

    Combustion temperature is retained as metadata and never enters this metering
    equation. No phase-equilibrium, condensation, composition, calorific-value or
    EOS calculation is implied. Geometric results must travel with their state.
    """
    from core.units.precision import pad_ieee754_decimal
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError("quantity must be a Pint volume quantity.")
    validate_registry_quantity(quantity)
    if not isinstance(source, GasReferenceConditions) or not isinstance(target, GasReferenceConditions):
        raise TypeError("source and target must be GasReferenceConditions.")
    # Revalidate copies: Pint Quantity members themselves are mutable, even in
    # a frozen dataclass, and input quantities must never be modified in place.
    from dataclasses import replace
    source, target = replace(source), replace(target)
    amount = _as_decimal(quantity.magnitude)
    if not amount.is_finite() or amount < 0:
        raise ValueError("Gas volume must be finite and non-negative.")
    if source.compressibility_factor is None or target.compressibility_factor is None:
        raise ValueError("Supply compressibility_factor for both states; use 1 explicitly for an ideal-gas assumption.")
    if conserve not in {"dry_gas", "total_gas"}:
        raise ValueError("conserve must be 'dry_gas' or 'total_gas'.")
    source_water, target_water = source.water_mole_fraction, target.water_mole_fraction
    water_factor = Decimal(1)
    if conserve == "dry_gas":
        if source_water is None or target_water is None:
            raise ValueError("Dry-gas conservation requires water_mole_fraction for both states.")
        water_factor = (1 - source_water) / (1 - target_water)
    elif source_water is not None and target_water is not None and source_water != target_water:
        raise ValueError("Fixed total moles cannot also represent a wet/dry composition change.")
    destination = parse_unit(target_unit)
    source_base, source_geometry = _reference_geometric_unit(quantity.units, source)
    target_base, target_geometry = _reference_geometric_unit(destination, target)
    geometrical = ureg.Quantity(quantity.to(source_base).magnitude, source_geometry).to(target_geometry)
    factor = (source.absolute_pressure.magnitude / target.absolute_pressure.magnitude
              * target.temperature.magnitude / source.temperature.magnitude
              * target.compressibility_factor / source.compressibility_factor * water_factor)
    result = ureg.Quantity(geometrical.magnitude * _scale_factor(geometrical.magnitude, factor), target_base).to(destination)
    if not _as_decimal(result.magnitude).is_finite():
        raise ValueError("Converted gas volume exceeds the numeric representation.")
    if isinstance(quantity.magnitude, Decimal):
        converted_qty = ureg.Quantity(pad_ieee754_decimal(_as_decimal(result.magnitude)), result.units)
    else:
        converted_qty = ureg.Quantity(float(result.magnitude), result.units)

    if not return_metadata:
        return converted_qty

    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    ref = ReferenceConditions(
        temperature=target.temperature,
        pressure=target.absolute_pressure,
        compressibility_factor=target.compressibility_factor,
        moisture_fraction=target.water_mole_fraction,
        analytical_basis=conserve,
        classification=ConversionClassification.STATE_AND_MODEL_DEPENDENT.value,
        parameters={
            "source_temperature": source.temperature,
            "source_pressure": source.absolute_pressure,
            "source_compressibility_factor": source.compressibility_factor,
            "source_water_mole_fraction": source.water_mole_fraction,
            "target_temperature": target.temperature,
            "target_pressure": target.absolute_pressure,
            "target_compressibility_factor": target.compressibility_factor,
            "target_water_mole_fraction": target.water_mole_fraction,
            "conservation_mode": conserve,
        },
    )
    base_params = {
        "quantity": _as_decimal(quantity.magnitude),
        "source_pressure": _as_decimal(source.absolute_pressure.magnitude),
        "source_temperature": _as_decimal(source.temperature.magnitude),
        "target_pressure": _as_decimal(target.absolute_pressure.magnitude),
        "target_temperature": _as_decimal(target.temperature.magnitude),
        "z_source": _as_decimal(source.compressibility_factor),
        "z_target": _as_decimal(target.compressibility_factor),
    }

    def _eval_gas_ref(p: Mapping[str, Decimal]) -> Decimal:
        from core.units.registry import ureg
        src_p = replace(
            source,
            temperature=ureg.Quantity(p["source_temperature"], source.temperature.units),
            absolute_pressure=ureg.Quantity(p["source_pressure"], source.absolute_pressure.units),
            compressibility_factor=p["z_source"],
        )
        dst_p = replace(
            target,
            temperature=ureg.Quantity(p["target_temperature"], target.temperature.units),
            absolute_pressure=ureg.Quantity(p["target_pressure"], target.absolute_pressure.units),
            compressibility_factor=p["z_target"],
        )
        q_in = ureg.Quantity(p["quantity"], quantity.units)
        out = convert_gas_reference_volume(q_in, destination, source=src_p, target=dst_p, conserve=conserve)
        return _as_decimal(out.magnitude)

    return build_conversion_result(
        converted_qty,
        quantity,
        commodity_id="natural_gas",
        conversion_path=(str(quantity.units), f"gas_state:{conserve}", str(converted_qty.units)),
        reference_conditions=ref,
        method="convert_gas_reference_volume",
        standard="ISO 13443 / AGA Report No. 8",
        estimated=False,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=_eval_gas_ref,
        base_parameters=base_params,
    )





def convert_actual_to_standard_gas_volume(
    actual_volume: pint.Quantity,
    actual_pressure: pint.Quantity,
    actual_temperature: pint.Quantity,
    *,
    target_standard: str = "ferc",
    relative_density: float = 0.60,
    method: str = "auto",
    molar_composition: dict[str, float] | None = None,
    z_actual: float | None = None,
    z_standard: float | None = None,
    co2_mol_frac: float | None = None,
    h2s_mol_frac: float | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts actual flowing or storage gas volume to standard reference volume.

    Applies the real-gas law:
        V_std = V_act * (P_act / P_std) * (T_std / T_act) * (Z_std / Z_act)

    This compatibility API assumes dry gas at both endpoints. Use
    convert_gas_reference_volume for arbitrary/wet contractual states.

    Supported `target_standard` baselines:
        - "ferc"   : 60°F, 14.73 psia (NYMEX / US interstate pipelines)
        - "tx"     : 60°F, 14.65 psia (Texas RRC)
        - "iso15"  : 15°C, 101.325 kPa (ISO 13443 / UK NBP / LNG)
        - "din0"   :  0°C, 101.325 kPa (DIN 1343 / Continental Europe)
        - "gost20" : 20°C, 101.325 kPa (GOST 2939-63 / Eurasian pipelines)
    """
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.registry import ureg
    from core.units.thermodynamics import calculate_z_factor
    if strict and (z_actual is None or z_standard is None) and method.strip().lower() == "auto":
        raise MissingReferenceConditionsError(
            "strict=True requires explicit z_actual/z_standard or an explicit EOS method instead of method='auto'."
        )
    std_key = target_standard.strip().lower()
    gost_unit = ureg.m3_20C_101325kPa_dry
    std_configs = {
        key: (ureg.Unit(tagged), ureg.Unit(geometric), pressure, temperature)
        for key, (tagged, geometric, temperature, pressure) in _GAS_REFERENCE_PRESETS.items()
    }
    if std_key not in std_configs:
        raise ValueError(f"Unknown target_standard {target_standard!r}. Choose from {tuple(std_configs)}.")

    dst_unit, geom_unit, p_std_kpa, t_std_k = std_configs[std_key]

    actual_pressure_pa = convert_pressure(actual_pressure, ureg.pascal, strict=strict)
    p_act_kpa = float(actual_pressure_pa.to(ureg.kilopascal).magnitude)
    t_act_k = float(actual_temperature.to(ureg.kelvin).magnitude)
    if not math.isfinite(p_act_kpa) or p_act_kpa <= 0 or not math.isfinite(t_act_k) or t_act_k <= 0:
        raise ValueError("Actual gas pressure and absolute temperature must be finite and strictly positive.")

    z_act = z_actual if z_actual is not None else calculate_z_factor(
        actual_pressure_pa,
        actual_temperature,
        relative_density,
        method=method,
        molar_composition=molar_composition,
        co2_mol_frac=co2_mol_frac,
        h2s_mol_frac=h2s_mol_frac,
    )
    z_std = z_standard if z_standard is not None else calculate_z_factor(
        ureg.Quantity(p_std_kpa, ureg.kilopascal),
        ureg.Quantity(t_std_k, ureg.kelvin),
        relative_density,
        method=method,
        molar_composition=molar_composition,
        co2_mol_frac=co2_mol_frac,
        h2s_mol_frac=h2s_mol_frac,
    )

    if not all(math.isfinite(value) and value > 0 for value in (z_act, z_std)):
        raise ValueError("Gas compressibility factors must be finite and strictly positive.")
    converted_qty = convert_gas_reference_volume(
        actual_volume, dst_unit,
        source=GasReferenceConditions(actual_temperature, actual_pressure_pa, _as_decimal(z_act), Decimal(0)),
        target=GasReferenceConditions.from_standard(std_key, compressibility_factor=_as_decimal(z_std)),
    )
    if not return_metadata:
        return converted_qty

    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    tgt_ref = GasReferenceConditions.from_standard(std_key, compressibility_factor=_as_decimal(z_std))
    ref = ReferenceConditions(
        temperature=tgt_ref.temperature,
        pressure=tgt_ref.absolute_pressure,
        compressibility_factor=tgt_ref.compressibility_factor,
        analytical_basis="dry_gas",
        classification=ConversionClassification.STATE_AND_MODEL_DEPENDENT.value,
        parameters={
            "actual_temperature": actual_temperature,
            "actual_pressure": actual_pressure_pa,
            "z_actual": z_act,
            "z_standard": z_std,
            "eos_method": method,
            "target_standard": std_key,
            "relative_density": relative_density,
        },
    )
    base_params = {
        "quantity": _as_decimal(actual_volume.magnitude),
        "actual_pressure": _as_decimal(actual_pressure_pa.magnitude),
        "actual_temperature": _as_decimal(actual_temperature.to(ureg.kelvin).magnitude),
        "z_actual": _as_decimal(z_act),
        "z_standard": _as_decimal(z_std),
    }

    def _eval_act_std(p: Mapping[str, Decimal]) -> Decimal:
        from core.units.registry import ureg
        use_z_act = float(p["z_actual"]) if (z_actual is not None or "z_actual" in (uncertainties if isinstance(uncertainties, Mapping) else {})) else None
        use_z_std = float(p["z_standard"]) if (z_standard is not None or "z_standard" in (uncertainties if isinstance(uncertainties, Mapping) else {})) else None
        out = convert_actual_to_standard_gas_volume(
            ureg.Quantity(p["quantity"], actual_volume.units),
            ureg.Quantity(p["actual_pressure"], ureg.pascal),
            ureg.Quantity(p["actual_temperature"], ureg.kelvin),
            target_standard=target_standard,
            relative_density=relative_density,
            method=method,
            molar_composition=molar_composition,
            z_actual=use_z_act,
            z_standard=use_z_std,
            co2_mol_frac=co2_mol_frac,
            h2s_mol_frac=h2s_mol_frac,
        )
        return _as_decimal(out.magnitude)

    return build_conversion_result(
        converted_qty,
        actual_volume,
        commodity_id="natural_gas",
        conversion_path=(str(actual_volume.units), f"real_gas:{method}", str(converted_qty.units)),
        reference_conditions=ref,
        method=f"convert_actual_to_standard_gas_volume:{method}",
        standard="AGA Report No. 8 / ISO 13443",
        standard_version=std_key,
        estimated=(z_actual is None or z_standard is None),
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=_eval_act_std,
        base_parameters=base_params,
    )





def convert_standard_to_actual_gas_volume(
    standard_volume: pint.Quantity,
    actual_pressure: pint.Quantity,
    actual_temperature: pint.Quantity,
    *,
    source_standard: str = "ferc",
    target_unit: str | pint.Unit | None = None,
    relative_density: float = 0.60,
    method: str = "auto",
    molar_composition: dict[str, float] | None = None,
    z_actual: float | None = None,
    z_standard: float | None = None,
    co2_mol_frac: float | None = None,
    h2s_mol_frac: float | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts standard reference gas volume to actual geometric volume at (P_act, T_act).

    Applies the inverse real-gas law:
        V_act = V_std * (P_std / P_act) * (T_act / T_std) * (Z_act / Z_std)

    This compatibility API assumes dry gas at both endpoints. Use
    convert_gas_reference_volume for arbitrary/wet contractual states.

    Supported `source_standard` baselines:
        - "ferc"   : 60°F, 14.73 psia (NYMEX / US interstate pipelines)
        - "tx"     : 60°F, 14.65 psia (Texas RRC)
        - "iso15"  : 15°C, 101.325 kPa (ISO 13443 / UK NBP / LNG)
        - "din0"   :  0°C, 101.325 kPa (DIN 1343 / Continental Europe)
        - "gost20" : 20°C, 101.325 kPa (GOST 2939-63 / Eurasian pipelines)
    """
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.thermodynamics import calculate_z_factor
    if strict and (z_actual is None or z_standard is None) and method.strip().lower() == "auto":
        raise MissingReferenceConditionsError(
            "strict=True requires explicit z_actual/z_standard or an explicit EOS method instead of method='auto'."
        )
    std_key = source_standard.strip().lower()
    gost_unit = ureg.m3_20C_101325kPa_dry
    std_configs = {
        key: (ureg.Unit(tagged), ureg.Unit(geometric), pressure, temperature)
        for key, (tagged, geometric, temperature, pressure) in _GAS_REFERENCE_PRESETS.items()
    }
    if std_key not in std_configs:
        raise ValueError(f"Unknown source_standard {source_standard!r}. Choose from {tuple(std_configs)}.")

    src_std_unit, geom_unit, p_std_kpa, t_std_k = std_configs[std_key]

    actual_pressure_pa = convert_pressure(actual_pressure, ureg.pascal, strict=strict)
    p_act_kpa = float(actual_pressure_pa.to(ureg.kilopascal).magnitude)
    t_act_k = float(actual_temperature.to(ureg.kelvin).magnitude)
    if not math.isfinite(p_act_kpa) or p_act_kpa <= 0 or not math.isfinite(t_act_k) or t_act_k <= 0:
        raise ValueError("Actual gas pressure and absolute temperature must be finite and strictly positive.")
    if p_act_kpa <= 0:
        raise ValueError(f"actual_pressure must be strictly positive, got {p_act_kpa} kPa.")

    z_act = z_actual if z_actual is not None else calculate_z_factor(
        actual_pressure_pa,
        actual_temperature,
        relative_density,
        method=method,
        molar_composition=molar_composition,
        co2_mol_frac=co2_mol_frac,
        h2s_mol_frac=h2s_mol_frac,
    )
    z_std = z_standard if z_standard is not None else calculate_z_factor(
        ureg.Quantity(p_std_kpa, ureg.kilopascal),
        ureg.Quantity(t_std_k, ureg.kelvin),
        relative_density,
        method=method,
        molar_composition=molar_composition,
        co2_mol_frac=co2_mol_frac,
        h2s_mol_frac=h2s_mol_frac,
    )

    if not all(math.isfinite(value) and value > 0 for value in (z_act, z_std)):
        raise ValueError("Gas compressibility factors must be finite and strictly positive.")
    converted_qty = convert_gas_reference_volume(
        standard_volume, parse_unit(target_unit) if target_unit is not None else geom_unit,
        source=GasReferenceConditions.from_standard(std_key, compressibility_factor=_as_decimal(z_std)),
        target=GasReferenceConditions(actual_temperature, actual_pressure_pa, _as_decimal(z_act), Decimal(0)),
    )
    if not return_metadata:
        return converted_qty

    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    ref = ReferenceConditions(
        temperature=actual_temperature,
        pressure=actual_pressure_pa,
        compressibility_factor=_as_decimal(z_act),
        analytical_basis="dry_gas",
        classification=ConversionClassification.STATE_AND_MODEL_DEPENDENT.value,
        parameters={
            "source_standard": std_key,
            "z_actual": z_act,
            "z_standard": z_std,
            "eos_method": method,
            "relative_density": relative_density,
        },
    )
    return build_conversion_result(
        converted_qty,
        standard_volume,
        commodity_id="natural_gas",
        conversion_path=(str(standard_volume.units), f"real_gas:{method}", str(converted_qty.units)),
        reference_conditions=ref,
        method=f"convert_standard_to_actual_gas_volume:{method}",
        standard="AGA Report No. 8 / ISO 13443",
        standard_version=std_key,
        estimated=(z_actual is None or z_standard is None),
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
    )





# ---------------------------------------------------------------------------
# Thermal & Metallurgical Coal Calorific Specification
# ---------------------------------------------------------------------------

# ISO 1928 / ASTM D5865 Gross-to-Net latent heat constants (kcal/kg per 1 wt%)
_ISO1928_H_KCAL_PER_PCT: Final[float] = 50.68



_ISO1928_M_KCAL_PER_PCT: Final[float] = 5.835



_KCAL_KG_TO_GJ_MT: Final[float] = 0.0041868      # 1 kcal_IT/kg = 4.1868 kJ/kg = 0.0041868 GJ/mt



_KCAL_KG_TO_BTU_LB: Final[float] = 1.8           # Exact IT definition





# ---------------------------------------------------------------------------
# Clean Fuel Carriers (Hydrogen, Ammonia & Methanol LHV / HHV Specification)
# ---------------------------------------------------------------------------

_LH2_LIQUID_DENSITY_KG_M3: Final[float] = 70.85       # Liquid H2 at -252.87°C, 1 atm



_NH3_LIQUID_DENSITY_KG_M3: Final[float] = 682.0       # Refrigerated liquid NH3 at -33.3°C, 1 atm



_NH3_GAS_DENSITY_0C_KG_NM3: Final[float] = 0.7714     # Gaseous NH3 at 0°C, 101.325 kPa (DIN 1343)



_MEOH_LHV_MJ_PER_KG: Final[float] = 19.90             # Pure Methanol (CH3OH) LHV



_MEOH_HHV_MJ_PER_KG: Final[float] = 22.70             # Pure Methanol (CH3OH) HHV



_MEOH_LIQUID_DENSITY_KG_M3: Final[float] = 792.0      # Liquid Methanol at 15°C, 1 atm






# ---------------------------------------------------------------------------
# 3. Pint Transformation Context Generators
# ---------------------------------------------------------------------------

def create_hydrocarbon_context(
    assay: CrudeAssay | None = None,
    *,
    api_gravity: float | None = None,
    specific_gravity: float | None = None,
    bbl_per_mt: float | None = None,
    name: str | None = None,
) -> pint.Context:
    """Creates a Pint context bridging [volume], [mass], and currency quotes for liquids."""
    from core.units.assays.energy.petroleum import CrudeAssay
    from core.units.registry import CURRENCY_CODES
    resolved_assay = assay or CrudeAssay(
        name=name or "custom_hydrocarbon",
        api_gravity=api_gravity,
        specific_gravity=specific_gravity,
        bbl_per_mt=bbl_per_mt,
    )
    ctx_id = name or f"hydrocarbon_{resolved_assay.name.lower().replace(' ', '_')}"
    ctx = pint.Context(ctx_id)

    rho_kg_m3 = resolved_assay.density_kg_per_m3

    def vol_to_mass(ureg_inst: pint.UnitRegistry, vol: pint.Quantity) -> pint.Quantity:
        vol_m3 = vol.to(ureg_inst.meter ** 3)
        factor = _scale_factor(vol_m3.magnitude, rho_kg_m3)
        return ureg_inst.Quantity(vol_m3.magnitude * factor, ureg_inst.kilogram)

    def mass_to_vol(ureg_inst: pint.UnitRegistry, mass: pint.Quantity) -> pint.Quantity:
        mass_kg = mass.to(ureg_inst.kilogram)
        factor = _scale_factor(mass_kg.magnitude, rho_kg_m3)
        return ureg_inst.Quantity(mass_kg.magnitude / factor, ureg_inst.meter ** 3)

    ctx.add_transformation("[volume]", "[mass]", vol_to_mass)
    ctx.add_transformation("[mass]", "[volume]", mass_to_vol)

    # Currency quote transformations for all known currencies (e.g. EUR/bbl <-> EUR/mt).
    for currency in sorted({_normalize_currency(code) for code in CURRENCY_CODES}):
        cdim = f"[currency_{currency}]"

        def _make_price_vol_to_mass(density: float, cur: str = currency):
            def _transform(ureg_inst: pint.UnitRegistry, price: pint.Quantity) -> pint.Quantity:
                per_m3 = price.to(ureg_inst.Unit(cur) / (ureg_inst.meter ** 3))
                factor = _scale_factor(per_m3.magnitude, density)
                return ureg_inst.Quantity(
                    per_m3.magnitude / factor,
                    ureg_inst.Unit(cur) / ureg_inst.kilogram,
                )

            return _transform

        def _make_price_mass_to_vol(density: float, cur: str = currency):
            def _transform(ureg_inst: pint.UnitRegistry, price: pint.Quantity) -> pint.Quantity:
                per_kg = price.to(ureg_inst.Unit(cur) / ureg_inst.kilogram)
                factor = _scale_factor(per_kg.magnitude, density)
                return ureg_inst.Quantity(
                    per_kg.magnitude * factor,
                    ureg_inst.Unit(cur) / (ureg_inst.meter ** 3),
                )

            return _transform

        ctx.add_transformation(f"{cdim} / [volume]", f"{cdim} / [mass]", _make_price_vol_to_mass(rho_kg_m3))
        ctx.add_transformation(f"{cdim} / [mass]", f"{cdim} / [volume]", _make_price_mass_to_vol(rho_kg_m3))

    return ctx





def create_gas_heating_context(
    assay: GasAssay | None = None,
    *,
    btu_per_scf: float | None = None,
    basis: str = "gcv",
    name: str | None = None,
) -> pint.Context:
    """Creates a Pint context bridging all standard gas reference volumes and [energy]."""
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.registry import ureg
    if basis.lower() not in {"gcv", "hhv", "ncv", "lhv"}:
        raise ValueError(f"Calorific basis must be 'gcv'/'hhv' or 'ncv'/'lhv', got {basis!r}")

    resolved_assay = assay or GasAssay(
        name=name or "custom_gas",
        gcv_btu_per_scf=btu_per_scf if btu_per_scf is not None else 1025.0,
    )
    use_net = basis.lower() in {"ncv", "lhv"}
    base_btu_scf = resolved_assay.ncv_btu_per_scf if use_net else resolved_assay.gcv_btu_per_scf

    ctx_id = name or f"gas_{resolved_assay.name.lower().replace(' ', '_')}_{basis.lower()}"
    ctx = pint.Context(ctx_id)

    ft3_per_m3 = 1.0 / (0.3048 ** 3)
    mj_per_btu = 1055.05585262 / 1e6

    btu_per_scf_ferc = base_btu_scf
    btu_per_scf_tx = base_btu_scf * _DENSITY_RATIO_TX
    mj_per_sm3_iso15 = base_btu_scf * mj_per_btu * ft3_per_m3 * _DENSITY_RATIO_ISO15
    mj_per_nm3_din0 = base_btu_scf * mj_per_btu * ft3_per_m3 * _DENSITY_RATIO_DIN0
    groningen_ratio = (resolved_assay.ncv_btu_per_scf / resolved_assay.gcv_btu_per_scf) if use_net else 1.0
    mj_per_m3_gron = _GRONINGEN_MJ_PER_M3 * groningen_ratio

    gas_specs = (
        ("[standard_gas_volume_20C_101325kPa_dry]", ureg.m3_20C_101325kPa_dry, ureg.megajoule, base_btu_scf * mj_per_btu * ft3_per_m3 * _DENSITY_RATIO_GOST20),
        ("[standard_gas_volume_60F_1atm_dry]", ureg.scf_60F_14696psi_dry, ureg.btu_it, base_btu_scf * (_P_1ATM_KPA / _P_14_73_KPA)),
        ("[standard_gas_volume_60F_14_73psia]", ureg.scf_60F_14_73psia, ureg.btu_it, btu_per_scf_ferc),
        ("[standard_gas_volume_60F_14_65psia]", ureg.scf_60F_14_65psia, ureg.btu_it, btu_per_scf_tx),
        ("[standard_gas_volume_15C_101325kPa_dry]", ureg.Sm3_15C_101325kPa_dry, ureg.megajoule, mj_per_sm3_iso15),
        ("[normal_gas_volume_0C_101325kPa_dry]", ureg.Nm3_0C_101325kPa_dry, ureg.megajoule, mj_per_nm3_din0),
        ("[groningen_gas_volume]", ureg.m3_Groningen, ureg.megajoule, mj_per_m3_gron),
    )

    for dim_str, ref_vol_unit, ref_energy_unit, factor_val in gas_specs:
        def _make_vol_to_energy(v_unit: Unit, e_unit: Unit, f: float):
            def _transform(ureg_inst: pint.UnitRegistry, vol: pint.Quantity) -> pint.Quantity:
                v_base = vol.to(v_unit)
                scale = _scale_factor(v_base.magnitude, f)
                return ureg_inst.Quantity(v_base.magnitude * scale, e_unit)

            return _transform

        def _make_energy_to_vol(v_unit: Unit, e_unit: Unit, f: float):
            def _transform(ureg_inst: pint.UnitRegistry, energy: pint.Quantity) -> pint.Quantity:
                e_base = energy.to(e_unit)
                scale = _scale_factor(e_base.magnitude, f)
                return ureg_inst.Quantity(e_base.magnitude / scale, v_unit)

            return _transform

        ctx.add_transformation(dim_str, "[energy]", _make_vol_to_energy(ref_vol_unit, ref_energy_unit, factor_val))
        ctx.add_transformation("[energy]", dim_str, _make_energy_to_vol(ref_vol_unit, ref_energy_unit, factor_val))

    molar_ratios = {
        "[standard_gas_volume_20C_101325kPa_dry]": (ureg.m3_20C_101325kPa_dry, _DENSITY_RATIO_GOST20 * ft3_per_m3),
        "[standard_gas_volume_60F_1atm_dry]": (ureg.scf_60F_14696psi_dry, _P_1ATM_KPA / _P_14_73_KPA),
        "[standard_gas_volume_60F_14_73psia]": (ureg.scf_60F_14_73psia, _DENSITY_RATIO_FERC),
        "[standard_gas_volume_60F_14_65psia]": (ureg.scf_60F_14_65psia, _DENSITY_RATIO_TX),
        "[standard_gas_volume_15C_101325kPa_dry]": (ureg.Sm3_15C_101325kPa_dry, _DENSITY_RATIO_ISO15 * ft3_per_m3),
        "[normal_gas_volume_0C_101325kPa_dry]": (ureg.Nm3_0C_101325kPa_dry, _DENSITY_RATIO_DIN0 * ft3_per_m3),
    }

    for src_dim, (src_unit, src_moles_per_unit) in molar_ratios.items():
        for dst_dim, (dst_unit, dst_moles_per_unit) in molar_ratios.items():
            if src_dim == dst_dim:
                continue
            conv_ratio = src_moles_per_unit / dst_moles_per_unit

            def _make_vol_to_vol(u_src: Unit, u_dst: Unit, r: float):
                def _transform(ureg_inst: pint.UnitRegistry, q: pint.Quantity) -> pint.Quantity:
                    q_base = q.to(u_src)
                    scale = _scale_factor(q_base.magnitude, r)
                    return ureg_inst.Quantity(q_base.magnitude * scale, u_dst)

                return _transform

            ctx.add_transformation(src_dim, dst_dim, _make_vol_to_vol(src_unit, dst_unit, conv_ratio))

    return ctx





# ---------------------------------------------------------------------------
# Liquid Biofuels & Sustainable Aviation Fuel (SAF) Specification
# (Ethanol, FAME Biodiesel, HVO / Renewable Diesel, HEFA-SPK SAF)
# ---------------------------------------------------------------------------

_LITRES_PER_US_GAL: Final[float] = 3.785411784



_LITRES_PER_BBL: Final[float] = 158.987294928




# ---------------------------------------------------------------------------
# Natural Gas Liquids (NGLs) & LPG Specification (GPA 2145-16)
# (Ethane, Propane, Isobutane, n-Butane, Natural Gasoline / C5+)
# ---------------------------------------------------------------------------

_WATER_DENSITY_60F_KG_M3: Final[float] = 999.016   # Pure water at 60°F (15.5556°C)




# ---------------------------------------------------------------------------
# Battery Metals & Precursor Chemicals Specification
# (Spodumene SC6, Li2CO3, LiOH·H2O, NiSO4·6H2O, CoSO4·7H2O, MnSO4·H2O, MHP)
# ---------------------------------------------------------------------------

# Exact IUPAC stoichiometric conversion factors to Lithium Carbonate Equivalent (LCE, Li2CO3)
_LI_TO_LCE_FACTOR: Final[float] = 5.32279        # M(Li2CO3) / (2 * M(Li))   = 73.891 / 13.882



_LI2O_TO_LCE_FACTOR: Final[float] = 2.47284      # M(Li2CO3) / M(Li2O)       = 73.891 / 29.881



_LIOH_H2O_TO_LCE_FACTOR: Final[float] = 0.88041  # M(Li2CO3) / (2*M(LiOH·H2O)) = 73.891 / 83.928




# ---------------------------------------------------------------------------
# Nuclear Fuel Cycle Specification (Uranium — U3O8, UF6, EUP & SWU)
# ---------------------------------------------------------------------------

# Exact IUPAC / WNA stoichiometric uranium mass fractions
_U_FRAC_IN_U3O8: Final[float] = 0.84800    # kgU per kg U3O8



_U_FRAC_IN_UF6: Final[float] = 0.67618     # kgU per kg UF6



_U_FRAC_IN_UO2: Final[float] = 0.88150     # kgU per kg UO2



_NATURAL_U235_PCT: Final[float] = 0.711    # Natural uranium U-235 isotopic weight %





def _dirac_value_function(x_frac: float) -> float:
    """Computes the Dirac separation potential V(x) = (2x - 1) * ln(x / (1 - x))."""
    if not (0.0 < x_frac < 1.0):
        raise ValueError(f"Isotopic fraction must be strictly in (0, 1), got {x_frac}")
    return (2.0 * x_frac - 1.0) * math.log(x_frac / (1.0 - x_frac))




# ---------------------------------------------------------------------------
# Solid Biomass & Wood Pellets Specification (ISO 17225 / ISO 18125)
# (Industrial Pellets I2P, ENplus A1, Torrefied Pellets, Wood Chips, PKS)
# ---------------------------------------------------------------------------

# ISO 18125 latent heat constants for solid biofuels (GJ/mt per 1 wt%, i.e. MJ/kg per 1 wt%)
_ISO18125_H_GJ_PER_MT_PCT: Final[float] = 0.2122    # Latent heat of water from H combustion (dry basis)



_ISO18125_M_GJ_PER_MT_PCT: Final[float] = 0.02443   # Latent heat of vaporisation of moisture at 25°C




# ---------------------------------------------------------------------------
# Petrochemicals Specification (Olefins, Aromatics & Intermediates)
# (Ethylene, Propylene, Benzene, Toluene, Paraxylene, Styrene, MEG)
# ---------------------------------------------------------------------------

_LB_PER_MT: Final[float] = 1000.0 / 0.45359237    # Exact: 2204.62262185 lb/mt




# ---------------------------------------------------------------------------
# Precious Metals & Doré Specification (LBMA / LPPM / COMEX)
# (Gold, Silver, Platinum, Palladium, Mine-Site Doré)
# ---------------------------------------------------------------------------

_GRAMS_PER_TROY_OUNCE: Final[float] = 31.1034768              # Exact international troy ounce (g)



_TROY_OUNCES_PER_KG: Final[float] = 1000.0 / 31.1034768       # 32.15074657 ozt/kg



_GRAMS_PER_TOLA: Final[float] = 11.6638038                    # Exact India/Dubai bullion tola (3/8 ozt)




# ---------------------------------------------------------------------------
# Fertilisers & Nutrient Equivalents Specification
# (Urea, UAN 28%/32%, DAP, MAP, Potash MOP, Ammonium Sulphate)
# ---------------------------------------------------------------------------

# Exact IUPAC elemental-to-oxide nutrient stoichiometry
_P_IN_P2O5: Final[float] = 0.43642          # 2 * M(P) / M(P2O5) = 61.948 / 141.945



_K_IN_K2O: Final[float] = 0.83015           # 2 * M(K) / M(K2O)  = 78.197 / 94.196



_MT_PER_SHORT_TON: Final[float] = 0.90718474  # Exact: 2000 lb in metric tonnes




# ---------------------------------------------------------------------------
# Paving Bitumen (Asphalt) Specification (ASTM D4311 / EN 12591)
# (Pen 60/70, Pen 80/100, PG 64-22, Oxidised Roofing Bitumen)
# ---------------------------------------------------------------------------

_ASTM_D4311_GROUP0_ALPHA_PER_C: Final[float] = 0.00062  # Standard Group 0 bitumen thermal expansion (/°C)





# ---------------------------------------------------------------------------
# Hydrogen & Clean Energy Carriers Specification
# (Gaseous H2, Compressed 350/700 bar, Liquid LH2, Ammonia NH3, LOHC)
# ---------------------------------------------------------------------------

_NH3_H2_MASS_FRACTION: Final[float] = 0.17755  # Approximate stoichiometric fraction: 3 * M(H) / M(NH3) = 3.0238 / 17.0305



# CI screening only: 70% saving against 94 g CO2e/MJ (LHV).
# https://observatory.clean-hydrogen.europa.eu/index.php/eu-policy/renewable-energy-directive
_EU_RFNBO_MAX_CI_G_PER_MJ: Final[float] = 28.2





# ---------------------------------------------------------------------------
# Ferroalloys Specification (Bulk & Noble Steelmaking Alloys)
# (HC/LC FeCr, HC FeMn, SiMn, FeSi 75, FeMo 60, FeV 80, FeNb 65)
# ---------------------------------------------------------------------------

_LB_PER_KG: Final[float] = 1.0 / 0.45359237  # Exact: 2.20462262185 lb/kg





# ---------------------------------------------------------------------------
# Pint Context Factories for Commodity Assays
# ---------------------------------------------------------------------------

def _default_ctx_name(assay_name: str, suffix: str | None = None) -> str:
    slug = assay_name.strip().lower().replace(" ", "_").replace("/", "_")
    return f"{slug}_{suffix}" if suffix else slug





def create_coal_heating_context(
    assay: CoalAssay,
    *,
    basis: str = "nar",
    name: str | None = None,
) -> pint.Context:
    """Bridges as-received coal mass and energy using NAR or GAR heating value."""
    from core.units.registry import ureg
    key = basis.strip().lower()
    if key in {"nar", "lhv", "ncv"}:
        gj_per_mt = assay.gj_per_mt_nar
    elif key in {"gar", "hhv", "gcv"}:
        gj_per_mt = assay.gj_per_mt_gar
    else:
        raise ValueError("Coal context basis must be 'nar'/'lhv'/'ncv' or 'gar'/'hhv'/'gcv'.")
    ctx = pint.Context(name or _default_ctx_name(assay.name, key))

    def mass_to_energy(_u: pint.UnitRegistry, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        tonnes = mass.to(ureg.metric_ton)
        factor = _scale_factor(tonnes.magnitude, gj_per_mt)
        return ureg.Quantity(tonnes.magnitude * factor, ureg.gigajoule)

    def energy_to_mass(_u: pint.UnitRegistry, energy: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        gj = energy.to(ureg.gigajoule)
        factor = _scale_factor(gj.magnitude, gj_per_mt)
        return ureg.Quantity(gj.magnitude / factor, ureg.metric_ton)

    ctx.add_transformation("[mass]", "[energy]", mass_to_energy)
    ctx.add_transformation("[energy]", "[mass]", energy_to_mass)
    return ctx





def create_clean_fuel_context(
    assay: CleanFuelAssay,
    *,
    basis: str = "lhv",
    name: str | None = None,
) -> pint.Context:
    """Bridges mass, energy and the assay's available liquid/normal gas volumes."""
    from core.units.registry import ureg
    assay._heating_value_mj_per_kg(basis)
    ctx = pint.Context(name or _default_ctx_name(assay.name, basis.strip().lower()))
    ctx.add_transformation("[mass]", "[energy]", lambda _u, q: assay.mass_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[mass]", lambda _u, q: assay.energy_to_mass(q, basis=basis))
    if assay.liquid_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.liquid_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_liquid_volume(q))
    if assay.density_kg_per_nm3 is not None:
        normal_dim = ureg.Nm3_0C_101325kPa_dry.dimensionality
        ctx.add_transformation(normal_dim, "[mass]", lambda _u, q: assay.gas_volume_to_mass(q))
        ctx.add_transformation("[mass]", normal_dim, lambda _u, q: assay.mass_to_gas_volume(q))
    return ctx





def create_biofuel_context(
    assay: BiofuelAssay,
    *,
    basis: str = "lhv",
    name: str | None = None,
) -> pint.Context:
    """Bridges [volume] <-> [mass] <-> [energy] for a BiofuelAssay."""
    assay._heating_value_mj_per_kg(basis)  # Validate basis upfront
    ctx = pint.Context(name or _default_ctx_name(assay.name, basis.lower()))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    ctx.add_transformation("[mass]", "[energy]", lambda _u, q: assay.mass_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[mass]", lambda _u, q: assay.energy_to_mass(q, basis=basis))
    ctx.add_transformation("[volume]", "[energy]", lambda _u, q: assay.volume_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[volume]", lambda _u, q: assay.energy_to_volume(q, basis=basis))
    return ctx





def create_ngl_context(
    assay: NGLAssay,
    *,
    basis: str = "hhv",
    name: str | None = None,
) -> pint.Context:
    """Bridges [volume] <-> [mass] <-> [energy] for an NGLAssay (60°F pressurised basis)."""
    assay._btu_per_gal(basis)  # Validate basis upfront
    ctx = pint.Context(name or _default_ctx_name(assay.name, basis.lower()))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    ctx.add_transformation("[mass]", "[energy]", lambda _u, q: assay.mass_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[mass]", lambda _u, q: assay.energy_to_mass(q, basis=basis))
    ctx.add_transformation("[volume]", "[energy]", lambda _u, q: assay.volume_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[volume]", lambda _u, q: assay.energy_to_volume(q, basis=basis))
    return ctx





def create_ore_context(
    assay: OreAssay,
    *,
    basis: str = "contained",
    name: str | None = None,
) -> pint.Context:
    """Bridges bulk stowage [volume] <-> WMT [mass] for an OreAssay.

    Note: Same-dimension mass conversions (WMT <-> DMT <-> contained/payable metal)
    must use the explicit OreAssay methods.
    """
    ctx = pint.Context(name or _default_ctx_name(assay.name, basis.lower()))
    if assay.bulk_density_wmt_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_wmt(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.wmt_to_volume(q))
    return ctx





def create_battery_material_context(
    assay: BatteryMaterialAssay,
    *,
    basis: str = "contained",
    name: str | None = None,
) -> pint.Context:
    """Returns a named Pint context for a BatteryMaterialAssay.

    Note: All battery material conversions (WMT <-> DMT <-> Metal <-> LCE) are
    [mass] <-> [mass] and must use the explicit BatteryMaterialAssay methods.
    """
    return pint.Context(name or _default_ctx_name(assay.name, basis.lower()))





def create_uranium_context(
    assay: UraniumAssay,
    *,
    burnup_gwd_per_mtu: float = 45.0,
    thermal_efficiency: float | None = None,
    name: str | None = None,
) -> pint.Context:
    """Bridges EUP [mass] (kgU) <-> reactor [energy] (MWh) for a UraniumAssay.

    If thermal_efficiency is None, bridges to thermal energy (MWh_th); otherwise
    bridges to net electrical energy (MWh_e).
    """
    from core.units.registry import ureg
    if not math.isfinite(burnup_gwd_per_mtu) or burnup_gwd_per_mtu <= 0.0:
        raise ValueError(f"burnup_gwd_per_mtu must be strictly positive, got {burnup_gwd_per_mtu}")
    if thermal_efficiency is not None and (not math.isfinite(thermal_efficiency) or not (0.0 < thermal_efficiency <= 1.0)):
        raise ValueError(f"thermal_efficiency must be in (0, 1], got {thermal_efficiency}")

    def mass_to_energy(_u: pint.UnitRegistry, mass: pint.Quantity) -> pint.Quantity:
        if thermal_efficiency is None:
            return assay.eup_to_thermal_energy(mass, burnup_gwd_per_mtu=burnup_gwd_per_mtu)
        return assay.eup_to_electrical_energy(
            mass,
            burnup_gwd_per_mtu=burnup_gwd_per_mtu,
            thermal_efficiency=thermal_efficiency,
        )

    # Derive the inverse from the assay's own forward conversion, preserving
    # Decimal magnitudes through quantity arithmetic.
    reference_mass = ureg.Quantity(1, ureg.kilogram)
    reference_energy = mass_to_energy(ureg, reference_mass)

    def energy_to_mass(_u: pint.UnitRegistry, energy: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        energy_q = energy.to(reference_energy.units)
        scale = _scale_factor(energy_q.magnitude, float(reference_energy.magnitude))
        return ureg.Quantity(energy_q.magnitude / scale, reference_mass.units)

    ctx = pint.Context(name or _default_ctx_name(assay.name))
    ctx.add_transformation("[mass]", "[energy]", mass_to_energy)
    ctx.add_transformation("[energy]", "[mass]", energy_to_mass)
    return ctx





def create_biomass_context(
    assay: BiomassAssay,
    *,
    basis: str = "lhv",
    name: str | None = None,
) -> pint.Context:
    """Bridges bulk [volume] <-> WMT [mass] <-> [energy] for a BiomassAssay."""
    assay._gj_per_mt(basis)  # Validate basis upfront
    ctx = pint.Context(name or _default_ctx_name(assay.name, basis.lower()))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_wmt(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.wmt_to_bulk_volume(q))
    ctx.add_transformation("[mass]", "[energy]", lambda _u, q: assay.mass_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[mass]", lambda _u, q: assay.energy_to_mass(q, basis=basis))
    ctx.add_transformation("[volume]", "[energy]", lambda _u, q: assay.bulk_volume_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[volume]", lambda _u, q: assay.energy_to_bulk_volume(q, basis=basis))
    return ctx





def create_petrochem_context(
    assay: PetrochemAssay,
    *,
    active_only: bool = False,
    name: str | None = None,
) -> pint.Context:
    """Bridges [volume] <-> [mass] <-> [substance] (moles) for a PetrochemAssay."""
    ctx = pint.Context(name or _default_ctx_name(assay.name))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    ctx.add_transformation("[mass]", "[substance]", lambda _u, q: assay.mass_to_moles(q, active_only=active_only))
    ctx.add_transformation("[substance]", "[mass]", lambda _u, q: assay.moles_to_mass(q, active_only=active_only))
    ctx.add_transformation(
        "[volume]",
        "[substance]",
        lambda _u, q: assay.mass_to_moles(assay.volume_to_mass(q), active_only=active_only),
    )
    ctx.add_transformation(
        "[substance]",
        "[volume]",
        lambda _u, q: assay.mass_to_volume(assay.moles_to_mass(q, active_only=active_only)),
    )
    return ctx





def create_bullion_context(
    assay: BullionAssay,
    *,
    basis: str = "fine",
    name: str | None = None,
) -> pint.Context:
    """Returns a named Pint context for a BullionAssay.

    Note: Gross <-> fine <-> payable conversions are [mass] <-> [mass] and must
    use the explicit BullionAssay methods.
    """
    return pint.Context(name or _default_ctx_name(assay.name, basis.lower()))





def create_fertilizer_context(
    assay: FertilizerAssay,
    *,
    nutrient: str = "N",
    name: str | None = None,
) -> pint.Context:
    """Bridges liquid/bulk [volume] <-> [mass] for a FertilizerAssay.

    Note: Product <-> nutrient mass conversions are [mass] <-> [mass] and must
    use product_to_nutrient() / nutrient_to_product().
    """
    ctx = pint.Context(name or _default_ctx_name(assay.name, nutrient.lower()))
    if assay.density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    return ctx





def create_bitumen_context(
    assay: BitumenAssay,
    *,
    name: str | None = None,
) -> pint.Context:
    """Bridges standard-temperature [volume] <-> [mass] for a BitumenAssay."""
    ctx = pint.Context(name or _default_ctx_name(assay.name))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    return ctx





def create_petcoke_context(
    assay: PetcokeAssay,
    *,
    basis: str = "lhv",
    name: str | None = None,
) -> pint.Context:
    """Bridges bulk [volume] <-> WMT [mass] <-> [energy] for a PetcokeAssay."""
    if basis.strip().lower() not in {"lhv", "nar", "hhv", "gar"}:
        raise ValueError(f"basis must be 'lhv'/'nar' or 'hhv'/'gar', got {basis!r}")
    ctx = pint.Context(name or _default_ctx_name(assay.name, basis.lower()))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_wmt(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.wmt_to_bulk_volume(q))
    ctx.add_transformation("[mass]", "[energy]", lambda _u, q: assay.mass_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[mass]", lambda _u, q: assay.energy_to_mass(q, basis=basis))
    ctx.add_transformation(
        "[volume]",
        "[energy]",
        lambda _u, q: assay.mass_to_energy(assay.bulk_volume_to_wmt(q), basis=basis),
    )
    ctx.add_transformation(
        "[energy]",
        "[volume]",
        lambda _u, q: assay.wmt_to_bulk_volume(assay.energy_to_mass(q, basis=basis)),
    )
    return ctx





def create_chemical_context(
    assay: ChemicalAssay,
    *,
    active_only: bool = True,
    name: str | None = None,
) -> pint.Context:
    """Bridges commercial [volume] <-> [mass] <-> [substance] (moles) for a ChemicalAssay.

    Note: Solution WMT <-> 100% active DMT conversions are [mass] <-> [mass] and
    must use solution_to_active_mass() / active_to_solution_mass().
    """
    ctx = pint.Context(name or _default_ctx_name(assay.name))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    ctx.add_transformation("[mass]", "[substance]", lambda _u, q: assay.mass_to_moles(q, active_only=active_only))
    ctx.add_transformation("[substance]", "[mass]", lambda _u, q: assay.moles_to_mass(q, active_only=active_only))
    ctx.add_transformation(
        "[volume]",
        "[substance]",
        lambda _u, q: assay.mass_to_moles(assay.volume_to_mass(q), active_only=active_only),
    )
    ctx.add_transformation(
        "[substance]",
        "[volume]",
        lambda _u, q: assay.mass_to_volume(assay.moles_to_mass(q, active_only=active_only)),
    )
    return ctx





def create_hydrogen_context(assay: HydrogenAssay, *, basis: str = "lhv", name: str | None = None,) -> pint.Context:
    """Bridges storage volume, product/carrier mass, delivered H2 energy and Nm3.

    Nm3 denotes delivered H2 at NTP; m^3 denotes the assay's storage state.
    Carrier-to-H2 mass conversions require the explicit assay methods.
    """
    from core.units.registry import ureg
    assay._mj_per_kg_h2(basis)
    ctx = pint.Context(name or _default_ctx_name(assay.name, basis.strip().lower()))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    ctx.add_transformation("[mass]", "[energy]", lambda _u, q: assay.mass_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[mass]", lambda _u, q: assay.energy_to_mass(q, basis=basis))
    ctx.add_transformation("[volume]", "[energy]", lambda _u, q: assay.volume_to_energy(q, basis=basis))
    ctx.add_transformation("[energy]", "[volume]", lambda _u, q: assay.energy_to_volume(q, basis=basis))
    normal_dim = ureg.Nm3_H2.dimensionality
    ctx.add_transformation(normal_dim, "[mass]", lambda _u, q: assay.h2_to_carrier_mass(assay.nm3_to_h2_mass(q)))
    ctx.add_transformation("[mass]", normal_dim, lambda _u, q: assay.h2_mass_to_nm3(assay.carrier_to_h2_mass(q)))
    return ctx





def create_ferroalloy_context(assay: FerroalloyAssay, *, name: str | None = None) -> pint.Context:
    """Bridges bulk volume and gross mass when bulk density is set.

    Contained-metal conversions require explicit assay methods.
    """
    ctx = pint.Context(name or _default_ctx_name(assay.name))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx





def create_metal_alloy_context(assay: MetalAlloyAssay, *, name: str | None = None) -> pint.Context:
    """Bridges solid volume and gross mass; remelt/element masses use assay methods."""
    ctx = pint.Context(name or _default_ctx_name(assay.name))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.solid_volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_solid_volume(q))
    return ctx





# ---------------------------------------------------------------------------
# 4. High-Level Functional Converters
# ---------------------------------------------------------------------------


def convert_hydrocarbon(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: CrudeAssay | None = None,
    api_gravity: float | None = None,
    specific_gravity: float | None = None,
    bbl_per_mt: float | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts a liquid hydrocarbon quantity (volume, mass, flow rate, or price) to target_unit."""
    from core.units.assays.energy.petroleum import CrudeAssay
    from core.units.precision import ensure_quantity_ieee754
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)

    resolved_target = parse_unit(target_unit)
    eff_assay = assay
    if eff_assay is None and any(v is not None for v in (api_gravity, specific_gravity, bbl_per_mt)):
        eff_assay = CrudeAssay(api_gravity=api_gravity, specific_gravity=specific_gravity, bbl_per_mt=bbl_per_mt)

    if quantity.dimensionality == resolved_target.dimensionality:
        out_q = ensure_quantity_ieee754(quantity.to(resolved_target))
        if not return_metadata:
            return out_q
        from core.units.provenance import ConversionClassification
        from core.units.provenance import build_conversion_result
        return build_conversion_result(
            out_q, quantity, method="convert_hydrocarbon",
            classification=ConversionClassification.EXACT,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    ctx = create_hydrocarbon_context(
        assay=eff_assay,
        api_gravity=api_gravity,
        specific_gravity=specific_gravity,
        bbl_per_mt=bbl_per_mt,
    )
    with ureg.context(ctx):
        out_q = ensure_quantity_ieee754(_convert_contextual_quantity(quantity, resolved_target, [ctx]))
    if not return_metadata:
        return out_q
    return convert(
        quantity,
        resolved_target,
        assay=eff_assay,
        return_metadata=True,
        strict=strict,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
    )





def convert_gas_energy(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: GasAssay | None = None,
    btu_per_scf: float | None = None,
    basis: str = "gcv",
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts a natural gas quantity between reference volumes, energy, flow rates, or prices."""
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.precision import ensure_quantity_ieee754
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)

    resolved_target = parse_unit(target_unit)
    eff_assay = assay if assay is not None else (GasAssay(name="custom_gas", gcv_btu_per_scf=btu_per_scf) if btu_per_scf is not None else None)
    if quantity.dimensionality == resolved_target.dimensionality:
        out_q = ensure_quantity_ieee754(quantity.to(resolved_target))
        if not return_metadata:
            return out_q
        from core.units.provenance import ConversionClassification
        from core.units.provenance import build_conversion_result
        return build_conversion_result(
            out_q, quantity, method="convert_gas_energy",
            classification=ConversionClassification.EXACT,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    ctx = create_gas_heating_context(assay=eff_assay, btu_per_scf=btu_per_scf, basis=basis)
    with ureg.context(ctx):
        out_q = ensure_quantity_ieee754(_convert_contextual_quantity(quantity, resolved_target, [ctx]))
    if not return_metadata:
        return out_q
    from core.units.provenance import ConversionClassification
    from core.units.provenance import build_conversion_result
    mass_dim = ureg.kilogram.dimensionality
    is_lng = eff_assay is not None and (quantity.dimensionality == mass_dim or resolved_target.dimensionality == mass_dim)
    return build_conversion_result(
        out_q,
        quantity,
        assay=eff_assay,
        commodity_id="lng" if is_lng else "natural_gas",
        heating_basis=basis,
        conversion_path=(str(quantity.units), f"gas_heating:{basis}", str(out_q.units)),
        method="convert_gas_energy",
        classification=ConversionClassification.COMPOSITION_DEPENDENT,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
    )





def convert_hydrogen(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: HydrogenAssay,
    basis: str = "lhv",
    mode: str = "auto",
    cracking_efficiency: float | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts hydrogen or carrier quantities across volume, mass, Nm^3, and energy (MWh, GJ, MMBtu)."""
    from core.units.assays.energy.hydrogen import HydrogenAssay
    from core.units.exceptions import UnsupportedConversionError
    from core.units.precision import ensure_quantity_ieee754
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(assay, HydrogenAssay):
        raise TypeError(f"Expected HydrogenAssay, got {type(assay).__name__}")
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)

    resolved_target = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        ctx = assay.as_context(basis=basis)
        with ureg.context(ctx):
            out_q = ensure_quantity_ieee754(_convert_contextual_quantity(quantity, resolved_target, [ctx]))
        if not return_metadata:
            return out_q
        from core.units.provenance import ConversionClassification
        from core.units.provenance import build_conversion_result
        cls_tag = ConversionClassification.EXACT if quantity.dimensionality == resolved_target.dimensionality else ConversionClassification.COMPOSITION_DEPENDENT
        return build_conversion_result(
            out_q, quantity, assay=assay, heating_basis=basis,
            method="convert_hydrogen", classification=cls_tag,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "carrier_to_h2": lambda q: assay.carrier_to_h2_mass(q, cracking_efficiency=cracking_efficiency),
        "h2_to_carrier": lambda q: assay.h2_to_carrier_mass(q, cracking_efficiency=cracking_efficiency),
        "h2_to_nm3": lambda q: HydrogenAssay.h2_mass_to_nm3(q),
        "nm3_to_h2": lambda q: HydrogenAssay.nm3_to_h2_mass(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported hydrogen mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    out_q = ensure_quantity_ieee754(dispatch[m](quantity).to(resolved_target))
    if not return_metadata:
        return out_q
    from core.units.provenance import ConversionClassification
    from core.units.provenance import build_conversion_result
    cls_tag = ConversionClassification.YIELD_DEPENDENT if "carrier" in m else ConversionClassification.COMPOSITION_DEPENDENT
    return build_conversion_result(
        out_q, quantity, assay=assay, heating_basis=basis,
        conversion_path=(str(quantity.units), m, str(out_q.units)),
        method=f"convert_hydrogen:{m}", classification=cls_tag,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
    )





def convert_biofuel(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: BiofuelAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using BiofuelAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.energy.biofuels import BiofuelAssay
    if not isinstance(assay, BiofuelAssay):
        raise TypeError("assay must be a BiofuelAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_ngl(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: NGLAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using NGLAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.energy.ngl import NGLAssay
    if not isinstance(assay, NGLAssay):
        raise TypeError("assay must be a NGLAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_ore(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: OreAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using OreAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.metals.ores import OreAssay
    if not isinstance(assay, OreAssay):
        raise TypeError("assay must be a OreAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_battery_material(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: BatteryMaterialAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using BatteryMaterialAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.metals.battery_materials import BatteryMaterialAssay
    if not isinstance(assay, BatteryMaterialAssay):
        raise TypeError("assay must be a BatteryMaterialAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_uranium(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: UraniumAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using UraniumAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.energy.nuclear import UraniumAssay
    if not isinstance(assay, UraniumAssay):
        raise TypeError("assay must be a UraniumAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_biomass(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: BiomassAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using BiomassAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.energy.biomass import BiomassAssay
    if not isinstance(assay, BiomassAssay):
        raise TypeError("assay must be a BiomassAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_petrochem(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: PetrochemAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using PetrochemAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.industrial.petrochemicals import PetrochemAssay
    if not isinstance(assay, PetrochemAssay):
        raise TypeError("assay must be a PetrochemAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_bullion(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: BullionAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using BullionAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.metals.precious_metals import BullionAssay
    if not isinstance(assay, BullionAssay):
        raise TypeError("assay must be a BullionAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_fertilizer(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: FertilizerAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using FertilizerAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.industrial.fertilisers import FertilizerAssay
    if not isinstance(assay, FertilizerAssay):
        raise TypeError("assay must be a FertilizerAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_bitumen(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: BitumenAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using BitumenAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.energy.bitumen import BitumenAssay
    if not isinstance(assay, BitumenAssay):
        raise TypeError("assay must be a BitumenAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_petcoke(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: PetcokeAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using PetcokeAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.energy.petcoke import PetcokeAssay
    if not isinstance(assay, PetcokeAssay):
        raise TypeError("assay must be a PetcokeAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_chemical(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: ChemicalAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using ChemicalAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.industrial.chemicals import ChemicalAssay
    if not isinstance(assay, ChemicalAssay):
        raise TypeError("assay must be a ChemicalAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_ferroalloy(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: FerroalloyAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using FerroalloyAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.metals.ferroalloys import FerroalloyAssay
    if not isinstance(assay, FerroalloyAssay):
        raise TypeError("assay must be a FerroalloyAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_metal_alloy(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: MetalAlloyAssay,
    **kwargs: Any,
) -> pint.Quantity:
    """Converts units using MetalAlloyAssay through the master router.

    Same-dimension composition/recovery changes require explicit assay methods,
    except labelled ore moisture and contained-metal bases handled by convert().
    """
    from core.units.assays.metals.alloys import MetalAlloyAssay
    if not isinstance(assay, MetalAlloyAssay):
        raise TypeError("assay must be a MetalAlloyAssay.")
    return convert(quantity, target_unit, assay=assay, **kwargs)





def convert_commodity(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: Any = None,
    **kwargs: Any,
) -> pint.Quantity:
    """Generic dispatcher that routes supported assays (or raw hydrocarbon/gas kwargs) to converters.

    Adapted to the currently defined assay types and converter functions in this module.
    """
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)

    # Dispatch only for assay types and converter functions that exist in this module.
    assay_converter_map: tuple[tuple[str, str], ...] = (
        ("BiofuelAssay", "convert_biofuel"),
        ("NGLAssay", "convert_ngl"),
        ("OreAssay", "convert_ore"),
        ("BatteryMaterialAssay", "convert_battery_material"),
        ("UraniumAssay", "convert_uranium"),
        ("BiomassAssay", "convert_biomass"),
        ("PetrochemAssay", "convert_petrochem"),
        ("BullionAssay", "convert_bullion"),
        ("FertilizerAssay", "convert_fertilizer"),
        ("BitumenAssay", "convert_bitumen"),
        ("PetcokeAssay", "convert_petcoke"),
        ("ChemicalAssay", "convert_chemical"),
        ("HydrogenAssay", "convert_hydrogen"),
        ("FerroalloyAssay", "convert_ferroalloy"),
        ("MetalAlloyAssay", "convert_metal_alloy"),
        ("BulkOreAssay", "convert_bulk_ore"),
        ("GasAssay", "convert_gas_energy"),
        ("CrudeAssay", "convert_hydrocarbon"),
    )
    module_globals = globals()
    if assay is not None:
        for assay_name, converter_name in assay_converter_map:
            assay_type = module_globals.get(assay_name)
            converter = module_globals.get(converter_name)
            if assay_type is not None and converter is not None and isinstance(assay, assay_type):
                return converter(quantity, target_unit, assay=assay, **kwargs)

        # Generic fallback for known assays that provide Pint context transformations.
        if hasattr(assay, "as_context"):
            resolved_target = parse_unit(target_unit)
            method = assay.as_context
            parameters = inspect.signature(method).parameters
            accepts_kwargs = any(p.kind == p.VAR_KEYWORD for p in parameters.values())
            context_kwargs = {k: v for k, v in kwargs.items() if k in {"basis", "nutrient", "active_only"} and (accepts_kwargs or k in parameters)}
            ctx = method(**context_kwargs)
            if ctx is None:
                return quantity.to(resolved_target)
            with ureg.context(ctx):
                return _convert_contextual_quantity(quantity, resolved_target, [ctx])

    if assay is None and any(k in kwargs for k in ("btu_per_scf", "basis")) and "api_gravity" not in kwargs:
        return convert_gas_energy(quantity, target_unit, assay=None, **kwargs)
    return convert_hydrocarbon(quantity, target_unit, assay=assay, **kwargs)





# ---------------------------------------------------------------------------
# 5. Pressure Datum Transformations (Gauge <-> Absolute) & Barometric Helpers
# ---------------------------------------------------------------------------

# STANDARD_ATMOSPHERE_PA is defined with the physical reference constants above.
STANDARD_ATMOSPHERE_PSI: Final[float] = 14.695948775513449



STANDARD_ATMOSPHERE_BAR: Final[float] = 1.01325




_GAUGE_PRESSURE_UNITS: Final[frozenset[str]] = frozenset({
    "psig", "barg", "kpag", "mpag", "pag", "psi_g", "bar_g",
    "kpa_g", "mpa_g", "pa_g", "lb_per_sq_in_gauge", "bar_gauge",
    "pound_per_square_inch_gauge", "pound_force_per_square_inch_gauge",
    "kilopascal_gauge", "megapascal_gauge", "pascal_gauge",
})



_ABSOLUTE_PRESSURE_UNITS: Final[frozenset[str]] = frozenset({
    "psia", "bara", "kpaa", "mpaa", "paa", "psi_a", "bar_a",
    "kpa_a", "mpa_a", "pa_a", "lb_per_sq_in_abs", "bar_absolute",
    "pound_per_square_inch_absolute", "pound_force_per_square_inch_absolute",
    "kilopascal_absolute", "megapascal_absolute", "pascal_absolute",
})



_NEUTRAL_PRESSURE_UNITS: Final[frozenset[str]] = frozenset({
    "psi", "pound_force_per_square_inch", "bar", "kpa", "kilopascal",
    "mpa", "megapascal", "pa", "pascal", "torr", "atm", "atmosphere",
    "mmhg", "inhg", "millimeter_hg", "inch_hg",
})





def _unit_tokens(unit: pint.Unit) -> set[str]:
    """Reads canonical names independently of the registry's display format.

    Inspect the unit's structural names rather than expanding to base units,
    which would erase gauge/absolute labels. Aliases resolve via get_name().
    """
    from core.units.registry import ureg
    tokens = {str(unit).strip().lower(), f"{unit:~}".strip().lower()}
    for name in unit._units:
        tokens.add(ureg.get_name(name).strip().lower())
    return tokens





def _classify_pressure_datum(unit: pint.Unit) -> str:
    """Classifies labelled pressure units; unqualified pressure is neutral."""
    tokens = _unit_tokens(unit)
    gauge = bool(tokens & _GAUGE_PRESSURE_UNITS)
    absolute = bool(tokens & _ABSOLUTE_PRESSURE_UNITS)
    if gauge or absolute:
        if gauge and absolute or len(unit._units) != 1 or next(iter(unit._units.values())) != 1:
            raise ValueError("Datum-labelled pressure units must be standalone units.")
        return "gauge" if gauge else "absolute"
    return "neutral"





def atmospheric_pressure_at_elevation(
    elevation: float | pint.Quantity,
    *,
    elevation_unit: str | pint.Unit = "meter",
    sea_level_pa: float = STANDARD_ATMOSPHERE_PA,
    unit: str | pint.Unit | None = None,
) -> pint.Quantity:
    """Returns ISA tropospheric pressure (Pa), valid from -500 to 11000 m.

    Numeric elevations use elevation_unit (or its alias unit); quantities carry
    their own units. This standard-atmosphere estimate is not a weather reading.
    """
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    sea_pa = float(sea_level_pa)
    if not math.isfinite(sea_pa) or sea_pa <= 0.0:
        raise ValueError("sea_level_pa must be finite and strictly positive.")
    if isinstance(elevation, pint.Quantity):
        validate_registry_quantity(elevation)
        height = elevation
    else:
        height = ureg.Quantity(elevation, parse_unit(unit or elevation_unit))
    h_m = float(height.to(ureg.meter).magnitude)
    if not math.isfinite(h_m) or not (-500.0 <= h_m <= 11000.0):
        raise ValueError("Elevation in metres must be finite and in [-500, 11000].")
    air_specific_gas_constant = UNIVERSAL_GAS_CONSTANT_J_MOL_K / (DRY_AIR_MOLAR_MASS_KG_KMOL / 1000.0)
    exponent = STANDARD_GRAVITY_M_S2 / (air_specific_gas_constant * BAROMETRIC_LAPSE_RATE_K_M)
    pressure = sea_pa * (1.0 - BAROMETRIC_LAPSE_RATE_K_M * h_m / STANDARD_SEA_LEVEL_TEMP_K) ** exponent
    return ureg.Quantity(pressure, ureg.pascal)


from core.units.validation import authorized_conversion


@authorized_conversion
def convert_pressure(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    p_atm: pint.Quantity | float | Decimal | None = None,
    atmospheric_pressure: pint.Quantity | float | Decimal | None = None,
    neutral_is_absolute: bool = True,
    allow_negative_absolute: bool = False,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts pressure scales with explicit gauge/absolute datum shifts.

    Unqualified units (psi, bar, Pa, torr, etc.) mean absolute pressure by
    default. With neutral_is_absolute=False, either unqualified endpoint raises
    ValueError; use explicit datum labels instead. Negative gauge pressure is
    valid down to perfect vacuum. Absolute pressure is checked even for same-
    datum conversions unless allow_negative_absolute=True is explicitly set.

    p_atm is the positive local absolute atmosphere, in Pa for numeric inputs;
    it defaults to 101325 Pa. atmospheric_pressure is its legacy alias, with
    p_atm taking precedence when both are supplied. Scalar Decimal magnitudes
    retain Decimal atmospheric precision. Datum shifts are not applied by .to().
    """
    from core.units.exceptions import InvalidPhysicalStateError
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.precision import ensure_quantity_ieee754
    from core.units.precision import pad_ieee754_decimal
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)
    if isinstance(p_atm, pint.Quantity):
        validate_registry_quantity(p_atm)
    if isinstance(atmospheric_pressure, pint.Quantity):
        validate_registry_quantity(atmospheric_pressure)
    resolved_target = parse_unit(target_unit)
    pressure_dim = ureg.pascal.dimensionality
    if quantity.dimensionality != pressure_dim:
        raise pint.DimensionalityError(quantity.units, ureg.pascal, quantity.dimensionality, pressure_dim)
    if resolved_target.dimensionality != pressure_dim:
        raise pint.DimensionalityError(quantity.units, resolved_target, pressure_dim, resolved_target.dimensionality)

    raw_src_datum = _classify_pressure_datum(quantity.units)
    raw_dst_datum = _classify_pressure_datum(resolved_target)
    if not neutral_is_absolute and "neutral" in (raw_src_datum, raw_dst_datum):
        raise MissingReferenceConditionsError("Unqualified pressure units require neutral_is_absolute=True or explicit datum labels.")
    if strict and "neutral" in (raw_src_datum, raw_dst_datum) and raw_src_datum != raw_dst_datum and ("gauge" in (raw_src_datum, raw_dst_datum)):
        raise MissingReferenceConditionsError("strict=True requires explicit gauge and absolute pressure units rather than unqualified neutral pressure units.")
    src_datum = "absolute" if raw_src_datum == "neutral" else raw_src_datum
    dst_datum = "absolute" if raw_dst_datum == "neutral" else raw_dst_datum

    atm_input = p_atm if p_atm is not None else atmospheric_pressure
    if atm_input is None:
        atm_magnitude = STANDARD_ATMOSPHERE_PA
    elif isinstance(atm_input, pint.Quantity):
        if _classify_pressure_datum(atm_input.units) == "gauge":
            raise ValueError("p_atm must be an absolute barometric pressure, not gauge pressure.")
        atm_magnitude = atm_input.to(ureg.pascal).magnitude
    else:
        atm_magnitude = atm_input
    atm_decimal = Decimal(str(atm_magnitude))
    if not atm_decimal.is_finite() or atm_decimal <= 0:
        raise ValueError("Atmospheric pressure must be finite and strictly positive.")

    is_decimal = isinstance(quantity.magnitude, Decimal)
    src_pa = quantity.to(ureg.pascal)
    src_magnitude = Decimal(str(src_pa.magnitude)) if is_decimal else float(src_pa.magnitude)
    if not Decimal(str(src_magnitude)).is_finite():
        raise ValueError("Input pressure must be finite.")
    atm_mag = atm_decimal if is_decimal else float(atm_decimal)
    if not Decimal(str(atm_mag)).is_finite():
        raise ValueError("Atmospheric pressure exceeds the source magnitude's numeric range.")
    absolute_magnitude = src_magnitude + atm_mag if src_datum == "gauge" else src_magnitude
    if not Decimal(str(absolute_magnitude)).is_finite():
        raise ValueError("Absolute pressure exceeds the source magnitude's numeric range.")
    if not allow_negative_absolute and absolute_magnitude < 0:
        raise InvalidPhysicalStateError("Absolute pressure cannot be negative.")
    if src_datum == dst_datum:
        result = quantity.to(resolved_target)
    else:
        shifted = absolute_magnitude if dst_datum == "absolute" else absolute_magnitude - atm_mag
        result = ureg.Quantity(shifted, ureg.pascal).to(resolved_target)
    if is_decimal:
        out_q = ureg.Quantity(pad_ieee754_decimal(Decimal(str(result.magnitude))), result.units)
    else:
        out_q = ensure_quantity_ieee754(result)
    if not return_metadata:
        return out_q

    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    datum_shifted = src_datum != dst_datum
    cls_tag = ConversionClassification.STATE_AND_MODEL_DEPENDENT if datum_shifted else ConversionClassification.EXACT
    ref = ReferenceConditions(
        pressure=ureg.Quantity(atm_decimal, ureg.pascal),
        classification=cls_tag.value,
        parameters={"source_datum": src_datum, "target_datum": dst_datum, "atmospheric_pressure_pa": atm_decimal},
    )
    return build_conversion_result(
        out_q,
        quantity,
        conversion_path=(str(quantity.units), str(out_q.units)) if not datum_shifted else (str(quantity.units), "datum_shift:pascal", str(out_q.units)),
        reference_conditions=ref,
        method="convert_pressure",
        standard="ISO 2533" if atm_input is None and datum_shifted else "SI_pressure",
        estimated=(datum_shifted and atm_input is None),
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=lambda p: convert_pressure(
            ureg.Quantity(p["quantity"], quantity.units),
            resolved_target,
            p_atm=p["atmospheric_pressure_pa"],
            neutral_is_absolute=neutral_is_absolute,
            allow_negative_absolute=allow_negative_absolute,
        ),
        base_parameters={"quantity": _as_decimal(quantity.magnitude), "atmospheric_pressure_pa": atm_decimal},
    )





def convert_pressure_datum(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    atmospheric_pressure: pint.Quantity | float | Decimal | None = None,
    p_atm: pint.Quantity | float | Decimal | None = None,
    neutral_is_absolute: bool = True,
    allow_negative_absolute: bool = False,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Backward-compatible alias for convert_pressure(), including its options."""
    return convert_pressure(
        quantity, target_unit, p_atm=p_atm, atmospheric_pressure=atmospheric_pressure,
        neutral_is_absolute=neutral_is_absolute, allow_negative_absolute=allow_negative_absolute,
        return_metadata=return_metadata, strict=strict,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
    )





def differential_pressure(
    p1: pint.Quantity,
    p2: pint.Quantity,
    target_unit: str | pint.Unit = "psi",
    *,
    p_atm: pint.Quantity | float | Decimal | None = None,
) -> pint.Quantity:
    """Returns p1 - p2 after normalising both readings to absolute pressure.

    p_atm is shared by both readings, so it cancels for two gauge readings.
    A negative difference is valid. Target units set scale only, without a
    gauge/absolute shift. Either Decimal input gives a Decimal result.
    """
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(p1, pint.Quantity) or not isinstance(p2, pint.Quantity):
        raise TypeError("p1 and p2 must be Pint quantities.")
    validate_registry_quantity(p1)
    validate_registry_quantity(p2)
    if isinstance(p_atm, pint.Quantity):
        validate_registry_quantity(p_atm)
    if isinstance(p1.magnitude, Decimal) or isinstance(p2.magnitude, Decimal):
        p1 = ureg.Quantity(Decimal(str(p1.magnitude)), p1.units)
        p2 = ureg.Quantity(Decimal(str(p2.magnitude)), p2.units)
    p1_abs = convert_pressure(p1, "pascal", p_atm=p_atm)
    p2_abs = convert_pressure(p2, "pascal", p_atm=p_atm)
    resolved_target = parse_unit(target_unit)
    return (p1_abs - p2_abs).to(resolved_target)





# ---------------------------------------------------------------------------
# 6. Dynamic FX Triangulation & Cross-Currency Settlement
# ---------------------------------------------------------------------------


def _normalize_currency(code: str) -> str:
    """Resolves registered currency aliases to the canonical currency dimension."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    try:
        unit = parse_unit(str(code).strip())
    except (pint.UndefinedUnitError, ValueError) as error:
        raise ValueError(f"Unknown or unsupported currency code: {code!r}") from error
    dimensions = unit.dimensionality
    if len(dimensions) != 1:
        raise ValueError(f"Expected an isolated currency, got {code!r}")
    dimension, exponent = next(iter(dimensions.items()))
    if exponent != 1 or not dimension.startswith("[currency_"):
        raise ValueError(f"Expected an isolated currency, got {code!r}")
    canonical = dimension[len("[currency_"):-1]
    if unit != ureg.Unit(canonical):
        raise ValueError(f"Expected a major currency code, got minor unit {code!r}")
    return canonical





def _as_decimal(value: Any) -> Decimal:
    if isinstance(value, Decimal):
        return value
    return Decimal(str(value))





def create_fx_context(rates: dict[str, Decimal | float | int], *, base_currency: str = "USD", name: str = "spot_fx_matrix",) -> pint.Context:
    """Creates a multi-currency Pint transformation context from spot FX fixings.

    Triangulates all pairwise cross-rates across the supplied currencies using
    exact Decimal division prior to registering directional Pint transformations.

    Parameters
    ----------
    rates : dict[str, Decimal | float | int]
        Exchange rates expressed as units of foreign currency per 1 unit of base_currency
        (e.g. {"EUR": 0.92, "GBP": 0.79, "JPY": 154.50} when base_currency="USD").
    base_currency : str
        Anchor ISO-4217 currency code for the rate table (default: "USD").
    name : str
        Identifier for the generated Pint Context.
    """
    from core.units.registry import CURRENCY_CODES
    from core.units.registry import ureg
    base = _normalize_currency(base_currency)
    if base not in CURRENCY_CODES and base not in ureg:
        raise ValueError(f"Unknown base currency: {base_currency!r}")

    # Normalise all rates to Decimal relative to 1 unit of base_currency
    clean_rates: dict[str, Decimal] = {base: Decimal("1")}
    for raw_ccy, raw_rate in rates.items():
        ccy = _normalize_currency(raw_ccy)
        if ccy == base:
            continue
        if ccy not in CURRENCY_CODES and ccy not in ureg:
            raise ValueError(f"Unknown currency code in FX rates: {raw_ccy!r}")

        dec_rate = Decimal(str(raw_rate))
        if not dec_rate.is_finite() or dec_rate <= 0:
            raise ValueError(f"FX rate for {ccy} must be finite and strictly positive, got {raw_rate}")
        clean_rates[ccy] = dec_rate

    ctx = pint.Context(name)
    currencies = tuple(clean_rates.keys())

    def _make_fx_transform(src_code: str, dst_code: str, rate: Decimal):
        # Each call binds its own units and rates, independently of the loops.
        from core.units.registry import ureg
        source_unit = ureg.Unit(src_code)
        destination_unit = ureg.Unit(dst_code)
        float_rate = float(rate)

        def _transform(ureg_inst: pint.UnitRegistry, quantity: pint.Quantity) -> pint.Quantity:
            normalised = quantity.to(source_unit)
            factor = rate if isinstance(normalised.magnitude, Decimal) else float_rate
            return ureg_inst.Quantity(normalised.magnitude * factor, destination_unit)

        return _transform

    # Register all directed pairwise edges (src_ccy -> dst_ccy).
    for src_ccy in currencies:
        for dst_ccy in currencies:
            if src_ccy == dst_ccy:
                continue
            cross_rate_dec = clean_rates[dst_ccy] / clean_rates[src_ccy]
            ctx.add_transformation(
                f"[currency_{src_ccy}]",
                f"[currency_{dst_ccy}]",
                _make_fx_transform(src_ccy, dst_ccy, cross_rate_dec),
            )

    return ctx





def _numerator_currency_iso(quantity: pint.Quantity) -> str:
    """Requires one numerator currency dimension and physical denominators."""
    currencies = [(dim, exp) for dim, exp in quantity.dimensionality.items() if dim.startswith("[currency_")]
    if len(currencies) != 1 or currencies[0][1] != 1:
        raise ValueError("Quantity must have exactly one numerator currency dimension.")
    return currencies[0][0][len("[currency_"):-1]





def _replace_numerator_currency(
    quantity: pint.Quantity,
    src_iso: str,
    target_unit_obj: pint.Unit,
    cross_rate_dec: Decimal,
) -> pint.Quantity:
    """Swaps currency factors, scaling both minor units and preserving quotes.

    cross_rate_dec always relates major ISO currency units. Any physical unit
    factors remain intact; source and target minor-unit scales are explicit.
    """
    from core.units.registry import ureg
    src_iso = _normalize_currency(src_iso)
    expected_dim = {f"[currency_{src_iso}]": 1}
    matched: pint.Unit | None = None
    for unit_name, exponent in quantity.units._units.items():
        candidate = ureg.Unit(unit_name)
        if exponent == 1 and candidate.dimensionality == expected_dim:
            matched = candidate
            break
    if matched is None:
        raise ValueError("Quantity must contain a standalone numerator currency unit.")
    target_iso = _numerator_currency_iso(ureg.Quantity(1, target_unit_obj))
    source_scale = _as_decimal(ureg.Quantity(Decimal("1"), matched).to(src_iso).magnitude)
    target_scale = _as_decimal(ureg.Quantity(Decimal("1"), target_unit_obj).to(target_iso).magnitude)
    factor = cross_rate_dec * source_scale / target_scale
    magnitude = quantity.magnitude
    if not _as_decimal(magnitude).is_finite():
        raise ValueError(f"FX amount must be finite, got {magnitude}")
    converted = magnitude * (factor if isinstance(magnitude, Decimal) else float(factor))
    destination = (quantity.units / matched) * target_unit_obj
    return ureg.Quantity(converted, destination)





def _quantize_currency_settlement(quantity: pint.Quantity, currency_unit: pint.Unit) -> pint.Quantity:
    """Rounds a value in currency_unit to its major currency's minor quantum."""
    from core.units.registry import CURRENCY_MINOR_DIGITS
    from core.units.registry import ureg
    iso = _numerator_currency_iso(ureg.Quantity(1, currency_unit))
    if iso not in CURRENCY_MINOR_DIGITS:
        raise ValueError(f"Settlement precision for {iso} is unknown/undefined; supply an explicit settlement rule.")
    digits = CURRENCY_MINOR_DIGITS[iso]
    major_quantum = Decimal("1") if digits <= 0 else Decimal("1").scaleb(-digits)
    target_scale = _as_decimal(ureg.Quantity(Decimal("1"), currency_unit).to(iso).magnitude)
    quantum = (major_quantum / target_scale).normalize()
    raw = _as_decimal(quantity.magnitude)
    if not raw.is_finite():
        raise ValueError("Settlement amount must be finite.")
    rounded = (raw / quantum).quantize(Decimal("1"), rounding=ROUND_HALF_EVEN) * quantum
    return ureg.Quantity(rounded, quantity.units)





def convert_currency(
    quantity: pint.Quantity,
    target_currency: str,
    *,
    rates: dict[str, Decimal | float | int],
    base_currency: str = "USD",
    quantize_settlement: bool = False,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts a currency amount or cross-currency commodity quote to target_currency.

    Preserves physical denominators when converting commodity price quotes
    (e.g. converting "$75.50/bbl" directly into "EUR / oil_barrel" or "GBP_penny / therm_it"
    into "USD / therm_it").
    """
    from core.units.exceptions import MissingFXRateError
    from core.units.precision import ensure_quantity_ieee754
    from core.units.registry import CURRENCY_CODES
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)

    target_unit_obj = parse_unit(target_currency)
    target_dims = target_unit_obj.dimensionality
    if len(target_dims) != 1:
        raise ValueError(f"target_currency must be an isolated currency unit, got {target_currency!r}")

    target_dim_name, target_exp = next(iter(target_dims.items()))
    if target_exp != 1 or not target_dim_name.startswith("[currency_"):
        raise ValueError(f"target_currency must be a currency dimension, got {target_currency!r}")

    target_iso = target_dim_name[len("[currency_") : -1]

    src_iso = _numerator_currency_iso(quantity)

    # Normalise rates relative to 1 unit of base_currency and compute cross-rate
    base = _normalize_currency(base_currency)
    clean_rates: dict[str, Decimal] = {base: Decimal("1")}
    for raw_ccy, raw_rate in rates.items():
        ccy = _normalize_currency(raw_ccy)
        if ccy == base:
            continue
        if ccy not in CURRENCY_CODES and ccy not in ureg:
            raise ValueError(f"Unknown currency code in FX rates: {raw_ccy!r}")
        dec_rate = Decimal(str(raw_rate))
        if not dec_rate.is_finite() or dec_rate <= 0:
            raise ValueError(f"FX rate for {ccy} must be finite and strictly positive, got {raw_rate}")
        clean_rates[ccy] = dec_rate

    if src_iso not in clean_rates:
        raise MissingFXRateError(f"Missing FX rate for source currency {src_iso!r} relative to base {base!r}.")
    if target_iso not in clean_rates:
        raise MissingFXRateError(f"Missing FX rate for target currency {target_iso!r} relative to base {base!r}.")

    cross_rate_dec = clean_rates[target_iso] / clean_rates[src_iso]

    conversion_input = quantity
    if quantize_settlement and not isinstance(quantity.magnitude, Decimal):
        conversion_input = ureg.Quantity(_as_decimal(quantity.magnitude), quantity.units)
    converted = _replace_numerator_currency(conversion_input, src_iso, target_unit_obj, cross_rate_dec)

    if quantize_settlement:
        out_q = _quantize_currency_settlement(converted, target_unit_obj)
    else:
        out_q = ensure_quantity_ieee754(converted)

    if not return_metadata:
        return out_q

    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    same_ccy = src_iso == target_iso
    cls_tag = ConversionClassification.EXACT if same_ccy else ConversionClassification.MARKET_DEPENDENT
    ref = ReferenceConditions(
        classification=cls_tag.value,
        parameters={
            "source_currency": src_iso,
            "target_currency": target_iso,
            "base_currency": base,
            "fx_rate": cross_rate_dec,
            "quantize_settlement": quantize_settlement,
        },
    )
    return build_conversion_result(
        out_q,
        quantity,
        commodity_id=f"{src_iso}/{target_iso}",
        conversion_path=(str(quantity.units), str(out_q.units)),
        reference_conditions=ref,
        method="convert_currency",
        standard="ISO 4217",
        estimated=False,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=lambda p: _replace_numerator_currency(
            ureg.Quantity(p["quantity"], quantity.units), src_iso, target_unit_obj, p["fx_rate"]
        ),
        base_parameters={"quantity": _as_decimal(quantity.magnitude), "fx_rate": cross_rate_dec},
    )





def triangulate_fx_rate(rates: Mapping[tuple[str, str] | str, Any], base: str, quote: str) -> Decimal:
    """Triangulate an FX cross rate from a sparse rate map.
    
    Rates mean 1 BASE = rate QUOTE. Direct, inverse, and multi-hop paths are supported.
    """
    from core.units.exceptions import MissingFXRateError
    base_c = _normalize_currency(base)
    quote_c = _normalize_currency(quote)
    if base_c == quote_c:
        return Decimal("1")

    graph: dict[str, list[tuple[str, Decimal]]] = {}
    for key, raw_rate in rates.items():
        if isinstance(key, str):
            if key.count("/") != 1:
                raise ValueError(f"FX pair string must be in 'BASE/QUOTE' format, got {key!r}")
            b, q = key.split("/", 1)
        elif isinstance(key, tuple) and len(key) == 2:
            b, q = key
        else:
            raise ValueError(f"FX pair key must be a (BASE, QUOTE) tuple or 'BASE/QUOTE' string, got {key!r}")
        b = _normalize_currency(b)
        q = _normalize_currency(q)
        rate = _as_decimal(raw_rate)
        if not rate.is_finite() or rate <= 0:
            raise ValueError(f"FX rate for {b}/{q} must be finite and positive, got {rate}")
        graph.setdefault(b, []).append((q, rate))
        graph.setdefault(q, []).append((b, Decimal("1") / rate))

    from collections import deque

    queue = deque([(base_c, Decimal("1"))])
    visited = {base_c}
    while queue:
        cur, acc = queue.popleft()
        for nxt, edge_rate in graph.get(cur, []):
            new_rate = acc * edge_rate
            if nxt == quote_c:
                return new_rate
            if nxt not in visited:
                visited.add(nxt)
                queue.append((nxt, new_rate))

    raise MissingFXRateError(f"No FX triangulation path found from {base_c} to {quote_c}.")





def convert_fx(
    amount: pint.Quantity | float | Decimal,
    target_currency: str,
    rates: Mapping[tuple[str, str] | str, Any],
    *,
    source_currency: str | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts currency amounts or commodity quotes using sparse FX rates."""
    from core.units.precision import ensure_quantity_ieee754
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    tgt_unit = parse_unit(target_currency)
    target_dims = tgt_unit.dimensionality
    if len(target_dims) != 1:
        raise ValueError("target_currency must be an isolated currency unit.")
    tgt = _numerator_currency_iso(ureg.Quantity(1, tgt_unit))

    if isinstance(amount, pint.Quantity):
        validate_registry_quantity(amount)
        src_qty = amount
        src = _numerator_currency_iso(amount)
        if source_currency is not None and _normalize_currency(source_currency) != src:
            raise ValueError("source_currency does not match the quantity's currency dimension.")
        rate = triangulate_fx_rate(rates, src, tgt)
        out_q = ensure_quantity_ieee754(_replace_numerator_currency(amount, src, tgt_unit, rate))
    else:
        if source_currency is None:
            raise ValueError("source_currency is required when amount is not a pint.Quantity.")
        magnitude = amount if isinstance(amount, Decimal) else float(amount)
        if not _as_decimal(magnitude).is_finite():
            raise ValueError(f"FX amount must be finite, got {amount}")
        src = _normalize_currency(source_currency)
        src_qty = ureg.Quantity(magnitude, src)
        rate = triangulate_fx_rate(rates, src, tgt)
        out_q = ensure_quantity_ieee754(_replace_numerator_currency(src_qty, src, tgt_unit, rate))

    if not return_metadata:
        return out_q

    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    same_ccy = src == tgt
    cls_tag = ConversionClassification.EXACT if same_ccy else ConversionClassification.MARKET_DEPENDENT
    ref = ReferenceConditions(
        classification=cls_tag.value,
        parameters={"source_currency": src, "target_currency": tgt, "fx_rate": rate},
    )
    return build_conversion_result(
        out_q,
        src_qty,
        commodity_id=f"{src}/{tgt}",
        reference_conditions=ref,
        method="convert_fx",
        standard="ISO 4217",
        estimated=False,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=lambda p: _replace_numerator_currency(
            ureg.Quantity(p["quantity"], src_qty.units), src, tgt_unit, p["fx_rate"]
        ),
        base_parameters={"quantity": _as_decimal(src_qty.magnitude), "fx_rate": rate},
    )





# ---------------------------------------------------------------------------
# 7. Bulk Mining Moisture & Contained Metal Grade Transformations
# ---------------------------------------------------------------------------


def _validate_pct_in_range(
    value: float | Decimal, *, low_inclusive: float, high_exclusive: bool,
) -> Decimal:
    """Validates percentages without discarding Decimal precision."""
    percentage = _as_decimal(value)
    low = _as_decimal(low_inclusive)
    if not percentage.is_finite():
        raise ValueError(f"Percentage must be finite, got {value}")
    upper_ok = percentage < 100 if high_exclusive else percentage <= 100
    if percentage < low or not upper_ok:
        bracket = ")" if high_exclusive else "]"
        raise ValueError(f"Percentage must be in [{low}, 100{bracket}, got {value}")
    return percentage





def _ore_multiply(magnitude: Any, factor: Decimal) -> Any:
    if not _as_decimal(magnitude).is_finite():
        raise ValueError("Ore mass or price must be finite.")
    return magnitude * (factor if isinstance(magnitude, Decimal) else float(factor))





def _ore_mass_factor(unit: pint.Unit) -> tuple[pint.Unit, int, str | None]:
    """Finds one standalone mass factor, including in rates and price quotes."""
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import ureg
    factors = []
    for name, exponent in unit._units.items():
        candidate = ureg.Unit(name)
        if candidate.dimensionality in (ureg.kilogram.dimensionality, ureg.wmt.dimensionality, ureg.dmt.dimensionality, ureg.contained_metal_tonne.dimensionality, ureg.payable_metal_tonne.dimensionality):
            tokens = _unit_tokens(candidate)
            if tokens & {"wmt", "wet_metric_ton", "wet_metric_tonne", "wet_tonne"}:
                basis = "wet"
            elif tokens & {"dmt", "dry_metric_ton", "dry_metric_tonne", "dry_tonne"}:
                basis = "dry"
            elif tokens & {"payable_metal_tonne", "payable_mt"}:
                basis = "payable"
            elif tokens & {"contained_metal_tonne", "contained_mt", "dmtu", "dltu", "dry_metric_tonne_unit", "dry_long_ton_unit", "mtu", "stu"}:
                basis = "contained"
            else:
                basis = None
            factors.append((candidate, exponent, basis))
    if len(factors) != 1 or factors[0][1] not in (-1, 1):
        raise UnsupportedConversionError("Ore conversions require one standalone mass factor with exponent +1 or -1.")
    return factors[0]


from core.units.validation import authorized_conversion


@authorized_conversion
def convert_wet_to_dry_mass(wet_mass: pint.Quantity, moisture_pct: float | Decimal) -> pint.Quantity:
    """Converts WMT to DMT: dry mass = wet mass * (1 - moisture/100)."""
    from core.units.analytical_basis import wet_to_dry_mass
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(wet_mass, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
    validate_registry_quantity(wet_mass)
    moisture = _validate_pct_in_range(moisture_pct, low_inclusive=0, high_exclusive=True)
    mass = wet_mass.to(ureg.wmt)
    return wet_to_dry_mass(mass, moisture_pct=moisture, target_unit=ureg.dmt, allow_negative=True)


from core.units.validation import authorized_conversion


@authorized_conversion
def convert_dry_to_wet_mass(dry_mass: pint.Quantity, moisture_pct: float | Decimal) -> pint.Quantity:
    """Converts DMT to WMT: wet mass = dry mass / (1 - moisture/100)."""
    from core.units.analytical_basis import dry_to_wet_mass
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(dry_mass, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(dry_mass).__name__}")
    validate_registry_quantity(dry_mass)
    moisture = _validate_pct_in_range(moisture_pct, low_inclusive=0, high_exclusive=True)
    mass = dry_mass.to(ureg.dmt)
    return dry_to_wet_mass(mass, moisture_pct=moisture, target_unit=ureg.wmt, allow_negative=True)





def convert_ore_moisture_price(
    price_quantity: pint.Quantity,
    target_mass_unit: str | pint.Unit,
    moisture_pct: float | Decimal,
) -> pint.Quantity:
    """Changes the mass denominator of a currency quote, preserving its currency.

    Dry-to-wet prices multiply by the dry fraction; wet-to-dry prices divide.
    Unqualified mass units retain the other endpoint's basis without an
    implicit moisture shift. Source and destination must use ore mass, not
    contained-metal units; use convert_bulk_ore() for grade-based prices.
    """
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(price_quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(price_quantity).__name__}")
    validate_registry_quantity(price_quantity)
    moisture = _validate_pct_in_range(moisture_pct, low_inclusive=0, high_exclusive=True)
    _numerator_currency_iso(price_quantity)
    source_mass, exponent, source_basis = _ore_mass_factor(price_quantity.units)
    target = parse_unit(target_mass_unit)
    if exponent != -1 or target.dimensionality not in (ureg.kilogram.dimensionality, ureg.wmt.dimensionality, ureg.dmt.dimensionality):
        raise ValueError("Price must have a mass denominator and target_mass_unit must be a mass unit.")
    _, _, target_basis = _ore_mass_factor(target)
    if "contained" in (source_basis, target_basis):
        raise UnsupportedConversionError("Moisture prices require wet/dry ore mass units, not contained-metal units.")
    destination = price_quantity.units * source_mass / target
    return convert_bulk_ore(price_quantity, destination, moisture_fraction=moisture / 100)


from core.units.validation import authorized_conversion


@authorized_conversion
def extract_contained_metal(
    dry_ore_mass: pint.Quantity,
    grade_pct: float | Decimal,
    *,
    target_unit: str | pint.Unit = "contained_metal_tonne",
) -> pint.Quantity:
    """Returns contained metal from dry ore; 1 DMTU represents 10 kg metal."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(dry_ore_mass, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(dry_ore_mass).__name__}")
    validate_registry_quantity(dry_ore_mass)
    _, _, output_basis = _ore_mass_factor(parse_unit(target_unit))
    if output_basis not in (None, "contained"):
        raise ValueError("Contained-metal extraction requires a contained-metal output unit.")
    grade = _validate_pct_in_range(grade_pct, low_inclusive=0, high_exclusive=False) / 100
    dry = dry_ore_mass.to(ureg.dmt)
    metal = ureg.Quantity(_ore_multiply(dry.magnitude, grade), ureg.metric_ton)
    return metal.to(parse_unit(target_unit))


from core.units.validation import authorized_conversion


@authorized_conversion
def extract_dry_mass_from_contained(
    contained_metal: pint.Quantity,
    grade_pct: float | Decimal,
    *,
    target_unit: str | pint.Unit = "dmt",
) -> pint.Quantity:
    """Returns dry ore required to supply a contained metal mass at grade_pct."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(contained_metal, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(contained_metal).__name__}")
    validate_registry_quantity(contained_metal)
    grade = _validate_pct_in_range(grade_pct, low_inclusive=0, high_exclusive=False) / 100
    if grade == 0:
        raise ValueError("Cannot compute dry ore mass with grade_pct=0.")
    _, _, supplied_basis = _ore_mass_factor(contained_metal.units)
    if supplied_basis not in (None, "contained"):
        raise ValueError("Contained-metal input is required for inverse grade extraction.")
    metal = contained_metal.to(ureg.metric_ton)
    _ore_multiply(metal.magnitude, Decimal(1))
    scale = grade if isinstance(metal.magnitude, Decimal) else float(grade)
    _, _, output_basis = _ore_mass_factor(parse_unit(target_unit))
    if output_basis not in (None, "dry"):
        raise ValueError("Inverse grade extraction requires a dry-bulk output unit.")
    dry = ureg.Quantity(metal.magnitude / scale, ureg.dmt)
    return dry.to(parse_unit(target_unit))


from core.units.validation import authorized_conversion


@authorized_conversion
def convert_bulk_ore(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: BulkOreAssay | OreAssay | None = None,
    moisture_fraction: float | Decimal | None = None,
    grade: float | Decimal | None = None,
    component: str = "Fe",
) -> pint.Quantity:
    """Converts wet/dry/contained mass factors in amounts, rates, and quotes.

    Grade is a dry-basis fraction, not a percentage. Unqualified mass units
    inherit the other endpoint's basis, so wmt -> kg is simple unit scaling.
    Explicit wmt/dmt/contained labels request a basis change. Compound units
    must have one standalone mass factor with exponent +1 or -1 and otherwise
    compatible dimensions. Decimal magnitudes and assay fractions retain their
    precision. Explicit assay values take precedence over scalar parameters.
    """
    from core.units.assays.base import MaterialSpec
    from core.units.assays.metals.ores import BulkOreAssay
    from core.units.assays.metals.ores import OreAssay
    from core.units.composition import CompositionBasis
    from core.units.composition import MaterialComposition
    from core.units.composition import MaterialProperties
    from core.units.exceptions import IncompatibleCommodityError
    from core.units.exceptions import MissingAssayError
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)
    if isinstance(assay, OreAssay):
        if component != assay.primary_element:
            raise IncompatibleCommodityError(f"OreAssay defines {assay.primary_element!r}, not {component!r}.")
        resolved = BulkOreAssay(assay.name, _as_decimal(assay.moisture_pct) / 100, {component: _as_decimal(assay.grade_pct) / 100})
    elif isinstance(assay, (MaterialComposition, MaterialProperties, MaterialSpec)):
        comp = assay.composition if isinstance(assay, MaterialSpec) else (assay if isinstance(assay, MaterialComposition) else None)
        props = assay.properties if isinstance(assay, MaterialSpec) else (assay if isinstance(assay, MaterialProperties) else None)
        moist_frac = (
            comp.moisture_fraction
            if comp is not None and comp.moisture_fraction is not None
            else (props.moisture if props is not None and props.moisture is not None else moisture_fraction)
        )
        dry_grades: dict[str, Decimal] = {}
        if comp is not None and comp.basis in {
            CompositionBasis.AS_RECEIVED,
            CompositionBasis.AIR_DRIED,
            CompositionBasis.DRY,
            CompositionBasis.DRY_ASH_FREE,
        }:
            if comp.basis == CompositionBasis.DRY:
                dry_comp: MaterialComposition | None = comp
            elif comp.basis == CompositionBasis.DRY_ASH_FREE and comp.ash_fraction is not None:
                dry_comp = comp.to_basis(CompositionBasis.DRY)
            elif moist_frac is not None:
                dry_comp = comp.to_basis(CompositionBasis.DRY, target_moisture_fraction=moist_frac)
            else:
                dry_comp = None
            if dry_comp is not None:
                dry_grades.update(dry_comp.components)
                dry_grades.update(dry_comp.impurities)
        if grade is not None and component not in dry_grades:
            dry_grades[component] = _as_decimal(grade)
        resolved = BulkOreAssay(
            name=assay.name if isinstance(assay, MaterialSpec) else "material_building_block",
            moisture_fraction=moist_frac,
            grades=dry_grades,
        )
    elif assay is None:
        resolved = BulkOreAssay(
            moisture_fraction=moisture_fraction,
            grades={component: grade} if grade is not None else {},
        )
    elif isinstance(assay, BulkOreAssay):
        resolved = assay
    else:
        raise TypeError("assay must be a BulkOreAssay, OreAssay, MaterialComposition, MaterialProperties, or MaterialSpec.")

    target = parse_unit(target_unit)
    _, source_exp, source_basis = _ore_mass_factor(quantity.units)
    _, target_exp, target_basis = _ore_mass_factor(target)
    if source_exp != target_exp:
        raise UnsupportedConversionError("Source and target mass factors must have matching exponents.")
    scaled = ureg.Quantity(ureg._convert_mass_scale(quantity.magnitude, quantity.units._units, target._units), target)
    source_basis = source_basis or target_basis
    target_basis = target_basis or source_basis
    magnitude = scaled.magnitude
    _ore_multiply(magnitude, Decimal(1))
    if source_basis != target_basis:
        if "payable" in (source_basis, target_basis):
            if isinstance(assay, OreAssay):
                payable_fraction = _as_decimal(assay.payable_grade_pct) / 100
            elif (
                isinstance(assay, MaterialSpec)
                and assay.processing is not None
                and component in assay.processing.recovery_rates
                and component in resolved.grades
            ):
                payable_fraction = _as_decimal(resolved.grades[component]) * assay.processing.recovery_rates[component]
            else:
                raise MissingAssayError("Payable mass requires an OreAssay or MaterialSpec with explicit smelter/recovery terms.")
        else:
            payable_fraction = None
        dry_fraction = _as_decimal(resolved.dry_fraction) if "wet" in (source_basis, target_basis) else Decimal(1)
        fraction = None
        if "contained" in (source_basis, target_basis):
            if component not in resolved.grades:
                raise MissingAssayError(f"Missing dry-basis grade for {component!r}.")
            fraction = _as_decimal(resolved.grades[component])

        def apply_factor(value: Any, factor: Decimal, *, divide: bool) -> Any:
            if divide:
                if factor == 0:
                    raise ValueError("A positive grade is required for this contained-metal conversion.")
                return value / (factor if isinstance(value, Decimal) else float(factor))
            return _ore_multiply(value, factor)

        # Apply factors directly: dividing by a rounded reciprocal would lose
        # precision at exact Decimal price or settlement boundaries.
        if source_basis == "wet":
            magnitude = apply_factor(magnitude, dry_fraction, divide=source_exp == -1)
        elif source_basis in {"contained", "payable"}:
            magnitude = apply_factor(magnitude, payable_fraction if source_basis == "payable" else fraction, divide=source_exp == 1)
        if target_basis == "wet":
            magnitude = apply_factor(magnitude, dry_fraction, divide=source_exp == 1)
        elif target_basis in {"contained", "payable"}:
            magnitude = apply_factor(magnitude, payable_fraction if target_basis == "payable" else fraction, divide=source_exp == -1)
    return ureg.Quantity(magnitude, target)





# ---------------------------------------------------------------------------
# 8. Unified Cross-Commodity & Cross-Currency Quote Arbitrage Router
# ---------------------------------------------------------------------------


def _split_ratio_units(value: pint.Quantity | pint.Unit) -> tuple[pint.Unit, pint.Unit]:
    """Extracts the actual currency numerator and its physical denominator.

    Minor currencies retain their own unit. A plain currency amount has a
    dimensionless denominator. Compound currency factors are rejected clearly.
    """
    from core.units.registry import ureg
    unit = value.units if isinstance(value, pint.Quantity) else value
    iso = _numerator_currency_iso(ureg.Quantity(1, unit))
    for name, exponent in unit._units.items():
        candidate = ureg.Unit(name)
        if exponent == 1 and candidate.dimensionality == {f"[currency_{iso}]": 1}:
            return candidate, candidate / unit
    raise ValueError("Quote must contain a standalone numerator currency unit.")


from core.units.validation import authorized_conversion


@authorized_conversion
def _convert_contextual_quantity(
    quantity: pint.Quantity,
    target: pint.Unit,
    contexts: list[pint.Context],
) -> pint.Quantity:
    """Converts physical quantities and matching time-based rates in active contexts.

    Pint applies context edges to full dimensions, so a volume-to-mass edge
    needs the shared time factor removed for volume/day -> mass/hour. Retain
    original time units in the physical quantities to preserve their scaling.
    """
    from core.units.registry import ureg
    try:
        return quantity.to(target)
    except pint.DimensionalityError as original_error:
        semantic = {"wet", "dry", "contained", "payable"}
        source_tags = set(ureg._semantic_tags(quantity.units._units)) & semantic
        target_tags = set(ureg._semantic_tags(target._units)) & semantic
        if source_tags and target_tags and source_tags != target_tags:
            raise original_error
        if source_tags or target_tags:
            # Retain the explicit assay context while normalizing semantic roots
            # for its physical/time edges. This supplies scaling, not assay factors.
            source_factor, source_units = ureg._project_mass_units(quantity.units._units)
            target_factor, target_units = ureg._project_mass_units(target._units)
            physical = ureg.Quantity(quantity.magnitude * _scale_factor(quantity.magnitude, source_factor), ureg.Unit(source_units))
            result = _convert_contextual_quantity(physical, ureg.Unit(target_units), contexts)
            return ureg.Quantity(result.magnitude / _scale_factor(result.magnitude, target_factor), target)
        # Runtime currencies may not be present in a context's pre-registered
        # quote edges. Convert their physical denominator independently.
        if any(dim.startswith("[currency_") for dim in quantity.dimensionality):
            return _convert_price_quantity(quantity, target, contexts=contexts)
        dimensions = {
            ureg.get_dimensionality(dim)
            for context in contexts for edge in context.funcs for dim in edge
        }
        for source_dimension in dimensions:
            remainder = quantity.dimensionality / source_dimension
            if not remainder or set(remainder) != {"[time]"}:
                continue
            for target_dimension in dimensions:
                if target.dimensionality / target_dimension != remainder:
                    continue
                time_unit = ureg.second ** remainder["[time]"]
                physical = ureg.Quantity(quantity.magnitude, quantity.units / time_unit)
                try:
                    result = physical.to(target / time_unit)
                except pint.DimensionalityError:
                    continue
                return ureg.Quantity(result.magnitude, target)
        raise original_error





def _convert_price_quantity(
    quantity: pint.Quantity,
    target_unit: pint.Unit,
    *,
    ctx: pint.Context | None = None,
    contexts: list[pint.Context] | None = None,
) -> pint.Quantity:
    """Converts a quote by converting one source denominator to target units.

    Uses the active physical contexts (or an explicitly supplied context).
    Handles zero and negative quotes without taking the reciprocal of a price.
    """
    from core.units.registry import ureg
    source_num, source_den = _split_ratio_units(quantity)
    target_num, target_den = _split_ratio_units(target_unit)
    if source_num.dimensionality != target_num.dimensionality:
        raise ValueError("Convert quote currencies before the physical denominator.")
    normalised = quantity.to(target_num / source_den)
    reference = ureg.Quantity(Decimal("1") if isinstance(quantity.magnitude, Decimal) else 1.0, source_den)
    if ctx is not None:
        with ureg.context(ctx):
            physical = _convert_contextual_quantity(reference, target_den, [ctx])
    else:
        physical = _convert_contextual_quantity(reference, target_den, contexts or [])
    factor = physical.magnitude
    if isinstance(normalised.magnitude, Decimal):
        factor = _as_decimal(factor)
    result = normalised.magnitude / factor
    return ureg.Quantity(result, target_unit)





def _build_quote_contexts(
    *,
    assay: Any = None,
    crude_assay: CrudeAssay | None = None,
    gas_assay: GasAssay | None = None,
    gas_basis: str = "gcv",
    energy_basis: str = "lhv",
    active_only: bool | None = None,
    product: str | None = None,
) -> list[pint.Context]:
    """Collects physical contexts; later contexts have Pint's normal precedence.

    active_only=None uses each assay's default (True for chemicals, False for
    petrochemicals); an explicit False is honoured for both. product selects
    the appropriate density when an assay supplies more than one product.
    """
    from core.units.assays.base import MaterialSpec
    from core.units.assays.energy.biofuels import BiofuelAssay
    from core.units.assays.energy.biofuels import CleanFuelAssay
    from core.units.assays.energy.biomass import BiomassAssay
    from core.units.assays.energy.coal import CoalAssay
    from core.units.assays.energy.hydrogen import HydrogenAssay
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.assays.energy.ngl import NGLAssay
    from core.units.assays.energy.petcoke import PetcokeAssay
    from core.units.assays.energy.petroleum import CrudeAssay
    from core.units.assays.industrial.chemicals import ChemicalAssay
    from core.units.assays.industrial.petrochemicals import PetrochemAssay
    from core.units.composition import MaterialProperties
    from core.units.exceptions import MissingAssayError
    from core.units.heat_content import HeatContentAssay

    contexts = []
    if product is not None and assay is None:
        raise MissingAssayError("product selection requires an assay.")
    if assay is not None:
        if product is not None:
            method = getattr(assay, "as_context", None)
            if not callable(method):
                raise TypeError("product selection requires assay.as_context(product=...).")
            context = method(product=product)
        elif isinstance(assay, HeatContentAssay):
            context = assay.as_context(basis=energy_basis)
        elif isinstance(assay, CrudeAssay):
            context = create_hydrocarbon_context(assay)
        elif isinstance(assay, GasAssay):
            context = create_gas_heating_context(assay, basis=gas_basis)
        elif isinstance(assay, (CoalAssay, CleanFuelAssay, BiofuelAssay, BiomassAssay, PetcokeAssay, HydrogenAssay)):
            context = assay.as_context(basis=energy_basis)
        elif isinstance(assay, NGLAssay):
            context = assay.as_context(basis=gas_basis)
        elif isinstance(assay, PetrochemAssay):
            context = create_petrochem_context(assay, active_only=False if active_only is None else active_only)
        elif isinstance(assay, ChemicalAssay):
            context = assay.as_context(active_only=True if active_only is None else active_only)
        elif isinstance(assay, (MaterialProperties, MaterialSpec)):
            context = assay.as_context(active_only=False if active_only is None else active_only)
        elif callable(getattr(assay, "as_context", None)):
            context = assay.as_context()
        else:
            raise TypeError(f"Unsupported assay type for quote conversion: {type(assay).__name__}")
        if not isinstance(context, pint.Context):
            raise TypeError("assay.as_context() must return a Pint Context.")
        contexts.append(context)
    if crude_assay is not None:
        contexts.append(create_hydrocarbon_context(crude_assay))
    if gas_assay is not None:
        contexts.append(create_gas_heating_context(gas_assay, basis=gas_basis))
    return contexts





def convert_commodity_quote(
    quote_quantity: pint.Quantity,
    target_quote_unit: str | pint.Unit,
    *,
    rates: Mapping[str, Decimal | float | int] | None = None,
    base_currency: str = "USD",
    assay: Any = None,
    crude_assay: CrudeAssay | None = None,
    gas_assay: GasAssay | None = None,
    gas_basis: str = "gcv",
    energy_basis: str = "lhv",
    active_only: bool | None = None,
    product: str | None = None,
    quantize_settlement: bool = False,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts quote currency and physical denominator, then optionally settles."""
    from core.units.exceptions import MissingAssayError
    from core.units.exceptions import MissingFXRateError
    from core.units.precision import ensure_quantity_ieee754
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quote_quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quote_quantity).__name__}")
    validate_registry_quantity(quote_quantity)
    target = parse_unit(target_quote_unit)
    _, source_den = _split_ratio_units(quote_quantity)
    target_num, target_den = _split_ratio_units(target)
    source_iso = _numerator_currency_iso(quote_quantity)
    target_iso = _numerator_currency_iso(ureg.Quantity(1, target_num))
    if not _as_decimal(quote_quantity.magnitude).is_finite():
        raise ValueError("Quote amount must be finite.")
    working = quote_quantity
    if quantize_settlement and not isinstance(working.magnitude, Decimal):
        working = ureg.Quantity(_as_decimal(working.magnitude), working.units)
    if source_iso != target_iso:
        if rates is None:
            raise MissingFXRateError(f"FX rates are required to convert {source_iso} to {target_iso}.")
        working = convert_currency(working, target_num, rates=dict(rates), base_currency=base_currency)

    contexts = _build_quote_contexts(
        assay=assay, crude_assay=crude_assay, gas_assay=gas_assay,
        gas_basis=gas_basis, energy_basis=energy_basis, active_only=active_only, product=product,
    )
    cross_dimension = source_den.dimensionality != target_den.dimensionality
    if cross_dimension and not contexts:
        raise MissingAssayError("Changing quote denominator dimensionality requires an assay, crude_assay, or gas_assay.")
    if contexts:
        with ureg.context(*contexts):
            converted = _convert_price_quantity(working, target, contexts=contexts) if cross_dimension else working.to(target)
    else:
        converted = working.to(target)
    if quantize_settlement:
        out_q = _quantize_currency_settlement(converted, target_num)
    else:
        out_q = ensure_quantity_ieee754(converted)

    if not return_metadata:
        return out_q

    from core.units.provenance import ConversionClassification
    from core.units.provenance import build_conversion_result
    eff_assay = assay or crude_assay or gas_assay
    cls_tag = (
        ConversionClassification.ASSAY_DEPENDENT
        if cross_dimension and eff_assay is not None
        else (ConversionClassification.MARKET_DEPENDENT if source_iso != target_iso else ConversionClassification.EXACT)
    )
    path = (
        (str(quote_quantity.units), str(working.units), str(out_q.units))
        if source_iso != target_iso and str(working.units) != str(out_q.units)
        else (str(quote_quantity.units), str(out_q.units))
    )
    return build_conversion_result(
        out_q,
        quote_quantity,
        assay=eff_assay,
        conversion_path=path,
        method="convert_commodity_quote",
        classification=cls_tag,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
    )





# ---------------------------------------------------------------------------
# 9. Master Conversion Router
# ---------------------------------------------------------------------------


def _is_energy_unit(unit_like: str | pint.Unit) -> bool:
    """Identify energy independently of its market alias or quoting scale."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.dimensions import roles
    from core.units.dimensions import physical_unit
    unit = parse_unit(unit_like)
    signature = roles(unit)
    if len(signature) == 1:
        unit = physical_unit(unit, role=next(iter(signature)))
    return unit.dimensionality == ureg.joule.dimensionality





def _is_thermal_endpoint(unit: pint.Unit) -> bool:
    """Recognise energy amounts, rates, calorific densities and quote legs."""
    from core.units.registry import ureg
    from core.units.dimensions import roles
    from core.units.dimensions import physical_unit
    if any(dimension.startswith("[currency_") for dimension in unit.dimensionality):
        _, unit = _split_ratio_units(unit)
    signature = roles(unit)
    if len(signature) == 1:
        unit = physical_unit(unit, role=next(iter(signature)))
    dimensions = unit.dimensionality
    energy = ureg.joule.dimensionality
    if _is_energy_unit(unit) or dimensions == ureg.watt.dimensionality:
        return True
    denominator = energy / dimensions
    return denominator in (
        ureg.kilogram.dimensionality,
        (ureg.meter ** 3).dimensionality,
        ureg.mole.dimensionality,
    ) or (
        len(denominator) == 1
        and "gas_volume" in next(iter(denominator))
        and next(iter(denominator.values())) == 1
    )





def _convert_with_heat_content(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: Any = None,
    heat_content: HeatContentAssay | None = None,
    basis: str = "lhv",
    strict: bool = False,
    **heat_overrides: Any,
) -> pint.Quantity:
    """Resolve a composable thermal spec without replacing the commodity assay."""
    from core.units.heat_content import HeatContentAssay
    from core.units.heat_content import convert_heat_content
    from core.units.heat_content import heat_content_from_assay

    if heat_content is not None and not isinstance(heat_content, HeatContentAssay):
        raise TypeError("heat_content must be a HeatContentAssay.")
    analytical_basis = heat_overrides.pop("analytical_basis", None)
    # An explicit spec is authoritative, including its density, molar mass and
    # reference state. Commodity adapters are used only when no spec is given.
    spec = heat_content_from_assay(heat_content if heat_content is not None else assay, **heat_overrides)
    return convert_heat_content(
        quantity, target_unit, assay=spec, basis=basis, analytical_basis=analytical_basis, strict=strict,
    )





def _perturb_assay_for_uncertainty(assay: Any, params: Mapping[str, Decimal], base_params: Mapping[str, Decimal]) -> Any:
    """Construct a perturbed copy of ``assay`` for first-order sensitivity evaluation."""
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.assays.energy.petroleum import CrudeAssay
    from core.units.registry import ureg
    from dataclasses import fields as dc_fields, is_dataclass, replace

    if assay is None:
        return None
    changed = {k for k in params if k in base_params and params[k] != base_params[k]}
    if not changed:
        return assay

    if isinstance(assay, CrudeAssay):
        kw = {"name": assay.name, "product_group": assay.product_group, "sulfur_percentage": assay.sulfur_percentage, "tan": assay.tan}
        if "density" in changed or "density_kg_per_m3" in changed:
            key = "density" if "density" in changed else "density_kg_per_m3"
            return CrudeAssay(**kw, density_kg_per_m3=float(params[key]))
        if "api_gravity" in changed:
            return CrudeAssay(**kw, api_gravity=float(params["api_gravity"]))
        if "specific_gravity" in changed:
            return CrudeAssay(**kw, specific_gravity=float(params["specific_gravity"]))
        if "bbl_per_mt" in changed:
            return CrudeAssay(**kw, bbl_per_mt=float(params["bbl_per_mt"]))
        return assay

    if isinstance(assay, GasAssay):
        gcv = float(params["gcv_btu_per_scf"]) if "gcv_btu_per_scf" in changed else (float(params["btu_per_scf"]) if "btu_per_scf" in changed else assay.gcv_btu_per_scf)
        ncv = float(params["ncv_btu_per_scf"]) if "ncv_btu_per_scf" in changed else assay.ncv_btu_per_scf
        rel = float(params["relative_density"]) if "relative_density" in changed else assay.relative_density
        lng = float(params["mmbtu_per_mt_lng"]) if "mmbtu_per_mt_lng" in changed else assay.mmbtu_per_mt_lng
        return GasAssay(
            name=assay.name,
            gcv_btu_per_scf=gcv,
            ncv_btu_per_scf=ncv,
            relative_density=rel,
            mmbtu_per_mt_lng=lng,
        )

    if is_dataclass(assay):
        updates: dict[str, Any] = {}
        fmap = {f.name: getattr(assay, f.name) for f in dc_fields(assay)}
        alias_targets = {
            "density": ("density", "density_kg_per_m3", "bulk_density_kg_per_m3", "oil_density_kg_per_m3", "milk_density_kg_per_m3"),
            "moisture_fraction": ("moisture_fraction",),
            "moisture_pct": ("moisture_pct", "moisture_percentage", "actual_moisture_pct"),
            "lhv": ("lhv_mj_per_kg", "ncv", "lhv_btu_per_gal"),
            "ncv": ("ncv", "lhv_mj_per_kg"),
            "hhv": ("hhv_mj_per_kg", "gcv", "hhv_btu_per_gal"),
            "gcv": ("gcv", "hhv_mj_per_kg"),
            "grade": ("grade_pct",),
        }
        for key in changed:
            if key in fmap:
                cur = fmap[key]
                if isinstance(cur, pint.Quantity):
                    updates[key] = ureg.Quantity(params[key], cur.units)
                elif isinstance(cur, Decimal):
                    updates[key] = params[key]
                elif isinstance(cur, (int, float)):
                    updates[key] = float(params[key])
            elif key in alias_targets:
                for cand in alias_targets[key]:
                    if cand in fmap and fmap[cand] is not None:
                        cur = fmap[cand]
                        if isinstance(cur, pint.Quantity):
                            updates[cand] = ureg.Quantity(params[key], cur.units)
                        elif isinstance(cur, Decimal):
                            updates[cand] = params[key]
                        else:
                            updates[cand] = float(params[key])
                        break
            elif key == "moisture_fraction" and ("moisture_pct" in fmap or "moisture_percentage" in fmap):
                m_field = "moisture_pct" if "moisture_pct" in fmap else "moisture_percentage"
                updates[m_field] = float(params[key] * Decimal(100))
        if updates:
            return replace(assay, **updates)
    return assay





def _extract_assay_base_parameters(
    quantity: pint.Quantity,
    assay: Any,
    *,
    fx_cross_rate: Decimal | None = None,
    moisture_fraction: Decimal | float | None = None,
    grade: Decimal | float | None = None,
) -> dict[str, Decimal]:
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.assays.energy.petroleum import CrudeAssay
    from dataclasses import fields as dc_fields, is_dataclass

    base: dict[str, Decimal] = {"quantity": _as_decimal(quantity.magnitude)}
    if fx_cross_rate is not None:
        base["fx_rate"] = fx_cross_rate
    if moisture_fraction is not None:
        base["moisture_fraction"] = _as_decimal(moisture_fraction)
    if grade is not None:
        base["grade"] = _as_decimal(grade)

    if isinstance(assay, CrudeAssay):
        base["density"] = _as_decimal(assay.density_kg_per_m3)
        base["density_kg_per_m3"] = _as_decimal(assay.density_kg_per_m3)
        base["api_gravity"] = _as_decimal(assay.api_gravity)
        base["specific_gravity"] = _as_decimal(assay.specific_gravity)
        base["bbl_per_mt"] = _as_decimal(assay.bbl_per_mt)
    elif isinstance(assay, GasAssay):
        base["btu_per_scf"] = _as_decimal(assay.gcv_btu_per_scf)
        base["gcv_btu_per_scf"] = _as_decimal(assay.gcv_btu_per_scf)
        base["ncv_btu_per_scf"] = _as_decimal(assay.ncv_btu_per_scf)
        base["relative_density"] = _as_decimal(assay.relative_density)
        base["mmbtu_per_mt_lng"] = _as_decimal(assay.mmbtu_per_mt_lng)
    elif assay is not None and is_dataclass(assay):
        for f in dc_fields(assay):
            val = getattr(assay, f.name)
            if isinstance(val, (int, float, Decimal)) and not isinstance(val, bool):
                dec_v = _as_decimal(val)
                if dec_v.is_finite():
                    base[f.name] = dec_v
            elif isinstance(val, pint.Quantity):
                dec_v = _as_decimal(val.magnitude)
                if dec_v.is_finite():
                    base[f.name] = dec_v
        for alias, targets in (
            ("density", ("density", "density_kg_per_m3", "bulk_density_kg_per_m3", "oil_density_kg_per_m3", "milk_density_kg_per_m3")),
            ("lhv", ("lhv_mj_per_kg", "ncv", "lhv_btu_per_gal")),
            ("ncv", ("ncv", "lhv_mj_per_kg")),
            ("hhv", ("hhv_mj_per_kg", "gcv", "hhv_btu_per_gal")),
            ("gcv", ("gcv", "hhv_mj_per_kg")),
            ("moisture_pct", ("moisture_pct", "moisture_percentage", "actual_moisture_pct")),
        ):
            if alias not in base:
                for cand in targets:
                    if cand in base:
                        base[alias] = base[cand]
                        break
        if "moisture_fraction" not in base and "moisture_pct" in base:
            base["moisture_fraction"] = base["moisture_pct"] / Decimal(100)
    return base





def convert(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: Any = None,
    heat_content: HeatContentAssay | None = None,
    source_reference: GasReferenceConditions | None = None,
    energy_process: Any = None,
    target_reference: GasReferenceConditions | None = None,
    gas_conservation: str = "dry_gas",
    api_gravity: float | None = None,
    specific_gravity: float | None = None,
    bbl_per_mt: float | None = None,
    btu_per_scf: float | None = None,
    basis: str | None = None,
    p_atm: pint.Quantity | float | Decimal | None = None,
    atmospheric_pressure: pint.Quantity | float | Decimal | None = None,
    rates: Mapping[tuple[str, str] | str, Any] | None = None,
    base_currency: str = "USD",
    target_currency: str | None = None,
    source_currency: str | None = None,
    moisture_fraction: float | Decimal | None = None,
    grade: float | Decimal | None = None,
    component: str = "Fe",
    observed_temp: float | pint.Quantity | None = None,
    liquid_method: LiquidCorrectionMethod | None = None,
    liquid_observed_state: LiquidState | None = None,
    quantize_settlement: bool = False,
    active_only: bool | None = None,
    product: str | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
    **heat_overrides: Any,
) -> Any:
    """Dispatches pressure, FX, ore bases and physical assay conversions.

    When ``return_metadata=True``, returns an auditable ``ConversionResult``
    reconstructing the exact conversion path, assay identity/version, reference
    conditions (density, temperature, pressure, moisture, Z), registry version
    and SHA-256 definition hash, reliability classification, estimation flag,
    and optional propagated standard uncertainty (with covariance support).
    """
    from core.units.assays.base import MaterialSpec
    from core.units.assays.energy.biofuels import CleanFuelAssay
    from core.units.assays.energy.hydrogen import HydrogenAssay
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.assays.energy.ngl import NGLAssay
    from core.units.assays.energy.nuclear import UraniumAssay
    from core.units.assays.energy.petroleum import CrudeAssay
    from core.units.assays.metals.ores import BulkOreAssay
    from core.units.assays.metals.ores import OreAssay
    from core.units.composition import MaterialComposition
    from core.units.composition import MaterialProperties
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    from core.units.exceptions import InvalidPhysicalStateError
    from core.units.exceptions import MissingFXRateError
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.exceptions import UnsupportedConversionError
    from core.units.precision import ensure_quantity_ieee754
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    validate_registry_quantity(quantity)
    if isinstance(observed_temp, pint.Quantity):
        validate_registry_quantity(observed_temp)
    if isinstance(p_atm, pint.Quantity):
        validate_registry_quantity(p_atm)
    if isinstance(atmospheric_pressure, pint.Quantity):
        validate_registry_quantity(atmospheric_pressure)
    if (liquid_method is None) != (liquid_observed_state is None):
        raise MissingReferenceConditionsError("Supply both liquid_method and liquid_observed_state.")
    if liquid_method is not None and observed_temp is not None:
        raise UnsupportedConversionError("Use either legacy observed_temp or explicit liquid CTPL correction, not both.")
    if strict and observed_temp is not None:
        raise StandardsApplicabilityError(
            "strict=True prohibits legacy temperature-only observed_temp VCF estimates; supply liquid_method and liquid_observed_state."
        )
    magnitude = _as_decimal(quantity.magnitude)
    if not magnitude.is_finite():
        raise ValueError("Conversion input must be finite.")
    physical_amount = quantity.dimensionality in (
        ureg.kilogram.dimensionality, ureg.wmt.dimensionality, ureg.dmt.dimensionality,
        ureg.contained_metal_tonne.dimensionality, ureg.payable_metal_tonne.dimensionality,
        (ureg.meter ** 3).dimensionality, ureg.mole.dimensionality,
    ) or (
        len(quantity.dimensionality) == 1
        and "gas_volume" in next(iter(quantity.dimensionality))
        and next(iter(quantity.dimensionality.values())) == 1
    )
    if physical_amount and magnitude < 0:
        raise InvalidPhysicalStateError("Physical mass, volume and molar amounts must be non-negative.")
    target = parse_unit(target_unit)
    if source_reference is not None or target_reference is not None:
        if source_reference is None or target_reference is None:
            raise MissingReferenceConditionsError("Both source_reference and target_reference are required.")
        if any(value is not None for value in (assay, heat_content, api_gravity, specific_gravity, bbl_per_mt, btu_per_scf, basis, p_atm, atmospheric_pressure, rates, source_currency, target_currency, moisture_fraction, grade, observed_temp, liquid_method, liquid_observed_state, active_only, product, energy_process)) or quantize_settlement or heat_overrides:
            raise UnsupportedConversionError("Gas reference-state routing cannot be combined with assay, thermal, FX, datum or processing overrides.")
        return convert_gas_reference_volume(
            quantity, target, source=source_reference, target=target_reference,
            conserve=gas_conservation, return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )
    if gas_conservation != "dry_gas":
        raise MissingReferenceConditionsError("gas_conservation requires source_reference and target_reference.")
    pressure_dim = ureg.pascal.dimensionality
    if quantity.dimensionality == pressure_dim and target.dimensionality == pressure_dim:
        if liquid_method is not None:
            raise ValueError("Liquid CTPL correction requires a volume endpoint, not a pressure conversion.")
        return convert_pressure(
            quantity, target, p_atm=p_atm, atmospheric_pressure=atmospheric_pressure,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    resolved_assay = assay
    crude_parameters = any(value is not None for value in (api_gravity, specific_gravity, bbl_per_mt))
    if resolved_assay is None and crude_parameters:
        resolved_assay = CrudeAssay(api_gravity=api_gravity, specific_gravity=specific_gravity, bbl_per_mt=bbl_per_mt)
    if resolved_assay is None and btu_per_scf is not None:
        resolved_assay = GasAssay(name="custom_gas", gcv_btu_per_scf=btu_per_scf)

    # Detect currencies structurally. Non-currency physical amounts skip FX.
    def currency_iso(unit: pint.Unit) -> str | None:
        from core.units.registry import ureg
        if any(dim.startswith("[currency_") for dim in unit.dimensionality):
            return _numerator_currency_iso(ureg.Quantity(1, unit))
        return None

    source_iso = currency_iso(quantity.units)
    target_iso = currency_iso(target)
    target_num = None
    source_den = None
    target_den = None
    if source_iso is not None or target_iso is not None or target_currency is not None:
        if source_iso is None or target_iso is None:
            raise ValueError("Financial conversions require currency units on both endpoints.")
        _, source_den = _split_ratio_units(quantity)
        target_num, target_den = _split_ratio_units(target)
        if target_currency is not None:
            override = parse_unit(target_currency)
            if len(override.dimensionality) != 1:
                raise ValueError("target_currency must be an isolated currency unit.")
            _numerator_currency_iso(ureg.Quantity(1, override))
            target_num = override
            target = target_num / target_den
        if target_den.dimensionless and not source_den.dimensionless:
            target_den = source_den
            target = target_num / target_den
        target_iso = _numerator_currency_iso(ureg.Quantity(1, target_num))
        if source_currency is not None and _normalize_currency(source_currency) != source_iso:
            raise ValueError("source_currency does not match the quantity's currency dimension.")

    working = quantity
    fx_changed = source_iso is not None and source_iso != target_iso
    fx_cross_rate: Decimal | None = None
    if quantize_settlement and source_iso is not None and not isinstance(working.magnitude, Decimal):
        working = ureg.Quantity(_as_decimal(working.magnitude), working.units)
    if fx_changed:
        if rates is None:
            raise MissingFXRateError(f"FX rates are required to convert {source_iso} to {target_iso}.")
        is_fixing_map = all(isinstance(key, str) and "/" not in key for key in rates)
        if is_fixing_map:
            fx_audit = convert_currency(working, target_num, rates=dict(rates), base_currency=base_currency, return_metadata=True)
            working = fx_audit.quantity
            fx_cross_rate = fx_audit.reference_conditions.parameters.get("fx_rate") if fx_audit.reference_conditions else None
        else:
            fx_audit = convert_fx(working, target_num, rates, source_currency=source_currency, return_metadata=True)
            working = fx_audit.quantity
            fx_cross_rate = fx_audit.reference_conditions.parameters.get("fx_rate") if fx_audit.reference_conditions else None

    # Temperature correction changes volume basis before the standard assay leg.
    vcf = None
    audit_warnings: list[str] = []
    volume_dim = (ureg.meter ** 3).dimensionality
    def is_volume_endpoint(dimensionality: Any) -> bool:
        remainder = dimensionality / volume_dim
        return not remainder or set(remainder) == {"[time]"}

    source_volume = is_volume_endpoint(source_den.dimensionality if source_den is not None else quantity.dimensionality)
    target_volume = is_volume_endpoint(target_den.dimensionality if target_den is not None else target.dimensionality)
    if liquid_method is not None:
        from core.units.assays.energy.petroleum import LiquidCorrectionMethod
        from core.units.assays.energy.petroleum import LiquidProductSpec
        if not isinstance(liquid_method, LiquidCorrectionMethod):
            raise TypeError("liquid_method must be a LiquidCorrectionMethod.")
        if not source_volume and not target_volume:
            raise ValueError("Liquid CTPL correction requires at least one volume endpoint.")
        if isinstance(resolved_assay, CrudeAssay):
            vcf = resolved_assay.correction_factors(liquid_observed_state, method=liquid_method).ctpl
        elif isinstance(resolved_assay, LiquidProductSpec):
            vcf = liquid_method.factors(resolved_assay, liquid_observed_state).ctpl
        else:
            raise TypeError("Liquid correction requires a CrudeAssay or explicit LiquidProductSpec with a matching density basis.")
        if source_volume:
            factor = 1 / vcf if source_iso is not None else vcf
            working = ureg.Quantity(_ore_multiply(working.magnitude, factor), working.units)
    elif observed_temp is not None:
        method_vcf = getattr(resolved_assay, "volume_correction_factor", None)
        if not callable(method_vcf):
            raise MissingReferenceConditionsError("observed_temp requires an assay with a volume-correction method.")
        vcf = _as_decimal(method_vcf(observed_temp))
        if not vcf.is_finite() or vcf <= 0:
            raise InvalidPhysicalStateError("Volume correction factor must be finite and positive.")
        audit_warnings.append("Legacy temperature-only observed_temp VCF approximation applied.")
        if source_volume:
            factor = 1 / vcf if source_iso is not None else vcf
            working = ureg.Quantity(_ore_multiply(working.magnitude, factor), working.units)

    thermal_endpoint = _is_thermal_endpoint(working.units) or _is_thermal_endpoint(target)
    thermal_requested = heat_content is not None or bool(heat_overrides) or (
        resolved_assay is not None
        # Nuclear burnup is an existing assay conversion, not calorific heat.
        and not isinstance(resolved_assay, UraniumAssay)
        and working.dimensionality != target.dimensionality
    )
    if heat_overrides and not thermal_endpoint:
        raise UnsupportedConversionError("Heat-content overrides require an energy conversion.")

    def _is_bulk_basis_route(src_u: pint.Unit, dst_u: pint.Unit) -> bool:
        try:
            _, src_exp, src_basis = _ore_mass_factor(src_u)
            _, dst_exp, dst_basis = _ore_mass_factor(dst_u)
        except ValueError:
            return False
        return src_exp == dst_exp and bool(src_basis and dst_basis and src_basis != dst_basis)

    bulk_basis_changed = _is_bulk_basis_route(working.units, target)
    bulk_route = (
        isinstance(resolved_assay, (BulkOreAssay, OreAssay, MaterialComposition))
        or (isinstance(resolved_assay, (MaterialProperties, MaterialSpec)) and bulk_basis_changed)
        or moisture_fraction is not None
        or grade is not None
    )
    specification_converter = getattr(resolved_assay, "convert_quantity", None)
    route_method = "pint_unit_registry"
    route_heating_basis: str | None = basis
    physical_cross_dim = working.dimensionality != target.dimensionality

    if energy_process is not None:
        from core.units.processing import GenerationProcess
        from core.units.heat_content import _kind
        if not isinstance(energy_process, GenerationProcess):
            raise TypeError("energy_process must be GenerationProcess.")
        if resolved_assay is not None or heat_content is not None or heat_overrides or bulk_route or vcf is not None:
            raise UnsupportedConversionError("Convert the fuel to thermal energy first, then apply the generation process.")
        if basis is not None and _kind(basis) != _kind(energy_process.fuel_basis):
            raise IncompatibleAnalyticalBasisError("Generation process fuel HHV/LHV basis disagrees with basis.")
        converted = energy_process.convert(working, target)
        route_method = "GenerationProcess.convert"
        route_heating_basis = energy_process.fuel_basis
    elif callable(specification_converter) and heat_content is None and not heat_overrides and not bulk_route:
        if active_only is not None or product is not None:
            raise UnsupportedConversionError("Use explicit specification processing methods or the Chemical/Petrochem adapter for active/product bases.")
        converted = specification_converter(working, target, basis=basis)
        if not isinstance(converted, pint.Quantity):
            raise TypeError("assay.convert_quantity() must return a Pint quantity.")
        route_method = f"{type(resolved_assay).__name__}.convert_quantity"
    elif thermal_endpoint and thermal_requested:
        # Preserve existing fuel defaults; explicit shared thermal specs use LHV.
        thermal_basis = basis or (
            "gcv" if heat_content is None and isinstance(resolved_assay, GasAssay)
            else "hhv" if heat_content is None and isinstance(resolved_assay, NGLAssay)
            else "lhv"
        )
        route_heating_basis = thermal_basis
        converted = _convert_with_heat_content(
            working, target, assay=resolved_assay, heat_content=heat_content,
            basis=thermal_basis, strict=strict, **heat_overrides,
        )
        route_method = "convert_heat_content"
    elif bulk_route:
        ore_assay = (
            resolved_assay
            if isinstance(resolved_assay, (BulkOreAssay, OreAssay, MaterialComposition, MaterialProperties, MaterialSpec))
            else None
        )
        if isinstance(resolved_assay, OreAssay) and component == "Fe":
            ore_component = resolved_assay.primary_element
        elif isinstance(resolved_assay, (MaterialComposition, MaterialSpec)) and component == "Fe":
            comp_ref = resolved_assay.composition if isinstance(resolved_assay, MaterialSpec) else resolved_assay
            if comp_ref is not None and "Fe" not in comp_ref.components and "Fe" not in comp_ref.impurities and len(comp_ref.components) == 1:
                ore_component = next(iter(comp_ref.components))
            else:
                ore_component = component
        else:
            ore_component = component
        converted = convert_bulk_ore(working, target, assay=ore_assay, moisture_fraction=moisture_fraction, grade=grade, component=ore_component)
        route_method = "convert_bulk_ore" if bulk_basis_changed else "pint_unit_registry"
    elif source_iso is not None:
        converted = convert_commodity_quote(
            working, target, assay=resolved_assay,
            gas_basis=basis or "gcv", energy_basis=basis or "lhv",
            active_only=active_only, product=product,
        )
        route_method = "convert_commodity_quote" if (physical_cross_dim or fx_changed) else "pint_unit_registry"
    elif resolved_assay is not None:
        contexts = _build_quote_contexts(
            assay=resolved_assay, gas_basis=basis or "gcv",
            energy_basis=basis or "lhv", active_only=active_only, product=product,
        )
        with ureg.context(*contexts):
            converted = _convert_contextual_quantity(working, target, contexts)
        route_method = f"{type(resolved_assay).__name__}.as_context" if physical_cross_dim else "pint_unit_registry"
    else:
        converted = working.to(target)
        route_method = "pint_unit_registry"

    if vcf is not None and target_volume:
        factor = vcf if source_iso is not None else 1 / vcf
        converted = ureg.Quantity(_ore_multiply(converted.magnitude, factor), converted.units)

    if quantize_settlement and target_num is not None:
        final_qty = _quantize_currency_settlement(converted, target_num)
    else:
        final_qty = ensure_quantity_ieee754(converted)

    if not return_metadata:
        return final_qty

    from core.units.provenance import ConversionClassification
    from core.units.provenance import build_conversion_result

    # Classify reliability & build conversion path
    path_steps: list[str] = [str(quantity.units)]
    if fx_changed and str(working.units) != str(quantity.units) and str(working.units) != str(final_qty.units):
        path_steps.append(str(working.units))
    if vcf is not None:
        path_steps.append("ctpl_correction" if liquid_method is not None else "legacy_vcf_correction")
    if route_method not in {"pint_unit_registry", "convert_commodity_quote"} and physical_cross_dim:
        path_steps.append(route_method)
    if path_steps[-1] != str(final_qty.units):
        path_steps.append(str(final_qty.units))

    extra_params: dict[str, Any] = {}
    if fx_changed:
        extra_params["source_currency"] = source_iso
        extra_params["target_currency"] = target_iso
        extra_params["fx_rate"] = fx_cross_rate
        extra_params["base_currency"] = base_currency
    if vcf is not None:
        extra_params["volume_correction_factor"] = vcf
        if observed_temp is not None:
            extra_params["observed_temperature"] = observed_temp
        if liquid_observed_state is not None:
            extra_params["observed_temperature"] = liquid_observed_state.temperature
            extra_params["observed_pressure"] = liquid_observed_state.absolute_pressure
    if moisture_fraction is not None:
        extra_params["moisture_fraction_override"] = _as_decimal(moisture_fraction)
    if grade is not None:
        extra_params["grade_override"] = _as_decimal(grade)

    audit_assay = heat_content if (thermal_endpoint and thermal_requested and heat_content is not None and resolved_assay is None) else resolved_assay
    commodity_override: str | None = None
    std_override: str | None = None
    estimated_override: bool | None = None

    if energy_process is not None:
        cls_tag = ConversionClassification.STATE_AND_MODEL_DEPENDENT
        commodity_override = "electricity_generation"
        extra_params["thermal_efficiency"] = energy_process.thermal_efficiency
        extra_params["output_basis"] = energy_process.output_basis
        estimated_override = False
    elif not physical_cross_dim and not bulk_basis_changed and vcf is None:
        if fx_changed:
            cls_tag = ConversionClassification.MARKET_DEPENDENT
            commodity_override = f"{source_iso}/{target_iso}"
            std_override = "ISO 4217"
            estimated_override = False
            audit_assay = None
        else:
            cls_tag = ConversionClassification.EXACT
            std_override = "SI / NIST SP 811"
            estimated_override = False
            audit_assay = None
    elif thermal_endpoint and thermal_requested:
        from core.units.heat_content import HeatContentAssay as _HCA
        from core.units.heat_content import heat_content_from_assay as _hc_from
        eff_hc = _hc_from(heat_content if heat_content is not None else resolved_assay, **heat_overrides)
        extra_params["calorific_method"] = str(eff_hc.method.value if hasattr(eff_hc.method, "value") else eff_hc.method)
        if isinstance(resolved_assay, GasAssay):
            cls_tag = ConversionClassification.COMPOSITION_DEPENDENT
            if quantity.dimensionality == ureg.kilogram.dimensionality or target.dimensionality == ureg.kilogram.dimensionality:
                commodity_override = "lng"
        elif isinstance(resolved_assay, (NGLAssay, CleanFuelAssay, HydrogenAssay)) or eff_hc.basis in {"DRY_GAS", "WET_GAS"}:
            cls_tag = ConversionClassification.COMPOSITION_DEPENDENT
        else:
            cls_tag = ConversionClassification.MOISTURE_AND_CALORIFIC_DEPENDENT
        if observed_temp is not None:
            estimated_override = True
    elif bulk_basis_changed:
        _, _, src_b = _ore_mass_factor(working.units)
        _, _, dst_b = _ore_mass_factor(target)
        cls_tag = (
            ConversionClassification.YIELD_DEPENDENT
            if "payable" in (src_b, dst_b)
            else ConversionClassification.ASSAY_DEPENDENT
        )
    elif type(resolved_assay).__name__ == "EthanolFeedstockAssay":
        cls_tag = ConversionClassification.YIELD_DEPENDENT
    elif isinstance(resolved_assay, GasAssay):
        cls_tag = ConversionClassification.COMPOSITION_DEPENDENT
    else:
        cls_tag = ConversionClassification.ASSAY_DEPENDENT
        if observed_temp is not None:
            estimated_override = True

    base_params = _extract_assay_base_parameters(
        quantity,
        audit_assay,
        fx_cross_rate=fx_cross_rate,
        moisture_fraction=moisture_fraction,
        grade=grade,
    )
    heat_ovr_copy = dict(heat_overrides)

    def _eval_convert(p: Mapping[str, Decimal]) -> Decimal:
        from core.units.registry import ureg
        q_in = ureg.Quantity(p["quantity"], quantity.units)
        p_assay = _perturb_assay_for_uncertainty(resolved_assay, p, base_params)
        p_hc = _perturb_assay_for_uncertainty(heat_content, p, base_params) if heat_content is not None else None
        p_rates = rates
        if fx_changed and fx_cross_rate is not None and "fx_rate" in p and p["fx_rate"] != fx_cross_rate:
            p_rates = {f"{source_iso}/{target_iso}": p["fx_rate"]}
        p_moist = p["moisture_fraction"] if "moisture_fraction" in p and moisture_fraction is not None else moisture_fraction
        p_grade = p["grade"] if "grade" in p and grade is not None else grade
        res_q = convert(
            q_in,
            target_unit,
            assay=p_assay,
            heat_content=p_hc,
            energy_process=energy_process,
            basis=basis,
            rates=p_rates,
            base_currency=base_currency if p_rates is rates else source_iso,
            target_currency=target_currency,
            source_currency=source_currency,
            moisture_fraction=p_moist,
            grade=p_grade,
            component=component,
            observed_temp=observed_temp,
            liquid_method=liquid_method,
            liquid_observed_state=liquid_observed_state,
            quantize_settlement=False,
            active_only=active_only,
            product=product,
            **heat_ovr_copy,
        )
        return _as_decimal(res_q.magnitude)

    return build_conversion_result(
        final_qty,
        quantity,
        assay=audit_assay,
        commodity_id=commodity_override,
        conversion_path=tuple(path_steps),
        method=route_method,
        standard=std_override,
        classification=cls_tag,
        heating_basis=route_heating_basis,
        extra_parameters=extra_params,
        estimated=estimated_override,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=_eval_convert,
        base_parameters=base_params,
        warnings=tuple(audit_warnings),
    )




_converters_viscosity_DATA_FILE: Final[Path] = (Path(__file__).resolve().parents[0] / "reference_data") / "units" / "astm_d2161.json"



_converters_viscosity_ABSOLUTE_ZERO_F: Final[float] = -459.67


from core.units.exceptions import StandardsApplicabilityError


class DomainError(StandardsApplicabilityError):
    """Raised when a physical reading lies outside the valid empirical domain of ASTM D2161."""





def _converters_viscosity_load_astm_d2161(path: Path = _converters_viscosity_DATA_FILE) -> tuple[
    MonotonicCubicTable, float, float, float, float,
    MonotonicCubicTable, float, float, float,
]:
    from core.units.thermodynamics import MonotonicCubicTable
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported astm_d2161.json schema_version: {payload.get('schema_version')!r}")

    sus_cfg = payload["sus_100f"]
    sus_table = MonotonicCubicTable(sus_cfg["table_cst_sus"])
    sus_thresh = float(sus_cfg["linear_threshold_cst"])
    sus_mult = float(sus_cfg["linear_multiplier"])
    sus_alpha = float(sus_cfg["temp_expansion_per_f"])
    sus_ref_t = float(sus_cfg["reference_temp_f"])

    sfs_cfg = payload["sfs_122f"]
    sfs_table = MonotonicCubicTable(sfs_cfg["table_cst_sfs"])
    sfs_thresh = float(sfs_cfg["linear_threshold_cst"])
    sfs_mult = float(sfs_cfg["linear_multiplier"])
    sfs_ref_t = float(sfs_cfg["reference_temp_f"])

    return (
        sus_table, sus_thresh, sus_mult, sus_alpha, sus_ref_t,
        sfs_table, sfs_thresh, sfs_mult, sfs_ref_t,
    )





(
    _converters_viscosity_SUS_TABLE_100F,
    _converters_viscosity_SUS_LINEAR_THRESHOLD_CST,
    ASTM_D2161_SUS_LINEAR_FACTOR,
    ASTM_D2161_SUS_TEMP_COEFF,
    _converters_viscosity_SUS_REF_TEMP_F,
    _converters_viscosity_SFS_TABLE_122F,
    _converters_viscosity_SFS_LINEAR_THRESHOLD_CST,
    ASTM_D2161_SFS_LINEAR_FACTOR,
    _converters_viscosity_SFS_REF_TEMP_F,
) = _converters_viscosity_load_astm_d2161()


from core.units.registry import ureg

_converters_viscosity_DYNAMIC_DIM: Final = ureg.centipoise.dimensionality


from core.units.registry import ureg
_converters_viscosity_KINEMATIC_DIM: Final = ureg.centistokes.dimensionality


from core.units.registry import ureg
_converters_viscosity_DENSITY_DIM: Final = (ureg.kilogram / (ureg.meter ** 3)).dimensionality


from core.units.registry import ureg
_converters_viscosity_SUS_DIM: Final = ureg.SUS.dimensionality


from core.units.registry import ureg
_converters_viscosity_SFS_DIM: Final = ureg.SFS.dimensionality





def _converters_viscosity_resolve_temp_f(temp_f: pint.Quantity | float | int | Decimal) -> float:
    """Resolve an affine temperature (float °F or Pint Quantity) to °F and check > 0 R."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if isinstance(temp_f, pint.Quantity):
        validate_registry_quantity(temp_f)
        _converters__interpolation_validate_real_scalar(temp_f.magnitude, "temperature")
        if temp_f.units in {ureg.delta_degC, ureg.delta_degF}:
            raise ValueError("Viscosity temperature must be an absolute/state temperature, not a delta interval.")
        t_val = float(temp_f.to(ureg.degF).magnitude)
    else:
        t_val = float(_converters__interpolation_validate_real_scalar(temp_f, "temp_f"))

    if t_val <= _converters_viscosity_ABSOLUTE_ZERO_F:
        raise DomainError(f"Temperature {t_val} °F is at or below absolute zero ({_converters_viscosity_ABSOLUTE_ZERO_F} °F).")
    return t_val





def _converters_viscosity_sus_temp_factor(t_f: float) -> float:
    factor = 1.0 + ASTM_D2161_SUS_TEMP_COEFF * (t_f - _converters_viscosity_SUS_REF_TEMP_F)
    if factor <= 0.0:
        raise DomainError(f"ASTM D2161 thermal expansion factor is non-positive at {t_f} °F.")
    return factor


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# Stage 1: Dynamic (cP / Pa*s) <-> Kinematic (cSt / mm^2/s) via Density
# ---------------------------------------------------------------------------


@authorized_conversion
def dynamic_to_kinematic_viscosity(
    dynamic_viscosity: pint.Quantity,
    density: pint.Quantity,
    *,
    target_unit: str | pint.Unit = "centistokes",
) -> pint.Quantity:
    """Convert dynamic viscosity (e.g. cP, Pa*s) to kinematic viscosity (nu = mu / rho)."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if not isinstance(dynamic_viscosity, pint.Quantity) or not isinstance(density, pint.Quantity):
        raise TypeError("dynamic_viscosity and density must be pint.Quantity instances.")
    validate_registry_quantity(dynamic_viscosity)
    validate_registry_quantity(density)
    if dynamic_viscosity.dimensionality != _converters_viscosity_DYNAMIC_DIM:
        raise pint.DimensionalityError(
            dynamic_viscosity.units, ureg.centipoise, dynamic_viscosity.dimensionality, _converters_viscosity_DYNAMIC_DIM
        )
    if density.dimensionality != _converters_viscosity_DENSITY_DIM:
        raise pint.DimensionalityError(
            density.units, ureg.kilogram / (ureg.meter ** 3), density.dimensionality, _converters_viscosity_DENSITY_DIM
        )

    mu_dec = _converters__interpolation_validate_real_scalar(dynamic_viscosity.magnitude, "dynamic_viscosity")
    rho_dec = _converters__interpolation_validate_real_scalar(density.magnitude, "density")
    if mu_dec <= 0:
        raise DomainError(f"Dynamic viscosity must be strictly positive, got {dynamic_viscosity}.")
    if rho_dec <= 0:
        raise ValueError(f"Fluid density must be strictly positive, got {density}.")

    dst = parse_unit(target_unit)
    if dst.dimensionality != _converters_viscosity_KINEMATIC_DIM:
        raise pint.DimensionalityError(ureg.centistokes, dst, _converters_viscosity_KINEMATIC_DIM, dst.dimensionality)

    mu_q = ureg.Quantity(_as_decimal(dynamic_viscosity.magnitude), dynamic_viscosity.units)
    rho_q = ureg.Quantity(_as_decimal(density.magnitude), density.units)
    nu_q = (mu_q / rho_q).to(dst)
    out_mag = _converters__interpolation_match_scalar_type(nu_q.magnitude, dynamic_viscosity, density)
    return ureg.Quantity(out_mag, dst)


from core.units.validation import authorized_conversion


@authorized_conversion
def kinematic_to_dynamic_viscosity(
    kinematic_viscosity: pint.Quantity,
    density: pint.Quantity,
    *,
    target_unit: str | pint.Unit = "centipoise",
) -> pint.Quantity:
    """Convert kinematic viscosity (e.g. cSt, mm^2/s) to dynamic viscosity (mu = nu * rho)."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if not isinstance(kinematic_viscosity, pint.Quantity) or not isinstance(density, pint.Quantity):
        raise TypeError("kinematic_viscosity and density must be pint.Quantity instances.")
    validate_registry_quantity(kinematic_viscosity)
    validate_registry_quantity(density)
    if kinematic_viscosity.dimensionality != _converters_viscosity_KINEMATIC_DIM:
        raise pint.DimensionalityError(
            kinematic_viscosity.units, ureg.centistokes, kinematic_viscosity.dimensionality, _converters_viscosity_KINEMATIC_DIM
        )
    if density.dimensionality != _converters_viscosity_DENSITY_DIM:
        raise pint.DimensionalityError(
            density.units, ureg.kilogram / (ureg.meter ** 3), density.dimensionality, _converters_viscosity_DENSITY_DIM
        )

    nu_dec = _converters__interpolation_validate_real_scalar(kinematic_viscosity.magnitude, "kinematic_viscosity")
    rho_dec = _converters__interpolation_validate_real_scalar(density.magnitude, "density")
    if nu_dec <= 0:
        raise DomainError(f"Kinematic viscosity must be strictly positive, got {kinematic_viscosity}.")
    if rho_dec <= 0:
        raise ValueError(f"Fluid density must be strictly positive, got {density}.")

    dst = parse_unit(target_unit)
    if dst.dimensionality != _converters_viscosity_DYNAMIC_DIM:
        raise pint.DimensionalityError(ureg.centipoise, dst, _converters_viscosity_DYNAMIC_DIM, dst.dimensionality)

    nu_q = ureg.Quantity(_as_decimal(kinematic_viscosity.magnitude), kinematic_viscosity.units)
    rho_q = ureg.Quantity(_as_decimal(density.magnitude), density.units)
    mu_q = (nu_q * rho_q).to(dst)
    out_mag = _converters__interpolation_match_scalar_type(mu_q.magnitude, kinematic_viscosity, density)
    return ureg.Quantity(out_mag, dst)


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# Stage 2: Kinematic (cSt) <-> Saybolt Universal (SUS) & Saybolt Furol (SFS)
# ---------------------------------------------------------------------------


@authorized_conversion
def cst_to_sus(
    viscosity_cst: pint.Quantity | float | int | Decimal,
    *,
    temp_f: pint.Quantity | float | int | Decimal = 100.0,
) -> pint.Quantity | float | Decimal:
    """Convert kinematic viscosity (cSt) to Saybolt Universal Seconds (SUS) per ASTM D2161."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    t_f = _converters_viscosity_resolve_temp_f(temp_f)
    is_qty = isinstance(viscosity_cst, pint.Quantity)
    if is_qty:
        validate_registry_quantity(viscosity_cst)
        if viscosity_cst.dimensionality != _converters_viscosity_KINEMATIC_DIM:
            raise pint.DimensionalityError(
                viscosity_cst.units, ureg.centistokes, viscosity_cst.dimensionality, _converters_viscosity_KINEMATIC_DIM
            )
        _converters__interpolation_validate_real_scalar(viscosity_cst.magnitude, "viscosity_cst")
        nu = float(viscosity_cst.to(ureg.centistokes).magnitude)
    else:
        nu = float(_converters__interpolation_validate_real_scalar(viscosity_cst, "viscosity_cst"))

    if nu < _converters_viscosity_SUS_TABLE_100F.min_x:
        raise DomainError(
            f"Kinematic viscosity {nu} cSt is below the ASTM D2161 SUS minimum domain "
            f"({_converters_viscosity_SUS_TABLE_100F.min_x} cSt / ~32.6 SUS at 100 °F)."
        )

    if nu < _converters_viscosity_SUS_LINEAR_THRESHOLD_CST:
        sus_100 = _converters_viscosity_SUS_TABLE_100F.evaluate(nu)
    else:
        sus_100 = ASTM_D2161_SUS_LINEAR_FACTOR * nu

    sus_t = sus_100 * _converters_viscosity_sus_temp_factor(t_f)
    out = _converters__interpolation_match_scalar_type(sus_t, viscosity_cst)
    return ureg.Quantity(out, ureg.SUS) if is_qty else out


from core.units.validation import authorized_conversion


@authorized_conversion
def sus_to_cst(
    sus: pint.Quantity | float | int | Decimal,
    *,
    temp_f: pint.Quantity | float | int | Decimal = 100.0,
    target_unit: str | pint.Unit = "centistokes",
) -> pint.Quantity | float | Decimal:
    """Invert Saybolt Universal Seconds (SUS) at ``temp_f`` back to kinematic viscosity (cSt)."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    t_f = _converters_viscosity_resolve_temp_f(temp_f)
    is_qty = isinstance(sus, pint.Quantity)
    if is_qty:
        validate_registry_quantity(sus)
        if sus.dimensionality != _converters_viscosity_SUS_DIM:
            raise pint.DimensionalityError(sus.units, ureg.SUS, sus.dimensionality, _converters_viscosity_SUS_DIM)
        _converters__interpolation_validate_real_scalar(sus.magnitude, "sus")
        sus_val = float(sus.to(ureg.SUS).magnitude)
    else:
        sus_val = float(_converters__interpolation_validate_real_scalar(sus, "sus"))

    sus_100 = sus_val / _converters_viscosity_sus_temp_factor(t_f)
    if sus_100 < _converters_viscosity_SUS_TABLE_100F.min_y - 1e-12:
        raise DomainError(
            f"Temperature-corrected SUS ({sus_100:.4f} s at 100 °F) is below the ASTM D2161 "
            f"minimum domain ({_converters_viscosity_SUS_TABLE_100F.min_y} s / {_converters_viscosity_SUS_TABLE_100F.min_x} cSt)."
        )
    sus_100 = max(sus_100, _converters_viscosity_SUS_TABLE_100F.min_y)

    if sus_100 < _converters_viscosity_SUS_TABLE_100F.max_y:
        nu_cst = _converters_viscosity_SUS_TABLE_100F.invert(sus_100)
    else:
        nu_cst = sus_100 / ASTM_D2161_SUS_LINEAR_FACTOR

    if is_qty:
        dst = parse_unit(target_unit)
        if dst.dimensionality != _converters_viscosity_KINEMATIC_DIM:
            raise pint.DimensionalityError(ureg.centistokes, dst, _converters_viscosity_KINEMATIC_DIM, dst.dimensionality)
        converted = ureg.Quantity(nu_cst, ureg.centistokes).to(dst)
        return ureg.Quantity(_converters__interpolation_match_scalar_type(converted.magnitude, sus), dst)
    return _converters__interpolation_match_scalar_type(nu_cst, sus)


from core.units.validation import authorized_conversion


@authorized_conversion
def cst_to_sfs(
    viscosity_cst: pint.Quantity | float | int | Decimal,
    *,
    temp_f: pint.Quantity | float | int | Decimal = 122.0,
) -> pint.Quantity | float | Decimal:
    """Convert kinematic viscosity (cSt) to Saybolt Furol Seconds (SFS) at 122 °F per ASTM D2161."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    t_f = _converters_viscosity_resolve_temp_f(temp_f)
    if abs(t_f - _converters_viscosity_SFS_REF_TEMP_F) > 1e-6:
        raise DomainError(
            f"ASTM D2161 Table 3 Saybolt Furol conversion requires 122 °F (50 °C), got {t_f} °F."
        )

    is_qty = isinstance(viscosity_cst, pint.Quantity)
    if is_qty:
        validate_registry_quantity(viscosity_cst)
        if viscosity_cst.dimensionality != _converters_viscosity_KINEMATIC_DIM:
            raise pint.DimensionalityError(
                viscosity_cst.units, ureg.centistokes, viscosity_cst.dimensionality, _converters_viscosity_KINEMATIC_DIM
            )
        _converters__interpolation_validate_real_scalar(viscosity_cst.magnitude, "viscosity_cst")
        nu = float(viscosity_cst.to(ureg.centistokes).magnitude)
    else:
        nu = float(_converters__interpolation_validate_real_scalar(viscosity_cst, "viscosity_cst"))

    if nu < _converters_viscosity_SFS_TABLE_122F.min_x:
        raise DomainError(
            f"Kinematic viscosity {nu} cSt is below the ASTM D2161 SFS minimum domain "
            f"({_converters_viscosity_SFS_TABLE_122F.min_x} cSt / ~25.1 SFS at 122 °F)."
        )

    if nu < _converters_viscosity_SFS_LINEAR_THRESHOLD_CST:
        sfs_val = _converters_viscosity_SFS_TABLE_122F.evaluate(nu)
    else:
        sfs_val = ASTM_D2161_SFS_LINEAR_FACTOR * nu

    out = _converters__interpolation_match_scalar_type(sfs_val, viscosity_cst)
    return ureg.Quantity(out, ureg.SFS) if is_qty else out


from core.units.validation import authorized_conversion


@authorized_conversion
def sfs_to_cst(
    sfs: pint.Quantity | float | int | Decimal,
    *,
    temp_f: pint.Quantity | float | int | Decimal = 122.0,
    target_unit: str | pint.Unit = "centistokes",
) -> pint.Quantity | float | Decimal:
    """Invert Saybolt Furol Seconds (SFS) at 122 °F back to kinematic viscosity (cSt)."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    t_f = _converters_viscosity_resolve_temp_f(temp_f)
    if abs(t_f - _converters_viscosity_SFS_REF_TEMP_F) > 1e-6:
        raise DomainError(
            f"ASTM D2161 Table 3 Saybolt Furol inversion requires 122 °F (50 °C), got {t_f} °F."
        )

    is_qty = isinstance(sfs, pint.Quantity)
    if is_qty:
        validate_registry_quantity(sfs)
        if sfs.dimensionality != _converters_viscosity_SFS_DIM:
            raise pint.DimensionalityError(sfs.units, ureg.SFS, sfs.dimensionality, _converters_viscosity_SFS_DIM)
        _converters__interpolation_validate_real_scalar(sfs.magnitude, "sfs")
        sfs_val = float(sfs.to(ureg.SFS).magnitude)
    else:
        sfs_val = float(_converters__interpolation_validate_real_scalar(sfs, "sfs"))

    if sfs_val < _converters_viscosity_SFS_TABLE_122F.min_y - 1e-12:
        raise DomainError(
            f"Saybolt Furol reading {sfs_val} s is below the ASTM D2161 minimum domain "
            f"({_converters_viscosity_SFS_TABLE_122F.min_y} s / {_converters_viscosity_SFS_TABLE_122F.min_x} cSt at 122 °F)."
        )
    sfs_val = max(sfs_val, _converters_viscosity_SFS_TABLE_122F.min_y)

    if sfs_val < _converters_viscosity_SFS_TABLE_122F.max_y:
        nu_cst = _converters_viscosity_SFS_TABLE_122F.invert(sfs_val)
    else:
        nu_cst = sfs_val / ASTM_D2161_SFS_LINEAR_FACTOR

    if is_qty:
        dst = parse_unit(target_unit)
        if dst.dimensionality != _converters_viscosity_KINEMATIC_DIM:
            raise pint.DimensionalityError(ureg.centistokes, dst, _converters_viscosity_KINEMATIC_DIM, dst.dimensionality)
        converted = ureg.Quantity(nu_cst, ureg.centistokes).to(dst)
        return ureg.Quantity(_converters__interpolation_match_scalar_type(converted.magnitude, sfs), dst)
    return _converters__interpolation_match_scalar_type(nu_cst, sfs)


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# Unified Two-Stage Viscosity Pipeline Router
# ---------------------------------------------------------------------------


@authorized_conversion
def convert_viscosity(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    density: pint.Quantity | None = None,
    temperature: pint.Quantity | float | int | Decimal | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Two-stage viscosity router: Dynamic (cP/Pa*s) <-> Kinematic (cSt) <-> SUS / SFS.

    - ``Dynamic <-> Kinematic`` requires ``density``.
    - ``Kinematic <-> SUS / SFS`` requires ``temperature`` (in °F or Pint temperature Quantity).
    - ``Dynamic <-> SUS / SFS`` requires both ``density`` and ``temperature`` and routes
      strictly through kinematic viscosity (cSt).
    """
    from core.units.exceptions import MissingAssayError
    from core.units.exceptions import MissingReferenceConditionsError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    from core.units.precision import ensure_quantity_ieee754

    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"quantity must be a pint.Quantity, got {type(quantity).__name__}.")
    validate_registry_quantity(quantity)
    if isinstance(density, pint.Quantity):
        validate_registry_quantity(density)
    if isinstance(temperature, pint.Quantity):
        validate_registry_quantity(temperature)
    mag_dec = _converters__interpolation_validate_real_scalar(quantity.magnitude, "viscosity")
    if mag_dec <= 0:
        raise DomainError(f"Viscosity must be strictly positive, got {quantity}.")

    src_dim = quantity.dimensionality
    dst = parse_unit(target_unit)
    dst_dim = dst.dimensionality

    valid_dims = {_converters_viscosity_DYNAMIC_DIM, _converters_viscosity_KINEMATIC_DIM, _converters_viscosity_SUS_DIM, _converters_viscosity_SFS_DIM}
    if src_dim not in valid_dims:
        raise pint.DimensionalityError(quantity.units, ureg.centistokes, src_dim, _converters_viscosity_KINEMATIC_DIM)
    if dst_dim not in valid_dims:
        raise pint.DimensionalityError(quantity.units, dst, src_dim, dst_dim)

    if src_dim == dst_dim:
        out_qty = ensure_quantity_ieee754(quantity.to(dst))
        if not return_metadata:
            return out_qty
        from core.units.provenance import ConversionClassification
        from core.units.provenance import build_conversion_result
        return build_conversion_result(
            out_qty,
            quantity,
            method="convert_viscosity",
            classification=ConversionClassification.EXACT,
            uncertainties=uncertainties,
            covariance=covariance,
            correlations=correlations,
        )

    # Step A: Reduce source to kinematic viscosity (cSt)
    path_steps: list[str] = [str(quantity.units)]
    if src_dim == _converters_viscosity_DYNAMIC_DIM:
        if density is None:
            raise MissingAssayError("Converting dynamic viscosity requires fluid density.")
        kin_cst = dynamic_to_kinematic_viscosity(quantity, density, target_unit="centistokes")
        path_steps.append("centistokes")
    elif src_dim == _converters_viscosity_SUS_DIM:
        if temperature is None:
            raise MissingReferenceConditionsError("Converting Saybolt Universal Seconds (SUS) requires temperature.")
        kin_cst = sus_to_cst(quantity, temp_f=temperature, target_unit="centistokes")  # type: ignore[assignment]
        path_steps.append("centistokes")
    elif src_dim == _converters_viscosity_SFS_DIM:
        if temperature is None:
            raise MissingReferenceConditionsError("Converting Saybolt Furol Seconds (SFS) requires temperature.")
        kin_cst = sfs_to_cst(quantity, temp_f=temperature, target_unit="centistokes")  # type: ignore[assignment]
        path_steps.append("centistokes")
    else:
        kin_cst = quantity.to(ureg.centistokes)

    # Step B: Map kinematic viscosity (cSt) to target regime
    if dst_dim == _converters_viscosity_KINEMATIC_DIM:
        converted = kin_cst.to(dst)
        out_qty = ensure_quantity_ieee754(ureg.Quantity(_converters__interpolation_match_scalar_type(converted.magnitude, quantity, density), dst))
    elif dst_dim == _converters_viscosity_DYNAMIC_DIM:
        if density is None:
            raise MissingAssayError("Converting to dynamic viscosity requires fluid density.")
        out_qty = ensure_quantity_ieee754(kinematic_to_dynamic_viscosity(kin_cst, density, target_unit=dst))
    elif dst_dim == _converters_viscosity_SUS_DIM:
        if temperature is None:
            raise MissingReferenceConditionsError("Converting to Saybolt Universal Seconds (SUS) requires temperature.")
        res_sus = cst_to_sus(kin_cst, temp_f=temperature)
        assert isinstance(res_sus, pint.Quantity)
        out_qty = ensure_quantity_ieee754(res_sus.to(dst))
    else:
        if temperature is None:
            raise MissingReferenceConditionsError("Converting to Saybolt Furol Seconds (SFS) requires temperature.")
        res_sfs = cst_to_sfs(kin_cst, temp_f=temperature)
        assert isinstance(res_sfs, pint.Quantity)
        out_qty = ensure_quantity_ieee754(res_sfs.to(dst))

    if path_steps[-1] != str(out_qty.units):
        path_steps.append(str(out_qty.units))
    if not return_metadata:
        return out_qty

    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    uses_saybolt = bool({src_dim, dst_dim} & {_converters_viscosity_SUS_DIM, _converters_viscosity_SFS_DIM})
    temp_qty = (
        temperature
        if isinstance(temperature, pint.Quantity)
        else (ureg.Quantity(float(temperature), "degF") if temperature is not None else None)
    )
    cls_tag = (
        ConversionClassification.STATE_AND_MODEL_DEPENDENT
        if uses_saybolt
        else ConversionClassification.ASSAY_DEPENDENT
    )
    ref = ReferenceConditions(
        temperature=temp_qty,
        density=density,
        classification=cls_tag.value,
        parameters={"saybolt_empirical": uses_saybolt},
    )
    return build_conversion_result(
        out_qty,
        quantity,
        conversion_path=tuple(path_steps),
        reference_conditions=ref,
        method="convert_viscosity",
        standard="ASTM D2161" if uses_saybolt else "SI_kinematic_viscosity",
        estimated=uses_saybolt,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
    )


__all__ = ('ASTM_D2161_SFS_LINEAR_FACTOR', 'ASTM_D2161_SUS_LINEAR_FACTOR', 'ASTM_D2161_SUS_TEMP_COEFF', 'BARREL_TO_M3', 'BatteryMaterialAssay', 'BiofuelAssay', 'BiomassAssay', 'BitumenAssay', 'BulkOreAssay', 'BullionAssay', 'ChemicalAssay', 'CleanFuelAssay', 'CoalAssay', 'CrudeAssay', 'DistillationCut', 'DomainError', 'FerroalloyAssay', 'FertilizerAssay', 'GasAssay', 'GasReferenceConditions', 'HydrogenAssay', 'MT_PER_BBL_WATER', 'MetalAlloyAssay', 'NGLAssay', 'OreAssay', 'PetcokeAssay', 'PetrochemAssay', 'RefineryCrudeAssay', 'STANDARD_ATMOSPHERE_BAR', 'STANDARD_ATMOSPHERE_PA', 'STANDARD_ATMOSPHERE_PSI', 'SolverBackendError', 'SolverBackendUnavailableError', 'SolverConvergenceError', 'SolverDomainError', 'SolverError', 'SolverResult', 'UraniumAssay', 'WATER_DENSITY_60F_KG_M3', 'WATER_DENSITY_60F_LB_GAL', 'atmospheric_pressure_at_elevation', 'calculate_z_factor', 'convert', 'convert_actual_to_standard_gas_volume', 'convert_battery_material', 'convert_biofuel', 'convert_biomass', 'convert_bitumen', 'convert_bulk_ore', 'convert_bullion', 'convert_chemical', 'convert_commodity', 'convert_commodity_quote', 'convert_currency', 'convert_dry_to_wet_mass', 'convert_ferroalloy', 'convert_fertilizer', 'convert_fx', 'convert_gas_energy', 'convert_gas_reference_volume', 'convert_hydrocarbon', 'convert_hydrogen', 'convert_metal_alloy', 'convert_ngl', 'convert_ore', 'convert_ore_moisture_price', 'convert_petcoke', 'convert_petrochem', 'convert_pressure', 'convert_pressure_datum', 'convert_standard_to_actual_gas_volume', 'convert_uranium', 'convert_viscosity', 'convert_wet_to_dry_mass', 'create_battery_material_context', 'create_biofuel_context', 'create_biomass_context', 'create_bitumen_context', 'create_bullion_context', 'create_chemical_context', 'create_clean_fuel_context', 'create_coal_heating_context', 'create_ferroalloy_context', 'create_fertilizer_context', 'create_fx_context', 'create_gas_heating_context', 'create_hydrocarbon_context', 'create_hydrogen_context', 'create_metal_alloy_context', 'create_ngl_context', 'create_ore_context', 'create_petcoke_context', 'create_petrochem_context', 'create_uranium_context', 'cst_to_sfs', 'cst_to_sus', 'differential_pressure', 'dynamic_to_kinematic_viscosity', 'extract_contained_metal', 'extract_dry_mass_from_contained', 'kinematic_to_dynamic_viscosity', 'sfs_to_cst', 'solve_z_factor', 'sus_to_cst', 'triangulate_fx_rate')
_PUBLIC_ALIASES = {'authorized_conversion': ('validation', 'authorized_conversion'), 'BiofuelAssay': ('assays.energy.biofuels', 'BiofuelAssay'), 'GAS_MOLAR_MASSES': ('composition', 'GAS_MOLAR_MASSES'), 'CURRENCY_MINOR_DIGITS': ('registry', 'CURRENCY_MINOR_DIGITS'), 'NGLAssay': ('assays.energy.ngl', 'NGLAssay'), 'validate_registry_quantity': ('registry', 'validate_registry_quantity'), 'MissingAssayError': ('exceptions', 'MissingAssayError'), '_solve_z_dak': ('thermodynamics', '_solve_z_dak'), 'SolverConvergenceError': ('thermodynamics', 'SolverConvergenceError'), 'SolverResult': ('thermodynamics', 'SolverResult'), '_QuantitySnapshot': ('validation', '_safety_QuantitySnapshot'), 'decompose_quote': ('registry', 'decompose_quote'), 'FerroalloyAssay': ('assays.metals.ferroalloys', 'FerroalloyAssay'), 'GasAssay': ('assays.energy.natural_gas', 'GasAssay'), 'ensure_quantity_ieee754': ('precision', 'ensure_quantity_ieee754'), 'ProcessingSpecification': ('processing', 'ProcessingSpecification'), 'PetcokeAssay': ('assays.energy.petcoke', 'PetcokeAssay'), 'HydrogenAssay': ('assays.energy.hydrogen', 'HydrogenAssay'), 'calculate_z_factor': ('thermodynamics', 'calculate_z_factor'), 'OreAssay': ('assays.metals.ores', 'OreAssay'), 'SolverBackendUnavailableError': ('thermodynamics', 'SolverBackendUnavailableError'), 'MissingReferenceConditionsError': ('exceptions', 'MissingReferenceConditionsError'), 'validate_dimensionality': ('registry', 'validate_dimensionality'), '_universal_wet_to_dry_mass': ('analytical_basis', 'wet_to_dry_mass'), 'sanitize_unit_string': ('registry', 'sanitize_unit_string'), 'CURRENCY_CODES': ('registry', 'CURRENCY_CODES'), 'CommodityUnitError': ('exceptions', 'CommodityUnitError'), '_universal_dry_matter_fraction': ('analytical_basis', 'dry_matter_fraction'), 'MissingFXRateError': ('exceptions', 'MissingFXRateError'), 'MaterialSpec': ('assays.base', 'MaterialSpec'), 'BiomassAssay': ('assays.energy.biomass', 'BiomassAssay'), 'ConversionConvergenceError': ('exceptions', 'ConversionConvergenceError'), 'InvalidPhysicalStateError': ('exceptions', 'InvalidPhysicalStateError'), 'BitumenAssay': ('assays.energy.bitumen', 'BitumenAssay'), 'MaterialComposition': ('composition', 'MaterialComposition'), '_universal_dry_to_wet_mass': ('analytical_basis', 'dry_to_wet_mass'), 'SolverError': ('thermodynamics', 'SolverError'), 'pad_ieee754_decimal': ('precision', 'pad_ieee754_decimal'), 'RefineryCrudeAssay': ('assays.energy.refining', 'RefineryCrudeAssay'), 'IncompatibleCommodityError': ('exceptions', 'IncompatibleCommodityError'), '_universal_convert_moisture_mass': ('analytical_basis', 'convert_moisture_mass'), 'SolverDomainError': ('thermodynamics', 'SolverDomainError'), 'validated_assay': ('validation', 'validated_assay'), 'BullionAssay': ('assays.metals.precious_metals', 'BullionAssay'), '_solve_z_coolprop': ('thermodynamics', '_solve_z_coolprop'), 'FertilizerAssay': ('assays.industrial.fertilisers', 'FertilizerAssay'), 'BatteryMaterialAssay': ('assays.metals.battery_materials', 'BatteryMaterialAssay'), 'CrudeAssay': ('assays.energy.petroleum', 'CrudeAssay'), 'StandardsApplicabilityError': ('exceptions', 'StandardsApplicabilityError'), 'PetrochemAssay': ('assays.industrial.petrochemicals', 'PetrochemAssay'), 'MetalAlloyAssay': ('assays.metals.alloys', 'MetalAlloyAssay'), 'Quantity': ('registry', 'Quantity'), 'CleanFuelAssay': ('assays.energy.biofuels', 'CleanFuelAssay'), 'GasReferenceConditions': ('reference_conditions', 'GasReferenceConditions'), 'ChemicalAssay': ('assays.industrial.chemicals', 'ChemicalAssay'), 'SolverBackendError': ('thermodynamics', 'SolverBackendError'), 'material_spec_from_assay': ('assays.base', 'material_spec_from_assay'), 'IncompatibleAnalyticalBasisError': ('exceptions', 'IncompatibleAnalyticalBasisError'), '_universal_rebase_moisture_intensive': ('analytical_basis', 'rebase_moisture_intensive'), 'CoalAssay': ('assays.energy.coal', 'CoalAssay'), 'UraniumAssay': ('assays.energy.nuclear', 'UraniumAssay'), 'Unit': ('registry', 'Unit'), 'ureg': ('registry', 'ureg'), 'solve_z_factor': ('thermodynamics', 'solve_z_factor'), 'BulkOreAssay': ('assays.metals.ores', 'BulkOreAssay'), '_reference_document': ('standards', '_reference_reference_document'), 'DistillationCut': ('assays.energy.refining', 'DistillationCut'), 'parse_unit': ('registry', 'parse_unit'), 'MaterialProperties': ('composition', 'MaterialProperties'), 'create_material_context': ('context', 'create_material_context'), 'UnsupportedConversionError': ('exceptions', 'UnsupportedConversionError'), 'CompositionBasis': ('composition', 'CompositionBasis'), '_solve_z_aga8_detail': ('thermodynamics', '_solve_z_aga8_detail')}

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
    from core.units.assays.metals.battery_materials import BatteryMaterialAssay
    from core.units.assays.energy.biofuels import BiofuelAssay
    from core.units.assays.energy.biomass import BiomassAssay
    from core.units.assays.energy.bitumen import BitumenAssay
    from core.units.assays.metals.ores import BulkOreAssay
    from core.units.assays.metals.precious_metals import BullionAssay
    from core.units.assays.industrial.chemicals import ChemicalAssay
    from core.units.assays.energy.biofuels import CleanFuelAssay
    from core.units.assays.energy.coal import CoalAssay
    from core.units.assays.energy.petroleum import CrudeAssay
    from core.units.assays.energy.refining import DistillationCut
    from core.units.assays.metals.ferroalloys import FerroalloyAssay
    from core.units.assays.industrial.fertilisers import FertilizerAssay
    from core.units.assays.energy.natural_gas import GasAssay
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.heat_content import HeatContentAssay
    from core.units.assays.energy.hydrogen import HydrogenAssay
    from core.units.assays.energy.petroleum import LiquidCorrectionMethod
    from core.units.assays.energy.petroleum import LiquidState
    from core.units.assays.metals.alloys import MetalAlloyAssay
    from core.units.thermodynamics import MonotonicCubicTable
    from core.units.assays.energy.ngl import NGLAssay
    from core.units.assays.metals.ores import OreAssay
    from core.units.assays.energy.petcoke import PetcokeAssay
    from core.units.assays.industrial.petrochemicals import PetrochemAssay
    from core.units.assays.energy.refining import RefineryCrudeAssay
    from core.units.thermodynamics import SolverBackendError
    from core.units.thermodynamics import SolverBackendUnavailableError
    from core.units.thermodynamics import SolverConvergenceError
    from core.units.thermodynamics import SolverDomainError
    from core.units.thermodynamics import SolverError
    from core.units.thermodynamics import SolverResult
    from core.units.registry import Unit
    from core.units.assays.energy.nuclear import UraniumAssay
    from core.units.thermodynamics import calculate_z_factor
    from core.units.thermodynamics import solve_z_factor
