"""commodities.core.units.conversions.

Physical, thermodynamic, and financial transformation engines.
Bridges isolated dimensional baselines (API gravity, heating values,
FX spot matrices, pressure datum shifts, and bulk ore grades).
"""

from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
import math
import inspect
import warnings
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final

import pint
from ._safety import authorized_conversion, validated_assay

from .registry import (CURRENCY_CODES, CURRENCY_MINOR_DIGITS, Quantity, Unit, decompose_quote, parse_unit, sanitize_unit_string, ureg, validate_dimensionality,)

if TYPE_CHECKING:
    from .heat_content import HeatContentAssay

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

# API MPMS Table 6A/6B Thermal Expansion Empirical Constants (K0, K1, K2)
# Used for Volume Correction Factor (VCF) temperature adjustments: alpha_60 = K0/rho^2 + K1/rho + K2
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
MOLAR_VOLUME_0C_1ATM_M3_KMOL: Final[float] = 22.41396954
# Ideal gas standard molar volume at 15°C (288.15 K) and 1 atm (101.325 kPa) in m³/kmol
MOLAR_VOLUME_15C_1ATM_M3_KMOL: Final[float] = 23.64480021
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
# 2. Hydrocarbon Physical Assay Specification (Crude Oil & Refined Products)
# ---------------------------------------------------------------------------

@validated_assay
@dataclass(frozen=True, slots=True)
class CrudeAssay:
    """Immutable physical specification for a crude oil or refined liquid hydrocarbon.

    Solves API Gravity, Specific Gravity (60°F/60°F), volumetric densities,
    barrel-to-tonne conversion factors, and API MPMS Ch. 11.1 thermal Volume
    Correction Factors (VCF) from any single density parameter.
    """

    name: str
    api_gravity: float
    specific_gravity: float
    density_kg_per_m3: float
    density_lb_per_gal: float
    bbl_per_mt: float
    mt_per_bbl: float
    product_group: str = "crude"
    sulfur_percentage: float | None = None
    tan: float | None = None

    def __init__(
        self,
        name: str = "custom_assay",
        *,
        api_gravity: float | None = None,
        specific_gravity: float | None = None,
        bbl_per_mt: float | None = None,
        density_kg_per_m3: float | None = None,
        product_group: str = "crude",
        sulfur_percentage: float | None = None,
        tan: float | None = None,
    ) -> None:
        """Constructs an immutable assay from one density/gravity parameter."""
        group = product_group.strip().lower()
        if group not in {"crude", "gasoline", "distillate"}:
            raise ValueError(f"product_group must be 'crude', 'gasoline', or 'distillate', got {product_group!r}")

        provided = sum(p is not None for p in (api_gravity, specific_gravity, density_kg_per_m3, bbl_per_mt))
        if provided != 1:
            raise ValueError(
                "Provide exactly one density parameter: 'api_gravity', 'specific_gravity', 'density_kg_per_m3', or 'bbl_per_mt'."
            )

        resolved_sg: float
        if specific_gravity is not None:
            resolved_sg = float(specific_gravity)
        elif api_gravity is not None:
            api_value = float(api_gravity)
            if not math.isfinite(api_value) or api_value <= -131.5:
                raise ValueError(f"API gravity must be finite and greater than -131.5, got {api_value}")
            resolved_sg = 141.5 / (api_value + 131.5)
        elif density_kg_per_m3 is not None:
            resolved_sg = float(density_kg_per_m3) / WATER_DENSITY_60F_KG_M3
        elif bbl_per_mt is not None:
            bbl_mt = float(bbl_per_mt)
            if not math.isfinite(bbl_mt) or bbl_mt <= 0:
                raise ValueError(f"bbl_per_mt must be finite and strictly positive, got {bbl_mt}")
            resolved_sg = (1.0 / bbl_mt) / MT_PER_BBL_WATER
        else:
            raise ValueError(
                "Must provide at least one density parameter: "
                "'api_gravity', 'specific_gravity', 'density_kg_per_m3', or 'bbl_per_mt'."
            )

        if not math.isfinite(resolved_sg) or resolved_sg <= 0:
            raise ValueError(f"Calculated Specific Gravity must be finite and strictly positive, got {resolved_sg}")

        calc_api = (141.5 / resolved_sg) - 131.5
        calc_rho_kg_m3 = resolved_sg * WATER_DENSITY_60F_KG_M3
        calc_rho_lb_gal = resolved_sg * WATER_DENSITY_60F_LB_GAL
        calc_mt_per_bbl = resolved_sg * MT_PER_BBL_WATER
        calc_bbl_per_mt = 1.0 / calc_mt_per_bbl

        if sulfur_percentage is not None and (sulfur_percentage < 0 or sulfur_percentage > 100):
            raise ValueError(f"Sulfur percentage must be between 0 and 100, got {sulfur_percentage}")
        if tan is not None and tan < 0:
            raise ValueError(f"Total Acid Number (TAN) cannot be negative, got {tan}")

        object.__setattr__(self, "name", str(name))
        object.__setattr__(self, "specific_gravity", (resolved_sg))
        object.__setattr__(self, "api_gravity", (calc_api))
        object.__setattr__(self, "density_kg_per_m3", (calc_rho_kg_m3))
        object.__setattr__(self, "density_lb_per_gal", (calc_rho_lb_gal))
        object.__setattr__(self, "bbl_per_mt", (calc_bbl_per_mt))
        object.__setattr__(self, "mt_per_bbl", (calc_mt_per_bbl))
        object.__setattr__(self, "product_group", group)
        object.__setattr__(self, "sulfur_percentage", sulfur_percentage)
        object.__setattr__(self, "tan", tan)

    @property
    def density_in_air_kg_per_m3(self) -> float:
        """Apparent density in air at 60°F (ASTM D1250 Table 56 buoyancy deduction ~1.1 kg/m³)."""
        return (self.density_kg_per_m3 - 1.1)

    @property
    def gravity_class(self) -> str:
        """Market API gravity classification ('light', 'medium', 'heavy', or 'extra_heavy')."""
        if self.api_gravity > 31.1:
            return "light"
        if self.api_gravity >= 22.3:
            return "medium"
        if self.api_gravity >= 10.0:
            return "heavy"
        return "extra_heavy"

    @property
    def is_sweet(self) -> bool | None:
        """Returns True if sulfur_percentage <= 0.5%, False if > 0.5%, or None if unspecified."""
        if self.sulfur_percentage is None:
            return None
        return self.sulfur_percentage <= 0.5

    def as_context(self, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume] and [mass] for this assay."""
        return create_hydrocarbon_context(self, name=context_name)

    def bbl_to_mt(self, barrels: float | Decimal, *, in_air: bool = False) -> float:
        """Converts volumetric barrels directly to metric tonnes (in vacuum or in air)."""
        if in_air:
            mt_per_bbl_air = (self.density_in_air_kg_per_m3 * BARREL_TO_M3) / 1000.0
            return float(barrels) * mt_per_bbl_air
        return float(barrels) * self.mt_per_bbl

    def mt_to_bbl(self, metric_tonnes: float | Decimal, *, in_air: bool = False) -> float:
        """Converts metric tonnes (in vacuum or in air) directly to volumetric barrels."""
        if in_air:
            mt_per_bbl_air = (self.density_in_air_kg_per_m3 * BARREL_TO_M3) / 1000.0
            return float(metric_tonnes) / mt_per_bbl_air
        return float(metric_tonnes) * self.bbl_per_mt

    def volume_correction_factor(self, observed_temp_f: float | pint.Quantity) -> float:
        """Computes API MPMS Table 6A/6B Volume Correction Factor (VCF) to 60°F.

        VCF = exp(-alpha_60 * delta_T * (1 + 0.8 * alpha_60 * delta_T))
        Multiply observed gross volume at `observed_temp_f` (float °F or temperature Quantity)
        by VCF to obtain net standard volume at 60°F.
        """
        if isinstance(observed_temp_f, pint.Quantity):
            temp_f = float(observed_temp_f.to("degF").magnitude)
        else:
            temp_f = float(observed_temp_f)

        if not math.isfinite(temp_f):
            raise ValueError(f"Observed temperature must be finite, got {observed_temp_f!r}")

        rho = self.density_kg_per_m3
        if self.product_group == "gasoline":
            k0, k1, k2 = API_VCF_GASOLINE_K0, API_VCF_GASOLINE_K1, 0.0
        elif self.product_group == "distillate":
            k0, k1, k2 = API_VCF_DISTILLATE_K0, API_VCF_DISTILLATE_K1, 0.0
        else:
            k0, k1, k2 = API_VCF_CRUDE_K0, API_VCF_CRUDE_K1, API_VCF_CRUDE_K2

        alpha_60 = (k0 / (rho ** 2)) + (k1 / rho) + k2
        delta_t = temp_f - 60.0
        return (math.exp(-alpha_60 * delta_t * (1.0 + 0.8 * alpha_60 * delta_t)))

    def correct_volume_to_60f(self, observed_volume: pint.Quantity, observed_temp_f: float | pint.Quantity, *, bsw_percentage: float = 0.0,) -> pint.Quantity:
        """Corrects an observed gross hydrocarbon volume at `observed_temp_f` to 60°F GSV/NSV."""
        if not isinstance(observed_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(observed_volume).__name__}")
        if not validate_dimensionality(observed_volume, "barrel"):
            raise pint.DimensionalityError(observed_volume.units, ureg.oil_barrel, observed_volume.dimensionality, ureg.oil_barrel.dimensionality,)
        if not (0.0 <= bsw_percentage < 100.0):
            raise ValueError(f"bsw_percentage must be in [0, 100), got {bsw_percentage}")

        vcf = self.volume_correction_factor(observed_temp_f)
        net_factor = vcf * (1.0 - (bsw_percentage / 100.0))
        scale = _scale_factor(observed_volume.magnitude, net_factor)
        return ureg.Quantity(observed_volume.magnitude * scale, observed_volume.units)

    def net_standard_volume(self, observed_volume: pint.Quantity, observed_temp_f: float | pint.Quantity, bsw_pct: float = 0.0,) -> pint.Quantity:
        """Computes Net Standard Volume (NSV) at 60°F after VCF and BS&W deduction."""
        return self.correct_volume_to_60f(observed_volume, observed_temp_f, bsw_percentage=bsw_pct,)

    def blend_with(self, other: CrudeAssay, *, self_volume_fraction: float, blend_name: str | None = None,) -> CrudeAssay:
        """Computes ideal volumetric blend assay (SG, API, mass-weighted sulfur & TAN)."""
        if not isinstance(other, CrudeAssay):
            raise TypeError(f"Expected CrudeAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")
        w_other = 1.0 - self_volume_fraction
        blended_sg = (self.specific_gravity * self_volume_fraction) + (other.specific_gravity * w_other)

        # Sulfur and TAN are mass-fraction properties
        mass_self = self.specific_gravity * self_volume_fraction
        mass_other = other.specific_gravity * w_other
        total_mass = mass_self + mass_other

        blended_sulfur: float | None = None
        if self.sulfur_percentage is not None and other.sulfur_percentage is not None:
            blended_sulfur = (((self.sulfur_percentage * mass_self) + (other.sulfur_percentage * mass_other)) / total_mass)

        blended_tan: float | None = None
        if self.tan is not None and other.tan is not None:
            blended_tan = (((self.tan * mass_self) + (other.tan * mass_other)) / total_mass)

        return CrudeAssay(name=blend_name or f"{self.name}/{other.name} Blend", specific_gravity=blended_sg, product_group=self.product_group, sulfur_percentage=blended_sulfur, tan=blended_tan,)

    # --- Crude Oil Benchmark Presets ---
    @classmethod
    def brent(cls) -> CrudeAssay:
        return cls("Brent", api_gravity=38.3, sulfur_percentage=0.40)

    @classmethod
    def wti(cls) -> CrudeAssay:
        return cls("WTI", api_gravity=40.8, sulfur_percentage=0.24)

    @classmethod
    def dubai(cls) -> CrudeAssay:
        return cls("Dubai", api_gravity=31.0, sulfur_percentage=2.13)

    @classmethod
    def arab_light(cls) -> CrudeAssay:
        return cls("Arab Light", api_gravity=33.4, sulfur_percentage=1.77)

    @classmethod
    def urals(cls) -> CrudeAssay:
        return cls("Urals", api_gravity=31.7, sulfur_percentage=1.70)

    @classmethod
    def maya(cls) -> CrudeAssay:
        return cls("Maya", api_gravity=21.8, sulfur_percentage=3.40)

    @classmethod
    def wcs(cls) -> CrudeAssay:
        return cls("WCS", api_gravity=20.5, sulfur_percentage=3.50)

    @classmethod
    def bonny_light(cls) -> CrudeAssay:
        return cls("Bonny Light", api_gravity=35.3, sulfur_percentage=0.15)

    # --- Refined Products & Marine Bunker Presets ---
    @classmethod
    def rbob_gasoline(cls) -> CrudeAssay:
        """NYMEX RBOB Reformulated Blendstock (~60.0 API, ~8.45 bbl/mt)."""
        return cls("RBOB Gasoline", api_gravity=60.0, product_group="gasoline", sulfur_percentage=0.001)

    @classmethod
    def ulsd_diesel(cls) -> CrudeAssay:
        """Ultra-Low Sulfur Diesel / ICE Gasoil (~36.5 API, ~7.45 bbl/mt)."""
        return cls("ULSD Diesel", api_gravity=36.5, product_group="distillate", sulfur_percentage=0.001)

    @classmethod
    def jet_a1(cls) -> CrudeAssay:
        """Aviation Turbine Fuel Jet A-1 (~43.5 API, ~7.88 bbl/mt)."""
        return cls("Jet A-1", api_gravity=43.5, product_group="distillate", sulfur_percentage=0.05)

    @classmethod
    def vlsfo_bunker(cls) -> CrudeAssay:
        """IMO 2020 Very Low Sulfur Fuel Oil 0.5% S (~18.5 API, ~6.35 bbl/mt)."""
        return cls("VLSFO 0.5%", api_gravity=18.5, product_group="crude", sulfur_percentage=0.50)


# ---------------------------------------------------------------------------
# Refutas Viscosity Blending Index (Standard Refinery LP Approximation)
# ---------------------------------------------------------------------------

def _viscosity_blending_number(viscosity_cst: float) -> float:
    """Computes Refutas Viscosity Blending Number (VBN) from kinematic viscosity (cSt)."""
    if not math.isfinite(viscosity_cst) or viscosity_cst <= 0.0:
        raise ValueError(f"Kinematic viscosity must be strictly positive, got {viscosity_cst}")
    return 14.534 * math.log(math.log(viscosity_cst + 0.8)) + 10.975


def _viscosity_from_vbn(vbn: float) -> float:
    """Inverts Refutas VBN back to kinematic viscosity (cSt)."""
    return math.exp(math.exp((vbn - 10.975) / 14.534)) - 0.8


# ---------------------------------------------------------------------------
# Distillation Cut & Full Refinery Crude Assay
# ---------------------------------------------------------------------------

@validated_assay
@dataclass(frozen=True, slots=True)
class DistillationCut:
    """Single boiling-point fraction from a True Boiling Point (TBP) assay."""

    name: str                  # e.g. "lpg", "naphtha", "jet", "gasoil", "vgo", "vac_resid"
    tbp_start_c: float         # Initial boiling point (°C)
    tbp_end_c: float           # Final boiling point (°C); use math.inf for Vacuum Residue
    yield_vol_pct: float       # Volumetric yield (% of whole crude)
    yield_wt_pct: float        # Mass yield (% of whole crude)
    specific_gravity: float    # Cut SG at 60°F/60°F
    sulfur_wt_pct: float | None = None
    viscosity_cst_50c: float | None = None
    nickel_ppm: float | None = None
    vanadium_ppm: float | None = None
    ccr_wt_pct: float | None = None  # Conradson Carbon Residue (wt%)

    def __post_init__(self) -> None:
        if self.tbp_end_c <= self.tbp_start_c:
            raise ValueError(f"Cut {self.name!r}: tbp_end_c must exceed tbp_start_c.")
        if not (0.0 <= self.yield_vol_pct <= 100.0) or not (0.0 <= self.yield_wt_pct <= 100.0):
            raise ValueError(f"Cut {self.name!r}: yields must be in [0, 100].")
        if not math.isfinite(self.specific_gravity) or self.specific_gravity <= 0.0:
            raise ValueError(f"Cut {self.name!r}: specific_gravity must be strictly positive.")


@validated_assay
@dataclass(frozen=True, slots=True)
class RefineryCrudeAssay:
    """Downstream refining specification composing a bulk CrudeAssay with TBP cuts.

    Drives Gross Product Worth (GPW), custom crack spreads, refinery LP blending,
    and contaminant constraints (Ni/V catalyst poisoning, CCR, viscosity).
    """

    bulk: CrudeAssay
    cuts: tuple[DistillationCut, ...]
    viscosity_cst_50c: float | None = None
    pour_point_c: float | None = None
    nickel_ppm: float | None = None
    vanadium_ppm: float | None = None
    ccr_wt_pct: float | None = None
    rvp_kpa: float | None = None  # Reid Vapour Pressure

    def __post_init__(self) -> None:
        if not isinstance(self.bulk, CrudeAssay):
            raise TypeError(f"Expected CrudeAssay for 'bulk', got {type(self.bulk).__name__}")
        if not self.cuts:
            raise ValueError("RefineryCrudeAssay requires at least one DistillationCut.")

        vol_sum = sum(c.yield_vol_pct for c in self.cuts)
        wt_sum = sum(c.yield_wt_pct for c in self.cuts)
        if not (99.5 <= vol_sum <= 100.5):
            raise ValueError(f"Volumetric cut yields must sum to ~100% (got {vol_sum:.2f}%).")
        if not (99.5 <= wt_sum <= 100.5):
            raise ValueError(f"Mass cut yields must sum to ~100% (got {wt_sum:.2f}%).")

        # Verify boiling-point ordering
        for prev, curr in zip(self.cuts, self.cuts[1:]):
            if curr.tbp_start_c < prev.tbp_end_c:
                raise ValueError(
                    f"Overlapping TBP cuts: {prev.name!r} ends at {prev.tbp_end_c}°C, "
                    f"but {curr.name!r} starts at {curr.tbp_start_c}°C."
                )

    def cut_yields_vol(self) -> dict[str, float]:
        """Returns a mapping of cut name -> volumetric yield fraction (0.0 to 1.0)."""
        return {c.name: c.yield_vol_pct / 100.0 for c in self.cuts}

    def gross_product_worth(
        self,
        product_prices_per_bbl: Mapping[str, float],
        *,
        refinery_fuel_loss_pct: float = 0.0,
    ) -> float:
        """Computes Gross Product Worth (GPW) in $/bbl from market product prices.

        GPW = (1 - fuel_loss) * sum(yield_vol_fraction_i * price_i)
        """
        if not (0.0 <= refinery_fuel_loss_pct < 25.0):
            raise ValueError(f"refinery_fuel_loss_pct must be in [0, 25), got {refinery_fuel_loss_pct}")

        gpw = 0.0
        for cut in self.cuts:
            if cut.name not in product_prices_per_bbl:
                raise KeyError(f"Missing market price for distillation cut {cut.name!r}.")
            gpw += (cut.yield_vol_pct / 100.0) * float(product_prices_per_bbl[cut.name])

        return (gpw * (1.0 - refinery_fuel_loss_pct / 100.0))

    def netback_margin(
        self,
        product_prices_per_bbl: Mapping[str, float],
        crude_delivered_price_per_bbl: float,
        *,
        variable_processing_cost_per_bbl: float = 0.0,
        refinery_fuel_loss_pct: float = 0.0,
    ) -> float:
        """Computes refinery netback / gross refining margin ($/bbl)."""
        gpw = self.gross_product_worth(
            product_prices_per_bbl,
            refinery_fuel_loss_pct=refinery_fuel_loss_pct,
        )
        return (gpw - crude_delivered_price_per_bbl - variable_processing_cost_per_bbl)

    def blend_with(
        self,
        other: RefineryCrudeAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> RefineryCrudeAssay:
        """Blends two refinery assays with identical cut structures.

        - Bulk properties: delegated to CrudeAssay.blend_with()
        - Cut yields: volume-weighted (vol%) and mass-weighted (wt%)
        - Contaminants (Ni, V, CCR): mass-weighted
        - Viscosity: Refutas Viscosity Blending Number (non-linear)
        """
        if not isinstance(other, RefineryCrudeAssay):
            raise TypeError(f"Expected RefineryCrudeAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be in [0.0, 1.0].")
        if len(self.cuts) != len(other.cuts) or any(a.name != b.name for a, b in zip(self.cuts, other.cuts)):
            raise ValueError("Both assays must share identical TBP cut names to blend directly.")

        w_self_vol = self_volume_fraction
        w_other_vol = 1.0 - self_volume_fraction

        mass_self = self.bulk.specific_gravity * w_self_vol
        mass_other = other.bulk.specific_gravity * w_other_vol
        total_mass = mass_self + mass_other
        w_self_wt = mass_self / total_mass
        w_other_wt = mass_other / total_mass

        blended_bulk = self.bulk.blend_with(
            other.bulk,
            self_volume_fraction=self_volume_fraction,
            blend_name=blend_name,
        )

        blended_cuts: list[DistillationCut] = []
        for c1, c2 in zip(self.cuts, other.cuts):
            b_vol = (c1.yield_vol_pct * w_self_vol) + (c2.yield_vol_pct * w_other_vol)
            b_wt = (c1.yield_wt_pct * w_self_wt) + (c2.yield_wt_pct * w_other_wt)
            b_sg = (
                (c1.specific_gravity * c1.yield_vol_pct * w_self_vol + c2.specific_gravity * c2.yield_vol_pct * w_other_vol)
                / b_vol
                if b_vol > 0
                else c1.specific_gravity
            )
            blended_cuts.append(
                DistillationCut(
                    name=c1.name,
                    tbp_start_c=c1.tbp_start_c,
                    tbp_end_c=c1.tbp_end_c,
                    yield_vol_pct=(b_vol),
                    yield_wt_pct=(b_wt),
                    specific_gravity=(b_sg),
                )
            )

        # Refutas non-linear viscosity blending (mass-weighted VBN)
        blended_visc: float | None = None
        if self.viscosity_cst_50c is not None and other.viscosity_cst_50c is not None:
            vbn_blend = (
                w_self_wt * _viscosity_blending_number(self.viscosity_cst_50c)
                + w_other_wt * _viscosity_blending_number(other.viscosity_cst_50c)
            )
            blended_visc = (_viscosity_from_vbn(vbn_blend))

        def _mass_blend(v1: float | None, v2: float | None) -> float | None:
            if v1 is None or v2 is None:
                return None
            return ((v1 * w_self_wt) + (v2 * w_other_wt))

        return RefineryCrudeAssay(
            bulk=blended_bulk,
            cuts=tuple(blended_cuts),
            viscosity_cst_50c=blended_visc,
            nickel_ppm=_mass_blend(self.nickel_ppm, other.nickel_ppm),
            vanadium_ppm=_mass_blend(self.vanadium_ppm, other.vanadium_ppm),
            ccr_wt_pct=_mass_blend(self.ccr_wt_pct, other.ccr_wt_pct),
        )


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


def _scale_factor(magnitude: Any, factor: float | Decimal | int) -> Decimal | float:
    """Returns a Decimal or float multiplier matching the incoming Quantity magnitude."""
    if isinstance(magnitude, Decimal):
        return Decimal(str(factor))
    return float(factor)


# ---------------------------------------------------------------------------
# Real-Gas Compressibility (Z-Factor) & Actual-to-Standard Volume Engine
# (Dranchuk–Abou-Kassem EoS + Optional AGA-8 Detail / GERG-2008 Wrapper)
# ---------------------------------------------------------------------------

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

# Canonical 21-component AGA-8 Detail molar ordering
_AGA8_COMPONENTS: Final[tuple[str, ...]] = (
    "CH4", "N2", "CO2", "C2H6", "C3H8", "iC4H10", "nC4H10",
    "iC5H12", "nC5H12", "C6H14", "C7H16", "C8H18", "C9H20", "C10H22",
    "H2", "O2", "CO", "H2O", "H2S", "He", "Ar",
)


def _validated_composition(composition):
    if not composition:
        raise ValueError("A non-empty molar composition is required.")
    unknown = set(composition) - set(_AGA8_COMPONENTS)
    if unknown:
        raise ValueError(f"Unsupported gas components: {sorted(unknown)}")
    values = {name: float(value) for name, value in composition.items()}
    if any(not math.isfinite(value) or value < 0 for value in values.values()):
        raise ValueError("Composition fractions must be finite and non-negative.")
    total = sum(values.values())
    if not math.isclose(total, 1., rel_tol=0., abs_tol=1e-6):
        raise ValueError("Composition fractions must sum to one.")
    return {name: value / total for name, value in values.items() if value > 0}


def _checked_z(value):
    value = float(value)
    if not math.isfinite(value) or value <= 0:
        raise ArithmeticError("Equation-of-state backend returned an invalid Z-factor.")
    return value


def _calculate_z_aga8_detail(p_kpa, t_k, molar_composition):
    """Strict AGA-8 Detail; never substitute a different equation of state."""
    composition = _validated_composition(molar_composition)
    try:
        import pyaga8
    except ImportError as error:
        raise ImportError("method='aga8' requires pyaga8; choose 'heos' or 'gerg2008' explicitly for CoolProp.") from error
    aga = pyaga8.Detail()
    aga.set_composition([composition.get(component, 0.) for component in _AGA8_COMPONENTS])
    aga.p, aga.t = p_kpa, t_k
    aga.density()
    aga.properties()
    return _checked_z(aga.z)


def _calculate_z_coolprop(p_kpa, t_k, molar_composition, backend):
    """Use the explicitly named CoolProp model without backend substitution."""
    norm_comp = _validated_composition(molar_composition)
    from CoolProp.CoolProp import PropsSI
    cp_map = {
        "CH4": "Methane", "N2": "Nitrogen", "CO2": "CarbonDioxide",
        "C2H6": "Ethane", "C3H8": "Propane", "iC4H10": "IsoButane",
        "nC4H10": "n-Butane", "iC5H12": "IsoPentane", "nC5H12": "n-Pentane",
        "C6H14": "n-Hexane", "C7H16": "n-Heptane", "C8H18": "n-Octane",
        "C9H20": "n-Nonane", "C10H22": "n-Decane", "H2": "Hydrogen",
        "O2": "Oxygen", "CO": "CarbonMonoxide", "H2O": "Water",
        "H2S": "HydrogenSulfide", "He": "Helium", "Ar": "Argon",
    }
    mixture = "&".join(f"{cp_map[name]}[{fraction:.17g}]" for name, fraction in norm_comp.items())
    return _checked_z(PropsSI("Z", "P", p_kpa * 1000., "T", t_k, backend + "::" + mixture))


def calculate_z_factor(
    pressure: pint.Quantity,
    temperature: pint.Quantity,
    relative_density: float = 0.60,
    *,
    method: str = "auto",
    molar_composition: Mapping[str, float] | None = None,
    co2_mol_frac: float = 0.0,
    h2s_mol_frac: float = 0.0,
    tol: float = 1e-10,
    max_iter: int = 50,
) -> float:
    """Computes real-gas compressibility factor (Z).

    Methods
    -------
    - "dak"  : Built-in Dranchuk–Abou-Kassem (1975) 11-constant EoS with Sutton (1985)
               pseudocritical properties and Wichert–Aziz (1972) sour-gas corrections.
    - "aga8" : Strict 21-component AGA-8 Detail via optional pyaga8.
    - "heos" / "gerg2008": Explicit CoolProp models; HEOS is not GERG-2008.
               Requires `molar_composition` (e.g. {"CH4": 0.94, "C2H6": 0.04, "N2": 0.02}).
    - "auto" : Uses AGA-8 Detail if `molar_composition` is supplied and a backend is
               installed; otherwise falls back to built-in "dak".
    """
    if not isinstance(pressure, pint.Quantity) or not isinstance(temperature, pint.Quantity):
        raise TypeError("pressure and temperature must be pint.Quantity instances.")

    method_key = method.strip().lower()
    if method_key not in {"auto", "dak", "aga8", "heos", "gerg2008"}:
        raise ValueError(f"method must be 'auto', 'dak', or 'aga8', got {method!r}")

    if not math.isfinite(tol) or tol <= 0 or isinstance(max_iter, bool) or not isinstance(max_iter, int) or max_iter <= 0:
        raise ValueError("tol must be finite/positive and max_iter a positive integer.")
    if molar_composition is not None:
        molar_composition = _validated_composition(molar_composition)
    if method_key in {"heos", "gerg2008"}:
        if molar_composition is None:
            raise ValueError("An explicit composition is required for CoolProp models.")
    pressure = convert_pressure(pressure, ureg.pascal)
    p_kpa = float(pressure.to(ureg.kilopascal).magnitude)
    t_k = float(temperature.to(ureg.kelvin).magnitude)

    if not math.isfinite(p_kpa) or not math.isfinite(t_k) or p_kpa < 0 or t_k <= 0:
        raise ValueError(f"Invalid state: pressure={p_kpa} kPa, temperature={t_k} K.")
    if p_kpa == 0.0:
        return 1.0

    if method_key in {"heos", "gerg2008"}:
        return _calculate_z_coolprop(p_kpa, t_k, molar_composition, "HEOS" if method_key == "heos" else "GERG2008")
    # Route to AGA-8 Detail wrapper if explicitly requested or available in auto mode
    if method_key == "aga8":
        if not molar_composition:
            raise ValueError("method='aga8' requires a non-empty `molar_composition` dictionary.")
        return _calculate_z_aga8_detail(p_kpa, t_k, molar_composition)

    if method_key == "auto" and molar_composition:
        try:
            return _calculate_z_aga8_detail(p_kpa, t_k, molar_composition)
        except ImportError:
            warnings.warn("AGA-8 backend unavailable; auto selected DAK, not AGA-8/GERG.", RuntimeWarning, stacklevel=2)
            # Gracefully fall back to DAK using CO2/H2S from composition if not overridden
            co2_mol_frac = co2_mol_frac or molar_composition.get("CO2", 0.0)
            h2s_mol_frac = h2s_mol_frac or molar_composition.get("H2S", 0.0)

    # Built-in Dranchuk–Abou-Kassem (1975) solver
    if not (0.50 <= relative_density <= 1.50):
        raise ValueError(f"relative_density out of valid hydrocarbon range [0.50, 1.50]: {relative_density}")
    if not (0.0 <= co2_mol_frac < 1.0) or not (0.0 <= h2s_mol_frac < 1.0) or co2_mol_frac + h2s_mol_frac >= 1:
        raise ValueError("CO2 and H2S mole fractions must be in [0, 1).")

    p_psia = float(pressure.to(ureg.psia).magnitude)
    t_rankine = float(temperature.to(ureg.degR).magnitude)

    # 1. Sutton (1985) pseudocritical pressure (psia) and temperature (°R)
    sg = float(relative_density)
    p_pc = 756.8 - 131.0 * sg - 3.6 * (sg ** 2)
    t_pc = 169.2 + 349.5 * sg - 74.0 * (sg ** 2)

    # 2. Wichert–Aziz (1972) sour-gas correction
    a_sour = co2_mol_frac + h2s_mol_frac
    if a_sour > 0.0:
        b_sour = h2s_mol_frac
        eps = 120.0 * (a_sour ** 0.9 - a_sour ** 1.6) + 15.0 * (b_sour ** 0.5 - b_sour ** 4)
        t_pc_adj = t_pc - eps
        p_pc_adj = (p_pc * t_pc_adj) / (t_pc + b_sour * (1.0 - b_sour) * eps)
        p_pc, t_pc = p_pc_adj, t_pc_adj

    if not math.isfinite(p_pc) or not math.isfinite(t_pc) or p_pc <= 0 or t_pc <= 0:
        raise ValueError("Invalid pseudocritical state after sour-gas correction.")

    # 3. Pseudoreduced pressure and temperature
    p_pr = p_psia / p_pc
    t_pr = t_rankine / t_pc

    a1, a2, a3, a4, a5, a6, a7, a8, a9, a10, a11 = _DAK_A
    c1 = a1 + a2 / t_pr + a3 / (t_pr ** 3) + a4 / (t_pr ** 4) + a5 / (t_pr ** 5)
    c2 = a6 + a7 / t_pr + a8 / (t_pr ** 2)
    c3 = a9 * (a7 / t_pr + a8 / (t_pr ** 2))
    c4 = a10 / (t_pr ** 3)
    rhs = 0.27 * p_pr / t_pr

    # 4. Newton–Raphson solver for reduced density (rho_r)
    rho_r = rhs
    for _ in range(max_iter):
        rho2 = rho_r * rho_r
        exp_term = math.exp(-a11 * rho2)

        f_val = (
            rho_r
            + c1 * rho2
            + c2 * (rho_r ** 3)
            - c3 * (rho_r ** 6)
            + c4 * (rho_r ** 3) * (1.0 + a11 * rho2) * exp_term
            - rhs
        )
        df_val = (
            1.0
            + 2.0 * c1 * rho_r
            + 3.0 * c2 * rho2
            - 6.0 * c3 * (rho_r ** 5)
            + c4 * rho2 * (3.0 + a11 * rho2 * (3.0 - 2.0 * a11 * rho2)) * exp_term
        )

        if math.isfinite(f_val) and abs(f_val) <= tol * max(1., abs(rhs)):
            return _checked_z(rhs / rho_r)
        if not math.isfinite(df_val) or abs(df_val) < 1e-15:
            raise ArithmeticError("DAK solver has a singular/nonfinite derivative.")
        step = f_val / df_val
        rho_r = max(1e-12, rho_r - step)
        if not math.isfinite(rho_r):
            raise ArithmeticError("DAK solver produced nonfinite reduced density.")

    raise ArithmeticError(f"DAK solver failed to establish residual convergence in {max_iter} iterations.")


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
    co2_mol_frac: float = 0.0,
    h2s_mol_frac: float = 0.0,
) -> pint.Quantity:
    """Converts actual flowing or storage gas volume to standard reference volume.

    Applies the real-gas law:
        V_std = V_act * (P_act / P_std) * (T_std / T_act) * (Z_std / Z_act)

    Supported `target_standard` baselines:
        - "ferc"   : 60°F, 14.73 psia (NYMEX / US interstate pipelines)
        - "tx"     : 60°F, 14.65 psia (Texas RRC)
        - "iso15"  : 15°C, 101.325 kPa (ISO 13443 / UK NBP / LNG)
        - "din0"   :  0°C, 101.325 kPa (DIN 1343 / Continental Europe)
        - "gost20" : 20°C, 101.325 kPa (GOST 2939-63 / Eurasian pipelines)
    """
    std_key = target_standard.strip().lower()
    gost_unit = ureg.m3_20C_101325kPa_dry
    std_configs: dict[str, tuple[Unit, Unit, float, float]] = {
        "ferc": (ureg.scf_60F_14_73psia, ureg.foot ** 3, _P_14_73_KPA, _T_60F_K),
        "tx": (ureg.scf_60F_14_65psia, ureg.foot ** 3, _P_14_65_KPA, _T_60F_K),
        "iso15": (ureg.Sm3_15C_101325kPa_dry, ureg.meter ** 3, _P_1ATM_KPA, _T_15C_K),
        "din0": (ureg.Nm3_0C_101325kPa_dry, ureg.meter ** 3, _P_1ATM_KPA, _T_0C_K),
        "gost20": (gost_unit, ureg.meter ** 3, _P_1ATM_KPA, _T_20C_K),
    }
    if std_key not in std_configs:
        raise ValueError(f"Unknown target_standard {target_standard!r}. Choose from {tuple(std_configs)}.")

    dst_unit, geom_unit, p_std_kpa, t_std_k = std_configs[std_key]

    actual_pressure = convert_pressure(actual_pressure, ureg.pascal)
    p_act_kpa = float(actual_pressure.to(ureg.kilopascal).magnitude)
    t_act_k = float(actual_temperature.to(ureg.kelvin).magnitude)
    if not math.isfinite(p_act_kpa) or p_act_kpa <= 0 or not math.isfinite(t_act_k) or t_act_k <= 0:
        raise ValueError("Actual gas pressure and absolute temperature must be finite and strictly positive.")

    z_act = z_actual if z_actual is not None else calculate_z_factor(
        actual_pressure,
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
    expansion_factor = (p_act_kpa / p_std_kpa) * (t_std_k / t_act_k) * (z_std / z_act)
    v_geom = actual_volume.to(geom_unit)
    scale = _scale_factor(v_geom.magnitude, expansion_factor)
    return ureg.Quantity(v_geom.magnitude * scale, dst_unit)


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
    co2_mol_frac: float = 0.0,
    h2s_mol_frac: float = 0.0,
) -> pint.Quantity:
    """Converts standard reference gas volume to actual geometric volume at (P_act, T_act).

    Applies the inverse real-gas law:
        V_act = V_std * (P_std / P_act) * (T_act / T_std) * (Z_act / Z_std)

    Supported `source_standard` baselines:
        - "ferc"   : 60°F, 14.73 psia (NYMEX / US interstate pipelines)
        - "tx"     : 60°F, 14.65 psia (Texas RRC)
        - "iso15"  : 15°C, 101.325 kPa (ISO 13443 / UK NBP / LNG)
        - "din0"   :  0°C, 101.325 kPa (DIN 1343 / Continental Europe)
        - "gost20" : 20°C, 101.325 kPa (GOST 2939-63 / Eurasian pipelines)
    """
    std_key = source_standard.strip().lower()
    gost_unit = ureg.m3_20C_101325kPa_dry
    std_configs: dict[str, tuple[Unit, Unit, float, float]] = {
        "ferc": (ureg.scf_60F_14_73psia, ureg.foot ** 3, _P_14_73_KPA, _T_60F_K),
        "tx": (ureg.scf_60F_14_65psia, ureg.foot ** 3, _P_14_65_KPA, _T_60F_K),
        "iso15": (ureg.Sm3_15C_101325kPa_dry, ureg.meter ** 3, _P_1ATM_KPA, _T_15C_K),
        "din0": (ureg.Nm3_0C_101325kPa_dry, ureg.meter ** 3, _P_1ATM_KPA, _T_0C_K),
        "gost20": (gost_unit, ureg.meter ** 3, _P_1ATM_KPA, _T_20C_K),
    }
    if std_key not in std_configs:
        raise ValueError(f"Unknown source_standard {source_standard!r}. Choose from {tuple(std_configs)}.")

    src_std_unit, geom_unit, p_std_kpa, t_std_k = std_configs[std_key]

    actual_pressure = convert_pressure(actual_pressure, ureg.pascal)
    p_act_kpa = float(actual_pressure.to(ureg.kilopascal).magnitude)
    t_act_k = float(actual_temperature.to(ureg.kelvin).magnitude)
    if not math.isfinite(p_act_kpa) or p_act_kpa <= 0 or not math.isfinite(t_act_k) or t_act_k <= 0:
        raise ValueError("Actual gas pressure and absolute temperature must be finite and strictly positive.")
    if p_act_kpa <= 0:
        raise ValueError(f"actual_pressure must be strictly positive, got {p_act_kpa} kPa.")

    z_act = z_actual if z_actual is not None else calculate_z_factor(
        actual_pressure,
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
    compression_factor = (p_std_kpa / p_act_kpa) * (t_act_k / t_std_k) * (z_act / z_std)
    v_std_base = standard_volume.to(src_std_unit)
    scale = _scale_factor(v_std_base.magnitude, compression_factor)
    actual_geom = ureg.Quantity(v_std_base.magnitude * scale, geom_unit)
    if target_unit is not None:
        return actual_geom.to(parse_unit(target_unit))
    return actual_geom

#-------------------------------------------------------------------------------
# Gas Assay
#--------------------------------------------------------------------------------

@validated_assay
@dataclass(frozen=True, slots=True)
class GasAssay:
    """Immutable thermodynamic specification for natural gas and LNG streams.

    Resolves calorific, Wobbe Index, and cryogenic LNG cargo conversion factors
    across international reference standards (FERC, Texas RRC, ISO 13443, DIN 1343, Groningen).
    """

    name: str
    gcv_btu_per_scf: float
    ncv_btu_per_scf: float
    gcv_mj_per_sm3: float
    gcv_mj_per_nm3: float
    mmbtu_per_mt_lng: float | None = None
    relative_density: float | None = None
    wobbe_index_mj_per_sm3: float | None = None
    molar_composition: Mapping[str, float] | None = None

    def __init__(
        self,
        name: str = "custom_gas",
        *,
        gcv_btu_per_scf: float | None = None,
        ncv_btu_per_scf: float | None = None,
        gcv_mj_per_sm3: float | None = None,
        gcv_mj_per_nm3: float | None = None,
        ncv_to_gcv_ratio: float = 0.9009,
        relative_density: float | None = None,
        wobbe_index: float | None = None,
        molar_composition: Mapping[str, float] | None = None,
    ) -> None:
        """Constructs an immutable gas assay from any single standard calorific benchmark."""
        provided = sum(p is not None for p in (gcv_btu_per_scf, gcv_mj_per_sm3, gcv_mj_per_nm3))
        if provided != 1:
            raise ValueError(
                "Provide exactly one calorific parameter: "
                "'gcv_btu_per_scf', 'gcv_mj_per_sm3', or 'gcv_mj_per_nm3'."
            )

        resolved_gcv_scf: float
        if gcv_btu_per_scf is not None:
            resolved_gcv_scf = float(gcv_btu_per_scf)
        elif gcv_mj_per_sm3 is not None:
            val_sm3 = float(gcv_mj_per_sm3)
            if not math.isfinite(val_sm3):
                raise ValueError(f"GCV must be finite, got {gcv_mj_per_sm3}")
            resolved_gcv_scf = val_sm3 / (MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_ISO15)
        else:
            val_nm3 = float(gcv_mj_per_nm3)  # type: ignore[arg-type]
            if not math.isfinite(val_nm3):
                raise ValueError(f"GCV must be finite, got {gcv_mj_per_nm3}")
            resolved_gcv_scf = val_nm3 / (MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_DIN0)

        if not math.isfinite(resolved_gcv_scf) or resolved_gcv_scf <= 0:
            raise ValueError(f"GCV must be finite and strictly positive, got {resolved_gcv_scf}")

        resolved_ncv_scf: float
        if ncv_btu_per_scf is not None:
            resolved_ncv_scf = float(ncv_btu_per_scf)
            if not math.isfinite(resolved_ncv_scf) or resolved_ncv_scf <= 0 or resolved_ncv_scf > resolved_gcv_scf:
                raise ValueError("NCV must be finite, strictly positive, and cannot exceed GCV.")
        else:
            if not math.isfinite(ncv_to_gcv_ratio) or not (0.5 < ncv_to_gcv_ratio <= 1.0):
                raise ValueError(f"ncv_to_gcv_ratio must be finite and in (0.5, 1.0], got {ncv_to_gcv_ratio}")
            resolved_ncv_scf = resolved_gcv_scf * ncv_to_gcv_ratio

        calc_mj_sm3 = resolved_gcv_scf * MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_ISO15
        calc_mj_nm3 = resolved_gcv_scf * MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_DIN0

        calc_wobbe: float | None = float(wobbe_index) if wobbe_index is not None else None
        calc_rel_dens: float | None = float(relative_density) if relative_density is not None else None
        if calc_rel_dens is not None:
            if not math.isfinite(calc_rel_dens) or calc_rel_dens <= 0:
                raise ValueError("Relative density must be finite and strictly positive.")
            if calc_wobbe is not None and (not math.isfinite(calc_wobbe) or calc_wobbe <= 0):
                raise ValueError("Wobbe index must be finite and strictly positive.")
            if calc_wobbe is None:
                calc_wobbe = calc_mj_sm3 / math.sqrt(calc_rel_dens)
        elif calc_wobbe is not None:
            if not math.isfinite(calc_wobbe) or calc_wobbe <= 0:
                raise ValueError("Wobbe index must be finite and strictly positive.")
            calc_rel_dens = (calc_mj_sm3 / calc_wobbe) ** 2

        clean_comp: Mapping[str, float] | None = None
        if molar_composition is not None:
            _validated_composition(molar_composition)
            if any(not math.isfinite(float(v)) or float(v) < 0.0 for v in molar_composition.values()):
                raise ValueError("Molar composition fractions must be finite and non-negative.")
            total_mol = sum(float(v) for v in molar_composition.values())
            if not (0.999 <= total_mol <= 1.001):
                raise ValueError(f"Molar composition fractions must sum to 1.0 (got {total_mol:.6f}).")
            clean_comp = MappingProxyType(
                {str(k): float(v) / total_mol for k, v in molar_composition.items() if float(v) > 0.0}
            )

        # Derive LNG gravimetric heating value (MMBtu/mt) if relative density is known
        calc_mmbtu_mt: float | None = None
        if calc_rel_dens is not None:
            gas_density_15c = calc_rel_dens * DRY_AIR_DENSITY_15C_KG_M3
            mj_per_kg = calc_mj_sm3 / gas_density_15c
            calc_mmbtu_mt = (mj_per_kg * 1000.0) / (MJ_PER_BTU_IT * 1e6)

        object.__setattr__(self, "name", str(name))
        object.__setattr__(self, "gcv_btu_per_scf", (resolved_gcv_scf))
        object.__setattr__(self, "ncv_btu_per_scf", (resolved_ncv_scf))
        object.__setattr__(self, "gcv_mj_per_sm3", (calc_mj_sm3))
        object.__setattr__(self, "gcv_mj_per_nm3", (calc_mj_nm3))
        object.__setattr__(self, "mmbtu_per_mt_lng", (calc_mmbtu_mt) if calc_mmbtu_mt is not None else None)
        object.__setattr__(self, "relative_density", (calc_rel_dens) if calc_rel_dens is not None else None)
        object.__setattr__(
            self,
            "wobbe_index_mj_per_sm3",
            (calc_wobbe) if calc_wobbe is not None else None,
        )
        object.__setattr__(self, "molar_composition", clean_comp)

    def __hash__(self) -> int:
        comp_tuple = (
            tuple(sorted(self.molar_composition.items()))
            if self.molar_composition is not None
            else None
        )
        return hash(
            (
                self.name,
                self.gcv_btu_per_scf,
                self.ncv_btu_per_scf,
                self.gcv_mj_per_sm3,
                self.gcv_mj_per_nm3,
                self.mmbtu_per_mt_lng,
                self.relative_density,
                self.wobbe_index_mj_per_sm3,
                comp_tuple,
            )
        )

    @property
    def wobbe_index_btu_per_scf(self) -> float | None:
        """US FERC Wobbe Index in Btu/scf (GCV_Btu_per_scf / sqrt(SG))."""
        if self.relative_density is None:
            return None
        return (self.gcv_btu_per_scf / math.sqrt(self.relative_density))

    @property
    def ncv_mj_per_sm3(self) -> float:
        """Net Calorific Value (LHV) in MJ/Sm³ at 15°C, 101.325 kPa."""
        return (self.gcv_mj_per_sm3 * (self.ncv_btu_per_scf / self.gcv_btu_per_scf))

    @property
    def ncv_mj_per_nm3(self) -> float:
        """Net Calorific Value (LHV) in MJ/Nm³ at 0°C, 101.325 kPa."""
        return (self.gcv_mj_per_nm3 * (self.ncv_btu_per_scf / self.gcv_btu_per_scf))

    @property
    def is_h_gas(self) -> bool | None:
        """Returns True for High-Calorific H-Gas (Wobbe >= 45.0 MJ/Sm³), False for L-Gas, or None."""
        if self.wobbe_index_mj_per_sm3 is None:
            return None
        return self.wobbe_index_mj_per_sm3 >= 45.0

    @property
    def gas_group(self) -> str | None:
        """European EN 16726 pipeline interchangeability group ('H-Gas', 'L-Gas', or None)."""
        if self.is_h_gas is None:
            return None
        return "H-Gas" if self.is_h_gas else "L-Gas"

    def _resolve_basis_ratio(self, basis: str) -> float:
        b = basis.strip().lower()
        if b not in {"gcv", "hhv", "ncv", "lhv"}:
            raise ValueError(f"Calorific basis must be 'gcv'/'hhv' or 'ncv'/'lhv', got {basis!r}")
        return (self.ncv_btu_per_scf / self.gcv_btu_per_scf) if b in {"ncv", "lhv"} else 1.0

    def as_context(self, *, basis: str = "gcv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging gas reference volumes and energy."""
        return create_gas_heating_context(self, basis=basis, name=context_name)

    def z_factor(
        self,
        pressure: pint.Quantity,
        temperature: pint.Quantity,
        *,
        method: str = "auto",
        co2_mol_frac: float = 0.0,
        h2s_mol_frac: float = 0.0,
    ) -> float:
        """Computes real-gas compressibility factor (Z) at (pressure, temperature) for this assay."""
        rel_dens = self.relative_density if self.relative_density is not None else 0.60
        return calculate_z_factor(
            pressure,
            temperature,
            relative_density=rel_dens,
            method=method,
            molar_composition=self.molar_composition,
            co2_mol_frac=co2_mol_frac,
            h2s_mol_frac=h2s_mol_frac,
        )

    def lng_liquid_to_energy(
        self,
        liquid_volume: pint.Quantity,
        *,
        liquid_density_kg_m3: float = LNG_LIQUID_DENSITY_KG_M3,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Converts cryogenic liquid LNG volume (e.g. cbm_lng) to total cargo energy (MMBtu_it)."""
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        if not math.isfinite(liquid_density_kg_m3) or liquid_density_kg_m3 <= 0:
            raise ValueError(f"liquid_density_kg_m3 must be finite and strictly positive, got {liquid_density_kg_m3}")
        ratio = self._resolve_basis_ratio(basis)
        cbm_q = liquid_volume.to(ureg.cbm_lng)
        mt_per_cbm = liquid_density_kg_m3 / 1000.0
        mmbtu_per_cbm = mt_per_cbm * self.mmbtu_per_mt_lng * ratio
        scale = _scale_factor(cbm_q.magnitude, mmbtu_per_cbm)
        return ureg.Quantity(cbm_q.magnitude * scale, ureg.MMBtu_it)

    def energy_to_lng_liquid(
        self,
        energy: pint.Quantity,
        *,
        liquid_density_kg_m3: float = LNG_LIQUID_DENSITY_KG_M3,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Converts cargo energy (e.g. MMBtu_it, MWh, GJ) to cryogenic liquid LNG volume (cbm_lng)."""
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        if not math.isfinite(liquid_density_kg_m3) or liquid_density_kg_m3 <= 0:
            raise ValueError(f"liquid_density_kg_m3 must be finite and strictly positive, got {liquid_density_kg_m3}")
        ratio = self._resolve_basis_ratio(basis)
        mmbtu_q = energy.to(ureg.MMBtu_it)
        mt_per_cbm = liquid_density_kg_m3 / 1000.0
        mmbtu_per_cbm = mt_per_cbm * self.mmbtu_per_mt_lng * ratio
        scale = _scale_factor(mmbtu_q.magnitude, mmbtu_per_cbm)
        return ureg.Quantity(mmbtu_q.magnitude / scale, ureg.cbm_lng)

    lng_volume_to_energy = lng_liquid_to_energy
    energy_to_lng_volume = energy_to_lng_liquid

    def lng_mass_to_energy(
        self,
        mass: pint.Quantity,
        *,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Converts LNG cargo mass (e.g. mt, tonne) directly to cargo energy (MMBtu_it)."""
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        ratio = self._resolve_basis_ratio(basis)
        mt_q = mass.to(ureg.metric_tonne)
        mmbtu_per_mt = self.mmbtu_per_mt_lng * ratio
        scale = _scale_factor(mt_q.magnitude, mmbtu_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.MMBtu_it)

    def energy_to_lng_mass(
        self,
        energy: pint.Quantity,
        *,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Converts cargo energy (e.g. MMBtu_it, MWh, GJ) directly to LNG cargo mass (metric_tonne)."""
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        ratio = self._resolve_basis_ratio(basis)
        mmbtu_q = energy.to(ureg.MMBtu_it)
        mmbtu_per_mt = self.mmbtu_per_mt_lng * ratio
        scale = _scale_factor(mmbtu_q.magnitude, mmbtu_per_mt)
        return ureg.Quantity(mmbtu_q.magnitude / scale, ureg.metric_tonne)

    def lng_liquid_to_gas_volume(
        self,
        liquid_volume: pint.Quantity,
        *,
        target_standard: str = "iso15",
        liquid_density_kg_m3: float = LNG_LIQUID_DENSITY_KG_M3,
    ) -> pint.Quantity:
        """Converts cryogenic liquid LNG volume (cbm_lng) to regasified pipeline gas volume."""
        if self.relative_density is None:
            raise ValueError("GasAssay requires relative_density to compute regasified gas expansion.")
        if not math.isfinite(liquid_density_kg_m3) or liquid_density_kg_m3 <= 0:
            raise ValueError(f"liquid_density_kg_m3 must be finite and strictly positive, got {liquid_density_kg_m3}")

        std_key = target_standard.strip().lower()
        gost_unit = ureg.m3_20C_101325kPa_dry
        gas_density_15c = self.relative_density * DRY_AIR_DENSITY_15C_KG_M3
        sm3_per_cbm = liquid_density_kg_m3 / gas_density_15c

        std_targets: dict[str, tuple[Unit, float]] = {
            "iso15": (ureg.Sm3_15C_101325kPa_dry, sm3_per_cbm),
            "din0": (ureg.Nm3_0C_101325kPa_dry, sm3_per_cbm * (_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_DIN0)),
            "gost20": (gost_unit, sm3_per_cbm * (_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_GOST20)),
            "ferc": (ureg.scf_60F_14_73psia, sm3_per_cbm * FT3_PER_M3 * (_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_FERC)),
            "tx": (ureg.scf_60F_14_65psia, sm3_per_cbm * FT3_PER_M3 * (_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_TX)),
        }
        if std_key not in std_targets:
            raise ValueError(f"Unknown target_standard {target_standard!r}. Choose from {tuple(std_targets)}.")

        dst_unit, factor = std_targets[std_key]
        cbm_q = liquid_volume.to(ureg.cbm_lng)
        scale = _scale_factor(cbm_q.magnitude, factor)
        return ureg.Quantity(cbm_q.magnitude * scale, dst_unit)

    def blend_with(
        self,
        other: GasAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> GasAssay:
        """Computes ideal volumetric blend of two gas streams (GCV, NCV, SG, Wobbe, and molar composition)."""
        if not isinstance(other, GasAssay):
            raise TypeError(f"Expected GasAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        w_other = 1.0 - self_volume_fraction
        blended_gcv = (self.gcv_btu_per_scf * self_volume_fraction) + (other.gcv_btu_per_scf * w_other)
        blended_ncv = (self.ncv_btu_per_scf * self_volume_fraction) + (other.ncv_btu_per_scf * w_other)

        blended_sg: float | None = None
        if self.relative_density is not None and other.relative_density is not None:
            blended_sg = (self.relative_density * self_volume_fraction) + (other.relative_density * w_other)

        blended_comp: dict[str, float] | None = None
        if self.molar_composition is not None and other.molar_composition is not None:
            all_keys = set(self.molar_composition) | set(other.molar_composition)
            blended_comp = {k: (self.molar_composition.get(k, 0.0) * self_volume_fraction) + (other.molar_composition.get(k, 0.0) * w_other) for k in all_keys}

        return GasAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            gcv_btu_per_scf=blended_gcv,
            ncv_btu_per_scf=blended_ncv,
            relative_density=blended_sg,
            molar_composition=blended_comp,
        )

    @classmethod
    def henry_hub(cls) -> GasAssay:
        return cls("Henry Hub", gcv_btu_per_scf=1025.0, relative_density=0.585)

    @classmethod
    def uk_nbp(cls) -> GasAssay:
        return cls("UK NBP", gcv_mj_per_sm3=39.50, relative_density=0.600)

    @classmethod
    def dutch_ttf(cls) -> GasAssay:
        return cls("Dutch TTF H-Gas", gcv_mj_per_nm3=40.00, relative_density=0.610)

    @classmethod
    def groningen_lgas(cls) -> GasAssay:
        return cls("Groningen L-Gas", gcv_mj_per_nm3=_GRONINGEN_MJ_PER_M3, relative_density=0.644)

    @classmethod
    def lng_lean(cls) -> GasAssay:
        return cls("Lean LNG", gcv_btu_per_scf=1035.0, relative_density=0.580)

    @classmethod
    def lng_rich(cls) -> GasAssay:
        return cls("Rich LNG", gcv_btu_per_scf=1120.0, relative_density=0.645)

    @classmethod
    def jkm_lng(cls) -> GasAssay:
        """Platts JKM (Japan/Korea Marker) rich seaborne LNG benchmark (~1,100 Btu/scf)."""
        return cls("JKM LNG", gcv_btu_per_scf=1100.0, relative_density=0.630)

    @classmethod
    def gost_russian(cls) -> GasAssay:
        """Russian/CIS export pipeline gas benchmark under GOST 2939-63 (20°C, 101.325 kPa)."""
        gost_mj_m3 = 38.20
        sm3_mj_m3 = gost_mj_m3 * (_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_GOST20)
        return cls("GOST Russian Pipeline", gcv_mj_per_sm3=sm3_mj_m3, relative_density=0.575)

    @classmethod
    def pure_methane(cls) -> GasAssay:
        return cls(
            "Pure Methane",
            gcv_btu_per_scf=1010.0,
            ncv_btu_per_scf=909.4,
            relative_density=0.5539,
            molar_composition={"CH4": 1.0},
        )


# ---------------------------------------------------------------------------
# Thermal & Metallurgical Coal Calorific Specification
# ---------------------------------------------------------------------------

# ISO 1928 / ASTM D5865 Gross-to-Net latent heat constants (kcal/kg per 1 wt%)
_ISO1928_H_KCAL_PER_PCT: Final[float] = 50.68
_ISO1928_M_KCAL_PER_PCT: Final[float] = 5.835
_KCAL_KG_TO_GJ_MT: Final[float] = 0.0041868      # 1 kcal_IT/kg = 4.1868 kJ/kg = 0.0041868 GJ/mt
_KCAL_KG_TO_BTU_LB: Final[float] = 1.8           # Exact IT definition


@validated_assay
@dataclass(frozen=True, slots=True)
class CoalAssay:
    """Immutable calorific and moisture specification for thermal and metallurgical coal.

    Resolves Net As Received (NAR), Gross As Received (GAR), Gross Air Dried (GAD),
    Dry Basis (DB), US Btu/lb, GJ/mt, MWh/mt, MMBtu/mt, and TCE/mt from any single
    calorific benchmark using ISO 1928 / ASTM D5865 moisture and hydrogen corrections.
    """

    name: str
    nar_kcal_per_kg: float
    gar_kcal_per_kg: float
    gad_kcal_per_kg: float | None
    gar_btu_per_lb: float
    gj_per_mt_nar: float
    mwh_per_mt_nar: float
    mmbtu_per_mt_nar: float
    tce_per_mt: float
    total_moisture_pct: float = 10.0
    inherent_moisture_pct: float | None = None
    hydrogen_ar_pct: float = 4.0
    ash_pct: float | None = None
    sulfur_pct: float | None = None

    def __init__(
        self,
        name: str = "custom_coal",
        *,
        nar_kcal_per_kg: float | None = None,
        gar_kcal_per_kg: float | None = None,
        gad_kcal_per_kg: float | None = None,
        gar_btu_per_lb: float | None = None,
        total_moisture_pct: float = 10.0,
        inherent_moisture_pct: float | None = None,
        hydrogen_ar_pct: float = 4.0,
        ash_pct: float | None = None,
        sulfur_pct: float | None = None,
    ) -> None:
        """Constructs an immutable coal assay from exactly one calorific benchmark."""
        provided = sum(p is not None for p in (nar_kcal_per_kg, gar_kcal_per_kg, gad_kcal_per_kg, gar_btu_per_lb))
        if provided != 1:
            raise ValueError(
                "Provide exactly one calorific parameter: "
                "'nar_kcal_per_kg', 'gar_kcal_per_kg', 'gad_kcal_per_kg', or 'gar_btu_per_lb'."
            )

        tm = float(total_moisture_pct)
        h_ar = float(hydrogen_ar_pct)
        if not math.isfinite(tm) or not (0.0 <= tm < 65.0):
            raise ValueError(f"total_moisture_pct must be in [0, 65), got {tm}")
        if not math.isfinite(h_ar) or not (0.0 <= h_ar <= 10.0):
            raise ValueError(f"hydrogen_ar_pct must be in [0, 10], got {h_ar}")

        im: float | None = float(inherent_moisture_pct) if inherent_moisture_pct is not None else None
        if im is not None:
            if not math.isfinite(im) or not (0.0 <= im <= tm):
                raise ValueError(f"inherent_moisture_pct must be in [0, total_moisture_pct], got {im}")

        if ash_pct is not None and (not math.isfinite(ash_pct) or not (0.0 <= ash_pct < 100.0)):
            raise ValueError(f"ash_pct must be in [0, 100), got {ash_pct}")
        if ash_pct is not None and (tm + float(ash_pct)) >= 100.0:
            raise ValueError(
                f"total_moisture_pct ({tm}) + ash_pct ({ash_pct}) must be strictly less than 100%."
            )
        if sulfur_pct is not None and (not math.isfinite(sulfur_pct) or not (0.0 <= sulfur_pct <= 15.0)):
            raise ValueError(f"sulfur_pct must be in [0, 15], got {sulfur_pct}")

        # ISO 1928 gross-to-net latent heat deduction (kcal/kg)
        latent_deduction = (_ISO1928_H_KCAL_PER_PCT * h_ar) + (_ISO1928_M_KCAL_PER_PCT * tm)

        resolved_gad: float | None = None
        if nar_kcal_per_kg is not None:
            resolved_nar = float(nar_kcal_per_kg)
            resolved_gar = resolved_nar + latent_deduction
        elif gar_kcal_per_kg is not None:
            resolved_gar = float(gar_kcal_per_kg)
            resolved_nar = resolved_gar - latent_deduction
        elif gar_btu_per_lb is not None:
            resolved_gar = float(gar_btu_per_lb) / _KCAL_KG_TO_BTU_LB
            resolved_nar = resolved_gar - latent_deduction
        else:
            assert gad_kcal_per_kg is not None
            if im is None:
                raise ValueError("Providing 'gad_kcal_per_kg' requires 'inherent_moisture_pct'.")
            resolved_gad = float(gad_kcal_per_kg)
            resolved_gar = resolved_gad * ((100.0 - tm) / (100.0 - im))
            resolved_nar = resolved_gar - latent_deduction

        if not math.isfinite(resolved_nar) or not math.isfinite(resolved_gar) or resolved_nar <= 0 or resolved_gar <= 0:
            raise ValueError(
                f"Resolved calorific values must be finite and positive (NAR={resolved_nar}, GAR={resolved_gar})."
            )

        if resolved_gad is None and im is not None:
            resolved_gad = resolved_gar * ((100.0 - im) / (100.0 - tm))

        gj_mt = resolved_nar * _KCAL_KG_TO_GJ_MT
        mwh_mt = gj_mt / 3.6
        # 1 MMBtu_IT = 1.05505585262 GJ = MJ_PER_BTU_IT * 1000 GJ
        mmbtu_mt = gj_mt / (MJ_PER_BTU_IT * 1000.0)
        tce_mt = resolved_nar / 7000.0

        object.__setattr__(self, "name", str(name))
        object.__setattr__(self, "nar_kcal_per_kg", (resolved_nar))
        object.__setattr__(self, "gar_kcal_per_kg", (resolved_gar))
        object.__setattr__(self, "gad_kcal_per_kg", (resolved_gad) if resolved_gad is not None else None)
        object.__setattr__(self, "gar_btu_per_lb", (resolved_gar * _KCAL_KG_TO_BTU_LB))
        object.__setattr__(self, "gj_per_mt_nar", (gj_mt))
        object.__setattr__(self, "mwh_per_mt_nar", (mwh_mt))
        object.__setattr__(self, "mmbtu_per_mt_nar", (mmbtu_mt))
        object.__setattr__(self, "tce_per_mt", (tce_mt))
        object.__setattr__(self, "total_moisture_pct", tm)
        object.__setattr__(self, "inherent_moisture_pct", im)
        object.__setattr__(self, "hydrogen_ar_pct", h_ar)
        object.__setattr__(self, "ash_pct", float(ash_pct) if ash_pct is not None else None)
        object.__setattr__(self, "sulfur_pct", float(sulfur_pct) if sulfur_pct is not None else None)

    @property
    def dry_basis_kcal_per_kg(self) -> float:
        """Gross Calorific Value on a Dry Basis (DB): GAR / (1 - TM/100)."""
        return (self.gar_kcal_per_kg * (100.0 / (100.0 - self.total_moisture_pct)))

    @property
    def nar_btu_per_lb(self) -> float:
        """Net Calorific Value (NAR) in Btu/lb."""
        return (self.nar_kcal_per_kg * _KCAL_KG_TO_BTU_LB)

    @property
    def gj_per_mt_gar(self) -> float:
        """Gross Calorific Value (GAR) in GJ/mt."""
        return (self.gar_kcal_per_kg * _KCAL_KG_TO_GJ_MT)

    @property
    def mwh_per_mt_gar(self) -> float:
        """Gross Calorific Value (GAR) in MWh/mt."""
        return (self.gj_per_mt_gar / 3.6)

    @property
    def mmbtu_per_mt_gar(self) -> float:
        """Gross Calorific Value (GAR) in MMBtu_IT/mt."""
        return (self.gj_per_mt_gar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def daf_kcal_per_kg(self) -> float | None:
        """Gross Calorific Value on a Dry Ash-Free (DAF) basis, or None if ash_pct is unset."""
        if self.ash_pct is None:
            return None
        return (self.gar_kcal_per_kg * (100.0 / (100.0 - self.total_moisture_pct - self.ash_pct)))

    def wmt_to_dmt(self, wet_metric_tonnes: float | Decimal) -> float:
        """Converts Wet Metric Tonnes (WMT) to Dry Metric Tonnes (DMT) using total moisture."""
        wmt = float(wet_metric_tonnes)
        if not math.isfinite(wmt):
            raise ValueError(f"wet_metric_tonnes must be finite, got {wet_metric_tonnes}")
        return wmt * (1.0 - self.total_moisture_pct / 100.0)

    def dmt_to_wmt(self, dry_metric_tonnes: float | Decimal) -> float:
        """Converts Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT) using total moisture."""
        dmt = float(dry_metric_tonnes)
        if not math.isfinite(dmt):
            raise ValueError(f"dry_metric_tonnes must be finite, got {dry_metric_tonnes}")
        return dmt / (1.0 - self.total_moisture_pct / 100.0)

    def as_context(self, *, basis: str = "nar", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [mass] and [energy] for this coal."""
        return create_coal_heating_context(self, basis=basis, name=context_name)

    def pro_rata_price(
        self,
        benchmark_price: float,
        benchmark_kcal_per_kg: float = 6000.0,
        *,
        basis: str = "nar",
    ) -> float:
        """Computes SCoTA / globalCOAL / ICI pro-rata calorific settlement price."""
        price = float(benchmark_price)
        bench_cv = float(benchmark_kcal_per_kg)
        if not math.isfinite(price) or not math.isfinite(bench_cv) or bench_cv <= 0.0:
            raise ValueError("benchmark_price must be finite and benchmark_kcal_per_kg strictly positive.")
        basis_key = basis.strip().lower()
        if basis_key == "nar":
            actual_cv = self.nar_kcal_per_kg
        elif basis_key == "gar":
            actual_cv = self.gar_kcal_per_kg
        else:
            raise ValueError(f"basis must be 'nar' or 'gar', got {basis!r}")
        return (price * (actual_cv / bench_cv))

    def blend_with(
        self,
        other: CoalAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> CoalAssay:
        """Computes mass-weighted coal blend assay (NAR, moisture, H, ash, sulfur)."""
        if not isinstance(other, CoalAssay):
            raise TypeError(f"Expected CoalAssay, got {type(other).__name__}")
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1

        def _mix(a: float | None, b: float | None) -> float | None:
            if a is None or b is None:
                return None
            return (a * w1) + (b * w2)

        mixed_tm = (self.total_moisture_pct * w1) + (other.total_moisture_pct * w2)
        mixed_im = _mix(self.inherent_moisture_pct, other.inherent_moisture_pct)
        if mixed_im is not None:
            mixed_im = min(mixed_im, mixed_tm)

        return CoalAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            nar_kcal_per_kg=(self.nar_kcal_per_kg * w1) + (other.nar_kcal_per_kg * w2),
            total_moisture_pct=mixed_tm,
            inherent_moisture_pct=mixed_im,
            hydrogen_ar_pct=(self.hydrogen_ar_pct * w1) + (other.hydrogen_ar_pct * w2),
            ash_pct=_mix(self.ash_pct, other.ash_pct),
            sulfur_pct=_mix(self.sulfur_pct, other.sulfur_pct),
        )

    # --- Global Thermal Coal Benchmark Presets ---
    @classmethod
    def newcastle_6000(cls) -> CoalAssay:
        """gC NEWC / ICE Newcastle FOB 6,000 kcal/kg NAR benchmark."""
        return cls("Newcastle 6000 NAR", nar_kcal_per_kg=6000.0, total_moisture_pct=10.0, ash_pct=14.0, sulfur_pct=0.75)

    @classmethod
    def api2_cif_ara(cls) -> CoalAssay:
        """Argus/McCloskey API 2 CIF ARA 6,000 kcal/kg NAR European benchmark."""
        return cls("API 2 CIF ARA", nar_kcal_per_kg=6000.0, total_moisture_pct=11.5, ash_pct=14.0, sulfur_pct=0.85)

    @classmethod
    def api4_fob_rbct(cls) -> CoalAssay:
        """Argus/McCloskey API 4 FOB Richards Bay 6,000 kcal/kg NAR benchmark."""
        return cls("API 4 FOB RBCT", nar_kcal_per_kg=6000.0, total_moisture_pct=10.0, ash_pct=15.0, sulfur_pct=0.70)

    @classmethod
    def indonesian_4200(cls) -> CoalAssay:
        """ICI 4 Indonesian sub-bituminous 4,200 kcal/kg GAR thermal coal."""
        return cls("Indonesian 4200 GAR", gar_kcal_per_kg=4200.0, total_moisture_pct=34.0, hydrogen_ar_pct=3.2, ash_pct=6.0, sulfur_pct=0.40)

    @classmethod
    def prb_8800(cls) -> CoalAssay:
        """US Powder River Basin 8,800 Btu/lb GAR sub-bituminous benchmark."""
        return cls("US PRB 8800", gar_btu_per_lb=8800.0, total_moisture_pct=27.0, hydrogen_ar_pct=3.5, ash_pct=5.5, sulfur_pct=0.35)

    @classmethod
    def capp_12500(cls) -> CoalAssay:
        """US Central Appalachian 12,500 Btu/lb GAR bituminous benchmark."""
        return cls("US CAPP 12500", gar_btu_per_lb=12500.0, total_moisture_pct=6.5, hydrogen_ar_pct=4.5, ash_pct=11.0, sulfur_pct=1.00)

    @classmethod
    def plv_hard_coking(cls) -> CoalAssay:
        """Platts Premium Low-Vol (PLV) Hard Coking Coal FOB Australia benchmark."""
        return cls(
            "PLV Hard Coking",
            nar_kcal_per_kg=7100.0,
            total_moisture_pct=10.0,
            hydrogen_ar_pct=4.3,
            ash_pct=9.5,
            sulfur_pct=0.60,
        )


# ---------------------------------------------------------------------------
# Clean Fuel Carriers (Hydrogen, Ammonia & Methanol LHV / HHV Specification)
# ---------------------------------------------------------------------------

_LH2_LIQUID_DENSITY_KG_M3: Final[float] = 70.85       # Liquid H2 at -252.87°C, 1 atm
_NH3_LIQUID_DENSITY_KG_M3: Final[float] = 682.0       # Refrigerated liquid NH3 at -33.3°C, 1 atm
_NH3_GAS_DENSITY_0C_KG_NM3: Final[float] = 0.7714     # Gaseous NH3 at 0°C, 101.325 kPa (DIN 1343)
_MEOH_LHV_MJ_PER_KG: Final[float] = 19.90             # Pure Methanol (CH3OH) LHV
_MEOH_HHV_MJ_PER_KG: Final[float] = 22.70             # Pure Methanol (CH3OH) HHV
_MEOH_LIQUID_DENSITY_KG_M3: Final[float] = 792.0      # Liquid Methanol at 15°C, 1 atm


@validated_assay
@dataclass(frozen=True, slots=True)
class CleanFuelAssay:
    """Immutable thermodynamic specification for zero/low-carbon energy carriers.

    Supports gravimetric (LHV/HHV), gaseous volumetric (Nm3), cryogenic/liquid
    shipping parity (m3), and hydrogen-carrier cracking equivalency for H2, NH3,
    and Methanol.
    """

    name: str
    lhv_mj_per_kg: float
    hhv_mj_per_kg: float
    density_kg_per_nm3: float | None = None
    liquid_density_kg_per_m3: float | None = None
    hydrogen_mass_fraction: float | None = None
    carbon_intensity_kg_co2e_per_kg: float | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.lhv_mj_per_kg) or self.lhv_mj_per_kg <= 0.0:
            raise ValueError(f"lhv_mj_per_kg must be finite and strictly positive, got {self.lhv_mj_per_kg}")
        if not math.isfinite(self.hhv_mj_per_kg) or self.hhv_mj_per_kg < self.lhv_mj_per_kg:
            raise ValueError(
                f"hhv_mj_per_kg ({self.hhv_mj_per_kg}) must be finite and >= lhv_mj_per_kg ({self.lhv_mj_per_kg})."
            )
        if self.density_kg_per_nm3 is not None:
            if not math.isfinite(self.density_kg_per_nm3) or self.density_kg_per_nm3 <= 0.0:
                raise ValueError(f"density_kg_per_nm3 must be strictly positive, got {self.density_kg_per_nm3}")
        if self.liquid_density_kg_per_m3 is not None:
            if not math.isfinite(self.liquid_density_kg_per_m3) or self.liquid_density_kg_per_m3 <= 0.0:
                raise ValueError(
                    f"liquid_density_kg_per_m3 must be strictly positive, got {self.liquid_density_kg_per_m3}"
                )
        if self.hydrogen_mass_fraction is not None:
            if not math.isfinite(self.hydrogen_mass_fraction) or not (0.0 < self.hydrogen_mass_fraction <= 1.0):
                raise ValueError(
                    f"hydrogen_mass_fraction must be in (0.0, 1.0], got {self.hydrogen_mass_fraction}"
                )
        if self.carbon_intensity_kg_co2e_per_kg is not None:
            if not math.isfinite(self.carbon_intensity_kg_co2e_per_kg) or self.carbon_intensity_kg_co2e_per_kg < 0.0:
                raise ValueError(
                    f"carbon_intensity_kg_co2e_per_kg cannot be negative, got {self.carbon_intensity_kg_co2e_per_kg}"
                )

    # --- Gravimetric & Volumetric Energy Density Properties ---
    @property
    def lhv_kwh_per_kg(self) -> float:
        return (self.lhv_mj_per_kg / 3.6)

    @property
    def hhv_kwh_per_kg(self) -> float:
        return (self.hhv_mj_per_kg / 3.6)

    @property
    def lhv_mwh_per_mt(self) -> float:
        return (self.lhv_mj_per_kg / 3.6)

    @property
    def hhv_mwh_per_mt(self) -> float:
        return (self.hhv_mj_per_kg / 3.6)

    @property
    def lhv_mmbtu_per_mt(self) -> float:
        return (self.lhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mmbtu_per_mt(self) -> float:
        return (self.hhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    @property
    def lhv_mj_per_nm3(self) -> float | None:
        if self.density_kg_per_nm3 is None:
            return None
        return (self.lhv_mj_per_kg * self.density_kg_per_nm3)

    @property
    def hhv_mj_per_nm3(self) -> float | None:
        if self.density_kg_per_nm3 is None:
            return None
        return (self.hhv_mj_per_kg * self.density_kg_per_nm3)

    @property
    def lhv_gj_per_m3_liquid(self) -> float | None:
        if self.liquid_density_kg_per_m3 is None:
            return None
        return ((self.lhv_mj_per_kg * self.liquid_density_kg_per_m3) / 1000.0)

    @property
    def carbon_intensity_g_co2e_per_mj_lhv(self) -> float | None:
        """Well-to-gate carbon intensity in gCO2e/MJ (LHV) for EU RFNBO / FuelEU Maritime."""
        if self.carbon_intensity_kg_co2e_per_kg is None:
            return None
        # (kg CO2e / kg fuel) * 1000 g/kg / (MJ / kg fuel) = gCO2e / MJ
        return ((self.carbon_intensity_kg_co2e_per_kg * 1000.0) / self.lhv_mj_per_kg)

    def _heating_value_mj_per_kg(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"lhv", "ncv"}:
            return self.lhv_mj_per_kg
        if key in {"hhv", "gcv"}:
            return self.hhv_mj_per_kg
        raise ValueError(f"basis must be 'lhv'/'ncv' or 'hhv'/'gcv', got {basis!r}")

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [mass] and [energy] for this carrier."""
        return create_clean_fuel_context(self, basis=basis, name=context_name)

    # --- Mass <-> Energy & Liquid Shipping Conversions ---
    def mass_to_energy(self, mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts carrier mass (e.g. kg, mt) to total energy (MWh)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        hv_mwh_per_mt = self._heating_value_mj_per_kg(basis) / 3.6
        mt_q = mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, hv_mwh_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.MWh)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. MWh, GJ, MMBtu) to required carrier mass (metric_ton)."""
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        hv_mwh_per_mt = self._heating_value_mj_per_kg(basis) / 3.6
        mwh_q = energy.to(ureg.MWh)
        scale = _scale_factor(mwh_q.magnitude, hv_mwh_per_mt)
        return ureg.Quantity(mwh_q.magnitude / scale, ureg.metric_ton)

    def gas_volume_to_mass(self, gas_volume: pint.Quantity) -> pint.Quantity:
        """Converts normal gaseous volume (DIN 1343 Nm3) to mass (kg)."""
        if not isinstance(gas_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(gas_volume).__name__}")
        if self.density_kg_per_nm3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_nm3.")
        nm3_q = gas_volume.to(ureg.Nm3_0C_101325kPa_dry)
        scale = _scale_factor(nm3_q.magnitude, self.density_kg_per_nm3)
        return ureg.Quantity(nm3_q.magnitude * scale, ureg.kilogram)

    def mass_to_gas_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts carrier mass (e.g. kg, mt) to normal gaseous volume (DIN 1343 Nm3)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.density_kg_per_nm3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_nm3.")
        kg_q = mass.to(ureg.kilogram)
        scale = _scale_factor(kg_q.magnitude, self.density_kg_per_nm3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.Nm3_0C_101325kPa_dry)

    def liquid_volume_to_mass(self, liquid_volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid/cryogenic cargo volume (e.g. m^3, L) to mass (metric_ton)."""
        if not isinstance(liquid_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(liquid_volume).__name__}")
        if self.liquid_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define liquid_density_kg_per_m3.")
        m3_q = liquid_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.liquid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_liquid_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts carrier mass (e.g. mt, kg) to liquid/cryogenic cargo volume (m^3)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.liquid_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define liquid_density_kg_per_m3.")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.liquid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def liquid_volume_to_energy(self, liquid_volume: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts liquid/cryogenic cargo volume (e.g. m^3) directly to energy (MWh)."""
        return self.mass_to_energy(self.liquid_volume_to_mass(liquid_volume), basis=basis)

    def energy_to_liquid_volume(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. MWh, GJ, MMBtu) directly to liquid/cryogenic cargo volume (m^3)."""
        return self.mass_to_liquid_volume(self.energy_to_mass(energy, basis=basis))

    # --- Hydrogen-Carrier Equivalency (Ammonia / Methanol Cracking) ---
    def contained_hydrogen_mass(
        self,
        carrier_mass: pint.Quantity,
        *,
        cracking_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes recoverable H2 mass from a carrier cargo given cracking efficiency."""
        if not isinstance(carrier_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(carrier_mass).__name__}")
        if not validate_dimensionality(carrier_mass, "kilogram"):
            raise pint.DimensionalityError(
                carrier_mass.units,
                ureg.kilogram,
                carrier_mass.dimensionality,
                ureg.kilogram.dimensionality,
            )
        if self.hydrogen_mass_fraction is None:
            raise ValueError(f"{self.name!r} does not define hydrogen_mass_fraction.")
        if not math.isfinite(cracking_efficiency) or not (0.0 < cracking_efficiency <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0.0, 1.0], got {cracking_efficiency}")
        factor = self.hydrogen_mass_fraction * cracking_efficiency
        scale = _scale_factor(carrier_mass.magnitude, factor)
        return ureg.Quantity(carrier_mass.magnitude * scale, carrier_mass.units)

    def carrier_mass_for_hydrogen(
        self,
        h2_mass: pint.Quantity,
        *,
        cracking_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes required carrier mass to deliver a target H2 mass after cracking."""
        if not isinstance(h2_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(h2_mass).__name__}")
        if not validate_dimensionality(h2_mass, "kilogram"):
            raise pint.DimensionalityError(
                h2_mass.units,
                ureg.kilogram,
                h2_mass.dimensionality,
                ureg.kilogram.dimensionality,
            )
        if self.hydrogen_mass_fraction is None:
            raise ValueError(f"{self.name!r} does not define hydrogen_mass_fraction.")
        if not math.isfinite(cracking_efficiency) or not (0.0 < cracking_efficiency <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0.0, 1.0], got {cracking_efficiency}")
        factor = self.hydrogen_mass_fraction * cracking_efficiency
        scale = _scale_factor(h2_mass.magnitude, factor)
        return ureg.Quantity(h2_mass.magnitude / scale, h2_mass.units)

    # --- Clean Fuel Carrier Presets ---
    @classmethod
    def hydrogen(cls) -> CleanFuelAssay:
        """Pure Hydrogen (H2) ISO 6976 / NIST specification."""
        return cls(
            name="Hydrogen (H2)",
            lhv_mj_per_kg=H2_LHV_MJ_PER_KG,
            hhv_mj_per_kg=H2_HHV_MJ_PER_KG,
            density_kg_per_nm3=H2_DENSITY_0C_KG_NM3,
            liquid_density_kg_per_m3=_LH2_LIQUID_DENSITY_KG_M3,
            hydrogen_mass_fraction=1.0,
        )

    @classmethod
    def ammonia(cls) -> CleanFuelAssay:
        """Anhydrous Ammonia (NH3) maritime & power carrier specification."""
        return cls(
            name="Ammonia (NH3)",
            lhv_mj_per_kg=NH3_LHV_MJ_PER_KG,
            hhv_mj_per_kg=NH3_HHV_MJ_PER_KG,
            density_kg_per_nm3=_NH3_GAS_DENSITY_0C_KG_NM3,
            liquid_density_kg_per_m3=_NH3_LIQUID_DENSITY_KG_M3,
            hydrogen_mass_fraction=0.1776,
        )

    @classmethod
    def methanol(cls) -> CleanFuelAssay:
        """Pure Methanol (CH3OH) maritime & chemical carrier specification."""
        return cls(name="Methanol (CH3OH)", lhv_mj_per_kg=_MEOH_LHV_MJ_PER_KG, hhv_mj_per_kg=_MEOH_HHV_MJ_PER_KG, density_kg_per_nm3=None, liquid_density_kg_per_m3=_MEOH_LIQUID_DENSITY_KG_M3,
            hydrogen_mass_fraction=0.1258,)



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


@validated_assay
@dataclass(frozen=True, slots=True)
class BiofuelAssay:
    """Immutable physical and regulatory specification for liquid biofuels and SAF.

    Bridges volumetric (gal, m^3, bbl), gravimetric (mt, kg), and energy (MJ, GJ,
    MWh, MMBtu) conversions alongside US EPA RIN equivalence values and LCFS/RED II
    carbon-intensity (gCO2e/MJ) credit calculations.
    """

    name: str
    density_kg_per_m3: float
    lhv_mj_per_kg: float
    hhv_mj_per_kg: float
    epa_rin_ev: float | None = None
    rin_d_code: str | None = None
    carbon_intensity_g_co2e_per_mj: float | None = None
    thermal_expansion_coeff_per_c: float | None = None
    red_energy_multiplier: float = 1.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.density_kg_per_m3) or not (500.0 <= self.density_kg_per_m3 <= 1200.0):
            raise ValueError(f"density_kg_per_m3 must be finite and in [500, 1200] kg/m^3, got {self.density_kg_per_m3}")
        if not math.isfinite(self.lhv_mj_per_kg) or self.lhv_mj_per_kg <= 0.0:
            raise ValueError(f"lhv_mj_per_kg must be finite and strictly positive, got {self.lhv_mj_per_kg}")
        if not math.isfinite(self.hhv_mj_per_kg) or self.hhv_mj_per_kg < self.lhv_mj_per_kg:
            raise ValueError(f"hhv_mj_per_kg ({self.hhv_mj_per_kg}) must be finite and >= lhv_mj_per_kg ({self.lhv_mj_per_kg}).")
        if self.epa_rin_ev is not None:
            if not math.isfinite(self.epa_rin_ev) or not (0.0 < self.epa_rin_ev <= 3.0):
                raise ValueError(f"epa_rin_ev must be in (0.0, 3.0], got {self.epa_rin_ev}")
        if self.carbon_intensity_g_co2e_per_mj is not None:
            if not math.isfinite(self.carbon_intensity_g_co2e_per_mj):
                raise ValueError(
                    f"carbon_intensity_g_co2e_per_mj must be finite, got {self.carbon_intensity_g_co2e_per_mj}"
                )
        if self.thermal_expansion_coeff_per_c is not None:
            if not math.isfinite(self.thermal_expansion_coeff_per_c) or self.thermal_expansion_coeff_per_c <= 0.0:
                raise ValueError(
                    f"thermal_expansion_coeff_per_c must be finite and strictly positive, got {self.thermal_expansion_coeff_per_c}"
                )
        if not math.isfinite(self.red_energy_multiplier) or self.red_energy_multiplier < 1.0:
            raise ValueError(f"red_energy_multiplier must be finite and >= 1.0, got {self.red_energy_multiplier}")

    # --- Volumetric & Gravimetric Parity Properties ---
    @property
    def specific_gravity(self) -> float:
        """Specific gravity relative to water at 15°C (999.102 kg/m^3)."""
        return (self.density_kg_per_m3 / 999.102)
    @property
    def bbl_per_mt(self) -> float:
        """Barrels per metric tonne at 15°C / 60°F."""
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_BBL / 1000.0)))

    @property
    def gal_per_mt(self) -> float:
        """US gallons per metric tonne at 15°C / 60°F."""
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def lhv_mj_per_l(self) -> float:
        """EU RED II volumetric Lower Heating Value (MJ/L)."""
        return ((self.lhv_mj_per_kg * self.density_kg_per_m3) / 1000.0)

    @property
    def hhv_mj_per_l(self) -> float:
        """Volumetric Higher Heating Value (MJ/L)."""
        return ((self.hhv_mj_per_kg * self.density_kg_per_m3) / 1000.0)

    @property
    def lhv_btu_per_gal(self) -> float:
        """US volumetric Lower Heating Value (Btu_IT/gal)."""
        return ((self.lhv_mj_per_l * _LITRES_PER_US_GAL) / MJ_PER_BTU_IT)

    @property
    def lhv_mwh_per_mt(self) -> float:
        return (self.lhv_mj_per_kg / 3.6)

    @property
    def lhv_mmbtu_per_mt(self) -> float:
        return (self.lhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mwh_per_mt(self) -> float:
        return (self.hhv_mj_per_kg / 3.6)

    @property
    def hhv_mmbtu_per_mt(self) -> float:
        return (self.hhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float = 15.0,
    ) -> float:
        """Computes thermal Volume Correction Factor (VCF) from observed temperature to reference.

        Accepts observed_temp as float (°C) or pint.Quantity (e.g. Q_(68, 'degF'), Q_(20, 'degC')).
        Use ref_temp_c=15.0 for ISO/ARA, 15.5556 (60°F) for US, or 20.0 for Brazilian ANP ethanol.
        """
        if self.thermal_expansion_coeff_per_c is None:
            raise ValueError(f"{self.name!r} does not define thermal_expansion_coeff_per_c.")
        if isinstance(observed_temp, pint.Quantity):
            t_obs_c = float(observed_temp.to(ureg.degC).magnitude)
        else:
            t_obs_c = float(observed_temp)
        if not math.isfinite(t_obs_c) or not math.isfinite(ref_temp_c):
            raise ValueError("Temperatures must be finite.")

        alpha = self.thermal_expansion_coeff_per_c
        dt = t_obs_c - ref_temp_c
        vcf = math.exp(-alpha * dt * (1.0 + 0.8 * alpha * dt))
        return (vcf)

    def ghg_savings_pct(self, fossil_comparator_g_co2e_per_mj: float = 94.0) -> float:
        """Computes EU RED II/III (94.0) or ICAO CORSIA (89.0) GHG emissions savings (%)."""
        if self.carbon_intensity_g_co2e_per_mj is None:
            raise ValueError(f"{self.name!r} does not define carbon_intensity_g_co2e_per_mj.")
        if not math.isfinite(fossil_comparator_g_co2e_per_mj) or fossil_comparator_g_co2e_per_mj <= 0.0:
            raise ValueError("fossil_comparator_g_co2e_per_mj must be strictly positive.")
        return ((1.0 - (self.carbon_intensity_g_co2e_per_mj / fossil_comparator_g_co2e_per_mj)) * 100.0)

    def _heating_value_mj_per_kg(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"lhv", "ncv"}:
            return self.lhv_mj_per_kg
        if key in {"hhv", "gcv"}:
            return self.hhv_mj_per_kg
        raise ValueError(f"basis must be 'lhv'/'ncv' or 'hhv'/'gcv', got {basis!r}")

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [energy]."""
        return create_biofuel_context(self, basis=basis, name=context_name)

    # --- Physical Volume <-> Mass <-> Energy Conversions ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid volume (e.g. gal, bbl, m^3) to mass (metric_ton)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts mass (e.g. mt, kg) to liquid volume (m^3)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def mass_to_energy(self, mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts mass (e.g. mt, kg) to total energy (GJ)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        # 1 MJ/kg = 1 GJ/mt
        hv_gj_per_mt = self._heating_value_mj_per_kg(basis)
        mt_q = mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, hv_gj_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.GJ)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. GJ, MWh, MMBtu) to mass (metric_ton)."""
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        hv_gj_per_mt = self._heating_value_mj_per_kg(basis)
        gj_q = energy.to(ureg.GJ)
        scale = _scale_factor(gj_q.magnitude, hv_gj_per_mt)
        return ureg.Quantity(gj_q.magnitude / scale, ureg.metric_ton)

    def volume_to_energy(self, volume: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts liquid volume (e.g. gal, m^3, bbl) directly to energy (GJ)."""
        return self.mass_to_energy(self.volume_to_mass(volume), basis=basis)

    def energy_to_volume(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. GJ, MWh, MMBtu) directly to liquid volume (m^3)."""
        return self.mass_to_volume(self.energy_to_mass(energy, basis=basis))

    # --- Regulatory Credit Calculations (EPA RFS RINs & LCFS / CORSIA) ---
    def rins_generated(self, volume: pint.Quantity) -> float:
        """Computes US EPA RFS Renewable Identification Numbers (RINs) generated by a fuel volume."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.epa_rin_ev is None:
            raise ValueError(f"{self.name!r} does not define an epa_rin_ev equivalence value.")
        gal = float(volume.to(ureg.gallon).magnitude)
        return (gal * self.epa_rin_ev)

    def lcfs_credits_mt(
        self,
        volume: pint.Quantity,
        benchmark_ci_g_co2e_per_mj: float,
        *,
        eer: float = 1.0,
    ) -> float:
        """Computes LCFS / Clean Fuel Standard credits (or deficits if negative) in mt CO2e.

        Formula:
            Credits (mt CO2e) = (CI_benchmark * EER - CI_fuel) * Energy_LHV_MJ / 1e6
        """
        if self.carbon_intensity_g_co2e_per_mj is None:
            raise ValueError(f"{self.name!r} does not define carbon_intensity_g_co2e_per_mj.")
        if not math.isfinite(benchmark_ci_g_co2e_per_mj) or benchmark_ci_g_co2e_per_mj < 0.0:
            raise ValueError(f"benchmark_ci_g_co2e_per_mj must be non-negative, got {benchmark_ci_g_co2e_per_mj}")
        if not math.isfinite(eer) or eer <= 0.0:
            raise ValueError(f"eer (Energy Economy Ratio) must be strictly positive, got {eer}")

        energy_mj = float(self.volume_to_energy(volume, basis="lhv").to(ureg.MJ).magnitude)
        ci_delta = (benchmark_ci_g_co2e_per_mj * eer) - self.carbon_intensity_g_co2e_per_mj
        return ((ci_delta * energy_mj) / 1e6)

    def red_compliance_energy(self, volume: pint.Quantity) -> pint.Quantity:
        """Computes EU RED II/III / UK RTFO compliance energy (GJ LHV) including double-counting."""
        phys_gj = self.volume_to_energy(volume, basis="lhv")
        scale = _scale_factor(phys_gj.magnitude, self.red_energy_multiplier)
        return ureg.Quantity(phys_gj.magnitude * scale, phys_gj.units)

    def blend_with(
        self,
        other: BiofuelAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> BiofuelAssay:
        """Computes a physically consistent liquid biofuel blend.

        - Density & RIN EV: volume-weighted
        - LHV & HHV (MJ/kg): mass-weighted
        - Carbon Intensity (gCO2e/MJ): LHV-energy-weighted
        """
        if not isinstance(other, BiofuelAssay):
            raise TypeError(f"Expected BiofuelAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        v1 = float(self_volume_fraction)
        v2 = 1.0 - v1

        blended_density = (self.density_kg_per_m3 * v1) + (other.density_kg_per_m3 * v2)
        m1 = self.density_kg_per_m3 * v1
        m2 = other.density_kg_per_m3 * v2
        total_m = m1 + m2
        w1, w2 = m1 / total_m, m2 / total_m

        blended_lhv = (self.lhv_mj_per_kg * w1) + (other.lhv_mj_per_kg * w2)
        blended_hhv = (self.hhv_mj_per_kg * w1) + (other.hhv_mj_per_kg * w2)

        blended_rin: float | None = None
        if self.epa_rin_ev is not None and other.epa_rin_ev is not None:
            blended_rin = ((self.epa_rin_ev * v1) + (other.epa_rin_ev * v2))

        e1 = m1 * self.lhv_mj_per_kg
        e2 = m2 * other.lhv_mj_per_kg

        blended_ci: float | None = None
        if self.carbon_intensity_g_co2e_per_mj is not None and other.carbon_intensity_g_co2e_per_mj is not None:
            blended_ci = ((self.carbon_intensity_g_co2e_per_mj * e1 + other.carbon_intensity_g_co2e_per_mj * e2) / (e1 + e2))

        blended_alpha: float | None = None
        if self.thermal_expansion_coeff_per_c is not None and other.thermal_expansion_coeff_per_c is not None:
            blended_alpha = ((self.thermal_expansion_coeff_per_c * v1) + (other.thermal_expansion_coeff_per_c * v2))

        blended_red_mult = ((self.red_energy_multiplier * e1 + other.red_energy_multiplier * e2) / (e1 + e2))

        return BiofuelAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            density_kg_per_m3=(blended_density),
            lhv_mj_per_kg=(blended_lhv),
            hhv_mj_per_kg=(blended_hhv),
            epa_rin_ev=blended_rin,
            rin_d_code=self.rin_d_code if self.rin_d_code == other.rin_d_code else None,
            carbon_intensity_g_co2e_per_mj=blended_ci,
            thermal_expansion_coeff_per_c=blended_alpha,
            red_energy_multiplier=blended_red_mult,
        )

    # --- Global Biofuel & SAF Benchmark Presets ---
    @classmethod
    def ethanol(cls) -> BiofuelAssay:
        """Anhydrous Fuel Ethanol (CBOT / Argus T2 ARA / US RFS D6)."""
        return cls(
            name="Fuel Ethanol",
            density_kg_per_m3=789.0,
            lhv_mj_per_kg=26.80,
            hhv_mj_per_kg=29.70,
            epa_rin_ev=1.0,
            rin_d_code="D6",
            thermal_expansion_coeff_per_c=0.00109,
        )

    @classmethod
    def hydrous_ethanol(cls) -> BiofuelAssay:
        """Brazilian B3 / ANP Hydrous Ethanol (E100, ~95.5% v/v ethanol)."""
        return cls(
            name="Hydrous Ethanol (E100)",
            density_kg_per_m3=809.0,
            lhv_mj_per_kg=25.10,
            hhv_mj_per_kg=27.80,
            thermal_expansion_coeff_per_c=0.00105,
        )

    @classmethod
    def fame_biodiesel(cls) -> BiofuelAssay:
        """FAME Biodiesel / B100 (ARA FAME 0 / US RFS D4)."""
        return cls(
            name="FAME Biodiesel (B100)",
            density_kg_per_m3=883.0,
            lhv_mj_per_kg=37.20,
            hhv_mj_per_kg=40.20,
            epa_rin_ev=1.5,
            rin_d_code="D4",
            thermal_expansion_coeff_per_c=0.00083,
        )

    @classmethod
    def ucome(cls) -> BiofuelAssay:
        """Used Cooking Oil Methyl Ester (ARA UCOME, EU RED II Annex IX-B double-counting)."""
        return cls(
            name="UCOME Biodiesel",
            density_kg_per_m3=883.0,
            lhv_mj_per_kg=37.20,
            hhv_mj_per_kg=40.20,
            epa_rin_ev=1.5,
            rin_d_code="D4",
            carbon_intensity_g_co2e_per_mj=14.0,
            thermal_expansion_coeff_per_c=0.00083,
            red_energy_multiplier=2.0,
        )

    @classmethod
    def hvo_renewable_diesel(cls) -> BiofuelAssay:
        """Hydrotreated Vegetable Oil / Renewable Diesel (HVO / RD100, US RFS D4)."""
        return cls(
            name="HVO / Renewable Diesel",
            density_kg_per_m3=780.0,
            lhv_mj_per_kg=44.00,
            hhv_mj_per_kg=47.00,
            epa_rin_ev=1.7,
            rin_d_code="D4",
            thermal_expansion_coeff_per_c=0.00085,
        )

    @classmethod
    def saf_hefa(cls) -> BiofuelAssay:
        """HEFA-SPK Sustainable Aviation Fuel (ASTM D7566 / CORSIA / US RFS D4)."""
        return cls(
            name="SAF (HEFA-SPK)",
            density_kg_per_m3=760.0,
            lhv_mj_per_kg=44.10,
            hhv_mj_per_kg=47.20,
            epa_rin_ev=1.6,
            rin_d_code="D4",
            thermal_expansion_coeff_per_c=0.00095,
        )

# ---------------------------------------------------------------------------
# Natural Gas Liquids (NGLs) & LPG Specification (GPA 2145-16)
# (Ethane, Propane, Isobutane, n-Butane, Natural Gasoline / C5+)
# ---------------------------------------------------------------------------

_WATER_DENSITY_60F_KG_M3: Final[float] = 999.016   # Pure water at 60°F (15.5556°C)


@validated_assay
@dataclass(frozen=True, slots=True)
class NGLAssay:
    """Immutable physical and thermodynamic specification for NGLs and LPG (GPA 2145).

    Bridges Mont Belvieu liquid volume (gal, bbl), international cargo mass (mt),
    energy (MMBtu, GJ), and gas-plant shrinkage volume (scf) for frac-spread and
    LPG arbitrage workflows.
    """

    name: str
    specific_gravity_60f: float
    hhv_btu_per_gal: float
    lhv_btu_per_gal: float
    gas_equivalent_scf_per_gal: float
    molecular_weight_g_per_mol: float | None = None
    vapor_pressure_psia_100f: float | None = None
    thermal_expansion_coeff_per_c: float | None = None
    refrigerated_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.specific_gravity_60f) or not (0.30 <= self.specific_gravity_60f <= 0.78):
            raise ValueError(
                f"specific_gravity_60f for NGL/LPG must be in [0.30, 0.78], got {self.specific_gravity_60f}"
            )
        if not math.isfinite(self.hhv_btu_per_gal) or self.hhv_btu_per_gal <= 0.0:
            raise ValueError(f"hhv_btu_per_gal must be finite and strictly positive, got {self.hhv_btu_per_gal}")
        if not math.isfinite(self.lhv_btu_per_gal) or not (0.0 < self.lhv_btu_per_gal <= self.hhv_btu_per_gal):
            raise ValueError(
                f"lhv_btu_per_gal ({self.lhv_btu_per_gal}) must be in (0, hhv_btu_per_gal ({self.hhv_btu_per_gal})]."
            )
        if not math.isfinite(self.gas_equivalent_scf_per_gal) or self.gas_equivalent_scf_per_gal <= 0.0:
            raise ValueError(
                f"gas_equivalent_scf_per_gal must be strictly positive, got {self.gas_equivalent_scf_per_gal}"
            )
        if self.molecular_weight_g_per_mol is not None:
            if not math.isfinite(self.molecular_weight_g_per_mol) or self.molecular_weight_g_per_mol <= 0.0:
                raise ValueError(f"molecular_weight_g_per_mol must be strictly positive, got {self.molecular_weight_g_per_mol}")
        if self.vapor_pressure_psia_100f is not None:
            if not math.isfinite(self.vapor_pressure_psia_100f) or self.vapor_pressure_psia_100f < 0.0:
                raise ValueError(f"vapor_pressure_psia_100f cannot be negative, got {self.vapor_pressure_psia_100f}")
        if self.thermal_expansion_coeff_per_c is not None:
            if not math.isfinite(self.thermal_expansion_coeff_per_c) or self.thermal_expansion_coeff_per_c <= 0.0:
                raise ValueError(
                    f"thermal_expansion_coeff_per_c must be strictly positive, got {self.thermal_expansion_coeff_per_c}"
                )
        if self.refrigerated_density_kg_per_m3 is not None:
            if not math.isfinite(self.refrigerated_density_kg_per_m3) or self.refrigerated_density_kg_per_m3 <= 0.0:
                raise ValueError(
                    f"refrigerated_density_kg_per_m3 must be strictly positive, got {self.refrigerated_density_kg_per_m3}"
                )

    # --- Mass, Volume & Energy Parity Properties ---
    @property
    def density_kg_per_m3(self) -> float:
        """Liquid density at 60°F (kg/m^3)."""
        return (self.specific_gravity_60f * _WATER_DENSITY_60F_KG_M3)

    @property
    def gal_per_mt(self) -> float:
        """US gallons per metric tonne at 60°F (Mont Belvieu <-> FEI/ARA factor)."""
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def bbl_per_mt(self) -> float:
        """Barrels per metric tonne at 60°F."""
        return (self.gal_per_mt / 42.0)

    @property
    def hhv_mmbtu_per_bbl(self) -> float:
        """Gross heating value per barrel (MMBtu_IT/bbl)."""
        return ((self.hhv_btu_per_gal * 42.0) / 1e6)

    @property
    def lhv_mmbtu_per_bbl(self) -> float:
        """Net heating value per barrel (MMBtu_IT/bbl)."""
        return ((self.lhv_btu_per_gal * 42.0) / 1e6)

    @property
    def hhv_mmbtu_per_mt(self) -> float:
        """Gross heating value per metric tonne (MMBtu_IT/mt)."""
        return ((self.hhv_btu_per_gal * self.gal_per_mt) / 1e6)

    @property
    def lhv_mmbtu_per_mt(self) -> float:
        """Net heating value per metric tonne (MMBtu_IT/mt)."""
        return ((self.lhv_btu_per_gal * self.gal_per_mt) / 1e6)

    @property
    def hhv_mj_per_kg(self) -> float:
        """Gross heating value (MJ/kg = GJ/mt)."""
        return (self.hhv_mmbtu_per_mt * MJ_PER_BTU_IT * 1000.0)

    @property
    def lhv_mj_per_kg(self) -> float:
        """Net heating value (MJ/kg = GJ/mt)."""
        return (self.lhv_mmbtu_per_mt * MJ_PER_BTU_IT * 1000.0)

    def _btu_per_gal(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"hhv", "gcv"}:
            return self.hhv_btu_per_gal
        if key in {"lhv", "ncv"}:
            return self.lhv_btu_per_gal
        raise ValueError(f"basis must be 'hhv'/'gcv' or 'lhv'/'ncv', got {basis!r}")

    def as_context(self, *, basis: str = "hhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [energy]."""
        return create_ngl_context(self, basis=basis, name=context_name)

    # --- Physical Volume <-> Mass <-> Energy <-> Shrinkage Conversions ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid NGL/LPG volume (e.g. gal, bbl, m^3) to mass (metric_ton)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        gal_q = volume.to(ureg.gallon)
        scale = _scale_factor(gal_q.magnitude, self.gal_per_mt)
        return ureg.Quantity(gal_q.magnitude / scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts NGL/LPG mass (e.g. mt, kg) to liquid volume (gallon at 60°F)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, self.gal_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.gallon)

    def volume_to_energy(self, volume: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts liquid NGL/LPG volume (e.g. gal, bbl) to energy (MMBtu_it)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        mmbtu_per_gal = self._btu_per_gal(basis) / 1e6
        gal_q = volume.to(ureg.gallon)
        scale = _scale_factor(gal_q.magnitude, mmbtu_per_gal)
        return ureg.Quantity(gal_q.magnitude * scale, ureg.MMBtu_it)

    def energy_to_volume(self, energy: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts energy (e.g. MMBtu_it, GJ) to liquid NGL/LPG volume (gallon)."""
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        mmbtu_per_gal = self._btu_per_gal(basis) / 1e6
        mmbtu_q = energy.to(ureg.MMBtu_it)
        scale = _scale_factor(mmbtu_q.magnitude, mmbtu_per_gal)
        return ureg.Quantity(mmbtu_q.magnitude / scale, ureg.gallon)

    def liquid_to_gas_shrinkage(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts extracted liquid NGL volume (gal/bbl) to gas-plant shrinkage volume (scf)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        gal_q = volume.to(ureg.gallon)
        scale = _scale_factor(gal_q.magnitude, self.gas_equivalent_scf_per_gal)
        return ureg.Quantity(gal_q.magnitude * scale, ureg.scf_60F_14696psi_dry)

    def mass_to_energy(self, mass: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts NGL/LPG mass (e.g. mt, kg) directly to energy (MMBtu_it)."""
        return self.volume_to_energy(self.mass_to_volume(mass), basis=basis)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts energy (e.g. MMBtu_it, GJ) directly to NGL/LPG mass (metric_ton)."""
        return self.volume_to_mass(self.energy_to_volume(energy, basis=basis))

    def gas_shrinkage_to_liquid(self, gas_volume: pint.Quantity) -> pint.Quantity:
        """Converts gas-plant shrinkage volume (e.g. scf, Mcf) to liquid NGL volume (gallon)."""
        if not isinstance(gas_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(gas_volume).__name__}")
        scf_q = gas_volume.to(ureg.scf_60F_14696psi_dry)
        scale = _scale_factor(scf_q.magnitude, self.gas_equivalent_scf_per_gal)
        return ureg.Quantity(scf_q.magnitude / scale, ureg.gallon)

    def refrigerated_volume_to_mass(self, liquid_volume: pint.Quantity) -> pint.Quantity:
        """Converts refrigerated VLGC cargo volume (m^3 at boiling point) to mass (metric_ton)."""
        if not isinstance(liquid_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(liquid_volume).__name__}")
        if self.refrigerated_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define refrigerated_density_kg_per_m3.")
        m3_q = liquid_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.refrigerated_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_refrigerated_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts LPG mass (e.g. mt) to refrigerated VLGC cargo volume (m^3 at boiling point)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.refrigerated_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define refrigerated_density_kg_per_m3.")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.refrigerated_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float = 15.5556,
    ) -> float:
        """Computes light-end thermal Volume Correction Factor (VCF) to 60°F (or 15°C)."""
        if self.thermal_expansion_coeff_per_c is None:
            raise ValueError(f"{self.name!r} does not define thermal_expansion_coeff_per_c.")
        if isinstance(observed_temp, pint.Quantity):
            t_obs_c = float(observed_temp.to(ureg.degC).magnitude)
        else:
            t_obs_c = float(observed_temp)
        if not math.isfinite(t_obs_c) or not math.isfinite(ref_temp_c):
            raise ValueError("Temperatures must be finite.")
        alpha = self.thermal_expansion_coeff_per_c
        dt = t_obs_c - ref_temp_c
        return (math.exp(-alpha * dt * (1.0 + 0.8 * alpha * dt)))

    def blend_with(
        self,
        other: NGLAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> NGLAssay:
        """Blends two NGL/LPG streams volumetrically (e.g. Propane/Butane autogas or Y-grade).

        - SG, HHV/LHV (Btu/gal), gas equivalent, thermal expansion: volume-weighted
        - Vapour pressure: mole-weighted (Raoult's Law) when molecular weights are known
        """
        if not isinstance(other, NGLAssay):
            raise TypeError(f"Expected NGLAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        v1 = float(self_volume_fraction)
        v2 = 1.0 - v1

        b_sg = (self.specific_gravity_60f * v1) + (other.specific_gravity_60f * v2)
        b_hhv = (self.hhv_btu_per_gal * v1) + (other.hhv_btu_per_gal * v2)
        b_lhv = (self.lhv_btu_per_gal * v1) + (other.lhv_btu_per_gal * v2)
        b_scf = (self.gas_equivalent_scf_per_gal * v1) + (other.gas_equivalent_scf_per_gal * v2)

        b_mw: float | None = None
        b_vp: float | None = None
        if self.molecular_weight_g_per_mol is not None and other.molecular_weight_g_per_mol is not None:
            n1 = (v1 * self.specific_gravity_60f) / self.molecular_weight_g_per_mol
            n2 = (v2 * other.specific_gravity_60f) / other.molecular_weight_g_per_mol
            x1, x2 = n1 / (n1 + n2), n2 / (n1 + n2)
            b_mw = ((x1 * self.molecular_weight_g_per_mol) + (x2 * other.molecular_weight_g_per_mol))
            if self.vapor_pressure_psia_100f is not None and other.vapor_pressure_psia_100f is not None:
                b_vp = ((x1 * self.vapor_pressure_psia_100f) + (x2 * other.vapor_pressure_psia_100f))
        elif self.vapor_pressure_psia_100f is not None and other.vapor_pressure_psia_100f is not None:
            b_vp = ((v1 * self.vapor_pressure_psia_100f) + (v2 * other.vapor_pressure_psia_100f))

        b_alpha: float | None = None
        if self.thermal_expansion_coeff_per_c is not None and other.thermal_expansion_coeff_per_c is not None:
            b_alpha = ((v1 * self.thermal_expansion_coeff_per_c) + (v2 * other.thermal_expansion_coeff_per_c))

        b_refrig_rho: float | None = None
        if self.refrigerated_density_kg_per_m3 is not None and other.refrigerated_density_kg_per_m3 is not None:
            b_refrig_rho = ((v1 * self.refrigerated_density_kg_per_m3) + (v2 * other.refrigerated_density_kg_per_m3))

        return NGLAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            specific_gravity_60f=(b_sg),
            hhv_btu_per_gal=(b_hhv),
            lhv_btu_per_gal=(b_lhv),
            gas_equivalent_scf_per_gal=(b_scf),
            molecular_weight_g_per_mol=b_mw,
            vapor_pressure_psia_100f=b_vp,
            thermal_expansion_coeff_per_c=b_alpha,
            refrigerated_density_kg_per_m3=b_refrig_rho,
        )

    # --- GPA 2145-16 Pure Component & Commercial LPG Presets ---
    @classmethod
    def ethane(cls) -> NGLAssay:
        """Pure Ethane (C2H6) — Mont Belvieu purity ethane."""
        return cls(
            name="Ethane (C2)",
            specific_gravity_60f=0.35628,
            hhv_btu_per_gal=66500.0,
            lhv_btu_per_gal=60910.0,
            gas_equivalent_scf_per_gal=37.498,
            molecular_weight_g_per_mol=30.069,
            vapor_pressure_psia_100f=None,  # Supercritical at 100°F (Tc = 89.9°F)
            thermal_expansion_coeff_per_c=0.00400,
        )

    @classmethod
    def propane(cls) -> NGLAssay:
        """Commercial Propane (C3H8) — Mont Belvieu / Argus FEI / Saudi CP."""
        return cls(
            name="Propane (C3)",
            specific_gravity_60f=0.50719,
            hhv_btu_per_gal=91420.0,
            lhv_btu_per_gal=84250.0,
            gas_equivalent_scf_per_gal=36.391,
            molecular_weight_g_per_mol=44.096,
            vapor_pressure_psia_100f=188.64,
            thermal_expansion_coeff_per_c=0.00295,
            refrigerated_density_kg_per_m3=582.0,
        )

    @classmethod
    def isobutane(cls) -> NGLAssay:
        """Isobutane (i-C4H10) — Alkylation feedstock."""
        return cls(
            name="Isobutane (iC4)",
            specific_gravity_60f=0.56283,
            hhv_btu_per_gal=98917.0,
            lhv_btu_per_gal=91340.0,
            gas_equivalent_scf_per_gal=30.639,
            molecular_weight_g_per_mol=58.122,
            vapor_pressure_psia_100f=72.58,
            thermal_expansion_coeff_per_c=0.00215,
            refrigerated_density_kg_per_m3=594.0,
        )

    @classmethod
    def n_butane(cls) -> NGLAssay:
        """Normal Butane (n-C4H10) — Gasoline blending & LPG benchmark."""
        return cls(
            name="Normal Butane (nC4)",
            specific_gravity_60f=0.58420,
            hhv_btu_per_gal=103769.0,
            lhv_btu_per_gal=95900.0,
            gas_equivalent_scf_per_gal=31.802,
            molecular_weight_g_per_mol=58.122,
            vapor_pressure_psia_100f=51.71,
            thermal_expansion_coeff_per_c=0.00200,
            refrigerated_density_kg_per_m3=601.0,
        )

    @classmethod
    def natural_gasoline(cls) -> NGLAssay:
        """Natural Gasoline (C5+ / Pentanes Plus) — Diluent & naphtha blendstock."""
        return cls(
            name="Natural Gasoline (C5+)",
            specific_gravity_60f=0.66500,
            hhv_btu_per_gal=115000.0,
            lhv_btu_per_gal=106500.0,
            gas_equivalent_scf_per_gal=26.970,
            molecular_weight_g_per_mol=78.000,
            vapor_pressure_psia_100f=12.00,
            thermal_expansion_coeff_per_c=0.00145,
        )

    @classmethod
    def y_grade(cls) -> NGLAssay:
        """US Y-Grade Raw NGL Mix (~42% C2, 28% C3, 7% iC4, 13% nC4, 10% C5+ by volume)."""
        return cls(
            name="Y-Grade Raw Mix",
            specific_gravity_60f=0.47350,
            hhv_btu_per_gal=85450.0,
            lhv_btu_per_gal=78750.0,
            gas_equivalent_scf_per_gal=34.920,
            molecular_weight_g_per_mol=43.800,
            thermal_expansion_coeff_per_c=0.00305,
        )


# ---------------------------------------------------------------------------
# Ferrous & Base Metal Ores / Concentrates Specification
# (Iron Ore, Copper/Zinc/Nickel Concentrates, Bauxite)
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class OreAssay:
    """Immutable physical and commercial specification for ores and metal concentrates.

    Bridges Wet Metric Tonnes (WMT), Dry Metric Tonnes (DMT), contained metal,
    smelter payable metal (percentage factor and minimum unit deduction),
    Dry Metric Tonne Unit (dmtu) pricing, and bulk stowage volume.
    """

    name: str
    primary_element: str
    grade_pct: float
    moisture_pct: float = 8.0
    payable_factor: float = 1.0
    min_deduction_pct: float = 0.0
    bulk_density_wmt_per_m3: float | None = None
    impurities_pct: Mapping[str, float] | None = None

    def __post_init__(self) -> None:
        if not self.primary_element or not self.primary_element.strip():
            raise ValueError("primary_element must be a non-empty string (e.g. 'Fe', 'Cu', 'Zn').")
        if not math.isfinite(self.grade_pct) or not (0.0 < self.grade_pct <= 100.0):
            raise ValueError(f"grade_pct must be finite and in (0, 100], got {self.grade_pct}")
        if not math.isfinite(self.moisture_pct) or not (0.0 <= self.moisture_pct < 50.0):
            raise ValueError(f"moisture_pct must be finite and in [0, 50), got {self.moisture_pct}")
        if not math.isfinite(self.payable_factor) or not (0.0 < self.payable_factor <= 1.0):
            raise ValueError(f"payable_factor must be finite and in (0, 1], got {self.payable_factor}")
        if not math.isfinite(self.min_deduction_pct) or not (0.0 <= self.min_deduction_pct < 100.0):
            raise ValueError(f"min_deduction_pct must be finite and in [0, 100), got {self.min_deduction_pct}")
        if self.bulk_density_wmt_per_m3 is not None:
            if not math.isfinite(self.bulk_density_wmt_per_m3) or self.bulk_density_wmt_per_m3 <= 0.0:
                raise ValueError(
                    f"bulk_density_wmt_per_m3 must be strictly positive, got {self.bulk_density_wmt_per_m3}"
                )
        if self.impurities_pct is not None:
            cleaned: dict[str, float] = {}
            for elem, pct in self.impurities_pct.items():
                val = float(pct)
                if not math.isfinite(val) or not (0.0 <= val <= 100.0):
                    raise ValueError(f"Impurity {elem!r} must be finite and in [0, 100], got {pct}")
                cleaned[str(elem)] = val
            if self.grade_pct + sum(cleaned.values()) > 100.0 + 1e-9:
                raise ValueError("grade_pct plus total impurities_pct cannot exceed 100%.")
            object.__setattr__(self, "impurities_pct", MappingProxyType(cleaned))

    def __hash__(self) -> int:
        impurities_key = tuple(sorted(self.impurities_pct.items())) if self.impurities_pct is not None else None
        return hash((
            self.name,
            self.primary_element,
            self.grade_pct,
            self.moisture_pct,
            self.payable_factor,
            self.min_deduction_pct,
            self.bulk_density_wmt_per_m3,
            impurities_key,
        ))

    # --- Core Grade & Stowage Properties ---
    @property
    def dry_mass_fraction(self) -> float:
        """Ratio of Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT)."""
        return 1.0 - (self.moisture_pct / 100.0)

    @property
    def payable_grade_pct(self) -> float:
        """Effective payable metal grade (%) after applying smelter payable factor and minimum deduction."""
        by_factor = self.grade_pct * self.payable_factor
        by_deduction = self.grade_pct - self.min_deduction_pct
        return (max(0.0, min(by_factor, by_deduction)))

    @property
    def stowage_factor_m3_per_wmt(self) -> float | None:
        """Vessel stowage factor (m^3 per WMT), or None if bulk_density_wmt_per_m3 is unset."""
        if self.bulk_density_wmt_per_m3 is None:
            return None
        return (1.0 / self.bulk_density_wmt_per_m3)

    def as_context(self, *, basis: str = "contained", context_name: str | None = None) -> pint.Context:
        """Bridges bulk volume and wet mass; use explicit methods for metal grades."""
        return create_ore_context(self, basis=basis, name=context_name)

    # --- Physical Mass Conversions (WMT <-> DMT <-> Contained / Payable Metal) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def wmt_to_dmt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to Dry Metric Tonnes (DMT)."""
        mt_q = self._validate_mass(wet_mass)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def dmt_to_wmt(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT)."""
        mt_q = self._validate_mass(dry_mass)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def dmt_to_contained_metal(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) of ore/concentrate to contained metal mass (metric_ton)."""
        mt_q = self._validate_mass(dry_mass)
        grade_frac = self.grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, grade_frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def contained_metal_to_dmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained metal mass (e.g. mt, lb) to required Dry Metric Tonnes (DMT)."""
        mt_q = self._validate_mass(metal_mass)
        grade_frac = self.grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, grade_frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def wmt_to_contained_metal(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to contained metal mass (metric_ton)."""
        return self.dmt_to_contained_metal(self.wmt_to_dmt(wet_mass))

    def contained_metal_to_wmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained metal mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.contained_metal_to_dmt(metal_mass))

    def dmt_to_payable_metal(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to smelter payable metal mass (metric_ton)."""
        mt_q = self._validate_mass(dry_mass)
        payable_frac = self.payable_grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, payable_frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def payable_metal_to_dmt(self, payable_mass: pint.Quantity) -> pint.Quantity:
        """Converts target smelter payable metal mass to required Dry Metric Tonnes (DMT)."""
        if self.payable_grade_pct <= 0.0:
            raise ValueError(f"{self.name!r} has zero payable grade.")
        mt_q = self._validate_mass(payable_mass)
        payable_frac = self.payable_grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, payable_frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def wmt_to_payable_metal(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to smelter payable metal mass (metric_ton)."""
        return self.dmt_to_payable_metal(self.wmt_to_dmt(wet_mass))

    def payable_metal_to_wmt(self, payable_mass: pint.Quantity) -> pint.Quantity:
        """Converts target payable metal mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.payable_metal_to_dmt(payable_mass))

    def wmt_to_volume(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to bulk vessel stowage volume (m^3)."""
        if self.bulk_density_wmt_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_wmt_per_m3.")
        mt_q = self._validate_mass(wet_mass)
        scale = _scale_factor(mt_q.magnitude, self.bulk_density_wmt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def volume_to_wmt(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts bulk vessel stowage volume (e.g. m^3) to Wet Metric Tonnes (WMT)."""
        if not isinstance(bulk_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bulk_volume).__name__}")
        if self.bulk_density_wmt_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_wmt_per_m3.")
        m3_q = bulk_volume.to(ureg.meter ** 3)
        scale = _scale_factor(m3_q.magnitude, self.bulk_density_wmt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def dmt_to_impurity_mass(self, dry_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Computes contained mass (metric_ton) of a specific impurity element from DMT."""
        if self.impurities_pct is None or element not in self.impurities_pct:
            raise KeyError(f"Impurity {element!r} is not defined on {self.name!r}.")
        mt_q = self._validate_mass(dry_mass)
        imp_frac = self.impurities_pct[element] / 100.0
        scale = _scale_factor(mt_q.magnitude, imp_frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    # --- Commercial Pricing Parity ($/dmtu <-> $/dmt <-> $/wmt) ---
    def dmtu_price_to_dmt_price(self, price_per_dmtu: float) -> float:
        """Converts a $/dmtu quote (1 dmtu = 1% grade per DMT) to $/dmt."""
        p = float(price_per_dmtu)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmtu must be finite, got {price_per_dmtu}")
        return (p * self.grade_pct)

    def dmt_price_to_dmtu_price(self, price_per_dmt: float) -> float:
        """Converts a $/dmt quote to $/dmtu."""
        p = float(price_per_dmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt must be finite, got {price_per_dmt}")
        return (p / self.grade_pct)

    def dmt_price_to_wmt_price(self, price_per_dmt: float) -> float:
        """Converts a $/dmt quote to $/wmt (wet basis)."""
        p = float(price_per_dmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt must be finite, got {price_per_dmt}")
        return (p * self.dry_mass_fraction)

    def wmt_price_to_dmt_price(self, price_per_wmt: float) -> float:
        """Converts a $/wmt quote to $/dmt (dry basis)."""
        p = float(price_per_wmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_wmt must be finite, got {price_per_wmt}")
        return (p / self.dry_mass_fraction)

    # --- Ore & Concentrate Blending (WMT Moisture, DMT Grade & Impurities) ---
    def blend_with(
        self,
        other: OreAssay,
        *,
        self_wet_mass_fraction: float,
        blend_name: str | None = None,
    ) -> OreAssay:
        """Computes a physically consistent ore blend from wet mass fractions.

        - Moisture blends linearly by Wet Metric Tonnes (WMT).
        - Dry-basis grade and impurities blend linearly by Dry Metric Tonnes (DMT).
        """
        if not isinstance(other, OreAssay):
            raise TypeError(f"Expected OreAssay, got {type(other).__name__}")
        if self.primary_element != other.primary_element:
            raise ValueError(
                f"Cannot blend ores with different primary elements: {self.primary_element!r} vs {other.primary_element!r}"
            )
        if not math.isfinite(self_wet_mass_fraction) or not (0.0 <= self_wet_mass_fraction <= 1.0):
            raise ValueError("self_wet_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_wet_mass_fraction)
        w2 = 1.0 - w1

        b_moisture = (self.moisture_pct * w1) + (other.moisture_pct * w2)

        d1 = w1 * self.dry_mass_fraction
        d2 = w2 * other.dry_mass_fraction
        total_d = d1 + d2
        dw1, dw2 = d1 / total_d, d2 / total_d

        b_grade = (self.grade_pct * dw1) + (other.grade_pct * dw2)

        b_impurities: dict[str, float] | None = None
        if self.impurities_pct is not None and other.impurities_pct is not None:
            all_keys = set(self.impurities_pct) | set(other.impurities_pct)
            b_impurities = {
                k: ((self.impurities_pct.get(k, 0.0) * dw1) + (other.impurities_pct.get(k, 0.0) * dw2))
                for k in sorted(all_keys)
            }

        b_density: float | None = None
        if self.bulk_density_wmt_per_m3 is not None and other.bulk_density_wmt_per_m3 is not None:
            # Harmonic mean of bulk densities (additive stowage volumes)
            vol_per_wmt = (w1 / self.bulk_density_wmt_per_m3) + (w2 / other.bulk_density_wmt_per_m3)
            b_density = (1.0 / vol_per_wmt)

        return OreAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            primary_element=self.primary_element,
            grade_pct=(b_grade),
            moisture_pct=(b_moisture),
            payable_factor=((self.payable_factor * dw1) + (other.payable_factor * dw2)),
            min_deduction_pct=((self.min_deduction_pct * dw1) + (other.min_deduction_pct * dw2)),
            bulk_density_wmt_per_m3=b_density,
            impurities_pct=b_impurities,
        )

    # --- Global Ferrous & Base Metal Benchmark Presets ---
    @classmethod
    def iron_ore_62fe(cls) -> OreAssay:
        """Platts IODEX 62% Fe CFR China fines benchmark."""
        return cls(
            name="Iron Ore 62% Fe Fines",
            primary_element="Fe",
            grade_pct=62.0,
            moisture_pct=8.0,
            bulk_density_wmt_per_m3=2.20,
            impurities_pct={"SiO2": 4.0, "Al2O3": 2.25, "P": 0.09, "S": 0.02},
        )

    @classmethod
    def iron_ore_65fe(cls) -> OreAssay:
        """Fastmarkets 65% Fe Brazil fines (Carajás) benchmark."""
        return cls(
            name="Iron Ore 65% Fe Fines",
            primary_element="Fe",
            grade_pct=65.0,
            moisture_pct=8.5,
            bulk_density_wmt_per_m3=2.30,
            impurities_pct={"SiO2": 2.0, "Al2O3": 1.30, "P": 0.05, "S": 0.01},
        )

    @classmethod
    def iron_ore_58fe(cls) -> OreAssay:
        """58% Fe low-grade iron ore fines benchmark."""
        return cls(
            name="Iron Ore 58% Fe Fines",
            primary_element="Fe",
            grade_pct=58.0,
            moisture_pct=9.0,
            bulk_density_wmt_per_m3=2.10,
            impurities_pct={"SiO2": 5.5, "Al2O3": 3.0, "P": 0.08, "S": 0.03},
        )

    @classmethod
    def copper_concentrate(cls) -> OreAssay:
        """Standard clean Copper Concentrate (28% Cu, 96.5% payable / 1.0% min deduction)."""
        return cls(
            name="Copper Concentrate (28% Cu)",
            primary_element="Cu",
            grade_pct=28.0,
            moisture_pct=9.0,
            payable_factor=0.965,
            min_deduction_pct=1.0,
            bulk_density_wmt_per_m3=1.80,
            impurities_pct={"As": 0.15, "Bi": 0.02, "Pb": 0.30, "Zn": 0.80},
        )

    @classmethod
    def zinc_concentrate(cls) -> OreAssay:
        """Standard Zinc Concentrate (52% Zn, 85% payable / 8.0% min deduction)."""
        return cls(
            name="Zinc Concentrate (52% Zn)",
            primary_element="Zn",
            grade_pct=52.0,
            moisture_pct=8.5,
            payable_factor=0.85,
            min_deduction_pct=8.0,
            bulk_density_wmt_per_m3=1.90,
            impurities_pct={"Fe": 7.5, "Pb": 1.2, "SiO2": 2.0},
        )

    @classmethod
    def bauxite_guinea(cls) -> OreAssay:
        """Guinea CBG metallurgical-grade Bauxite (45% Al2O3)."""
        return cls(
            name="Guinea Bauxite (45% Al2O3)",
            primary_element="Al2O3",
            grade_pct=45.0,
            moisture_pct=10.0,
            bulk_density_wmt_per_m3=1.40,
            impurities_pct={"SiO2": 3.0, "Fe2O3": 22.0},
        )

    @classmethod
    def nickel_concentrate(cls) -> OreAssay:
        """Standard Nickel Sulphide Concentrate (12% Ni, 75% smelter payable)."""
        return cls(
            name="Nickel Concentrate (12% Ni)",
            primary_element="Ni",
            grade_pct=12.0,
            moisture_pct=9.0,
            payable_factor=0.75,
            min_deduction_pct=0.0,
            bulk_density_wmt_per_m3=1.85,
            impurities_pct={"Fe": 35.0, "S": 30.0, "MgO": 4.5, "As": 0.05},
        )

# ---------------------------------------------------------------------------
# Battery Metals & Precursor Chemicals Specification
# (Spodumene SC6, Li2CO3, LiOH·H2O, NiSO4·6H2O, CoSO4·7H2O, MnSO4·H2O, MHP)
# ---------------------------------------------------------------------------

# Exact IUPAC stoichiometric conversion factors to Lithium Carbonate Equivalent (LCE, Li2CO3)
_LI_TO_LCE_FACTOR: Final[float] = 5.32279        # M(Li2CO3) / (2 * M(Li))   = 73.891 / 13.882
_LI2O_TO_LCE_FACTOR: Final[float] = 2.47284      # M(Li2CO3) / M(Li2O)       = 73.891 / 29.881
_LIOH_H2O_TO_LCE_FACTOR: Final[float] = 0.88041  # M(Li2CO3) / (2*M(LiOH·H2O)) = 73.891 / 83.928


@validated_assay
@dataclass(frozen=True, slots=True)
class BatteryMaterialAssay:
    """Immutable physical and stoichiometric specification for battery raw materials.

    Bridges Wet Metric Tonnes (WMT), Dry Metric Tonnes (DMT), contained battery
    metal (Li, Ni, Co, Mn), and Lithium Carbonate Equivalent (LCE) including
    refinery metallurgical recovery rates.
    """

    name: str
    active_metal: str
    active_metal_pct: float
    moisture_pct: float = 0.0
    purity_pct: float | None = None
    lce_factor_per_dmt: float | None = None
    default_recovery_rate: float = 1.0

    def __post_init__(self) -> None:
        if not self.active_metal or not self.active_metal.strip():
            raise ValueError("active_metal must be a non-empty string (e.g. 'Li', 'Ni', 'Co', 'Mn').")
        if not math.isfinite(self.active_metal_pct) or not (0.0 < self.active_metal_pct <= 100.0):
            raise ValueError(f"active_metal_pct must be finite and in (0, 100], got {self.active_metal_pct}")
        if not math.isfinite(self.moisture_pct) or not (0.0 <= self.moisture_pct < 60.0):
            raise ValueError(f"moisture_pct must be finite and in [0, 60), got {self.moisture_pct}")
        if self.purity_pct is not None:
            if not math.isfinite(self.purity_pct) or not (0.0 < self.purity_pct <= 100.0):
                raise ValueError(f"purity_pct must be finite and in (0, 100], got {self.purity_pct}")
        if self.lce_factor_per_dmt is not None:
            if not math.isfinite(self.lce_factor_per_dmt) or self.lce_factor_per_dmt <= 0.0:
                raise ValueError(f"lce_factor_per_dmt must be strictly positive, got {self.lce_factor_per_dmt}")
        if not math.isfinite(self.default_recovery_rate) or not (0.0 < self.default_recovery_rate <= 1.0):
            raise ValueError(f"default_recovery_rate must be in (0, 1], got {self.default_recovery_rate}")

    # --- Core Stoichiometric & Dry-Basis Properties ---
    @property
    def dry_mass_fraction(self) -> float:
        """Ratio of Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT)."""
        return 1.0 - (self.moisture_pct / 100.0)

    @property
    def dmt_per_mt_contained_metal(self) -> float:
        """Dry Metric Tonnes of material required per 1 mt of contained active metal."""
        return (100.0 / self.active_metal_pct)

    @property
    def dmt_per_mt_lce_stoichiometric(self) -> float | None:
        """Stoichiometric DMT required per 1 mt LCE (100% recovery), or None if non-lithium."""
        if self.lce_factor_per_dmt is None:
            return None
        return (1.0 / self.lce_factor_per_dmt)

    @property
    def dmt_per_mt_lce_recoverable(self) -> float | None:
        """Effective DMT required per 1 mt LCE at default_recovery_rate, or None if non-lithium."""
        if self.lce_factor_per_dmt is None:
            return None
        return (1.0 / (self.lce_factor_per_dmt * self.default_recovery_rate))

    @property
    def li2o_pct(self) -> float | None:
        """Equivalent Li2O grade (%) on a dry basis for lithium-bearing materials."""
        if self.lce_factor_per_dmt is None:
            return None
        return ((self.lce_factor_per_dmt / _LI2O_TO_LCE_FACTOR) * 100.0)

    def as_context(self, *, basis: str = "contained", context_name: str | None = None) -> pint.Context:
        """Returns a named context; use explicit methods for metal and LCE masses."""
        return create_battery_material_context(self, basis=basis, name=context_name)

    # --- Physical Mass Conversions (WMT <-> DMT <-> Contained Metal <-> LCE) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def wmt_to_dmt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to Dry Metric Tonnes (DMT)."""
        mt_q = self._validate_mass(wet_mass)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def dmt_to_wmt(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT)."""
        mt_q = self._validate_mass(dry_mass)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def dmt_to_contained_metal(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to contained active metal mass (metric_ton)."""
        mt_q = self._validate_mass(dry_mass)
        frac = self.active_metal_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def contained_metal_to_dmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained active metal mass to required Dry Metric Tonnes (DMT)."""
        mt_q = self._validate_mass(metal_mass)
        frac = self.active_metal_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def wmt_to_contained_metal(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to contained active metal mass (metric_ton)."""
        return self.dmt_to_contained_metal(self.wmt_to_dmt(wet_mass))

    def contained_metal_to_wmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained active metal mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.contained_metal_to_dmt(metal_mass))

    def dmt_to_lce(
        self,
        dry_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) of lithium feedstock/chemical to LCE (metric_ton)."""
        if self.lce_factor_per_dmt is None:
            raise ValueError(f"{self.name!r} is not a lithium-bearing material (lce_factor_per_dmt is None).")
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError(f"recovery_rate must be in (0, 1], got {rec}")
        mt_q = self._validate_mass(dry_mass)
        factor = self.lce_factor_per_dmt * rec
        scale = _scale_factor(mt_q.magnitude, factor)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def lce_to_dmt(
        self,
        lce_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts target LCE mass (metric_ton) to required Dry Metric Tonnes (DMT) of feedstock/chemical."""
        if self.lce_factor_per_dmt is None:
            raise ValueError(f"{self.name!r} is not a lithium-bearing material (lce_factor_per_dmt is None).")
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError(f"recovery_rate must be in (0, 1], got {rec}")
        mt_q = self._validate_mass(lce_mass)
        factor = self.lce_factor_per_dmt * rec
        scale = _scale_factor(mt_q.magnitude, factor)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def wmt_to_lce(
        self,
        wet_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to LCE (metric_ton)."""
        return self.dmt_to_lce(self.wmt_to_dmt(wet_mass), recovery_rate=recovery_rate)

    def lce_to_wmt(
        self,
        lce_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts target LCE mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.lce_to_dmt(lce_mass, recovery_rate=recovery_rate))

    # --- Commercial Pricing Parity ($/dmt Product <-> $/mt Contained Metal or LCE) ---
    def dmt_price_to_metal_price(self, price_per_dmt: float) -> float:
        """Converts a $/dmt product quote (e.g. NiSO4, CoSO4) to $/mt contained active metal."""
        p = float(price_per_dmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt must be finite, got {price_per_dmt}")
        return (p / (self.active_metal_pct / 100.0))

    def metal_price_to_dmt_price(self, price_per_mt_metal: float) -> float:
        """Converts a $/mt contained metal quote to $/dmt product equivalent."""
        p = float(price_per_mt_metal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_metal must be finite, got {price_per_mt_metal}")
        return (p * (self.active_metal_pct / 100.0))

    def dmt_price_to_lce_price(
        self,
        price_per_dmt: float,
        *,
        recovery_rate: float | None = None,
        conversion_cost_per_mt_lce: float = 0.0,
    ) -> float:
        """Converts feedstock/chemical $/dmt (e.g. SC6) to implied $/mt LCE parity."""
        if self.lce_factor_per_dmt is None:
            raise ValueError(f"{self.name!r} is not a lithium-bearing material.")
        p = float(price_per_dmt)
        cost = float(conversion_cost_per_mt_lce)
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(p) or not math.isfinite(cost) or not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError("Inputs must be finite and recovery_rate must be in (0, 1].")
        return ((p / (self.lce_factor_per_dmt * rec)) + cost)

    def lce_price_to_dmt_price(
        self,
        price_per_mt_lce: float,
        *,
        recovery_rate: float | None = None,
        conversion_cost_per_mt_lce: float = 0.0,
    ) -> float:
        """Computes implied feedstock/chemical $/dmt netback from spot $/mt LCE."""
        if self.lce_factor_per_dmt is None:
            raise ValueError(f"{self.name!r} is not a lithium-bearing material.")
        p_lce = float(price_per_mt_lce)
        cost = float(conversion_cost_per_mt_lce)
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(p_lce) or not math.isfinite(cost) or not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError("Inputs must be finite and recovery_rate must be in (0, 1].")
        return ((p_lce - cost) * (self.lce_factor_per_dmt * rec))

    def blend_with(
        self,
        other: BatteryMaterialAssay,
        *,
        self_wet_mass_fraction: float,
        blend_name: str | None = None,
    ) -> BatteryMaterialAssay:
        """Blends two battery material lots of the same active metal.

        - Moisture blends linearly by Wet Metric Tonnes (WMT).
        - Active metal %, LCE factor, purity, and recovery blend by Dry Metric Tonnes (DMT).
        """
        if not isinstance(other, BatteryMaterialAssay):
            raise TypeError(f"Expected BatteryMaterialAssay, got {type(other).__name__}")
        if self.active_metal != other.active_metal:
            raise ValueError(
                f"Cannot blend materials with different active metals: {self.active_metal!r} vs {other.active_metal!r}"
            )
        if not math.isfinite(self_wet_mass_fraction) or not (0.0 <= self_wet_mass_fraction <= 1.0):
            raise ValueError("self_wet_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_wet_mass_fraction)
        w2 = 1.0 - w1
        b_moisture = (self.moisture_pct * w1) + (other.moisture_pct * w2)

        d1 = w1 * self.dry_mass_fraction
        d2 = w2 * other.dry_mass_fraction
        total_d = d1 + d2
        dw1, dw2 = d1 / total_d, d2 / total_d

        b_metal = (self.active_metal_pct * dw1) + (other.active_metal_pct * dw2)
        b_rec = (self.default_recovery_rate * dw1) + (other.default_recovery_rate * dw2)

        b_purity: float | None = None
        if self.purity_pct is not None and other.purity_pct is not None:
            b_purity = ((self.purity_pct * dw1) + (other.purity_pct * dw2))

        b_lce: float | None = None
        if self.lce_factor_per_dmt is not None and other.lce_factor_per_dmt is not None:
            b_lce = ((self.lce_factor_per_dmt * dw1) + (other.lce_factor_per_dmt * dw2))

        return BatteryMaterialAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            active_metal=self.active_metal,
            active_metal_pct=(b_metal),
            moisture_pct=(b_moisture),
            purity_pct=b_purity,
            lce_factor_per_dmt=b_lce,
            default_recovery_rate=(b_rec),
        )

    # --- Global Battery Material Benchmark Presets & Factories ---
    @classmethod
    def from_spodumene_grade(
        cls,
        li2o_pct: float,
        *,
        moisture_pct: float = 8.0,
        recovery_rate: float = 0.85,
        name: str | None = None,
    ) -> BatteryMaterialAssay:
        """Constructs a custom Spodumene concentrate assay directly from % Li2O."""
        grade = float(li2o_pct)
        if not math.isfinite(grade) or not (0.5 <= grade <= 8.03):
            raise ValueError(f"Spodumene li2o_pct must be in [0.5, 8.03] (pure spodumene max), got {li2o_pct}")
        li_pct = grade * (_LI2O_TO_LCE_FACTOR / _LI_TO_LCE_FACTOR)
        return cls(
            name=name or f"Spodumene Concentrate ({grade:.2f}% Li2O)",
            active_metal="Li",
            active_metal_pct=(li_pct),
            moisture_pct=moisture_pct,
            lce_factor_per_dmt=((grade / 100.0) * _LI2O_TO_LCE_FACTOR),
            default_recovery_rate=recovery_rate,
        )
    @classmethod
    def spodumene_sc6(cls) -> BatteryMaterialAssay:
        """Spodumene Concentrate 6.0% Li2O (SC6 CIF China, ~85% typical refinery recovery)."""
        li_pct = 6.0 * (_LI2O_TO_LCE_FACTOR / _LI_TO_LCE_FACTOR)  # ~2.7875% Li
        return cls(
            name="Spodumene Concentrate (SC6, 6% Li2O)",
            active_metal="Li",
            active_metal_pct=(li_pct),
            moisture_pct=8.0,
            lce_factor_per_dmt=(0.06 * _LI2O_TO_LCE_FACTOR),
            default_recovery_rate=0.85,
        )

    @classmethod
    def lithium_carbonate_bg(cls) -> BatteryMaterialAssay:
        """Battery-Grade Lithium Carbonate (Li2CO3, 99.5% min purity)."""
        li_pct = 0.995 * (100.0 / _LI_TO_LCE_FACTOR)  # ~18.693% Li
        return cls(
            name="Lithium Carbonate (BG 99.5%)",
            active_metal="Li",
            active_metal_pct=(li_pct),
            moisture_pct=0.0,
            purity_pct=99.5,
            lce_factor_per_dmt=0.995,
            default_recovery_rate=1.0,
        )

    @classmethod
    def lithium_hydroxide_bg(cls) -> BatteryMaterialAssay:
        """Battery-Grade Lithium Hydroxide Monohydrate (LiOH·H2O, 56.5% min LiOH)."""
        li_pct = 0.995 * (_LIOH_H2O_TO_LCE_FACTOR / _LI_TO_LCE_FACTOR) * 100.0  # ~16.458% Li
        return cls(
            name="Lithium Hydroxide Monohydrate (BG)",
            active_metal="Li",
            active_metal_pct=(li_pct),
            moisture_pct=0.0,
            purity_pct=99.5,
            lce_factor_per_dmt=(0.995 * _LIOH_H2O_TO_LCE_FACTOR),
            default_recovery_rate=1.0,
        )

    @classmethod
    def nickel_sulphate_hexahydrate(cls) -> BatteryMaterialAssay:
        """Battery-Grade Nickel Sulphate Hexahydrate (NiSO4·6H2O, 22.2% Ni)."""
        return cls(
            name="Nickel Sulphate (NiSO4·6H2O)",
            active_metal="Ni",
            active_metal_pct=22.2,
            moisture_pct=0.0,
            purity_pct=99.5,
        )

    @classmethod
    def cobalt_sulphate_heptahydrate(cls) -> BatteryMaterialAssay:
        """Battery-Grade Cobalt Sulphate Heptahydrate (CoSO4·7H2O, 20.5% Co)."""
        return cls(
            name="Cobalt Sulphate (CoSO4·7H2O)",
            active_metal="Co",
            active_metal_pct=20.5,
            moisture_pct=0.0,
            purity_pct=99.5,
        )

    @classmethod
    def manganese_sulphate_monohydrate(cls) -> BatteryMaterialAssay:
        """Battery-Grade Manganese Sulphate Monohydrate (MnSO4·H2O, 32.0% Mn)."""
        return cls(
            name="Manganese Sulphate (MnSO4·H2O)",
            active_metal="Mn",
            active_metal_pct=32.0,
            moisture_pct=0.0,
            purity_pct=99.5,
        )

    @classmethod
    def mhp_nickel(cls) -> BatteryMaterialAssay:
        """Mixed Hydroxide Precipitate (MHP — 38.0% Ni dry basis, ~30% moisture)."""
        return cls(
            name="Mixed Hydroxide Precipitate (MHP)",
            active_metal="Ni",
            active_metal_pct=38.0,
            moisture_pct=30.0,
        )

    @classmethod
    def cobalt_hydroxide(cls) -> BatteryMaterialAssay:
        """DRC Cobalt Hydroxide intermediate (30.0% Co dry basis, ~35% moisture)."""
        return cls(
            name="Cobalt Hydroxide (30% Co)",
            active_metal="Co",
            active_metal_pct=30.0,
            moisture_pct=35.0,
        )

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


@validated_assay
@dataclass(frozen=True, slots=True)
class UraniumAssay:
    """Immutable isotopic and front-end fuel-cycle specification for Uranium.

    Bridges U3O8 yellowcake (lb), natural UF6 (kgU / kg UF6), Enriched Uranium
    Product (EUP, kgU), and Separative Work Units (SWU) via U-235 mass balance
    and the Dirac value function.
    """

    name: str
    product_u235_pct: float
    feed_u235_pct: float = _NATURAL_U235_PCT
    tails_u235_pct: float = 0.25
    conversion_loss_pct: float = 0.5

    def __post_init__(self) -> None:
        if not math.isfinite(self.product_u235_pct) or not (0.0 < self.product_u235_pct < 100.0):
            raise ValueError(f"product_u235_pct must be in (0, 100), got {self.product_u235_pct}")
        if not math.isfinite(self.feed_u235_pct) or not (0.0 < self.feed_u235_pct < 100.0):
            raise ValueError(f"feed_u235_pct must be in (0, 100), got {self.feed_u235_pct}")
        if not math.isfinite(self.tails_u235_pct) or not (0.0 < self.tails_u235_pct < 100.0):
            raise ValueError(f"tails_u235_pct must be in (0, 100), got {self.tails_u235_pct}")
        if not math.isfinite(self.conversion_loss_pct) or not (0.0 <= self.conversion_loss_pct < 10.0):
            raise ValueError(f"conversion_loss_pct must be in [0, 10), got {self.conversion_loss_pct}")

        # For enriched product (product > feed), tails must be strictly below feed
        if self.product_u235_pct > self.feed_u235_pct and self.tails_u235_pct >= self.feed_u235_pct:
            raise ValueError(
                f"tails_u235_pct ({self.tails_u235_pct}%) must be strictly less than "
                f"feed_u235_pct ({self.feed_u235_pct}%) for enrichment."
            )

    # --- Enrichment Mass Balance & SWU Ratios ---
    @property
    def is_enriched(self) -> bool:
        """True if product U-235 assay exceeds feed U-235 assay."""
        return self.product_u235_pct > self.feed_u235_pct

    @property
    def feed_to_product_ratio(self) -> float:
        """Feed uranium (kgU as UF6) required per 1 kgU of product (F/P)."""
        if not self.is_enriched:
            return 1.0
        xp = self.product_u235_pct / 100.0
        xf = self.feed_u235_pct / 100.0
        xw = self.tails_u235_pct / 100.0
        return ((xp - xw) / (xf - xw))

    @property
    def tails_to_product_ratio(self) -> float:
        """Depleted uranium tails (kgU) produced per 1 kgU of product (W/P)."""
        if not self.is_enriched:
            return 0.0
        return (self.feed_to_product_ratio - 1.0)

    @property
    def swu_per_kg_u(self) -> float:
        """Separative Work Units (kg-SWU) required per 1 kgU of enriched product."""
        if not self.is_enriched:
            return 0.0
        xp = self.product_u235_pct / 100.0
        xf = self.feed_u235_pct / 100.0
        xw = self.tails_u235_pct / 100.0
        fp = (xp - xw) / (xf - xw)
        wp = fp - 1.0
        swu = _dirac_value_function(xp) + (wp * _dirac_value_function(xw)) - (fp * _dirac_value_function(xf))
        return (swu)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Bridges EUP mass and thermal energy; use explicit methods for fuel-cycle masses."""
        return create_uranium_context(self, name=context_name)

    # --- Front-End Stoichiometric Conversions (U3O8 <-> kgU <-> UF6 <-> UO2) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.kilogram)

    def u3o8_to_natural_kgu(
        self,
        u3o8_mass: pint.Quantity,
        *,
        include_conversion_loss: bool = True,
    ) -> pint.Quantity:
        """Converts U3O8 yellowcake mass (e.g. lb, kg) to natural uranium (kgU as UF6)."""
        kg_q = self._validate_mass(u3o8_mass)
        loss_mult = (1.0 - self.conversion_loss_pct / 100.0) if include_conversion_loss else 1.0
        factor = _U_FRAC_IN_U3O8 * loss_mult
        scale = _scale_factor(kg_q.magnitude, factor)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def natural_kgu_to_u3o8(
        self,
        natural_u_mass: pint.Quantity,
        *,
        include_conversion_loss: bool = True,
    ) -> pint.Quantity:
        """Converts natural uranium (kgU as UF6) to required U3O8 yellowcake (lb)."""
        kg_q = self._validate_mass(natural_u_mass)
        loss_mult = (1.0 - self.conversion_loss_pct / 100.0) if include_conversion_loss else 1.0
        factor = _U_FRAC_IN_U3O8 * loss_mult
        scale = _scale_factor(kg_q.magnitude, factor)
        u3o8_kg = ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)
        return u3o8_kg.to(ureg.pound)

    def kgu_to_uf6_mass(self, u_mass: pint.Quantity) -> pint.Quantity:
        """Converts elemental uranium mass (kgU) to total UF6 compound mass (kg UF6)."""
        kg_q = self._validate_mass(u_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UF6)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def uf6_mass_to_kgu(self, uf6_mass: pint.Quantity) -> pint.Quantity:
        """Converts total UF6 compound mass (kg UF6) to elemental uranium mass (kgU)."""
        kg_q = self._validate_mass(uf6_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UF6)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def kgu_to_uo2_mass(self, u_mass: pint.Quantity) -> pint.Quantity:
        """Converts elemental uranium mass (kgU) to fabricated UO2 pellet mass (kg UO2)."""
        kg_q = self._validate_mass(u_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UO2)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def uo2_mass_to_kgu(self, uo2_mass: pint.Quantity) -> pint.Quantity:
        """Converts fabricated UO2 pellet mass (kg UO2) to elemental uranium mass (kgU)."""
        kg_q = self._validate_mass(uo2_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UO2)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    # --- Enrichment Conversions (EUP <-> Feed UF6 <-> U3O8 & SWU) ---
    def eup_to_feed_kgu(self, eup_u_mass: pint.Quantity) -> pint.Quantity:
        """Converts Enriched Uranium Product mass (kgU) to required natural feed (kgU as UF6)."""
        kg_q = self._validate_mass(eup_u_mass)
        scale = _scale_factor(kg_q.magnitude, self.feed_to_product_ratio)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def feed_kgu_to_eup(self, feed_u_mass: pint.Quantity) -> pint.Quantity:
        """Converts natural feed (kgU as UF6) to producible Enriched Uranium Product (kgU)."""
        kg_q = self._validate_mass(feed_u_mass)
        scale = _scale_factor(kg_q.magnitude, self.feed_to_product_ratio)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def eup_to_u3o8(
        self,
        eup_u_mass: pint.Quantity,
        *,
        include_conversion_loss: bool = True,
    ) -> pint.Quantity:
        """Converts Enriched Uranium Product mass (kgU) directly to required U3O8 (lb)."""
        return self.natural_kgu_to_u3o8(
            self.eup_to_feed_kgu(eup_u_mass),
            include_conversion_loss=include_conversion_loss,
        )

    def required_swu(self, eup_u_mass: pint.Quantity) -> float:
        """Computes total Separative Work Units (kg-SWU) required for an EUP quantity."""
        kg_u = float(self._validate_mass(eup_u_mass).magnitude)
        return (kg_u * self.swu_per_kg_u)

    def u3o8_to_eup(
        self,
        u3o8_mass: pint.Quantity,
        *,
        include_conversion_loss: bool = True,
    ) -> pint.Quantity:
        """Converts U3O8 yellowcake mass (e.g. lb) directly to producible EUP (kgU)."""
        return self.feed_kgu_to_eup(
            self.u3o8_to_natural_kgu(u3o8_mass, include_conversion_loss=include_conversion_loss)
        )

    def eup_to_tails_kgu(self, eup_u_mass: pint.Quantity) -> pint.Quantity:
        """Computes depleted uranium tails mass (kgU) generated when producing EUP."""
        kg_q = self._validate_mass(eup_u_mass)
        scale = _scale_factor(kg_q.magnitude, self.tails_to_product_ratio)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def eup_to_thermal_energy(
        self,
        eup_u_mass: pint.Quantity,
        *,
        burnup_gwd_per_mtu: float = 45.0,
    ) -> pint.Quantity:
        """Converts EUP mass (kgU) to reactor thermal energy output (MWh_th).

        1 GWd/mtU = 24 MWh_th / kgU.
        """
        if not math.isfinite(burnup_gwd_per_mtu) or burnup_gwd_per_mtu <= 0.0:
            raise ValueError(f"burnup_gwd_per_mtu must be strictly positive, got {burnup_gwd_per_mtu}")
        kg_q = self._validate_mass(eup_u_mass)
        mwh_th_per_kgu = burnup_gwd_per_mtu * 24.0
        scale = _scale_factor(kg_q.magnitude, mwh_th_per_kgu)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.MWh)

    def eup_to_electrical_energy(
        self,
        eup_u_mass: pint.Quantity,
        *,
        burnup_gwd_per_mtu: float = 45.0,
        thermal_efficiency: float = 0.33,
    ) -> pint.Quantity:
        """Converts EUP mass (kgU) to net electrical energy output (MWh_e)."""
        if not math.isfinite(thermal_efficiency) or not (0.0 < thermal_efficiency <= 1.0):
            raise ValueError(f"thermal_efficiency must be in (0, 1], got {thermal_efficiency}")
        th_q = self.eup_to_thermal_energy(eup_u_mass, burnup_gwd_per_mtu=burnup_gwd_per_mtu)
        scale = _scale_factor(th_q.magnitude, thermal_efficiency)
        return ureg.Quantity(th_q.magnitude * scale, ureg.MWh)

    # --- Commercial Fuel-Cycle Economics & Optimal Tails Assay ---
    def eup_cost_per_kgu(
        self,
        u3o8_price_per_lb: float,
        conversion_price_per_kgu: float,
        swu_price_per_swu: float,
        *,
        include_conversion_loss: bool = True,
    ) -> float:
        """Computes total Enriched Uranium Product cost ($/kgU as enriched UF6)."""
        p_u3o8 = float(u3o8_price_per_lb)
        p_conv = float(conversion_price_per_kgu)
        p_swu = float(swu_price_per_swu)
        if not all(math.isfinite(x) and x >= 0.0 for x in (p_u3o8, p_conv, p_swu)):
            raise ValueError("All fuel-cycle prices must be finite and non-negative.")

        lb_u3o8_per_feed_kgu = float(
            self.natural_kgu_to_u3o8(
                ureg.Quantity(1.0, ureg.kilogram),
                include_conversion_loss=include_conversion_loss,
            ).magnitude
        )
        feed_cost_per_kgu = (p_u3o8 * lb_u3o8_per_feed_kgu) + p_conv
        return ((self.feed_to_product_ratio * feed_cost_per_kgu) + (self.swu_per_kg_u * p_swu))

    @staticmethod
    def optimal_tails_u235_pct(
        feed_cost_per_kgu: float,
        swu_price_per_swu: float,
        *,
        feed_u235_pct: float = _NATURAL_U235_PCT,
    ) -> float:
        """Solves for the cost-minimising tails assay (% U-235) via bisection.

        At the optimum, C_feed / C_swu = (V(x_f) - V(x_w)) - (x_f - x_w) * V'(x_w).
        """
        cf = float(feed_cost_per_kgu)
        cs = float(swu_price_per_swu)
        xf = float(feed_u235_pct) / 100.0
        if not (math.isfinite(cf) and cf > 0.0 and math.isfinite(cs) and cs > 0.0):
            raise ValueError("feed_cost_per_kgu and swu_price_per_swu must be strictly positive.")
        if not (0.0 < xf < 1.0):
            raise ValueError("feed_u235_pct must be strictly in (0, 100).")

        target = cf / cs
        v_xf = _dirac_value_function(xf)

        def _rhs(xw: float) -> float:
            v_xw = _dirac_value_function(xw)
            dv_xw = 2.0 * math.log(xw / (1.0 - xw)) + (2.0 * xw - 1.0) / (xw * (1.0 - xw))
            return (v_xf - v_xw) - (dv_xw * (xf - xw))

        lo, hi = 1e-5, xf - 1e-6
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            if _rhs(mid) > target:
                lo = mid
            else:
                hi = mid
        return (0.5 * (lo + hi) * 100.0)

    def blend_with(
        self,
        other: UraniumAssay,
        *,
        self_u_mass_fraction: float,
        blend_name: str | None = None,
    ) -> UraniumAssay:
        """Computes isotopic U-235 mass-weighted downblending / blending of two uranium streams."""
        if not isinstance(other, UraniumAssay):
            raise TypeError(f"Expected UraniumAssay, got {type(other).__name__}")
        if not math.isfinite(self_u_mass_fraction) or not (0.0 <= self_u_mass_fraction <= 1.0):
            raise ValueError("self_u_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_u_mass_fraction)
        w2 = 1.0 - w1
        blended_u235 = (self.product_u235_pct * w1) + (other.product_u235_pct * w2)
        return UraniumAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            product_u235_pct=(blended_u235),
            feed_u235_pct=self.feed_u235_pct,
            tails_u235_pct=self.tails_u235_pct,
            conversion_loss_pct=self.conversion_loss_pct,
        )

    # --- Nuclear Fuel Cycle Benchmark Presets ---
    @classmethod
    def natural_uranium(cls) -> UraniumAssay:
        """Natural Uranium (0.711% U-235) — U3O8 yellowcake / natural UF6."""
        return cls(name="Natural Uranium (0.711% U-235)", product_u235_pct=_NATURAL_U235_PCT)

    @classmethod
    def leu_4_5pct(cls, *, tails_u235_pct: float = 0.25) -> UraniumAssay:
        """Standard Light-Water Reactor (LWR) Low-Enriched Uranium (4.5% U-235)."""
        return cls(
            name="LWR LEU (4.5% U-235)",
            product_u235_pct=4.5,
            tails_u235_pct=tails_u235_pct,
        )

    @classmethod
    def haleu_19_75pct(cls, *, tails_u235_pct: float = 0.20) -> UraniumAssay:
        """High-Assay Low-Enriched Uranium (HALEU, 19.75% U-235) for SMRs."""
        return cls(
            name="SMR HALEU (19.75% U-235)",
            product_u235_pct=19.75,
            tails_u235_pct=tails_u235_pct,
        )

    @classmethod
    def depleted_tails(cls, *, tails_u235_pct: float = 0.25) -> UraniumAssay:
        """Depleted Uranium (DU) enrichment tails."""
        return cls(
            name=f"Depleted Uranium ({tails_u235_pct:.2f}% U-235)",
            product_u235_pct=tails_u235_pct,
            tails_u235_pct=tails_u235_pct,
        )

# ---------------------------------------------------------------------------
# Solid Biomass & Wood Pellets Specification (ISO 17225 / ISO 18125)
# (Industrial Pellets I2P, ENplus A1, Torrefied Pellets, Wood Chips, PKS)
# ---------------------------------------------------------------------------

# ISO 18125 latent heat constants for solid biofuels (GJ/mt per 1 wt%, i.e. MJ/kg per 1 wt%)
_ISO18125_H_GJ_PER_MT_PCT: Final[float] = 0.2122    # Latent heat of water from H combustion (dry basis)
_ISO18125_M_GJ_PER_MT_PCT: Final[float] = 0.02443   # Latent heat of vaporisation of moisture at 25°C


@validated_assay
@dataclass(frozen=True, slots=True)
class BiomassAssay:
    """Immutable calorific, moisture, and stowage specification for solid biomass.

    Bridges Wet Metric Tonnes (WMT), Oven-Dry Tonnes (ODT/DMT), bulk stowage
    volume (m^3), and energy (GJ, MWh, MMBtu) on Net As Received (LHV/NAR),
    Gross As Received (HHV/GAR), and Dry Basis (DB) per ISO 18125.
    """

    name: str
    lhv_gj_per_mt_ar: float
    hhv_gj_per_mt_ar: float
    hhv_gj_per_mt_db: float
    moisture_pct: float = 8.0
    bulk_density_kg_per_m3: float = 650.0
    hydrogen_db_pct: float = 6.0
    ash_db_pct: float | None = None
    chlorine_db_pct: float | None = None

    def __init__(
        self,
        name: str = "custom_biomass",
        *,
        lhv_gj_per_mt_ar: float | None = None,
        hhv_gj_per_mt_ar: float | None = None,
        hhv_gj_per_mt_db: float | None = None,
        moisture_pct: float = 8.0,
        bulk_density_kg_per_m3: float = 650.0,
        hydrogen_db_pct: float = 6.0,
        ash_db_pct: float | None = None,
        chlorine_db_pct: float | None = None,
    ) -> None:
        """Constructs an immutable biomass assay from exactly one calorific benchmark."""
        provided = sum(p is not None for p in (lhv_gj_per_mt_ar, hhv_gj_per_mt_ar, hhv_gj_per_mt_db))
        if provided != 1:
            raise ValueError(
                "Provide exactly one calorific parameter: "
                "'lhv_gj_per_mt_ar', 'hhv_gj_per_mt_ar', or 'hhv_gj_per_mt_db'."
            )

        m = float(moisture_pct)
        h_db = float(hydrogen_db_pct)
        rho = float(bulk_density_kg_per_m3)
        if not math.isfinite(m) or not (0.0 <= m < 75.0):
            raise ValueError(f"moisture_pct must be finite and in [0, 75), got {m}")
        if not math.isfinite(h_db) or not (0.0 <= h_db <= 10.0):
            raise ValueError(f"hydrogen_db_pct must be finite and in [0, 10], got {h_db}")
        if not math.isfinite(rho) or not (50.0 <= rho <= 1200.0):
            raise ValueError(f"bulk_density_kg_per_m3 must be finite and in [50, 1200], got {rho}")

        ash = float(ash_db_pct) if ash_db_pct is not None else None
        if ash is not None and (not math.isfinite(ash) or not (0.0 <= ash < 50.0)):
            raise ValueError(f"ash_db_pct must be finite and in [0, 50), got {ash}")

        cl = float(chlorine_db_pct) if chlorine_db_pct is not None else None
        if cl is not None and (not math.isfinite(cl) or not (0.0 <= cl <= 5.0)):
            raise ValueError(f"chlorine_db_pct must be finite and in [0, 5], got {cl}")

        dry_frac = 1.0 - (m / 100.0)
        # ISO 18125 gross-to-net latent heat deduction (GJ/mt as-received)
        latent_deduction = (_ISO18125_H_GJ_PER_MT_PCT * h_db * dry_frac) + (_ISO18125_M_GJ_PER_MT_PCT * m)

        if lhv_gj_per_mt_ar is not None:
            res_lhv_ar = float(lhv_gj_per_mt_ar)
            res_hhv_ar = res_lhv_ar + latent_deduction
            res_hhv_db = res_hhv_ar / dry_frac
        elif hhv_gj_per_mt_ar is not None:
            res_hhv_ar = float(hhv_gj_per_mt_ar)
            res_lhv_ar = res_hhv_ar - latent_deduction
            res_hhv_db = res_hhv_ar / dry_frac
        else:
            assert hhv_gj_per_mt_db is not None
            res_hhv_db = float(hhv_gj_per_mt_db)
            res_hhv_ar = res_hhv_db * dry_frac
            res_lhv_ar = res_hhv_ar - latent_deduction

        if not all(math.isfinite(x) and x > 0.0 for x in (res_lhv_ar, res_hhv_ar, res_hhv_db)):
            raise ValueError(
                f"Resolved calorific values must be finite and positive "
                f"(LHV_ar={res_lhv_ar}, HHV_ar={res_hhv_ar}, HHV_db={res_hhv_db})."
            )

        object.__setattr__(self, "name", str(name))
        object.__setattr__(self, "lhv_gj_per_mt_ar", (res_lhv_ar))
        object.__setattr__(self, "hhv_gj_per_mt_ar", (res_hhv_ar))
        object.__setattr__(self, "hhv_gj_per_mt_db", (res_hhv_db))
        object.__setattr__(self, "moisture_pct", m)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)
        object.__setattr__(self, "hydrogen_db_pct", h_db)
        object.__setattr__(self, "ash_db_pct", ash)
        object.__setattr__(self, "chlorine_db_pct", cl)

    # --- Energy Density & Stowage Properties ---
    @property
    def dry_mass_fraction(self) -> float:
        """Ratio of Oven-Dry Tonnes (ODT/DMT) to Wet Metric Tonnes (WMT)."""
        return 1.0 - (self.moisture_pct / 100.0)

    @property
    def lhv_mwh_per_mt_ar(self) -> float:
        return (self.lhv_gj_per_mt_ar / 3.6)

    @property
    def hhv_mwh_per_mt_ar(self) -> float:
        return (self.hhv_gj_per_mt_ar / 3.6)

    @property
    def lhv_mmbtu_per_mt_ar(self) -> float:
        return (self.lhv_gj_per_mt_ar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mmbtu_per_mt_ar(self) -> float:
        return (self.hhv_gj_per_mt_ar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def stowage_factor_m3_per_wmt(self) -> float:
        """Vessel stowage factor (m^3 per WMT)."""
        return (1000.0 / self.bulk_density_kg_per_m3)

    @property
    def lhv_gj_per_m3_bulk(self) -> float:
        """Volumetric energy density (GJ LHV per bulk m^3)."""
        return (self.lhv_gj_per_mt_ar * (self.bulk_density_kg_per_m3 / 1000.0))

    @property
    def lhv_gj_per_mt_db(self) -> float:
        """Net Calorific Value on a Dry Basis (ISO 18125): HHV_db - 0.2122 * H_db."""
        return (self.hhv_gj_per_mt_db - (_ISO18125_H_GJ_PER_MT_PCT * self.hydrogen_db_pct))

    @property
    def lhv_kcal_per_kg_ar(self) -> float:
        """Net As Received calorific value in kcal_IT/kg (Asian PKS / Japanese FIT basis)."""
        return (self.lhv_gj_per_mt_ar / _KCAL_KG_TO_GJ_MT)

    @property
    def hhv_kcal_per_kg_ar(self) -> float:
        """Gross As Received calorific value in kcal_IT/kg."""
        return (self.hhv_gj_per_mt_ar / _KCAL_KG_TO_GJ_MT)

    @property
    def ash_ar_pct(self) -> float | None:
        """As-received ash content (wt%), or None if ash_db_pct is unset."""
        if self.ash_db_pct is None:
            return None
        return (self.ash_db_pct * self.dry_mass_fraction)

    def _gj_per_mt(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"lhv", "nar", "ncv"}:
            return self.lhv_gj_per_mt_ar
        if key in {"hhv", "gar", "gcv"}:
            return self.hhv_gj_per_mt_ar
        raise ValueError(f"basis must be 'lhv'/'nar' or 'hhv'/'gar', got {basis!r}")

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [mass], [volume], and [energy]."""
        return create_biomass_context(self, basis=basis, name=context_name)

    # --- Physical Conversions (WMT <-> ODT/DMT <-> Bulk m^3 <-> Energy) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def wmt_to_odt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to Oven-Dry Tonnes (ODT / DMT)."""
        mt_q = self._validate_mass(wet_mass)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def odt_to_wmt(self, oven_dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Oven-Dry Tonnes (ODT / DMT) to Wet Metric Tonnes (WMT)."""
        mt_q = self._validate_mass(oven_dry_mass)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def wmt_to_bulk_volume(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to bulk vessel/silo stowage volume (m^3)."""
        mt_q = self._validate_mass(wet_mass)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def bulk_volume_to_wmt(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts bulk stowage volume (e.g. m^3) to Wet Metric Tonnes (WMT)."""
        if not isinstance(bulk_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bulk_volume).__name__}")
        m3_q = bulk_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_energy(self, wet_mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts as-received mass (WMT) to total energy (GJ)."""
        mt_q = self._validate_mass(wet_mass)
        gj_per_mt = self._gj_per_mt(basis)
        scale = _scale_factor(mt_q.magnitude, gj_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.GJ)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. GJ, MWh, MMBtu) to required as-received mass (metric_ton)."""
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        gj_per_mt = self._gj_per_mt(basis)
        gj_q = energy.to(ureg.GJ)
        scale = _scale_factor(gj_q.magnitude, gj_per_mt)
        return ureg.Quantity(gj_q.magnitude / scale, ureg.metric_ton)

    def bulk_volume_to_energy(self, bulk_volume: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts bulk stowage volume (m^3) directly to energy (GJ)."""
        return self.mass_to_energy(self.bulk_volume_to_wmt(bulk_volume), basis=basis)

    def energy_to_bulk_volume(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. GJ, MWh) directly to bulk stowage volume (m^3)."""
        return self.wmt_to_bulk_volume(self.energy_to_mass(energy, basis=basis))

    def at_moisture(
        self,
        target_moisture_pct: float,
        *,
        new_bulk_density_kg_per_m3: float | None = None,
        name: str | None = None,
    ) -> BiomassAssay:
        """Projects this biomass assay to a new moisture content holding dry-matter properties constant.

        Useful for modelling wood-chip drying, pelletising, or open-air moisture uptake.
        """
        m_new = float(target_moisture_pct)
        if not math.isfinite(m_new) or not (0.0 <= m_new < 75.0):
            raise ValueError(f"target_moisture_pct must be finite and in [0, 75), got {target_moisture_pct}")
        return BiomassAssay(
            name=name or f"{self.name} ({m_new:.1f}% M)",
            hhv_gj_per_mt_db=self.hhv_gj_per_mt_db,
            moisture_pct=m_new,
            bulk_density_kg_per_m3=(
                self.bulk_density_kg_per_m3 if new_bulk_density_kg_per_m3 is None else new_bulk_density_kg_per_m3
            ),
            hydrogen_db_pct=self.hydrogen_db_pct,
            ash_db_pct=self.ash_db_pct,
            chlorine_db_pct=self.chlorine_db_pct,
        )

    def wmt_price_to_odt_price(self, price_per_wmt: float) -> float:
        """Converts a $/WMT quote to $/ODT (Oven-Dry / Bone-Dry Tonne)."""
        p = float(price_per_wmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_wmt must be finite, got {price_per_wmt}")
        return (p / self.dry_mass_fraction)

    def odt_price_to_wmt_price(self, price_per_odt: float) -> float:
        """Converts a $/ODT (Oven-Dry / Bone-Dry Tonne) quote to $/WMT."""
        p = float(price_per_odt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_odt must be finite, got {price_per_odt}")
        return (p * self.dry_mass_fraction)

    def blend_with(
        self,
        other: BiomassAssay,
        *,
        self_wet_mass_fraction: float,
        blend_name: str | None = None,
    ) -> BiomassAssay:
        """Computes a physically consistent solid biomass blend from wet mass fractions.

        - Moisture blends linearly by Wet Metric Tonnes (WMT).
        - Dry-basis HHV, hydrogen, ash, and chlorine blend by Oven-Dry Tonnes (ODT).
        - Bulk density blends via harmonic mean (additive stowage volumes).
        """
        if not isinstance(other, BiomassAssay):
            raise TypeError(f"Expected BiomassAssay, got {type(other).__name__}")
        if not math.isfinite(self_wet_mass_fraction) or not (0.0 <= self_wet_mass_fraction <= 1.0):
            raise ValueError("self_wet_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_wet_mass_fraction)
        w2 = 1.0 - w1
        b_moisture = (self.moisture_pct * w1) + (other.moisture_pct * w2)

        d1 = w1 * self.dry_mass_fraction
        d2 = w2 * other.dry_mass_fraction
        total_d = d1 + d2
        dw1, dw2 = d1 / total_d, d2 / total_d

        b_hhv_db = (self.hhv_gj_per_mt_db * dw1) + (other.hhv_gj_per_mt_db * dw2)
        b_h_db = (self.hydrogen_db_pct * dw1) + (other.hydrogen_db_pct * dw2)

        vol_per_mt = (w1 / self.bulk_density_kg_per_m3) + (w2 / other.bulk_density_kg_per_m3)
        b_rho = 1.0 / vol_per_mt

        def _mix_db(a: float | None, b: float | None) -> float | None:
            if a is None or b is None:
                return None
            return ((a * dw1) + (b * dw2))

        return BiomassAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            hhv_gj_per_mt_db=b_hhv_db,
            moisture_pct=(b_moisture),
            bulk_density_kg_per_m3=(b_rho),
            hydrogen_db_pct=(b_h_db),
            ash_db_pct=_mix_db(self.ash_db_pct, other.ash_db_pct),
            chlorine_db_pct=_mix_db(self.chlorine_db_pct, other.chlorine_db_pct),
        )

    # --- Global Solid Biomass & Wood Pellet Benchmark Presets ---
    @classmethod
    def industrial_pellets_i2(cls) -> BiomassAssay:
        """Argus CIF ARA I2 Industrial Wood Pellets (17.0 GJ/mt NAR, 8% moisture)."""
        return cls(
            name="Industrial Wood Pellets (I2)",
            lhv_gj_per_mt_ar=17.00,
            moisture_pct=8.0,
            bulk_density_kg_per_m3=650.0,
            hydrogen_db_pct=6.0,
            ash_db_pct=1.5,
            chlorine_db_pct=0.05,
        )

    @classmethod
    def enplus_a1_pellets(cls) -> BiomassAssay:
        """ENplus A1 Premium Wood Pellets (17.2 GJ/mt NAR, 7% moisture, <=0.7% ash DB)."""
        return cls(
            name="ENplus A1 Wood Pellets",
            lhv_gj_per_mt_ar=17.20,
            moisture_pct=7.0,
            bulk_density_kg_per_m3=660.0,
            hydrogen_db_pct=6.0,
            ash_db_pct=0.5,
            chlorine_db_pct=0.02,
        )

    @classmethod
    def torrefied_pellets(cls) -> BiomassAssay:
        """Torrefied / Black Biomass Pellets (21.0 GJ/mt NAR, 4% moisture, coal-replacement grade)."""
        return cls(
            name="Torrefied Black Pellets",
            lhv_gj_per_mt_ar=21.00,
            moisture_pct=4.0,
            bulk_density_kg_per_m3=750.0,
            hydrogen_db_pct=5.2,
            ash_db_pct=2.0,
            chlorine_db_pct=0.03,
        )

    @classmethod
    def forestry_wood_chips(cls) -> BiomassAssay:
        """Seasoned Forestry Wood Chips (35% moisture, 19.2 GJ/mt HHV dry basis)."""
        return cls(
            name="Forestry Wood Chips",
            hhv_gj_per_mt_db=19.20,
            moisture_pct=35.0,
            bulk_density_kg_per_m3=320.0,
            hydrogen_db_pct=6.0,
            ash_db_pct=2.0,
            chlorine_db_pct=0.03,
        )

    @classmethod
    def palm_kernel_shells(cls) -> BiomassAssay:
        """Southeast Asian Palm Kernel Shells (PKS — ~3,800 kcal/kg NAR, Japanese FIT biomass)."""
        return cls(
            name="Palm Kernel Shells (PKS)",
            lhv_gj_per_mt_ar=3800.0 * _KCAL_KG_TO_GJ_MT,
            moisture_pct=15.0,
            bulk_density_kg_per_m3=580.0,
            hydrogen_db_pct=5.8,
            ash_db_pct=3.5,
            chlorine_db_pct=0.08,
        )

# ---------------------------------------------------------------------------
# Petrochemicals Specification (Olefins, Aromatics & Intermediates)
# (Ethylene, Propylene, Benzene, Toluene, Paraxylene, Styrene, MEG)
# ---------------------------------------------------------------------------

_LB_PER_MT: Final[float] = 1000.0 / 0.45359237    # Exact: 2204.62262185 lb/mt


@validated_assay
@dataclass(frozen=True, slots=True)
class PetrochemAssay:
    """Immutable physical and stoichiometric specification for petrochemicals.

    Bridges liquid/cryogenic volume (gal, bbl, m^3), mass (mt, lb, kg), molar
    quantities (kmol), ASTM D1555 thermal volume correction, and stoichiometric
    feedstock-to-product conversion yields.
    """

    name: str
    formula: str
    molecular_weight_g_per_mol: float
    density_kg_per_m3: float
    reference_temp_c: float = 15.0
    purity_pct: float = 99.5
    thermal_expansion_coeff_per_c: float | None = None
    freezing_point_c: float | None = None

    def __post_init__(self) -> None:
        if not self.formula or not self.formula.strip():
            raise ValueError("formula must be a non-empty string (e.g. 'C2H4', 'C6H6').")
        if not math.isfinite(self.molecular_weight_g_per_mol) or self.molecular_weight_g_per_mol <= 0.0:
            raise ValueError(
                f"molecular_weight_g_per_mol must be strictly positive, got {self.molecular_weight_g_per_mol}"
            )
        if not math.isfinite(self.density_kg_per_m3) or not (300.0 <= self.density_kg_per_m3 <= 2000.0):
            raise ValueError(
                f"density_kg_per_m3 must be finite and in [300, 2000] kg/m^3, got {self.density_kg_per_m3}"
            )
        if not math.isfinite(self.reference_temp_c) or not (-200.0 <= self.reference_temp_c <= 100.0):
            raise ValueError(f"reference_temp_c must be in [-200, 100] °C, got {self.reference_temp_c}")
        if not math.isfinite(self.purity_pct) or not (0.0 < self.purity_pct <= 100.0):
            raise ValueError(f"purity_pct must be in (0, 100], got {self.purity_pct}")
        if self.thermal_expansion_coeff_per_c is not None:
            if not math.isfinite(self.thermal_expansion_coeff_per_c) or self.thermal_expansion_coeff_per_c <= 0.0:
                raise ValueError(
                    f"thermal_expansion_coeff_per_c must be strictly positive, got {self.thermal_expansion_coeff_per_c}"
                )
        if self.freezing_point_c is not None and not math.isfinite(self.freezing_point_c):
            raise ValueError(f"freezing_point_c must be finite, got {self.freezing_point_c}")

    # --- Volumetric & Gravimetric Parity Properties ---
    @property
    def specific_gravity(self) -> float:
        """Specific gravity relative to water at 15°C (999.102 kg/m^3)."""
        return (self.density_kg_per_m3 / 999.102)

    @property
    def gal_per_mt(self) -> float:
        """US gallons per metric tonne at reference_temp_c."""
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def bbl_per_mt(self) -> float:
        """Barrels per metric tonne at reference_temp_c."""
        return (self.gal_per_mt / 42.0)

    @property
    def lb_per_gal(self) -> float:
        """US pounds per US gallon at reference_temp_c (USGC aromatics standard)."""
        return (_LB_PER_MT / self.gal_per_mt)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [substance]."""
        return create_petrochem_context(self, name=context_name)

    # --- Physical Volume <-> Mass <-> Moles Conversions ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid/cryogenic volume (e.g. gal, bbl, m^3) to mass (metric_ton)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts mass (e.g. mt, lb, kg) to liquid/cryogenic volume (m^3)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def mass_to_moles(self, mass: pint.Quantity, *, active_only: bool = False) -> pint.Quantity:
        """Converts mass (e.g. mt, kg) to molar amount (kilomole).

        If active_only=True, scales mass by purity_pct before dividing by molecular weight.
        """
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        kg_q = mass.to(ureg.kilogram)
        purity_mult = (self.purity_pct / 100.0) if active_only else 1.0
        # 1 kg / (MW g/mol) = 1 / MW kmol
        kmol_per_kg = purity_mult / self.molecular_weight_g_per_mol
        scale = _scale_factor(kg_q.magnitude, kmol_per_kg)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilomole)

    def moles_to_mass(self, moles: pint.Quantity, *, active_only: bool = False) -> pint.Quantity:
        """Converts molar amount (e.g. kmol, mol) to mass (metric_ton)."""
        if not isinstance(moles, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(moles).__name__}")
        kmol_q = moles.to(ureg.kilomole)
        purity_mult = (self.purity_pct / 100.0) if active_only else 1.0
        # 1 kmol * MW (kg/kmol) / 1000 = MW / 1000 mt
        mt_per_kmol = (self.molecular_weight_g_per_mol / 1000.0) / purity_mult
        scale = _scale_factor(kmol_q.magnitude, mt_per_kmol)
        return ureg.Quantity(kmol_q.magnitude * scale, ureg.metric_ton)

    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float | None = None,
    ) -> float:
        """Computes ASTM D1555 thermal Volume Correction Factor (VCF) to reference temperature."""
        if self.thermal_expansion_coeff_per_c is None:
            raise ValueError(f"{self.name!r} does not define thermal_expansion_coeff_per_c.")
        t_ref = self.reference_temp_c if ref_temp_c is None else float(ref_temp_c)
        if isinstance(observed_temp, pint.Quantity):
            t_obs_c = float(observed_temp.to(ureg.degC).magnitude)
        else:
            t_obs_c = float(observed_temp)
        if not math.isfinite(t_obs_c) or not math.isfinite(t_ref):
            raise ValueError("Temperatures must be finite.")
        alpha = self.thermal_expansion_coeff_per_c
        dt = t_obs_c - t_ref
        return (math.exp(-alpha * dt * (1.0 + 0.8 * alpha * dt)))

    def density_at_temp(self, observed_temp: float | pint.Quantity) -> float:
        """Computes hot-loaded liquid density (kg/m^3) at observed_temp via ASTM D4311."""
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        return (self.density_kg_per_m3 * vcf)

    def observed_volume_to_mass(
        self,
        observed_volume: pint.Quantity,
        observed_temp: float | pint.Quantity,
    ) -> pint.Quantity:
        """Converts observed liquid volume at observed_temp directly to mass (metric_ton) via ASTM D1555 VCF."""
        if not isinstance(observed_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(observed_volume).__name__}")
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        scale = _scale_factor(observed_volume.magnitude, vcf)
        std_vol = ureg.Quantity(observed_volume.magnitude * scale, observed_volume.units)
        return self.volume_to_mass(std_vol)

    # --- Stoichiometric Feedstock-to-Product Yield Parity ---
    def stoichiometric_product_mass(
        self,
        feedstock_mass: pint.Quantity,
        product_assay: PetrochemAssay,
        *,
        molar_ratio: float = 1.0,
        yield_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes downstream product mass (metric_ton) from this feedstock mass.

        Args:
            feedstock_mass: Input mass of this feedstock (e.g. Q_(1000, 'metric_ton')).
            product_assay: Target downstream PetrochemAssay.
            molar_ratio: Moles of product formed per 1 mole of this feedstock.
            yield_efficiency: Plant stoichiometric yield efficiency in (0, 1].
        """
        if not isinstance(product_assay, PetrochemAssay):
            raise TypeError(f"Expected PetrochemAssay, got {type(product_assay).__name__}")
        if not math.isfinite(molar_ratio) or molar_ratio <= 0.0:
            raise ValueError(f"molar_ratio must be strictly positive, got {molar_ratio}")
        if not math.isfinite(yield_efficiency) or not (0.0 < yield_efficiency <= 1.0):
            raise ValueError(f"yield_efficiency must be in (0, 1], got {yield_efficiency}")

        feed_kmol = self.mass_to_moles(feedstock_mass, active_only=True)
        scale = _scale_factor(feed_kmol.magnitude, molar_ratio * yield_efficiency)
        prod_kmol = ureg.Quantity(feed_kmol.magnitude * scale, ureg.kilomole)
        return product_assay.moles_to_mass(prod_kmol, active_only=True)

    def stoichiometric_feedstock_mass(
        self,
        product_mass: pint.Quantity,
        feedstock_assay: PetrochemAssay,
        *,
        molar_ratio: float = 1.0,
        yield_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes required feedstock mass (metric_ton) to produce a target mass of this product.

        Args:
            product_mass: Target mass of this product (e.g. Q_(1000, 'metric_ton')).
            feedstock_assay: Upstream feedstock PetrochemAssay.
            molar_ratio: Moles of feedstock consumed per 1 mole of this product.
            yield_efficiency: Plant stoichiometric yield efficiency in (0, 1].
        """
        if not isinstance(feedstock_assay, PetrochemAssay):
            raise TypeError(f"Expected PetrochemAssay, got {type(feedstock_assay).__name__}")
        if not math.isfinite(molar_ratio) or molar_ratio <= 0.0:
            raise ValueError(f"molar_ratio must be strictly positive, got {molar_ratio}")
        if not math.isfinite(yield_efficiency) or not (0.0 < yield_efficiency <= 1.0):
            raise ValueError(f"yield_efficiency must be in (0, 1], got {yield_efficiency}")

        prod_kmol = self.mass_to_moles(product_mass, active_only=True)
        scale = _scale_factor(prod_kmol.magnitude, molar_ratio / yield_efficiency)
        feed_kmol = ureg.Quantity(prod_kmol.magnitude * scale, ureg.kilomole)
        return feedstock_assay.moles_to_mass(feed_kmol, active_only=True)

    # --- Commercial Regional Arbitrage Pricing ($/gal & $/lb <-> $/mt) ---
    def gal_price_to_mt_price(self, price_per_gal: float) -> float:
        """Converts USGC $/gal quote to Europe/Asia $/mt equivalent."""
        p = float(price_per_gal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gal must be finite, got {price_per_gal}")
        return (p * self.gal_per_mt)

    def mt_price_to_gal_price(self, price_per_mt: float) -> float:
        """Converts Europe/Asia $/mt quote to USGC $/gal equivalent."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.gal_per_mt)

    @staticmethod
    def lb_price_to_mt_price(price_per_lb: float) -> float:
        """Converts USGC $/lb quote (e.g. Ethylene/Propylene) to $/mt."""
        p = float(price_per_lb)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb must be finite, got {price_per_lb}")
        return (p * _LB_PER_MT)

    @staticmethod
    def mt_price_to_lb_price(price_per_mt: float) -> float:
        """Converts $/mt quote to USGC $/lb equivalent."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / _LB_PER_MT)

    def blend_with(
        self,
        other: PetrochemAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> PetrochemAssay:
        """Blends two miscible liquid petrochemical streams (e.g. BTX / Mixed Xylenes).

        - Density & thermal expansion: volume-weighted
        - Molecular weight: mole-weighted
        - Purity: mass-weighted
        """
        if not isinstance(other, PetrochemAssay):
            raise TypeError(f"Expected PetrochemAssay, got {type(other).__name__}")
        if abs(self.reference_temp_c - other.reference_temp_c) > 1e-6:
            raise ValueError(
                f"Cannot blend streams with different reference temperatures: "
                f"{self.reference_temp_c}°C vs {other.reference_temp_c}°C."
            )
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        v1 = float(self_volume_fraction)
        v2 = 1.0 - v1

        b_rho = (self.density_kg_per_m3 * v1) + (other.density_kg_per_m3 * v2)
        m1 = self.density_kg_per_m3 * v1
        m2 = other.density_kg_per_m3 * v2
        total_m = m1 + m2
        w1, w2 = m1 / total_m, m2 / total_m

        n1 = m1 / self.molecular_weight_g_per_mol
        n2 = m2 / other.molecular_weight_g_per_mol
        x1, x2 = n1 / (n1 + n2), n2 / (n1 + n2)
        b_mw = (self.molecular_weight_g_per_mol * x1) + (other.molecular_weight_g_per_mol * x2)
        b_purity = (self.purity_pct * w1) + (other.purity_pct * w2)

        b_alpha: float | None = None
        if self.thermal_expansion_coeff_per_c is not None and other.thermal_expansion_coeff_per_c is not None:
            b_alpha = ((self.thermal_expansion_coeff_per_c * v1) + (other.thermal_expansion_coeff_per_c * v2))

        return PetrochemAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            formula=self.formula if self.formula == other.formula else "Blend",
            molecular_weight_g_per_mol=(b_mw),
            density_kg_per_m3=(b_rho),
            reference_temp_c=self.reference_temp_c,
            purity_pct=(b_purity),
            thermal_expansion_coeff_per_c=b_alpha,
            freezing_point_c=None,
        )

    # --- Global Petrochemical Benchmark Presets ---
    @classmethod
    def ethylene(cls) -> PetrochemAssay:
        """Polymer-Grade Ethylene (C2H4 — cryogenic liquid at -103.7°C, 1 atm)."""
        return cls(
            name="Ethylene (C2H4)",
            formula="C2H4",
            molecular_weight_g_per_mol=28.054,
            density_kg_per_m3=568.0,
            reference_temp_c=-103.7,
            purity_pct=99.9,
            freezing_point_c=-169.2,
        )

    @classmethod
    def propylene_pgp(cls) -> PetrochemAssay:
        """Polymer-Grade Propylene (C3H6 — pressurised liquid at 15°C)."""
        return cls(
            name="Propylene (PGP, 99.5%)",
            formula="C3H6",
            molecular_weight_g_per_mol=42.081,
            density_kg_per_m3=522.0,
            reference_temp_c=15.0,
            purity_pct=99.5,
            thermal_expansion_coeff_per_c=0.00290,
            freezing_point_c=-185.2,
        )

    @classmethod
    def benzene(cls) -> PetrochemAssay:
        """Nitration-Grade Benzene (C6H6 — ASTM D2359, 15°C)."""
        return cls(
            name="Benzene (C6H6)",
            formula="C6H6",
            molecular_weight_g_per_mol=78.114,
            density_kg_per_m3=884.2,
            reference_temp_c=15.0,
            purity_pct=99.85,
            thermal_expansion_coeff_per_c=0.00124,
            freezing_point_c=5.53,
        )

    @classmethod
    def toluene(cls) -> PetrochemAssay:
        """Nitration-Grade Toluene (C7H8 — ASTM D841, 15°C)."""
        return cls(
            name="Toluene (C7H8)",
            formula="C7H8",
            molecular_weight_g_per_mol=92.141,
            density_kg_per_m3=871.5,
            reference_temp_c=15.0,
            purity_pct=99.5,
            thermal_expansion_coeff_per_c=0.00107,
            freezing_point_c=-95.0,
        )

    @classmethod
    def paraxylene(cls) -> PetrochemAssay:
        """Paraxylene (PX, C8H10 — polyester chain benchmark, 15°C)."""
        return cls(
            name="Paraxylene (PX)",
            formula="C8H10",
            molecular_weight_g_per_mol=106.167,
            density_kg_per_m3=865.2,
            reference_temp_c=15.0,
            purity_pct=99.7,
            thermal_expansion_coeff_per_c=0.00101,
            freezing_point_c=13.26,
        )

    @classmethod
    def styrene_monomer(cls) -> PetrochemAssay:
        """Styrene Monomer (SM, C8H8 — 15°C)."""
        return cls(
            name="Styrene Monomer (SM)",
            formula="C8H8",
            molecular_weight_g_per_mol=104.152,
            density_kg_per_m3=910.2,
            reference_temp_c=15.0,
            purity_pct=99.8,
            thermal_expansion_coeff_per_c=0.00097,
            freezing_point_c=-30.6,
        )

    @classmethod
    def meg(cls) -> PetrochemAssay:
        """Fibre-Grade Monoethylene Glycol (MEG, C2H6O2 — 15°C)."""
        return cls(
            name="Monoethylene Glycol (MEG)",
            formula="C2H6O2",
            molecular_weight_g_per_mol=62.068,
            density_kg_per_m3=1117.0,
            reference_temp_c=15.0,
            purity_pct=99.8,
            thermal_expansion_coeff_per_c=0.00065,
            freezing_point_c=-12.9,
        )

    @classmethod
    def mixed_xylenes(cls) -> PetrochemAssay:
        """Isomer-Grade Mixed Xylenes (MX, C8H10 — USGC / Asian aromatics benchmark, 15°C)."""
        return cls(
            name="Mixed Xylenes (MX)",
            formula="C8H10",
            molecular_weight_g_per_mol=106.167,
            density_kg_per_m3=870.0,
            reference_temp_c=15.0,
            purity_pct=99.0,
            thermal_expansion_coeff_per_c=0.00099,
            freezing_point_c=-47.4,
        )

    @classmethod
    def butadiene(cls) -> PetrochemAssay:
        """1,3-Butadiene (BD, C4H6 — pressurised liquid at 15°C)."""
        return cls(
            name="1,3-Butadiene (C4H6)",
            formula="C4H6",
            molecular_weight_g_per_mol=54.091,
            density_kg_per_m3=627.0,
            reference_temp_c=15.0,
            purity_pct=99.5,
            thermal_expansion_coeff_per_c=0.00180,
            freezing_point_c=-108.9,
        )

# ---------------------------------------------------------------------------
# Precious Metals & Doré Specification (LBMA / LPPM / COMEX)
# (Gold, Silver, Platinum, Palladium, Mine-Site Doré)
# ---------------------------------------------------------------------------

_GRAMS_PER_TROY_OUNCE: Final[float] = 31.1034768              # Exact international troy ounce (g)
_TROY_OUNCES_PER_KG: Final[float] = 1000.0 / 31.1034768       # 32.15074657 ozt/kg
_GRAMS_PER_TOLA: Final[float] = 11.6638038                    # Exact India/Dubai bullion tola (3/8 ozt)


@validated_assay
@dataclass(frozen=True, slots=True)
class BullionAssay:
    """Immutable fineness and custody specification for precious metals and doré.

    Bridges gross mass (g, kg, ozt), fine troy ounces (fine ozt), refinery payable
    ounces, doré byproduct metals (e.g. Ag in Au doré), standard bar counts, and
    regional $/ozt <-> $/g <-> $/kg pricing parity.
    """

    name: str
    primary_metal: str
    fineness: float
    payable_factor: float = 1.0
    standard_bar_mass_ozt: float | None = None
    byproduct_fineness: Mapping[str, float] | None = None

    def __post_init__(self) -> None:
        if not self.primary_metal or not self.primary_metal.strip():
            raise ValueError("primary_metal must be a non-empty string (e.g. 'Au', 'Ag', 'Pt', 'Pd').")
        if not math.isfinite(self.fineness) or not (0.0 < self.fineness <= 1.0):
            raise ValueError(f"fineness must be finite and in (0.0, 1.0], got {self.fineness}")
        if not math.isfinite(self.payable_factor) or not (0.0 < self.payable_factor <= 1.0):
            raise ValueError(f"payable_factor must be finite and in (0.0, 1.0], got {self.payable_factor}")
        if self.standard_bar_mass_ozt is not None:
            if not math.isfinite(self.standard_bar_mass_ozt) or self.standard_bar_mass_ozt <= 0.0:
                raise ValueError(f"standard_bar_mass_ozt must be strictly positive, got {self.standard_bar_mass_ozt}")
        if self.byproduct_fineness is not None:
            cleaned: dict[str, float] = {}
            for elem, fin in self.byproduct_fineness.items():
                val = float(fin)
                if not math.isfinite(val) or not (0.0 <= val <= 1.0):
                    raise ValueError(f"Byproduct fineness for {elem!r} must be in [0.0, 1.0], got {fin}")
                cleaned[str(elem)] = val
            if self.fineness + sum(cleaned.values()) > 1.0 + 1e-9:
                raise ValueError("fineness plus total byproduct_fineness cannot exceed 1.0.")
            object.__setattr__(self, "byproduct_fineness", MappingProxyType(cleaned))

    def __hash__(self) -> int:
        byprod_key = (
            tuple(sorted(self.byproduct_fineness.items())) if self.byproduct_fineness is not None else None
        )
        return hash((
            self.name,
            self.primary_metal,
            self.fineness,
            self.payable_factor,
            self.standard_bar_mass_ozt,
            byprod_key,
        ))

    # --- Fineness & Karat Properties ---
    @property
    def millesimal_fineness(self) -> float:
        """Fineness expressed in parts per thousand (e.g. 995.0, 999.9)."""
        return (self.fineness * 1000.0)

    @property
    def karat(self) -> float | None:
        """Karat purity (24K = 1.0) for gold assays; None for non-gold metals."""
        if self.primary_metal.strip().upper() not in {"AU", "GOLD"}:
            return None
        return (self.fineness * 24.0)

    @property
    def effective_payable_fineness(self) -> float:
        """Net payable fineness after applying refinery payable_factor."""
        return (self.fineness * self.payable_factor)

    def as_context(self, *, basis: str = "fine", context_name: str | None = None) -> pint.Context:
        """Returns a named context; use explicit methods for fine and payable masses."""
        return create_bullion_context(self, basis=basis, name=context_name)

    # --- Physical Mass Conversions (Gross <-> Fine <-> Payable <-> Bars) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.troy_ounce)

    def gross_to_fine(self, gross_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross bullion/doré mass (e.g. kg, g, ozt) to fine troy ounces."""
        ozt_q = self._validate_mass(gross_mass)
        scale = _scale_factor(ozt_q.magnitude, self.fineness)
        return ureg.Quantity(ozt_q.magnitude * scale, ureg.troy_ounce)

    def fine_to_gross(self, fine_mass: pint.Quantity) -> pint.Quantity:
        """Converts fine metal mass (e.g. ozt, kg, g) to required gross troy ounces."""
        ozt_q = self._validate_mass(fine_mass)
        scale = _scale_factor(ozt_q.magnitude, self.fineness)
        return ureg.Quantity(ozt_q.magnitude / scale, ureg.troy_ounce)

    def gross_to_payable_fine(self, gross_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross bullion/doré mass to refinery payable fine troy ounces."""
        ozt_q = self._validate_mass(gross_mass)
        factor = self.fineness * self.payable_factor
        scale = _scale_factor(ozt_q.magnitude, factor)
        return ureg.Quantity(ozt_q.magnitude * scale, ureg.troy_ounce)

    def payable_fine_to_gross(self, payable_mass: pint.Quantity) -> pint.Quantity:
        """Converts target payable fine metal mass to required gross troy ounces."""
        ozt_q = self._validate_mass(payable_mass)
        factor = self.fineness * self.payable_factor
        scale = _scale_factor(ozt_q.magnitude, factor)
        return ureg.Quantity(ozt_q.magnitude / scale, ureg.troy_ounce)

    def gross_to_byproduct_fine(self, gross_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Computes contained fine troy ounces of a doré byproduct metal (e.g. 'Ag')."""
        if self.byproduct_fineness is None or element not in self.byproduct_fineness:
            raise KeyError(f"Byproduct {element!r} is not defined on {self.name!r}.")
        ozt_q = self._validate_mass(gross_mass)
        scale = _scale_factor(ozt_q.magnitude, self.byproduct_fineness[element])
        return ureg.Quantity(ozt_q.magnitude * scale, ureg.troy_ounce)

    def gross_to_bar_count(self, gross_mass: pint.Quantity) -> float:
        """Converts gross mass to nominal standard bar count."""
        if self.standard_bar_mass_ozt is None:
            raise ValueError(f"{self.name!r} does not define standard_bar_mass_ozt.")
        ozt = float(self._validate_mass(gross_mass).magnitude)
        return (ozt / self.standard_bar_mass_ozt)

    def bar_count_to_gross(self, bar_count: float) -> pint.Quantity:
        """Converts a nominal standard bar count to gross troy ounces."""
        if self.standard_bar_mass_ozt is None:
            raise ValueError(f"{self.name!r} does not define standard_bar_mass_ozt.")
        n = float(bar_count)
        if not math.isfinite(n) or n < 0.0:
            raise ValueError(f"bar_count must be finite and non-negative, got {bar_count}")
        return ureg.Quantity(n * self.standard_bar_mass_ozt, ureg.troy_ounce)

    # --- Commercial Pricing Parity ($/ozt <-> $/g <-> $/kg & Fine <-> Gross) ---
    @staticmethod
    def ozt_price_to_gram_price(price_per_ozt: float) -> float:
        """Converts London/COMEX $/ozt quote to Shanghai/Dubai $/g equivalent."""
        p = float(price_per_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_ozt must be finite, got {price_per_ozt}")
        return (p / _GRAMS_PER_TROY_OUNCE)

    @staticmethod
    def gram_price_to_ozt_price(price_per_gram: float) -> float:
        """Converts $/g quote to $/ozt equivalent."""
        p = float(price_per_gram)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gram must be finite, got {price_per_gram}")
        return (p * _GRAMS_PER_TROY_OUNCE)

    @staticmethod
    def ozt_price_to_kg_price(price_per_ozt: float) -> float:
        """Converts $/ozt quote to $/kg kilobar equivalent."""
        p = float(price_per_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_ozt must be finite, got {price_per_ozt}")
        return (p * _TROY_OUNCES_PER_KG)

    @staticmethod
    def kg_price_to_ozt_price(price_per_kg: float) -> float:
        """Converts $/kg quote to $/ozt equivalent."""
        p = float(price_per_kg)
        if not math.isfinite(p):
            raise ValueError(f"price_per_kg must be finite, got {price_per_kg}")
        return (p / _TROY_OUNCES_PER_KG)

    def fine_ozt_price_to_gross_ozt_price(self, price_per_fine_ozt: float) -> float:
        """Converts a $/fine ozt spot quote to implied $/gross ozt value for this assay."""
        p = float(price_per_fine_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_fine_ozt must be finite, got {price_per_fine_ozt}")
        return (p * self.fineness * self.payable_factor)

    def gross_ozt_price_to_fine_ozt_price(self, price_per_gross_ozt: float) -> float:
        """Converts a $/gross ozt doré/scrap quote to implied $/fine ozt spot parity."""
        p = float(price_per_gross_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gross_ozt must be finite, got {price_per_gross_ozt}")
        return (p / (self.fineness * self.payable_factor))

    @staticmethod
    def ozt_price_to_tola_price(price_per_ozt: float) -> float:
        """Converts London/COMEX $/ozt quote to India/Dubai $/tola (11.6638038 g) equivalent."""
        p = float(price_per_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_ozt must be finite, got {price_per_ozt}")
        return (p * (_GRAMS_PER_TOLA / _GRAMS_PER_TROY_OUNCE))

    @staticmethod
    def tola_price_to_ozt_price(price_per_tola: float) -> float:
        """Converts India/Dubai $/tola quote to $/ozt equivalent."""
        p = float(price_per_tola)
        if not math.isfinite(p):
            raise ValueError(f"price_per_tola must be finite, got {price_per_tola}")
        return (p * (_GRAMS_PER_TROY_OUNCE / _GRAMS_PER_TOLA))

    def blend_with(
        self,
        other: BullionAssay,
        *,
        self_gross_mass_fraction: float,
        blend_name: str | None = None,
    ) -> BullionAssay:
        """Blends two bullion or doré melts of the same primary metal by gross mass."""
        if not isinstance(other, BullionAssay):
            raise TypeError(f"Expected BullionAssay, got {type(other).__name__}")
        if self.primary_metal != other.primary_metal:
            raise ValueError(
                f"Cannot blend bullion with different primary metals: {self.primary_metal!r} vs {other.primary_metal!r}"
            )
        if not math.isfinite(self_gross_mass_fraction) or not (0.0 <= self_gross_mass_fraction <= 1.0):
            raise ValueError("self_gross_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_gross_mass_fraction)
        w2 = 1.0 - w1

        f1 = w1 * self.fineness
        f2 = w2 * other.fineness
        total_f = f1 + f2
        blended_payable = ((self.payable_factor * f1) + (other.payable_factor * f2)) / total_f

        b_byprod: dict[str, float] | None = None
        if self.byproduct_fineness is not None or other.byproduct_fineness is not None:
            d1 = self.byproduct_fineness or {}
            d2 = other.byproduct_fineness or {}
            b_byprod = {
                k: ((d1.get(k, 0.0) * w1) + (d2.get(k, 0.0) * w2))
                for k in sorted(set(d1) | set(d2))
            }

        return BullionAssay(
            name=blend_name or f"{self.name}/{other.name} Melt",
            primary_metal=self.primary_metal,
            fineness=((self.fineness * w1) + (other.fineness * w2)),
            payable_factor=(blended_payable),
            standard_bar_mass_ozt=None,
            byproduct_fineness=b_byprod,
        )

    # --- Global Precious Metals & Doré Benchmark Presets & Factories ---
    @classmethod
    def from_karat(
        cls,
        karat: float,
        *,
        payable_factor: float = 1.0,
        name: str | None = None,
    ) -> BullionAssay:
        """Constructs a gold bullion/scrap assay directly from karat purity (e.g. 22, 18, 14)."""
        k = float(karat)
        if not math.isfinite(k) or not (0.0 < k <= 24.0):
            raise ValueError(f"karat must be finite and in (0, 24], got {karat}")
        return cls(
            name=name or f"Gold ({k:g}K)",
            primary_metal="Au",
            fineness=(k / 24.0),
            payable_factor=payable_factor,
        )

    @classmethod
    def lbma_gold_good_delivery(cls) -> BullionAssay:
        """LBMA Good Delivery Gold Bar (995.0 min fineness, 400 ozt nominal)."""
        return cls(
            name="LBMA Good Delivery Gold",
            primary_metal="Au",
            fineness=0.995,
            standard_bar_mass_ozt=400.0,
        )

    @classmethod
    def comex_gold_100oz(cls) -> BullionAssay:
        """COMEX Good Delivery Gold Bar (995.0 min fineness, 100 ozt nominal)."""
        return cls(
            name="COMEX Good Delivery Gold (100 ozt)",
            primary_metal="Au",
            fineness=0.995,
            standard_bar_mass_ozt=100.0,
        )

    @classmethod
    def gold_9999_kilobar(cls) -> BullionAssay:
        """Four-Nines Gold Kilobar (999.9 fineness, 1 kg / 32.1507 ozt — SGE / Dubai / Swiss)."""
        return cls(
            name="Gold 999.9 Kilobar",
            primary_metal="Au",
            fineness=0.9999,
            standard_bar_mass_ozt=(_TROY_OUNCES_PER_KG),
        )

    @classmethod
    def lbma_silver_good_delivery(cls) -> BullionAssay:
        """LBMA Good Delivery Silver Bar (999.0 min fineness, 1,000 ozt nominal)."""
        return cls(
            name="LBMA Good Delivery Silver",
            primary_metal="Ag",
            fineness=0.999,
            standard_bar_mass_ozt=1000.0,
        )

    @classmethod
    def lppm_platinum_good_delivery(cls) -> BullionAssay:
        """LPPM Good Delivery Platinum Plate/Ingot (999.5 min fineness)."""
        return cls(
            name="LPPM Good Delivery Platinum",
            primary_metal="Pt",
            fineness=0.9995,
        )

    @classmethod
    def lppm_palladium_good_delivery(cls) -> BullionAssay:
        """LPPM Good Delivery Palladium Plate/Ingot (999.5 min fineness)."""
        return cls(
            name="LPPM Good Delivery Palladium",
            primary_metal="Pd",
            fineness=0.9995,
        )

    @classmethod
    def gold_dore(cls) -> BullionAssay:
        """Typical Mine-Site Gold Doré Bar (85% Au, 10% Ag, 99.9% refinery payable)."""
        return cls(
            name="Gold Doré (85% Au / 10% Ag)",
            primary_metal="Au",
            fineness=0.85,
            payable_factor=0.999,
            byproduct_fineness={"Ag": 0.10},
        )

# ---------------------------------------------------------------------------
# Fertilisers & Nutrient Equivalents Specification
# (Urea, UAN 28%/32%, DAP, MAP, Potash MOP, Ammonium Sulphate)
# ---------------------------------------------------------------------------

# Exact IUPAC elemental-to-oxide nutrient stoichiometry
_P_IN_P2O5: Final[float] = 0.43642          # 2 * M(P) / M(P2O5) = 61.948 / 141.945
_K_IN_K2O: Final[float] = 0.83015           # 2 * M(K) / M(K2O)  = 78.197 / 94.196
_MT_PER_SHORT_TON: Final[float] = 0.90718474  # Exact: 2000 lb in metric tonnes


@validated_assay
@dataclass(frozen=True, slots=True)
class FertilizerAssay:
    """Immutable nutrient-grade and physical specification for fertilisers.

    Bridges physical product mass (mt, short_ton, lb), nutrient-equivalent mass
    (N, P2O5, P, K2O, K, S), liquid/bulk volume (gal, m^3), NOLA $/st <-> $/mt
    pricing parity, and NPK blending.
    """

    name: str
    n_pct: float = 0.0
    p2o5_pct: float = 0.0
    k2o_pct: float = 0.0
    s_pct: float = 0.0
    density_kg_per_m3: float | None = None
    is_liquid: bool = False

    def __post_init__(self) -> None:
        for field_name, val in (
            ("n_pct", self.n_pct),
            ("p2o5_pct", self.p2o5_pct),
            ("k2o_pct", self.k2o_pct),
            ("s_pct", self.s_pct),
        ):
            if not math.isfinite(val) or not (0.0 <= val <= 100.0):
                raise ValueError(f"{field_name} must be finite and in [0, 100], got {val}")
        if (self.n_pct + self.p2o5_pct + self.k2o_pct + self.s_pct) <= 0.0:
            raise ValueError("At least one nutrient percentage (N, P2O5, K2O, S) must be strictly positive.")
        if (self.n_pct + self.p2o5_pct + self.k2o_pct + self.s_pct) > 100.0 + 1e-9:
            raise ValueError("Sum of nutrient percentages cannot exceed 100%.")
        if self.density_kg_per_m3 is not None:
            if not math.isfinite(self.density_kg_per_m3) or not (300.0 <= self.density_kg_per_m3 <= 2500.0):
                raise ValueError(
                    f"density_kg_per_m3 must be finite and in [300, 2500] kg/m^3, got {self.density_kg_per_m3}"
                )

    # --- Nutrient Grade & Volumetric Parity Properties ---
    @property
    def npk_grade(self) -> str:
        """Standard N-P2O5-K2O (and optional -S) grade string (e.g. '46-0-0', '18-46-0')."""
        base = f"{self.n_pct:g}-{self.p2o5_pct:g}-{self.k2o_pct:g}"
        return f"{base}-{self.s_pct:g}S" if self.s_pct > 0.0 else base

    @property
    def elemental_p_pct(self) -> float:
        """Elemental phosphorus content (wt% P) converted from P2O5."""
        return (self.p2o5_pct * _P_IN_P2O5)

    @property
    def elemental_k_pct(self) -> float:
        """Elemental potassium content (wt% K) converted from K2O."""
        return (self.k2o_pct * _K_IN_K2O)

    @property
    def gal_per_mt(self) -> float | None:
        """US gallons per metric tonne (primarily for liquid UAN), or None if density unset."""
        if self.density_kg_per_m3 is None:
            return None
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def lb_per_gal(self) -> float | None:
        """US pounds per US gallon (US liquid fertiliser standard), or None if density unset."""
        if self.gal_per_mt is None:
            return None
        return (_LB_PER_MT / self.gal_per_mt)

    @property
    def lb_n_per_gal(self) -> float | None:
        """Pounds of elemental Nitrogen per US gallon (UAN application benchmark)."""
        if self.lb_per_gal is None:
            return None
        return (self.lb_per_gal * (self.n_pct / 100.0))

    @property
    def stowage_factor_m3_per_mt(self) -> float | None:
        """Bulk stowage factor (m^3 per mt), or None if density unset."""
        if self.density_kg_per_m3 is None:
            return None
        return (1000.0 / self.density_kg_per_m3)

    def _nutrient_fraction(self, nutrient: str) -> float:
        key = nutrient.strip().upper()
        mapping = {
            "N": self.n_pct / 100.0,
            "P2O5": self.p2o5_pct / 100.0,
            "P": (self.p2o5_pct * _P_IN_P2O5) / 100.0,
            "K2O": self.k2o_pct / 100.0,
            "K": (self.k2o_pct * _K_IN_K2O) / 100.0,
            "S": self.s_pct / 100.0,
        }
        if key not in mapping:
            raise ValueError(f"nutrient must be one of {sorted(mapping)}, got {nutrient!r}")
        frac = mapping[key]
        if frac <= 0.0:
            raise ValueError(f"{self.name!r} contains 0% {key}.")
        return frac

    def as_context(self, *, nutrient: str = "N", context_name: str | None = None) -> pint.Context:
        """Bridges volume and product mass; use explicit methods for nutrient masses."""
        return create_fertilizer_context(self, nutrient=nutrient, name=context_name)

    # --- Physical Mass & Volume Conversions ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def product_to_nutrient(self, product_mass: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts physical fertiliser mass (e.g. mt, short_ton, lb) to nutrient mass (metric_ton)."""
        mt_q = self._validate_mass(product_mass)
        frac = self._nutrient_fraction(nutrient)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def nutrient_to_product(self, nutrient_mass: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts target nutrient mass (e.g. mt N, mt P2O5, mt K2O) to required product mass (metric_ton)."""
        mt_q = self._validate_mass(nutrient_mass)
        frac = self._nutrient_fraction(nutrient)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid or bulk stowage volume (e.g. gal, m^3) to mass (metric_ton)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts mass (e.g. mt, short_ton) to liquid or bulk stowage volume (m^3)."""
        if self.density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def volume_to_nutrient(self, volume: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts liquid/bulk volume (e.g. gal, m^3) directly to nutrient mass (metric_ton)."""
        return self.product_to_nutrient(self.volume_to_mass(volume), nutrient)

    def nutrient_to_volume(self, nutrient_mass: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts target nutrient mass (e.g. mt N, lb N) directly to liquid/bulk volume (m^3)."""
        return self.mass_to_volume(self.nutrient_to_product(nutrient_mass, nutrient))

    def equivalent_product_mass(
        self,
        product_mass: pint.Quantity,
        other: FertilizerAssay,
        nutrient: str,
    ) -> pint.Quantity:
        """Computes the mass of `other` required to deliver the same `nutrient` tonnage as `product_mass`."""
        if not isinstance(other, FertilizerAssay):
            raise TypeError(f"Expected FertilizerAssay, got {type(other).__name__}")
        target_nutrient = self.product_to_nutrient(product_mass, nutrient)
        return other.nutrient_to_product(target_nutrient, nutrient)

    # --- Commercial Pricing Parity ($/st <-> $/mt, $/product <-> $/nutrient, $/gal <-> $/mt) ---
    @staticmethod
    def st_price_to_mt_price(price_per_st: float) -> float:
        """Converts US NOLA barge $/short_ton quote to international $/metric_ton."""
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / _MT_PER_SHORT_TON)

    @staticmethod
    def mt_price_to_st_price(price_per_mt: float) -> float:
        """Converts international $/metric_ton quote to US NOLA barge $/short_ton."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p * _MT_PER_SHORT_TON)

    def product_price_to_nutrient_price(self, price_per_mt_product: float, nutrient: str) -> float:
        """Converts a $/mt product quote to implied $/mt contained nutrient (e.g. $/mt N)."""
        p = float(price_per_mt_product)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_product must be finite, got {price_per_mt_product}")
        return (p / self._nutrient_fraction(nutrient))

    def nutrient_price_to_product_price(self, price_per_mt_nutrient: float, nutrient: str) -> float:
        """Converts a $/mt nutrient quote to implied $/mt product equivalent."""
        p = float(price_per_mt_nutrient)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_nutrient must be finite, got {price_per_mt_nutrient}")
        return (p * self._nutrient_fraction(nutrient))

    def gal_price_to_mt_price(self, price_per_gal: float) -> float:
        """Converts liquid fertiliser $/gal quote to $/mt equivalent."""
        if self.gal_per_mt is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        p = float(price_per_gal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gal must be finite, got {price_per_gal}")
        return (p * self.gal_per_mt)

    def mt_price_to_gal_price(self, price_per_mt: float) -> float:
        """Converts liquid fertiliser $/mt quote to $/gal equivalent."""
        if self.gal_per_mt is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.gal_per_mt)

    def st_price_to_lb_nutrient_price(self, price_per_st: float, nutrient: str) -> float:
        """Converts a US $/short_ton product quote to $/lb contained nutrient (e.g. $/lb N)."""
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / (2000.0 * self._nutrient_fraction(nutrient)))

    def st_price_to_unit_price(self, price_per_st: float, nutrient: str) -> float:
        """Converts a US $/short_ton product quote to $/short-ton unit (1 unit = 20 lb nutrient)."""
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / (self._nutrient_fraction(nutrient) * 100.0))

    def blend_with(
        self,
        other: FertilizerAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> FertilizerAssay:
        """Computes a mass-weighted NPK/S fertiliser blend (e.g. bulk-blended NPK)."""
        if not isinstance(other, FertilizerAssay):
            raise TypeError(f"Expected FertilizerAssay, got {type(other).__name__}")
        if self.is_liquid != other.is_liquid:
            raise ValueError("Cannot directly blend liquid and solid fertilisers.")
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1

        b_rho: float | None = None
        if self.density_kg_per_m3 is not None and other.density_kg_per_m3 is not None:
            vol_per_kg = (w1 / self.density_kg_per_m3) + (w2 / other.density_kg_per_m3)
            b_rho = (1.0 / vol_per_kg)

        return FertilizerAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            n_pct=((self.n_pct * w1) + (other.n_pct * w2)),
            p2o5_pct=((self.p2o5_pct * w1) + (other.p2o5_pct * w2)),
            k2o_pct=((self.k2o_pct * w1) + (other.k2o_pct * w2)),
            s_pct=((self.s_pct * w1) + (other.s_pct * w2)),
            density_kg_per_m3=b_rho,
            is_liquid=self.is_liquid,
        )

    # --- Global Fertiliser Benchmark Presets ---
    @classmethod
    def urea_granular(cls) -> FertilizerAssay:
        """Granular Urea (46-0-0 — Middle East FOB / NOLA Barge / Brazil CFR)."""
        return cls(
            name="Granular Urea (46-0-0)",
            n_pct=46.0,
            density_kg_per_m3=750.0,
            is_liquid=False,
        )

    @classmethod
    def uan_32(cls) -> FertilizerAssay:
        """Urea Ammonium Nitrate 32% N liquid solution (32-0-0, 11.06 lb/gal)."""
        return cls(
            name="UAN 32% Solution",
            n_pct=32.0,
            density_kg_per_m3=1325.0,
            is_liquid=True,
        )

    @classmethod
    def uan_28(cls) -> FertilizerAssay:
        """Urea Ammonium Nitrate 28% N liquid solution (28-0-0, 10.67 lb/gal)."""
        return cls(
            name="UAN 28% Solution",
            n_pct=28.0,
            density_kg_per_m3=1280.0,
            is_liquid=True,
        )

    @classmethod
    def dap(cls) -> FertilizerAssay:
        """Diammonium Phosphate (DAP, 18-46-0 — Tampa / Morocco / India CFR)."""
        return cls(
            name="DAP (18-46-0)",
            n_pct=18.0,
            p2o5_pct=46.0,
            density_kg_per_m3=960.0,
            is_liquid=False,
        )

    @classmethod
    def map_fertilizer(cls) -> FertilizerAssay:
        """Monoammonium Phosphate (MAP, 11-52-0 — Brazil CFR / NOLA)."""
        return cls(
            name="MAP (11-52-0)",
            n_pct=11.0,
            p2o5_pct=52.0,
            density_kg_per_m3=980.0,
            is_liquid=False,
        )

    @classmethod
    def potash_mop(cls) -> FertilizerAssay:
        """Muriate of Potash (MOP / KCl, 0-0-60 — Vancouver / Brazil CFR / SE Asia)."""
        return cls(
            name="Potash MOP (0-0-60)",
            k2o_pct=60.0,
            density_kg_per_m3=1050.0,
            is_liquid=False,
        )

    @classmethod
    def ammonium_sulphate(cls) -> FertilizerAssay:
        """Granular Ammonium Sulphate (21-0-0-24S)."""
        return cls(
            name="Ammonium Sulphate (21-0-0-24S)",
            n_pct=21.0,
            s_pct=24.0,
            density_kg_per_m3=960.0,
            is_liquid=False,
        )

    @classmethod
    def anhydrous_ammonia(cls) -> FertilizerAssay:
        """Anhydrous Ammonia agricultural grade (82-0-0, refrigerated liquid at -33°C)."""
        return cls(
            name="Anhydrous Ammonia (82-0-0)",
            n_pct=82.0,
            density_kg_per_m3=682.0,
            is_liquid=True,
        )

    @classmethod
    def tsp(cls) -> FertilizerAssay:
        """Triple Superphosphate (TSP, 0-46-0 — Morocco / Brazil CFR)."""
        return cls(
            name="TSP (0-46-0)",
            p2o5_pct=46.0,
            density_kg_per_m3=1040.0,
            is_liquid=False,
        )

    @classmethod
    def can_27(cls) -> FertilizerAssay:
        """Calcium Ammonium Nitrate (CAN 27% N — European benchmark)."""
        return cls(
            name="CAN (27-0-0)",
            n_pct=27.0,
            density_kg_per_m3=1000.0,
            is_liquid=False,
        )

    @classmethod
    def potash_sop(cls) -> FertilizerAssay:
        """Sulphate of Potash (SOP / K2SO4, 0-0-50-18S — specialty chloride-free potash)."""
        return cls(
            name="Potash SOP (0-0-50-18S)",
            k2o_pct=50.0,
            s_pct=18.0,
            density_kg_per_m3=1250.0,
            is_liquid=False,
        )

# ---------------------------------------------------------------------------
# Paving Bitumen (Asphalt) Specification (ASTM D4311 / EN 12591)
# (Pen 60/70, Pen 80/100, PG 64-22, Oxidised Roofing Bitumen)
# ---------------------------------------------------------------------------

_ASTM_D4311_GROUP0_ALPHA_PER_C: Final[float] = 0.00062  # Standard Group 0 bitumen thermal expansion (/°C)


@validated_assay
@dataclass(frozen=True, slots=True)
class BitumenAssay:
    """Immutable physical and rheological specification for paving bitumen (asphalt).

    Bridges standard-temperature and hot-loaded liquid volume (m^3, bbl, gal)
    with mass (mt, short_ton) via ASTM D4311 thermal volume correction, plus
    log-linear penetration blending and regional pricing parity.
    """

    name: str
    density_kg_per_m3: float
    reference_temp_c: float = 15.0
    thermal_expansion_coeff_per_c: float = _ASTM_D4311_GROUP0_ALPHA_PER_C
    penetration_dmm: float | None = None
    softening_point_c: float | None = None
    sulphur_pct: float | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.density_kg_per_m3) or not (950.0 <= self.density_kg_per_m3 <= 1150.0):
            raise ValueError(
                f"density_kg_per_m3 for bitumen must be in [950, 1150] kg/m^3, got {self.density_kg_per_m3}"
            )
        if not math.isfinite(self.reference_temp_c) or not (0.0 <= self.reference_temp_c <= 40.0):
            raise ValueError(f"reference_temp_c must be in [0, 40] °C, got {self.reference_temp_c}")
        if not math.isfinite(self.thermal_expansion_coeff_per_c) or self.thermal_expansion_coeff_per_c <= 0.0:
            raise ValueError(
                f"thermal_expansion_coeff_per_c must be strictly positive, got {self.thermal_expansion_coeff_per_c}"
            )
        if self.penetration_dmm is not None:
            if not math.isfinite(self.penetration_dmm) or self.penetration_dmm <= 0.0:
                raise ValueError(f"penetration_dmm must be strictly positive, got {self.penetration_dmm}")
        if self.softening_point_c is not None:
            if not math.isfinite(self.softening_point_c) or not (20.0 <= self.softening_point_c <= 180.0):
                raise ValueError(f"softening_point_c must be in [20, 180] °C, got {self.softening_point_c}")
        if self.sulphur_pct is not None:
            if not math.isfinite(self.sulphur_pct) or not (0.0 <= self.sulphur_pct <= 12.0):
                raise ValueError(f"sulphur_pct must be in [0, 12], got {self.sulphur_pct}")

    # --- Volumetric & Gravimetric Parity Properties ---
    @property
    def specific_gravity(self) -> float:
        """Specific gravity relative to water at 15°C (999.102 kg/m^3)."""
        return (self.density_kg_per_m3 / 999.102)

    @property
    def bbl_per_mt(self) -> float:
        """Standard barrels per metric tonne at reference_temp_c."""
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_BBL / 1000.0)))

    @property
    def gal_per_mt(self) -> float:
        """Standard US gallons per metric tonne at reference_temp_c."""
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging standard bitumen [volume] and [mass]."""
        return create_bitumen_context(self, name=context_name)

    # --- ASTM D4311 Thermal Volume Correction & Custody Transfer ---
    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float | None = None,
    ) -> float:
        """Computes ASTM D4311 thermal Volume Correction Factor (VCF) from hot observed temp to reference."""
        t_ref = self.reference_temp_c if ref_temp_c is None else float(ref_temp_c)
        if isinstance(observed_temp, pint.Quantity):
            t_obs_c = float(observed_temp.to(ureg.degC).magnitude)
        else:
            t_obs_c = float(observed_temp)
        if not math.isfinite(t_obs_c) or not math.isfinite(t_ref):
            raise ValueError("Temperatures must be finite.")
        alpha = self.thermal_expansion_coeff_per_c
        dt = t_obs_c - t_ref
        return (math.exp(-alpha * dt * (1.0 + 0.8 * alpha * dt)))

    def density_at_temp(self, observed_temp: float | pint.Quantity) -> float:
        """Computes hot-loaded liquid density (kg/m^3) at observed_temp via ASTM D4311."""
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        return (self.density_kg_per_m3 * vcf)

    def volume_to_mass(self, std_volume: pint.Quantity) -> pint.Quantity:
        """Converts standard-temperature bitumen volume (m^3, bbl, gal) to mass (metric_ton)."""
        if not isinstance(std_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(std_volume).__name__}")
        m3_q = std_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts bitumen mass (mt, short_ton) to standard-temperature volume (m^3)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def observed_volume_to_mass(
        self,
        observed_volume: pint.Quantity,
        observed_temp: float | pint.Quantity,
    ) -> pint.Quantity:
        """Converts hot-loaded observed volume (e.g. at 150°C) directly to mass (metric_ton)."""
        if not isinstance(observed_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(observed_volume).__name__}")
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        scale = _scale_factor(observed_volume.magnitude, vcf)
        std_vol = ureg.Quantity(observed_volume.magnitude * scale, observed_volume.units)
        return self.volume_to_mass(std_vol)

    def mass_to_observed_volume(
        self,
        mass: pint.Quantity,
        observed_temp: float | pint.Quantity,
    ) -> pint.Quantity:
        """Converts mass (mt) to hot-loaded tank/vessel volume (m^3) at observed_temp."""
        std_m3 = self.mass_to_volume(mass)
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        scale = _scale_factor(std_m3.magnitude, vcf)
        return ureg.Quantity(std_m3.magnitude / scale, ureg.meter ** 3)

    # --- Commercial Pricing Parity ($/st <-> $/mt, $/gal <-> $/mt) ---
    @staticmethod
    def st_price_to_mt_price(price_per_st: float) -> float:
        """Converts US $/short_ton asphalt quote to international $/metric_ton."""
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / _MT_PER_SHORT_TON)

    @staticmethod
    def mt_price_to_st_price(price_per_mt: float) -> float:
        """Converts international $/metric_ton quote to US $/short_ton."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p * _MT_PER_SHORT_TON)

    def gal_price_to_mt_price(self, price_per_gal: float) -> float:
        """Converts $/gal bitumen quote to $/mt equivalent at reference temperature."""
        p = float(price_per_gal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gal must be finite, got {price_per_gal}")
        return (p * self.gal_per_mt)

    def mt_price_to_gal_price(self, price_per_mt: float) -> float:
        """Converts $/mt bitumen quote to $/gal equivalent at reference temperature."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.gal_per_mt)

    def bbl_price_to_mt_price(self, price_per_bbl: float) -> float:
        """Converts $/bbl bitumen quote to $/mt equivalent (crude/HSFO crack-spread parity)."""
        p = float(price_per_bbl)
        if not math.isfinite(p):
            raise ValueError(f"price_per_bbl must be finite, got {price_per_bbl}")
        return (p * self.bbl_per_mt)

    def mt_price_to_bbl_price(self, price_per_mt: float) -> float:
        """Converts $/mt bitumen quote to $/bbl equivalent."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.bbl_per_mt)

    def blend_with(
        self,
        other: BitumenAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> BitumenAssay:
        """Blends two bitumen grades by mass fraction.

        - Density: harmonic mean (additive standard volumes)
        - Penetration (dmm): log-linear mass-weighted (ASTM/Chevron blending rule)
        - Softening point & sulphur: linear mass-weighted
        """
        if not isinstance(other, BitumenAssay):
            raise TypeError(f"Expected BitumenAssay, got {type(other).__name__}")
        if abs(self.reference_temp_c - other.reference_temp_c) > 1e-6:
            raise ValueError(
                f"Cannot blend bitumens with different reference temperatures: "
                f"{self.reference_temp_c}°C vs {other.reference_temp_c}°C."
            )
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1

        vol_per_kg = (w1 / self.density_kg_per_m3) + (w2 / other.density_kg_per_m3)
        b_rho = (1.0 / vol_per_kg)

        b_pen: float | None = None
        if self.penetration_dmm is not None and other.penetration_dmm is not None:
            log_pen = (w1 * math.log10(self.penetration_dmm)) + (w2 * math.log10(other.penetration_dmm))
            b_pen = (10.0 ** log_pen)

        b_sp = (
            ((self.softening_point_c * w1) + (other.softening_point_c * w2))
            if (self.softening_point_c is not None and other.softening_point_c is not None)
            else None
        )
        b_s = (
            ((self.sulphur_pct * w1) + (other.sulphur_pct * w2))
            if (self.sulphur_pct is not None and other.sulphur_pct is not None)
            else None
        )

        return BitumenAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            density_kg_per_m3=b_rho,
            reference_temp_c=self.reference_temp_c,
            thermal_expansion_coeff_per_c=((self.thermal_expansion_coeff_per_c * w1) + (other.thermal_expansion_coeff_per_c * w2)),
            penetration_dmm=b_pen,
            softening_point_c=b_sp,
            sulphur_pct=b_s,
        )

    # --- Global Bitumen Benchmark Presets ---
    @classmethod
    def pen_60_70(cls) -> BitumenAssay:
        """Standard Paving Bitumen 60/70 (Singapore FOB / Med / Middle East benchmark)."""
        return cls(
            name="Bitumen Pen 60/70",
            density_kg_per_m3=1030.0,
            penetration_dmm=65.0,
            softening_point_c=49.0,
            sulphur_pct=4.0,
        )

    @classmethod
    def pen_80_100(cls) -> BitumenAssay:
        """Softer Paving Bitumen 80/100 (cooler-climate road grade)."""
        return cls(
            name="Bitumen Pen 80/100",
            density_kg_per_m3=1025.0,
            penetration_dmm=90.0,
            softening_point_c=45.0,
            sulphur_pct=3.8,
        )

    @classmethod
    def pg_64_22(cls) -> BitumenAssay:
        """US Superpave Performance Grade PG 64-22 (primary North American paving binder)."""
        return cls(
            name="Asphalt PG 64-22",
            density_kg_per_m3=1032.0,
            penetration_dmm=60.0,
            softening_point_c=50.0,
            sulphur_pct=4.2,
        )

    @classmethod
    def oxidised_roofing(cls) -> BitumenAssay:
        """Oxidised / Blown Industrial Roofing Bitumen (95/25 grade)."""
        return cls(
            name="Oxidised Bitumen (95/25)",
            density_kg_per_m3=1045.0,
            penetration_dmm=25.0,
            softening_point_c=95.0,
            sulphur_pct=4.5,
        )


# ---------------------------------------------------------------------------
# Petroleum Coke (Petcoke) Calorific & Proximate Specification
# (USGC Fuel-Grade 6.5%S / 4.5%S, Calcined Anode-Grade, Shot Coke)
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class PetcokeAssay:
    """Immutable calorific, proximate, and bulk-stowage specification for petroleum coke.

    Covers fuel-grade green petcoke (sponge and shot coke traded for cement kilns
    and power generation on WMT / kcal/kg NAR/GAR and $/MMBtu parity) as well as
    calcined anode-grade petcoke (traded for aluminium smelting on DMT and dry-basis
    Fixed Carbon / Sulphur specs).
    """

    name: str
    lhv_kcal_per_kg_ar: float
    hhv_kcal_per_kg_ar: float
    moisture_pct: float = 8.0
    sulphur_db_pct: float = 5.5
    volatile_matter_db_pct: float = 11.0
    ash_db_pct: float = 0.5
    hgi: float | None = None
    bulk_density_kg_per_m3: float = 800.0
    is_calcined: bool = False

    def __post_init__(self) -> None:
        if not math.isfinite(self.lhv_kcal_per_kg_ar) or self.lhv_kcal_per_kg_ar <= 0.0:
            raise ValueError(f"lhv_kcal_per_kg_ar must be finite and positive, got {self.lhv_kcal_per_kg_ar}")
        if not math.isfinite(self.hhv_kcal_per_kg_ar) or self.hhv_kcal_per_kg_ar < self.lhv_kcal_per_kg_ar:
            raise ValueError(
                f"hhv_kcal_per_kg_ar ({self.hhv_kcal_per_kg_ar}) must be >= lhv_kcal_per_kg_ar ({self.lhv_kcal_per_kg_ar})"
            )
        if not math.isfinite(self.moisture_pct) or not (0.0 <= self.moisture_pct < 30.0):
            raise ValueError(f"moisture_pct must be finite and in [0, 30), got {self.moisture_pct}")
        if not math.isfinite(self.sulphur_db_pct) or not (0.0 <= self.sulphur_db_pct <= 15.0):
            raise ValueError(f"sulphur_db_pct must be finite and in [0, 15], got {self.sulphur_db_pct}")
        if not math.isfinite(self.volatile_matter_db_pct) or not (0.0 <= self.volatile_matter_db_pct <= 25.0):
            raise ValueError(
                f"volatile_matter_db_pct must be finite and in [0, 25], got {self.volatile_matter_db_pct}"
            )
        if not math.isfinite(self.ash_db_pct) or not (0.0 <= self.ash_db_pct <= 15.0):
            raise ValueError(f"ash_db_pct must be finite and in [0, 15], got {self.ash_db_pct}")
        if (self.volatile_matter_db_pct + self.ash_db_pct) >= 100.0:
            raise ValueError("Sum of volatile_matter_db_pct and ash_db_pct must be strictly less than 100%.")
        if self.hgi is not None:
            if not math.isfinite(self.hgi) or not (10.0 <= self.hgi <= 150.0):
                raise ValueError(f"Hardgrove Grindability Index (hgi) must be in [10, 150], got {self.hgi}")
        if not math.isfinite(self.bulk_density_kg_per_m3) or not (400.0 <= self.bulk_density_kg_per_m3 <= 1400.0):
            raise ValueError(
                f"bulk_density_kg_per_m3 must be in [400, 1400] kg/m^3, got {self.bulk_density_kg_per_m3}"
            )

    # --- Calorific, Proximate & Stowage Parity Properties ---
    @property
    def dry_mass_fraction(self) -> float:
        """Dry matter mass fraction (DMT per WMT)."""
        return (1.0 - (self.moisture_pct / 100.0))

    @property
    def fixed_carbon_db_pct(self) -> float:
        """Proximate-analysis Fixed Carbon on a dry basis: 100 - VM_db - Ash_db."""
        return (max(0.0, 100.0 - self.volatile_matter_db_pct - self.ash_db_pct))

    @property
    def lhv_gj_per_mt_ar(self) -> float:
        """As-received Net Calorific Value (LHV / NAR) in GJ/mt."""
        return (self.lhv_kcal_per_kg_ar * _KCAL_KG_TO_GJ_MT)

    @property
    def hhv_gj_per_mt_ar(self) -> float:
        """As-received Gross Calorific Value (HHV / GAR) in GJ/mt."""
        return (self.hhv_kcal_per_kg_ar * _KCAL_KG_TO_GJ_MT)

    @property
    def lhv_mmbtu_per_mt_ar(self) -> float:
        """As-received LHV in MMBtu_IT per metric tonne."""
        return (self.lhv_gj_per_mt_ar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mmbtu_per_mt_ar(self) -> float:
        """As-received HHV in MMBtu_IT per metric tonne."""
        return (self.hhv_gj_per_mt_ar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def lhv_mwh_per_mt_ar(self) -> float:
        return (self.lhv_gj_per_mt_ar / 3.6)

    @property
    def hhv_mwh_per_mt_ar(self) -> float:
        return (self.hhv_gj_per_mt_ar / 3.6)

    @property
    def stowage_factor_m3_per_wmt(self) -> float:
        """Bulk vessel stowage factor in cubic metres per wet metric tonne."""
        return (1000.0 / self.bulk_density_kg_per_m3)

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging petcoke [mass], [volume], and [energy]."""
        return create_petcoke_context(self, basis=basis, name=context_name)

    # --- Physical Mass, Stowage Volume & Energy Conversions ---
    def wmt_to_dmt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts as-received wet petcoke mass (WMT) to dry metric tonnes (DMT)."""
        if not isinstance(wet_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
        mt_q = wet_mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.dmt)

    def dmt_to_wmt(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts dry petcoke mass (DMT) to as-received wet metric tonnes (WMT)."""
        if not isinstance(dry_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(dry_mass).__name__}")
        mt_q = dry_mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, self.dry_mass_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.wmt)

    def wmt_to_bulk_volume(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts wet petcoke mass to vessel hold / stockpile bulk volume (m^3)."""
        if not isinstance(wet_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
        mt_q = wet_mass.to(ureg.metric_ton)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def bulk_volume_to_wmt(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts vessel hold / stockpile bulk volume (m^3) to wet petcoke mass (metric_ton)."""
        if not isinstance(bulk_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bulk_volume).__name__}")
        m3_q = bulk_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_energy(self, wet_mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts as-received petcoke mass to thermal energy (GJ)."""
        if not isinstance(wet_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
        mt_q = wet_mass.to(ureg.metric_ton)
        gj_per_mt = self.lhv_gj_per_mt_ar if basis.strip().lower() in {"lhv", "nar"} else self.hhv_gj_per_mt_ar
        scale = _scale_factor(mt_q.magnitude, gj_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.gigajoule)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts target thermal energy (GJ, MWh, MMBtu) to required petcoke mass (metric_ton)."""
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        gj_q = energy.to(ureg.gigajoule)
        gj_per_mt = self.lhv_gj_per_mt_ar if basis.strip().lower() in {"lhv", "nar"} else self.hhv_gj_per_mt_ar
        scale = _scale_factor(gj_q.magnitude, gj_per_mt)
        return ureg.Quantity(gj_q.magnitude / scale, ureg.metric_ton)

    def at_moisture(
        self,
        target_moisture_pct: float,
        *,
        name: str | None = None,
    ) -> PetcokeAssay:
        """Projects this petcoke assay to a new moisture content holding dry-basis specs constant."""
        m_new = float(target_moisture_pct)
        if not math.isfinite(m_new) or not (0.0 <= m_new < 30.0):
            raise ValueError(f"target_moisture_pct must be finite and in [0, 30), got {target_moisture_pct}")

        old_dry = self.dry_mass_fraction
        new_dry = 1.0 - (m_new / 100.0)

        # Latent heat of moisture vaporisation: 5.835 kcal_IT/kg per 1 wt% moisture
        hhv_db = self.hhv_kcal_per_kg_ar / old_dry
        latent_h_db = max(0.0, (self.hhv_kcal_per_kg_ar - self.lhv_kcal_per_kg_ar) - (5.835 * self.moisture_pct)) / old_dry

        new_hhv_ar = hhv_db * new_dry
        new_lhv_ar = new_hhv_ar - (latent_h_db * new_dry) - (5.835 * m_new)

        return PetcokeAssay(
            name=name or f"{self.name} ({m_new:.1f}% M)",
            lhv_kcal_per_kg_ar=(new_lhv_ar),
            hhv_kcal_per_kg_ar=(new_hhv_ar),
            moisture_pct=m_new,
            sulphur_db_pct=self.sulphur_db_pct,
            volatile_matter_db_pct=self.volatile_matter_db_pct,
            ash_db_pct=self.ash_db_pct,
            hgi=self.hgi,
            bulk_density_kg_per_m3=self.bulk_density_kg_per_m3,
            is_calcined=self.is_calcined,
        )

    # --- Commercial Pricing Parity ($/mt <-> $/MMBtu, $/WMT <-> $/DMT) ---
    def mt_price_to_mmbtu_price(self, price_per_mt: float, *, basis: str = "hhv") -> float:
        """Converts $/mt petcoke quote to $/MMBtu fuel-switching parity."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        mmbtu_per_mt = self.lhv_mmbtu_per_mt_ar if basis.strip().lower() in {"lhv", "nar"} else self.hhv_mmbtu_per_mt_ar
        return (p / mmbtu_per_mt)

    def mmbtu_price_to_mt_price(self, price_per_mmbtu: float, *, basis: str = "hhv") -> float:
        """Converts $/MMBtu energy quote to implied $/mt petcoke parity."""
        p = float(price_per_mmbtu)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mmbtu must be finite, got {price_per_mmbtu}")
        mmbtu_per_mt = self.lhv_mmbtu_per_mt_ar if basis.strip().lower() in {"lhv", "nar"} else self.hhv_mmbtu_per_mt_ar
        return (p * mmbtu_per_mt)

    def wmt_price_to_dmt_price(self, price_per_wmt: float) -> float:
        """Converts a $/WMT quote to $/DMT (dry basis)."""
        p = float(price_per_wmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_wmt must be finite, got {price_per_wmt}")
        return (p / self.dry_mass_fraction)

    def dmt_price_to_wmt_price(self, price_per_dmt: float) -> float:
        """Converts a $/DMT quote to $/WMT (wet basis)."""
        p = float(price_per_dmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt must be finite, got {price_per_dmt}")
        return (p * self.dry_mass_fraction)

    def blend_with(
        self,
        other: PetcokeAssay,
        *,
        self_wmt_fraction: float,
        blend_name: str | None = None,
    ) -> PetcokeAssay:
        """Blends two petcoke cargoes by wet mass (WMT) fraction.

        - AR calorific values & moisture: weighted by WMT
        - Dry-basis sulphur, VM, ash, and HGI: weighted by dry matter (DMT)
        - Bulk density: harmonic mean (additive bulk volumes)
        """
        if not isinstance(other, PetcokeAssay):
            raise TypeError(f"Expected PetcokeAssay, got {type(other).__name__}")
        if not math.isfinite(self_wmt_fraction) or not (0.0 <= self_wmt_fraction <= 1.0):
            raise ValueError("self_wmt_fraction must be between 0.0 and 1.0.")

        w1 = float(self_wmt_fraction)
        w2 = 1.0 - w1

        d1 = w1 * self.dry_mass_fraction
        d2 = w2 * other.dry_mass_fraction
        d_tot = d1 + d2
        dw1 = d1 / d_tot
        dw2 = d2 / d_tot

        vol_per_kg = (w1 / self.bulk_density_kg_per_m3) + (w2 / other.bulk_density_kg_per_m3)
        b_rho = (1.0 / vol_per_kg)

        b_hgi = (
            ((self.hgi * dw1) + (other.hgi * dw2))
            if (self.hgi is not None and other.hgi is not None)
            else None
        )

        return PetcokeAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            lhv_kcal_per_kg_ar=((self.lhv_kcal_per_kg_ar * w1) + (other.lhv_kcal_per_kg_ar * w2)),
            hhv_kcal_per_kg_ar=((self.hhv_kcal_per_kg_ar * w1) + (other.hhv_kcal_per_kg_ar * w2)),
            moisture_pct=((self.moisture_pct * w1) + (other.moisture_pct * w2)),
            sulphur_db_pct=((self.sulphur_db_pct * dw1) + (other.sulphur_db_pct * dw2)),
            volatile_matter_db_pct=((self.volatile_matter_db_pct * dw1) + (other.volatile_matter_db_pct * dw2)),
            ash_db_pct=((self.ash_db_pct * dw1) + (other.ash_db_pct * dw2)),
            hgi=b_hgi,
            bulk_density_kg_per_m3=b_rho,
            is_calcined=(self.is_calcined and other.is_calcined),
        )

    # --- Global Petcoke Benchmark Presets ---
    @classmethod
    def usgc_fuel_grade_65s(cls) -> PetcokeAssay:
        """US Gulf Coast High-Sulphur Fuel-Grade Sponge Coke (6.5% S, 7,500 kcal/kg NAR, HGI 45)."""
        return cls(
            name="USGC Fuel-Grade Petcoke (6.5% S)",
            lhv_kcal_per_kg_ar=7500.0,
            hhv_kcal_per_kg_ar=7750.0,
            moisture_pct=8.0,
            sulphur_db_pct=6.5,
            volatile_matter_db_pct=11.0,
            ash_db_pct=0.5,
            hgi=45.0,
            bulk_density_kg_per_m3=800.0,
            is_calcined=False,
        )

    @classmethod
    def usgc_fuel_grade_45s(cls) -> PetcokeAssay:
        """US Gulf Coast Mid-Sulphur Fuel-Grade Sponge Coke (4.5% S, 7,550 kcal/kg NAR, HGI 50)."""
        return cls(
            name="USGC Fuel-Grade Petcoke (4.5% S)",
            lhv_kcal_per_kg_ar=7550.0,
            hhv_kcal_per_kg_ar=7800.0,
            moisture_pct=8.0,
            sulphur_db_pct=4.5,
            volatile_matter_db_pct=11.5,
            ash_db_pct=0.5,
            hgi=50.0,
            bulk_density_kg_per_m3=800.0,
            is_calcined=False,
        )

    @classmethod
    def calcined_anode_grade(cls) -> PetcokeAssay:
        """Calcined Petroleum Coke (CPC) Anode Grade for aluminium smelting (>99% FC, 0.2% moisture)."""
        return cls(
            name="Calcined Anode-Grade Petcoke (CPC)",
            lhv_kcal_per_kg_ar=7950.0,
            hhv_kcal_per_kg_ar=7980.0,
            moisture_pct=0.2,
            sulphur_db_pct=2.5,
            volatile_matter_db_pct=0.3,
            ash_db_pct=0.2,
            hgi=None,
            bulk_density_kg_per_m3=880.0,
            is_calcined=True,
        )

    @classmethod
    def shot_coke(cls) -> PetcokeAssay:
        """Hard USGC Fuel-Grade Shot Coke (6.8% S, low HGI 35)."""
        return cls(
            name="USGC Shot Coke (6.8% S, HGI 35)",
            lhv_kcal_per_kg_ar=7450.0,
            hhv_kcal_per_kg_ar=7700.0,
            moisture_pct=8.0,
            sulphur_db_pct=6.8,
            volatile_matter_db_pct=10.0,
            ash_db_pct=0.6,
            hgi=35.0,
            bulk_density_kg_per_m3=850.0,
            is_calcined=False,
        )



# ---------------------------------------------------------------------------
# Industrial & Inorganic Chemicals Specification
# (Caustic Soda, Sulphuric/HCl Acids, Methanol, Soda Ash, Sulphur, Chlorine)
# ---------------------------------------------------------------------------

@validated_assay
@dataclass(frozen=True, slots=True)
class ChemicalAssay:
    """Immutable physical, concentration, and molar specification for industrial chemicals.

    Bridges commercial/solution mass (WMT), 100% active solute mass (DMT),
    liquid/bulk volume (m^3, gal), molar quantities (kmol), Baumé gravity,
    dilution water requirements, and regional pricing parity.
    """

    name: str
    formula: str
    molecular_weight_g_per_mol: float
    concentration_pct: float
    density_kg_per_m3: float
    reference_temp_c: float = 15.0
    is_liquid: bool = True

    def __post_init__(self) -> None:
        if not self.formula or not self.formula.strip():
            raise ValueError("formula must be a non-empty string (e.g. 'NaOH', 'H2SO4', 'CH3OH').")
        if not math.isfinite(self.molecular_weight_g_per_mol) or self.molecular_weight_g_per_mol <= 0.0:
            raise ValueError(
                f"molecular_weight_g_per_mol must be strictly positive, got {self.molecular_weight_g_per_mol}"
            )
        if not math.isfinite(self.concentration_pct) or not (0.0 < self.concentration_pct <= 100.0):
            raise ValueError(f"concentration_pct must be finite and in (0, 100], got {self.concentration_pct}")
        if not math.isfinite(self.density_kg_per_m3) or not (300.0 <= self.density_kg_per_m3 <= 3000.0):
            raise ValueError(
                f"density_kg_per_m3 must be finite and in [300, 3000] kg/m^3, got {self.density_kg_per_m3}"
            )
        if not math.isfinite(self.reference_temp_c) or not (-100.0 <= self.reference_temp_c <= 200.0):
            raise ValueError(f"reference_temp_c must be in [-100, 200] °C, got {self.reference_temp_c}")

    # --- Concentration, Gravity & Volumetric Properties ---
    @property
    def active_mass_fraction(self) -> float:
        """Mass fraction of active chemical solute (e.g. 0.50 for 50% Caustic Soda)."""
        return self.concentration_pct / 100.0

    @property
    def specific_gravity(self) -> float:
        """Specific gravity relative to water at 15°C (999.102 kg/m^3)."""
        return (self.density_kg_per_m3 / WATER_DENSITY_15C_KG_M3)

    @property
    def baume_degrees(self) -> float | None:
        """Baumé gravity (°Bé) for liquid chemicals; None for bulk solids.

        Uses the standard US/industrial scales:
        - Heavy liquids (SG >= 1.0): 145 - (145 / SG)
        - Light liquids (SG < 1.0):  (140 / SG) - 130
        """
        if not self.is_liquid:
            return None
        sg = self.density_kg_per_m3 / WATER_DENSITY_15C_KG_M3
        if sg >= 1.0:
            return (145.0 - (145.0 / sg))
        return ((140.0 / sg) - 130.0)

    @property
    def gal_per_mt(self) -> float:
        """US gallons per commercial metric tonne at reference_temp_c."""
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def lb_per_gal(self) -> float:
        """US pounds per US gallon at reference_temp_c."""
        return (_LB_PER_MT / self.gal_per_mt)

    @property
    def active_kg_per_m3(self) -> float:
        """Active solute mass (kg) contained in 1 m^3 of commercial solution/bulk."""
        return (self.density_kg_per_m3 * self.active_mass_fraction)

    def as_context(self, *, active_only: bool = True, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [substance]."""
        return create_chemical_context(self, active_only=active_only, name=context_name)

    # --- Physical Conversions (Volume <-> Solution Mass <-> Active Mass <-> Moles) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts commercial liquid/bulk volume (m^3, gal, L) to total mass (metric_ton)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts commercial total mass (mt, short_ton, lb) to volume (m^3)."""
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def solution_to_active_mass(self, solution_mass: pint.Quantity) -> pint.Quantity:
        """Converts commercial/wet solution mass (WMT) to 100% active solute mass (DMT)."""
        mt_q = self._validate_mass(solution_mass)
        scale = _scale_factor(mt_q.magnitude, self.active_mass_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def active_to_solution_mass(self, active_mass: pint.Quantity) -> pint.Quantity:
        """Converts 100% active solute mass (DMT) to required commercial solution mass (WMT)."""
        mt_q = self._validate_mass(active_mass)
        scale = _scale_factor(mt_q.magnitude, self.active_mass_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def volume_to_active_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts commercial volume (e.g. m^3, gal) directly to 100% active solute mass (metric_ton)."""
        return self.solution_to_active_mass(self.volume_to_mass(volume))

    def active_mass_to_volume(self, active_mass: pint.Quantity) -> pint.Quantity:
        """Converts target 100% active solute mass directly to required commercial volume (m^3)."""
        return self.mass_to_volume(self.active_to_solution_mass(active_mass))

    def mass_to_moles(self, mass: pint.Quantity, *, active_only: bool = True) -> pint.Quantity:
        """Converts mass to kilomoles; active_only applies the solution concentration."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        kg_q = mass.to(ureg.kilogram)
        conc_mult = self.active_mass_fraction if active_only else 1.0
        kmol_per_kg = conc_mult / self.molecular_weight_g_per_mol
        scale = _scale_factor(kg_q.magnitude, kmol_per_kg)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilomole)

    def moles_to_mass(self, moles: pint.Quantity, *, active_only: bool = True) -> pint.Quantity:
        """Converts moles to tonnes; active_only returns the required solution mass."""
        if not isinstance(moles, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(moles).__name__}")
        kmol_q = moles.to(ureg.kilomole)
        conc_mult = self.active_mass_fraction if active_only else 1.0
        mt_per_kmol = (self.molecular_weight_g_per_mol / 1000.0) / conc_mult
        scale = _scale_factor(kmol_q.magnitude, mt_per_kmol)
        return ureg.Quantity(kmol_q.magnitude * scale, ureg.metric_ton)

    # --- Dilution & Strength Blending ---
    def water_required_to_dilute(
        self,
        solution_mass: pint.Quantity,
        target_concentration_pct: float,
    ) -> pint.Quantity:
        """Computes mass of water (metric_ton) required to dilute solution_mass to target_concentration_pct."""
        c_target = float(target_concentration_pct)
        if not math.isfinite(c_target) or not (0.0 < c_target < self.concentration_pct):
            raise ValueError(
                f"target_concentration_pct must be strictly in (0, {self.concentration_pct}), got {target_concentration_pct}"
            )
        mt_q = self._validate_mass(solution_mass)
        water_ratio = (self.concentration_pct / c_target) - 1.0
        scale = _scale_factor(mt_q.magnitude, water_ratio)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def dilute_to(
        self,
        target_concentration_pct: float,
        *,
        new_density_kg_per_m3: float | None = None,
        name: str | None = None,
    ) -> ChemicalAssay:
        """Returns a diluted ChemicalAssay at target_concentration_pct.

        If new_density_kg_per_m3 is omitted, estimates density via ideal volume
        additivity with water at 15 degrees C. This approximation does not model
        thermal effects or volume contraction; supply a measured density when needed.
        """
        c_target = float(target_concentration_pct)
        if not math.isfinite(c_target) or not (0.0 < c_target < self.concentration_pct):
            raise ValueError(
                f"target_concentration_pct must be strictly in (0, {self.concentration_pct}), got {target_concentration_pct}"
            )
        if new_density_kg_per_m3 is None:
            sol_frac = c_target / self.concentration_pct
            water_frac = 1.0 - sol_frac
            vol_per_kg = (sol_frac / self.density_kg_per_m3) + (water_frac / WATER_DENSITY_15C_KG_M3)
            est_rho = (1.0 / vol_per_kg)
        else:
            est_rho = float(new_density_kg_per_m3)

        return ChemicalAssay(
            name=name or f"{self.formula} ({c_target:g}% Solution)",
            formula=self.formula,
            molecular_weight_g_per_mol=self.molecular_weight_g_per_mol,
            concentration_pct=c_target,
            density_kg_per_m3=est_rho,
            reference_temp_c=self.reference_temp_c,
            is_liquid=True,
        )

    def blend_with(
        self,
        other: ChemicalAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> ChemicalAssay:
        """Blends two lots/strengths of the same chemical by total mass fraction."""
        if not isinstance(other, ChemicalAssay):
            raise TypeError(f"Expected ChemicalAssay, got {type(other).__name__}")
        if self.formula != other.formula:
            raise ValueError(f"Cannot blend different chemicals: {self.formula!r} vs {other.formula!r}")
        if self.is_liquid != other.is_liquid:
            raise ValueError("Cannot directly blend a liquid solution with a solid grade.")
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1
        b_conc = (self.concentration_pct * w1) + (other.concentration_pct * w2)
        vol_per_kg = (w1 / self.density_kg_per_m3) + (w2 / other.density_kg_per_m3)

        return ChemicalAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            formula=self.formula,
            molecular_weight_g_per_mol=self.molecular_weight_g_per_mol,
            concentration_pct=b_conc,
            density_kg_per_m3=(1.0 / vol_per_kg),
            reference_temp_c=self.reference_temp_c,
            is_liquid=self.is_liquid,
        )

    # --- Commercial Pricing Parity ($/WMT <-> $/DMT 100%, $/gal <-> $/mt, $/st <-> $/mt) ---
    def solution_price_to_active_price(self, price_per_wmt: float) -> float:
        """Converts a $/WMT commercial solution quote to $/DMT (100% active basis)."""
        p = float(price_per_wmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_wmt must be finite, got {price_per_wmt}")
        return (p / self.active_mass_fraction)

    def active_price_to_solution_price(self, price_per_dmt_active: float) -> float:
        """Converts a $/DMT (100% active basis) quote to $/WMT commercial solution."""
        p = float(price_per_dmt_active)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt_active must be finite, got {price_per_dmt_active}")
        return (p * self.active_mass_fraction)

    def gal_price_to_mt_price(self, price_per_gal: float) -> float:
        """Converts USGC $/gal quote (e.g. Methanol) to international $/mt equivalent."""
        p = float(price_per_gal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gal must be finite, got {price_per_gal}")
        return (p * self.gal_per_mt)

    def mt_price_to_gal_price(self, price_per_mt: float) -> float:
        """Converts international $/mt quote to USGC $/gal equivalent."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.gal_per_mt)

    @staticmethod
    def st_price_to_mt_price(price_per_st: float) -> float:
        """Converts US $/short_ton quote (e.g. US Caustic/Soda Ash) to $/metric_ton."""
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / _MT_PER_SHORT_TON)

    @staticmethod
    def mt_price_to_st_price(price_per_mt: float) -> float:
        """Converts $/metric_ton quote to US $/short_ton."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p * _MT_PER_SHORT_TON)

    # --- Global Industrial Chemical Benchmark Presets ---
    @classmethod
    def caustic_soda_50pct(cls) -> ChemicalAssay:
        """Liquid Caustic Soda (50% NaOH aqueous — global chlor-alkali benchmark)."""
        return cls(
            name="Caustic Soda Liquid (50% NaOH)",
            formula="NaOH",
            molecular_weight_g_per_mol=39.997,
            concentration_pct=50.0,
            density_kg_per_m3=1525.0,
            reference_temp_c=20.0,
            is_liquid=True,
        )

    @classmethod
    def caustic_soda_flakes(cls) -> ChemicalAssay:
        """Solid Caustic Soda Flakes (99% NaOH)."""
        return cls(
            name="Caustic Soda Flakes (99% NaOH)",
            formula="NaOH",
            molecular_weight_g_per_mol=39.997,
            concentration_pct=99.0,
            density_kg_per_m3=1050.0,
            reference_temp_c=20.0,
            is_liquid=False,
        )

    @classmethod
    def sulphuric_acid_98pct(cls) -> ChemicalAssay:
        """Commercial Sulphuric Acid (98% H2SO4 — smelter/fertiliser benchmark)."""
        return cls(
            name="Sulphuric Acid (98% H2SO4)",
            formula="H2SO4",
            molecular_weight_g_per_mol=98.079,
            concentration_pct=98.0,
            density_kg_per_m3=1836.0,
            reference_temp_c=15.0,
            is_liquid=True,
        )

    @classmethod
    def sulphuric_acid_66be(cls) -> ChemicalAssay:
        """66°Bé Sulphuric Acid (93.19% H2SO4 — North American commercial grade)."""
        return cls(
            name="Sulphuric Acid (66°Bé / 93.2% H2SO4)",
            formula="H2SO4",
            molecular_weight_g_per_mol=98.079,
            concentration_pct=93.19,
            density_kg_per_m3=1834.0,
            reference_temp_c=15.0,
            is_liquid=True,
        )

    @classmethod
    def hydrochloric_acid_32pct(cls) -> ChemicalAssay:
        """Commercial Hydrochloric Acid (32% HCl, ~20°Bé)."""
        return cls(
            name="Hydrochloric Acid (32% HCl)",
            formula="HCl",
            molecular_weight_g_per_mol=36.461,
            concentration_pct=32.0,
            density_kg_per_m3=1160.0,
            reference_temp_c=15.0,
            is_liquid=True,
        )

    @classmethod
    def methanol(cls) -> ChemicalAssay:
        """IMPCA-Grade Methanol (99.85% CH3OH — USGC $/gal vs. Rotterdam/China $/mt)."""
        return cls(
            name="Methanol (IMPCA 99.85%)",
            formula="CH3OH",
            molecular_weight_g_per_mol=32.042,
            concentration_pct=99.85,
            density_kg_per_m3=796.0,
            reference_temp_c=15.0,
            is_liquid=True,
        )

    @classmethod
    def soda_ash_dense(cls) -> ChemicalAssay:
        """Dense Soda Ash (99.2% Na2CO3 — glass & lithium refining benchmark)."""
        return cls(
            name="Dense Soda Ash (Na2CO3)",
            formula="Na2CO3",
            molecular_weight_g_per_mol=105.988,
            concentration_pct=99.2,
            density_kg_per_m3=1050.0,
            reference_temp_c=20.0,
            is_liquid=False,
        )

    @classmethod
    def molten_sulphur(cls) -> ChemicalAssay:
        """Molten Elemental Sulphur (99.5% S — pumped at ~135°C)."""
        return cls(
            name="Molten Sulphur (99.5% S)",
            formula="S",
            molecular_weight_g_per_mol=32.06,
            concentration_pct=99.5,
            density_kg_per_m3=1790.0,
            reference_temp_c=135.0,
            is_liquid=True,
        )

    @classmethod
    def hydrogen_peroxide_50pct(cls) -> ChemicalAssay:
        """Industrial Hydrogen Peroxide (50% H2O2 aqueous)."""
        return cls(
            name="Hydrogen Peroxide (50% H2O2)",
            formula="H2O2",
            molecular_weight_g_per_mol=34.015,
            concentration_pct=50.0,
            density_kg_per_m3=1196.0,
            reference_temp_c=20.0,
            is_liquid=True,
        )

    @classmethod
    def liquid_chlorine(cls) -> ChemicalAssay:
        """Liquid Chlorine (99.5% Cl2 — pressurised at 0°C)."""
        return cls(
            name="Liquid Chlorine (Cl2)",
            formula="Cl2",
            molecular_weight_g_per_mol=70.906,
            concentration_pct=99.5,
            density_kg_per_m3=1468.0,
            reference_temp_c=0.0,
            is_liquid=True,
        )


# ---------------------------------------------------------------------------
# Hydrogen & Clean Energy Carriers Specification
# (Gaseous H2, Compressed 350/700 bar, Liquid LH2, Ammonia NH3, LOHC)
# ---------------------------------------------------------------------------

_NH3_H2_MASS_FRACTION: Final[float] = 0.17755  # Approximate stoichiometric fraction: 3 * M(H) / M(NH3) = 3.0238 / 17.0305
# CI screening only: 70% saving against 94 g CO2e/MJ (LHV).
# https://observatory.clean-hydrogen.europa.eu/index.php/eu-policy/renewable-energy-directive
_EU_RFNBO_MAX_CI_G_PER_MJ: Final[float] = 28.2


@validated_assay
@dataclass(frozen=True, slots=True)
class HydrogenAssay:
    """Immutable physical, energetic, and regulatory specification for hydrogen and carriers.

    Bridges storage/transport volume (m^3, Nm^3, L), carrier and contained H2 mass
    (kg, metric_ton), energy (MWh, GJ, MMBtu), electrolyser power demand, ammonia/LOHC
    cracking parity, and US 45V / EU RFNBO carbon-intensity screening.

    Density describes storage volume at the stated temperature and pressure,
    not an equation of state. Purity is metadata; h2_mass_fraction controls yield.
    Context masses always represent total product/carrier mass.
    """

    name: str
    state: str
    density_kg_per_m3: float
    h2_mass_fraction: float = 1.0
    reference_temp_c: float = 0.0
    reference_pressure_bar: float = 1.01325
    lhv_mj_per_kg_h2: float = H2_LHV_MJ_PER_KG
    hhv_mj_per_kg_h2: float = H2_HHV_MJ_PER_KG
    default_cracking_efficiency: float = 1.0
    purity_pct: float = 99.97
    carbon_intensity_kg_co2e_per_kg_h2: float | None = None

    def __post_init__(self) -> None:
        if self.state not in {"gas", "liquid", "carrier"}:
            raise ValueError(f"state must be 'gas', 'liquid', or 'carrier', got {self.state!r}")
        if not math.isfinite(self.density_kg_per_m3) or not (0.01 <= self.density_kg_per_m3 <= 2000.0):
            raise ValueError(
                f"density_kg_per_m3 must be finite and in [0.01, 2000] kg/m^3, got {self.density_kg_per_m3}"
            )
        if not math.isfinite(self.h2_mass_fraction) or not (0.0 < self.h2_mass_fraction <= 1.0):
            raise ValueError(f"h2_mass_fraction must be in (0, 1], got {self.h2_mass_fraction}")
        if not math.isfinite(self.reference_temp_c) or not (-273.15 < self.reference_temp_c <= 200.0):
            raise ValueError(f"reference_temp_c must be in (-273.15, 200] °C, got {self.reference_temp_c}")
        if not math.isfinite(self.reference_pressure_bar) or self.reference_pressure_bar <= 0.0:
            raise ValueError(f"reference_pressure_bar must be strictly positive, got {self.reference_pressure_bar}")
        if not math.isfinite(self.lhv_mj_per_kg_h2) or self.lhv_mj_per_kg_h2 <= 0.0:
            raise ValueError(f"lhv_mj_per_kg_h2 must be strictly positive, got {self.lhv_mj_per_kg_h2}")
        if not math.isfinite(self.hhv_mj_per_kg_h2) or self.hhv_mj_per_kg_h2 < self.lhv_mj_per_kg_h2:
            raise ValueError("hhv_mj_per_kg_h2 must be finite and >= lhv_mj_per_kg_h2.")
        if not math.isfinite(self.default_cracking_efficiency) or not (0.0 < self.default_cracking_efficiency <= 1.0):
            raise ValueError(
                f"default_cracking_efficiency must be in (0, 1], got {self.default_cracking_efficiency}"
            )
        if not math.isfinite(self.purity_pct) or not (0.0 < self.purity_pct <= 100.0):
            raise ValueError(f"purity_pct must be in (0, 100], got {self.purity_pct}")
        if self.carbon_intensity_kg_co2e_per_kg_h2 is not None:
            if not math.isfinite(self.carbon_intensity_kg_co2e_per_kg_h2) or self.carbon_intensity_kg_co2e_per_kg_h2 < 0.0:
                raise ValueError(
                    f"carbon_intensity_kg_co2e_per_kg_h2 must be non-negative, got {self.carbon_intensity_kg_co2e_per_kg_h2}"
                )

    # --- Energy Density, Volumetric Storage & Regulatory Properties ---
    @property
    def lhv_kwh_per_kg_h2(self) -> float:
        """Lower Heating Value in kWh per kg of contained H2 (~33.32 kWh/kg)."""
        return (self.lhv_mj_per_kg_h2 / 3.6)

    @property
    def hhv_kwh_per_kg_h2(self) -> float:
        """Higher Heating Value in kWh per kg of contained H2 (~39.41 kWh/kg)."""
        return (self.hhv_mj_per_kg_h2 / 3.6)

    @property
    def lhv_mmbtu_per_kg_h2(self) -> float:
        """Lower Heating Value in MMBtu per kg of contained H2."""
        return (self.lhv_mj_per_kg_h2 / (MJ_PER_BTU_IT * 1_000_000.0))

    @property
    def hhv_mmbtu_per_kg_h2(self) -> float:
        """Higher Heating Value in MMBtu per kg of contained H2."""
        return (self.hhv_mj_per_kg_h2 / (MJ_PER_BTU_IT * 1_000_000.0))

    @property
    def h2_kg_per_m3(self) -> float:
        """Volumetric H2 storage density (kg contained H2 per m^3 of gas/liquid/carrier)."""
        return (self.density_kg_per_m3 * self.h2_mass_fraction)

    @property
    def is_eu_rfnbo_compliant(self) -> bool | None:
        """CI-only RFNBO screen (28.2 g CO2e/MJ LHV); None if CI is unset.

        Assumes supplied CI includes the required lifecycle/delivery boundary.
        This does not verify renewable origin, additionality, or electricity
        temporal/geographical correlation and cannot certify full compliance.
        """
        if self.carbon_intensity_kg_co2e_per_kg_h2 is None:
            return None
        return self.carbon_intensity_kg_co2e_per_kg_h2 <= (_EU_RFNBO_MAX_CI_G_PER_MJ * self.lhv_mj_per_kg_h2 / 1000.0)

    @property
    def us_45v_credit_usd_per_kg_h2(self) -> float | None:
        """Unadjusted statutory 45V tier reference, assuming the 5x PWA rate.

        Uses a $0.60 base and 33.4% for the 0.45-to-1.5 tier ($1.002/kg).
        Not a current-year credit or eligibility determination. Use
        us_45v_credit_for_year() for supported inflation-adjusted amounts.
        https://uscode.house.gov/us/usc/t26/s45V/
        """
        ci = self.carbon_intensity_kg_co2e_per_kg_h2
        if ci is None:
            return None
        if ci < 0.45:
            return 3.00
        if ci < 1.50:
            return 1.002
        if ci < 2.50:
            return 0.75
        if ci <= 4.00:
            return 0.60
        return 0.0

    def us_45v_credit_for_year(
        self,
        production_year: int,
        *,
        prevailing_wage_and_apprenticeship: bool = False,
    ) -> float | None:
        """Returns the 2026 CI-tier rate, conditional on all other eligibility rules.

        PWA must be explicitly selected for the 5x increased credit. Other
        years require an updated rate table rather than silently reusing 2026.
        Source: IRS Notice 2026-41, section 3.02.
        https://www.irs.gov/irb/2026-29_irb
        """
        if production_year != 2026:
            raise ValueError("Only production_year=2026 has a verified rate table.")
        ci = self.carbon_intensity_kg_co2e_per_kg_h2
        if ci is None:
            return None
        if ci < 0.45:
            rate = 0.656
        elif ci < 1.5:
            rate = 0.219
        elif ci < 2.5:
            rate = 0.164
        elif ci <= 4.0:
            rate = 0.131
        else:
            rate = 0.0
        return float(Decimal(str(rate)) * (5 if prevailing_wage_and_apprenticeship else 1))

    def _mj_per_kg_h2(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"lhv", "ncv"}:
            return self.lhv_mj_per_kg_h2
        if key in {"hhv", "gcv"}:
            return self.hhv_mj_per_kg_h2
        raise ValueError(f"basis must be 'lhv'/'ncv' or 'hhv'/'gcv', got {basis!r}")

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [energy]."""
        return create_hydrogen_context(self, basis=basis, name=context_name)

    # --- Physical Conversions (Volume <-> Carrier Mass <-> H2 Mass <-> Nm^3 <-> Energy) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.kilogram)

    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts storage/transport volume (m^3, L) to total product/carrier mass (kg)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        scale = _scale_factor(m3_q.magnitude, self.density_kg_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.kilogram)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts total product/carrier mass (kg, metric_ton) to storage/transport volume (m^3)."""
        kg_q = self._validate_mass(mass)
        scale = _scale_factor(kg_q.magnitude, self.density_kg_per_m3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.meter ** 3)

    def carrier_to_h2_mass(
        self,
        carrier_mass: pint.Quantity,
        *,
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts carrier mass (e.g. mt NH3, mt LOHC) to delivered pure H2 mass (kg)."""
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0, 1], got {eff}")
        kg_q = self._validate_mass(carrier_mass)
        factor = self.h2_mass_fraction * eff
        scale = _scale_factor(kg_q.magnitude, factor)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def h2_to_carrier_mass(
        self,
        h2_mass: pint.Quantity,
        *,
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts target delivered H2 mass (kg, mt) to required carrier mass (metric_ton)."""
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0, 1], got {eff}")
        kg_q = self._validate_mass(h2_mass)
        factor = self.h2_mass_fraction * eff
        scale = _scale_factor(kg_q.magnitude, factor)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram).to(ureg.metric_ton)

    @staticmethod
    def h2_mass_to_nm3(h2_mass: pint.Quantity) -> pint.Quantity:
        """Converts pure H2 mass (kg, mt) to Normal cubic metres (Nm^3 at 0°C, 1 atm)."""
        if not isinstance(h2_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(h2_mass).__name__}")
        kg_q = h2_mass.to(ureg.kilogram)
        scale = _scale_factor(kg_q.magnitude, H2_DENSITY_0C_KG_NM3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.Nm3_H2)

    @staticmethod
    def nm3_to_h2_mass(nm3_volume: pint.Quantity) -> pint.Quantity:
        """Converts Normal cubic metres (Nm^3 at 0°C, 1 atm) to pure H2 mass (kg)."""
        if not isinstance(nm3_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(nm3_volume).__name__}")
        if nm3_volume.dimensionality == ureg.Nm3_H2.dimensionality:
            m3_q = nm3_volume.to(ureg.Nm3_H2)
        else:
            m3_q = nm3_volume.to(ureg.meter ** 3)
        scale = _scale_factor(m3_q.magnitude, H2_DENSITY_0C_KG_NM3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.kilogram)

    def mass_to_energy(
        self,
        mass: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts product/carrier mass to delivered H2 chemical energy (MWh)."""
        h2_kg = self.carrier_to_h2_mass(mass, cracking_efficiency=cracking_efficiency)
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        scale = _scale_factor(h2_kg.magnitude, mwh_per_kg)
        return ureg.Quantity(h2_kg.magnitude * scale, ureg.MWh)

    def energy_to_mass(
        self,
        energy: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts target delivered H2 energy (MWh, GJ, MMBtu) to required product/carrier mass (kg)."""
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        mwh_q = energy.to(ureg.MWh)
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        scale = _scale_factor(mwh_q.magnitude, mwh_per_kg)
        req_h2_kg = ureg.Quantity(mwh_q.magnitude / scale, ureg.kilogram)
        return self.h2_to_carrier_mass(req_h2_kg, cracking_efficiency=cracking_efficiency).to(ureg.kilogram)

    def volume_to_energy(
        self,
        volume: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts storage/transport volume (m^3) directly to delivered H2 energy (MWh)."""
        return self.mass_to_energy(
            self.volume_to_mass(volume),
            basis=basis,
            cracking_efficiency=cracking_efficiency,
        )

    def energy_to_volume(
        self,
        energy: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts target delivered H2 energy (MWh, GJ) directly to storage/transport volume (m^3)."""
        return self.mass_to_volume(
            self.energy_to_mass(energy, basis=basis, cracking_efficiency=cracking_efficiency)
        )

    def electrolysis_power_required(
        self,
        h2_mass: pint.Quantity,
        *,
        electrolyser_efficiency_lhv: float = 0.65,
    ) -> pint.Quantity:
        """Computes electrical energy input (MWh_e) required to produce h2_mass via electrolysis.

        At 65% LHV efficiency, 1 kg H2 requires ~51.26 kWh_e (0.05126 MWh_e).
        """
        if not math.isfinite(electrolyser_efficiency_lhv) or not (0.0 < electrolyser_efficiency_lhv <= 1.0):
            raise ValueError(
                f"electrolyser_efficiency_lhv must be in (0, 1], got {electrolyser_efficiency_lhv}"
            )
        kg_q = self._validate_mass(h2_mass)
        mwh_e_per_kg = (self.lhv_mj_per_kg_h2 / 3600.0) / electrolyser_efficiency_lhv
        scale = _scale_factor(kg_q.magnitude, mwh_e_per_kg)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.MWh)

    # --- Commercial Pricing Parity ($/kg H2 <-> $/MWh, $/MMBtu, $/mt Carrier) ---
    def kg_price_to_mwh_price(self, price_per_kg_h2: float, *, basis: str = "lhv") -> float:
        """Converts $/kg H2 quote to $/MWh energy parity."""
        p = float(price_per_kg_h2)
        if not math.isfinite(p):
            raise ValueError(f"price_per_kg_h2 must be finite, got {price_per_kg_h2}")
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        return (p / mwh_per_kg)

    def mwh_price_to_kg_price(self, price_per_mwh: float, *, basis: str = "lhv") -> float:
        """Converts $/MWh energy quote to $/kg H2 parity."""
        p = float(price_per_mwh)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mwh must be finite, got {price_per_mwh}")
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        return (p * mwh_per_kg)

    def kg_price_to_mmbtu_price(self, price_per_kg_h2: float, *, basis: str = "hhv") -> float:
        """Converts $/kg H2 quote to $/MMBtu natural gas switching parity."""
        p = float(price_per_kg_h2)
        if not math.isfinite(p):
            raise ValueError(f"price_per_kg_h2 must be finite, got {price_per_kg_h2}")
        mmbtu_per_kg = self._mj_per_kg_h2(basis) / (MJ_PER_BTU_IT * 1_000_000.0)
        return (p / mmbtu_per_kg)

    def mmbtu_price_to_kg_price(self, price_per_mmbtu: float, *, basis: str = "hhv") -> float:
        """Converts $/MMBtu quote to $/kg H2 parity."""
        p = float(price_per_mmbtu)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mmbtu must be finite, got {price_per_mmbtu}")
        mmbtu_per_kg = self._mj_per_kg_h2(basis) / (MJ_PER_BTU_IT * 1_000_000.0)
        return (p * mmbtu_per_kg)

    def carrier_mt_price_to_h2_kg_price(
        self,
        price_per_mt_carrier: float,
        *,
        cracking_efficiency: float | None = None,
        cracking_cost_per_kg_h2: float = 0.0,
    ) -> float:
        """Converts carrier $/mt (e.g. $/mt NH3) to implied delivered $/kg H2."""
        p = float(price_per_mt_carrier)
        cost = float(cracking_cost_per_kg_h2)
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(p) or not math.isfinite(cost) or not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError("Inputs must be finite and cracking_efficiency must be in (0, 1].")
        delivered_kg_h2_per_mt = 1000.0 * self.h2_mass_fraction * eff
        return ((p / delivered_kg_h2_per_mt) + cost)

    def h2_kg_price_to_carrier_mt_price(
        self,
        price_per_kg_h2: float,
        *,
        cracking_efficiency: float | None = None,
        cracking_cost_per_kg_h2: float = 0.0,
    ) -> float:
        """Computes implied carrier $/mt netback from delivered $/kg H2."""
        p = float(price_per_kg_h2)
        cost = float(cracking_cost_per_kg_h2)
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(p) or not math.isfinite(cost) or not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError("Inputs must be finite and cracking_efficiency must be in (0, 1].")
        delivered_kg_h2_per_mt = 1000.0 * self.h2_mass_fraction * eff
        return ((p - cost) * delivered_kg_h2_per_mt)

    def blend_with(
        self,
        other: HydrogenAssay,
        *,
        self_h2_mass_fraction: float,
        blend_name: str | None = None,
    ) -> HydrogenAssay:
        """Blends two compatible H2 streams (e.g. Green + Blue H2 mass-balancing).

        Requires identical physical state and carrier fraction; weights carbon
        intensity and purity linearly by H2 mass.
        """
        if not isinstance(other, HydrogenAssay):
            raise TypeError(f"Expected HydrogenAssay, got {type(other).__name__}")
        if self.state != other.state or abs(self.h2_mass_fraction - other.h2_mass_fraction) > 1e-6:
            raise ValueError("Can only blend hydrogen streams with matching state and h2_mass_fraction.")
        for attribute in (
            "density_kg_per_m3", "reference_temp_c", "reference_pressure_bar",
            "lhv_mj_per_kg_h2", "hhv_mj_per_kg_h2",
        ):
            if not math.isclose(getattr(self, attribute), getattr(other, attribute), rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"Cannot blend streams with different {attribute}; normalise their physical basis first.")
        if not math.isfinite(self_h2_mass_fraction) or not (0.0 <= self_h2_mass_fraction <= 1.0):
            raise ValueError("self_h2_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_h2_mass_fraction)
        w2 = 1.0 - w1
        if w1 == 1.0:
            b_ci = self.carbon_intensity_kg_co2e_per_kg_h2
        elif w1 == 0.0:
            b_ci = other.carbon_intensity_kg_co2e_per_kg_h2
        elif self.carbon_intensity_kg_co2e_per_kg_h2 is None or other.carbon_intensity_kg_co2e_per_kg_h2 is None:
            b_ci = None
        else:
            b_ci = self.carbon_intensity_kg_co2e_per_kg_h2 * w1 + other.carbon_intensity_kg_co2e_per_kg_h2 * w2
        return HydrogenAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            state=self.state,
            density_kg_per_m3=self.density_kg_per_m3,
            h2_mass_fraction=self.h2_mass_fraction,
            reference_temp_c=self.reference_temp_c,
            reference_pressure_bar=self.reference_pressure_bar,
            lhv_mj_per_kg_h2=self.lhv_mj_per_kg_h2,
            hhv_mj_per_kg_h2=self.hhv_mj_per_kg_h2,
            default_cracking_efficiency=(self.default_cracking_efficiency * w1) + (other.default_cracking_efficiency * w2),
            purity_pct=((self.purity_pct * w1) + (other.purity_pct * w2)),
            carbon_intensity_kg_co2e_per_kg_h2=b_ci,
        )

    # --- Global Hydrogen & Carrier Benchmark Presets ---
    @classmethod
    def gaseous_h2_ntp(cls) -> HydrogenAssay:
        """Pure Gaseous Hydrogen at NTP (0°C, 1 atm — Nm^3 metering basis)."""
        return cls(
            name="Gaseous H2 (NTP, 0°C / 1 atm)",
            state="gas",
            density_kg_per_m3=H2_DENSITY_0C_KG_NM3,
        )

    @classmethod
    def compressed_h2_350bar(cls) -> HydrogenAssay:
        """Compressed Hydrogen at 350 bar, 15°C (heavy-duty mobility & tube trailers)."""
        return cls(
            name="Compressed H2 (350 bar)",
            state="gas",
            density_kg_per_m3=23.7,
            reference_temp_c=15.0,
            reference_pressure_bar=350.0,
        )

    @classmethod
    def compressed_h2_700bar(cls) -> HydrogenAssay:
        """Compressed Hydrogen at 700 bar, 15°C (light-duty FCEV standard)."""
        return cls(
            name="Compressed H2 (700 bar)",
            state="gas",
            density_kg_per_m3=39.7,
            reference_temp_c=15.0,
            reference_pressure_bar=700.0,
        )

    @classmethod
    def liquid_h2(cls) -> HydrogenAssay:
        """Cryogenic Liquid Hydrogen (LH2 at -252.87°C, 1 atm)."""
        return cls(
            name="Liquid Hydrogen (LH2)",
            state="liquid",
            density_kg_per_m3=70.85,
            reference_temp_c=-252.87,
            reference_pressure_bar=1.01325,
        )

    @classmethod
    def ammonia_carrier(cls) -> HydrogenAssay:
        """Liquid Ammonia as H2 Carrier (NH3 at -33°C, 17.755 wt% H2, ~85% cracking efficiency)."""
        return cls(
            name="Ammonia H2 Carrier (Liquid NH3)",
            state="carrier",
            density_kg_per_m3=682.0,
            h2_mass_fraction=_NH3_H2_MASS_FRACTION,
            reference_temp_c=-33.3,
            reference_pressure_bar=1.01325,
            default_cracking_efficiency=0.85,
        )

    @classmethod
    def lohc_dbt(cls) -> HydrogenAssay:
        """Liquid Organic Hydrogen Carrier (H18-Dibenzyltoluene, 6.2 wt% H2, ~90% release)."""
        return cls(
            name="LOHC (H18-Dibenzyltoluene)",
            state="carrier",
            density_kg_per_m3=910.0,
            h2_mass_fraction=0.062,
            reference_temp_c=15.0,
            reference_pressure_bar=1.01325,
            default_cracking_efficiency=0.90,
        )

    @classmethod
    def green_h2_rfnbo(cls) -> HydrogenAssay:
        """Electrolytic Green Hydrogen (illustrative low-CI stream, 0.4 kg CO2e/kg H2)."""
        return cls(
            name="Green H2 (Low-CI Reference)",
            state="gas",
            density_kg_per_m3=H2_DENSITY_0C_KG_NM3,
            carbon_intensity_kg_co2e_per_kg_h2=0.40,
        )

    @classmethod
    def blue_h2_atr(cls) -> HydrogenAssay:
        """Autothermal Reforming Blue Hydrogen with CCS (~2.0 kg CO2e/kg H2)."""
        return cls(
            name="Blue H2 (ATR + CCS)",
            state="gas",
            density_kg_per_m3=H2_DENSITY_0C_KG_NM3,
            carbon_intensity_kg_co2e_per_kg_h2=2.00,
        )


# ---------------------------------------------------------------------------
# Ferroalloys Specification (Bulk & Noble Steelmaking Alloys)
# (HC/LC FeCr, HC FeMn, SiMn, FeSi 75, FeMo 60, FeV 80, FeNb 65)
# ---------------------------------------------------------------------------

_LB_PER_KG: Final[float] = 1.0 / 0.45359237  # Exact: 2.20462262185 lb/kg


@validated_assay
@dataclass(frozen=True, slots=True)
class FerroalloyAssay:
    """Immutable metallurgical and commercial specification for ferroalloys.

    Bridges gross alloy mass (mt, kg, lb), contained alloying element mass
    (primary and secondary elements), bulk stowage volume (m^3), ladle dosing
    requirements, pro-rata basis-grade settlement, and $/mt <-> $/kg <-> $/lb
    contained-metal pricing parity.
    """

    name: str
    primary_element: str
    primary_element_pct: float
    basis_grade_pct: float | None = None
    carbon_pct: float | None = None
    secondary_elements: Mapping[str, float] | None = None
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.primary_element, str) or not self.primary_element.strip():
            raise ValueError("primary_element must be a non-empty string (e.g. 'Cr', 'Mn', 'Si', 'Mo', 'V').")
        object.__setattr__(self, "primary_element", self.primary_element.strip())
        if not math.isfinite(self.primary_element_pct) or not (0.0 < self.primary_element_pct <= 100.0):
            raise ValueError(
                f"primary_element_pct must be finite and in (0, 100], got {self.primary_element_pct}"
            )
        if self.basis_grade_pct is not None:
            if not math.isfinite(self.basis_grade_pct) or not (0.0 < self.basis_grade_pct <= 100.0):
                raise ValueError(f"basis_grade_pct must be in (0, 100], got {self.basis_grade_pct}")
        if self.carbon_pct is not None:
            if not math.isfinite(self.carbon_pct) or not (0.0 <= self.carbon_pct <= 15.0):
                raise ValueError(f"carbon_pct must be in [0, 15], got {self.carbon_pct}")
        if self.bulk_density_kg_per_m3 is not None:
            if not math.isfinite(self.bulk_density_kg_per_m3) or not (500.0 <= self.bulk_density_kg_per_m3 <= 7000.0):
                raise ValueError(
                    f"bulk_density_kg_per_m3 must be in [500, 7000] kg/m^3, got {self.bulk_density_kg_per_m3}"
                )

        if self.primary_element == "C" and self.carbon_pct is not None:
            raise ValueError("carbon_pct duplicates primary_element='C'.")
        total_pct = self.primary_element_pct + (self.carbon_pct or 0.0)
        if self.secondary_elements is not None:
            cleaned: dict[str, float] = {}
            for elem, pct in self.secondary_elements.items():
                if not isinstance(elem, str) or not elem.strip():
                    raise ValueError("Secondary element names must be non-empty strings.")
                elem = elem.strip()
                if elem in cleaned:
                    raise ValueError(f"Duplicate secondary element {elem!r} after trimming whitespace.")
                if elem == "C" and self.carbon_pct is not None:
                    raise ValueError("Specify carbon using carbon_pct or secondary_elements['C'], not both.")
                val = float(pct)
                if not math.isfinite(val) or not (0.0 <= val <= 100.0):
                    raise ValueError(f"Secondary element {elem!r} must be in [0, 100], got {pct}")
                if str(elem) == self.primary_element:
                    raise ValueError(f"Secondary element {elem!r} duplicates primary_element.")
                cleaned[str(elem)] = val
            total_pct += sum(cleaned.values())
            object.__setattr__(self, "secondary_elements", MappingProxyType(cleaned))

        if total_pct > 100.0 + 1e-9:
            raise ValueError(f"Sum of specified element percentages ({total_pct:.2f}%) cannot exceed 100%.")

    def __hash__(self) -> int:
        sec_key = (
            tuple(sorted(self.secondary_elements.items())) if self.secondary_elements is not None else None
        )
        return hash((
            self.name,
            self.primary_element,
            self.primary_element_pct,
            self.basis_grade_pct,
            self.carbon_pct,
            sec_key,
            self.bulk_density_kg_per_m3,
        ))

    # --- Metallurgical & Stowage Properties ---
    @property
    def primary_mass_fraction(self) -> float:
        """Mass fraction of the primary alloying element (e.g. 0.60 for FeMo 60%)."""
        return self.primary_element_pct / 100.0

    @property
    def implied_iron_balance_pct(self) -> float:
        """Approximate iron (Fe) + trace balance (wt%)."""
        sec_sum = sum(self.secondary_elements.values()) if self.secondary_elements else 0.0
        return (max(0.0, 100.0 - self.primary_element_pct - sec_sum - (self.carbon_pct or 0.0)))

    @property
    def stowage_factor_m3_per_mt(self) -> float | None:
        """Bulk vessel/container stowage factor (m^3 per mt), or None if bulk density unset."""
        if self.bulk_density_kg_per_m3 is None:
            return None
        return (1000.0 / self.bulk_density_kg_per_m3)

    def _element_fraction(self, element: str | None) -> float:
        if element is None or element == self.primary_element:
            return self.primary_mass_fraction
        if self.secondary_elements is not None and element in self.secondary_elements:
            frac = self.secondary_elements[element] / 100.0
            if frac <= 0.0:
                raise ValueError(f"{self.name!r} contains 0% {element!r}.")
            return frac
        raise KeyError(f"Element {element!r} is not defined on {self.name!r}.")

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging bulk [volume] and [mass]."""
        return create_ferroalloy_context(self, name=context_name)

    # --- Physical Mass, Volume & Ladle Dosing Conversions ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def alloy_to_contained_metal(
        self,
        alloy_mass: pint.Quantity,
        *,
        element: str | None = None,
    ) -> pint.Quantity:
        """Converts gross ferroalloy mass (mt, kg, lb) to contained metal mass (metric_ton)."""
        mt_q = self._validate_mass(alloy_mass)
        frac = self._element_fraction(element)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def contained_metal_to_alloy(
        self,
        metal_mass: pint.Quantity,
        *,
        element: str | None = None,
    ) -> pint.Quantity:
        """Converts target contained metal mass to required gross ferroalloy mass (metric_ton)."""
        mt_q = self._validate_mass(metal_mass)
        frac = self._element_fraction(element)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def bulk_volume_to_mass(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts bulk stowage volume (m^3) to gross ferroalloy mass (metric_ton)."""
        if not isinstance(bulk_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bulk_volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = bulk_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_bulk_volume(self, alloy_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross ferroalloy mass (mt) to bulk stowage volume (m^3)."""
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        mt_q = self._validate_mass(alloy_mass)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def ladle_addition_required(
        self,
        steel_mass: pint.Quantity,
        target_pickup_pct: float,
        *,
        recovery_rate: float = 0.85,
        element: str | None = None,
    ) -> pint.Quantity:
        """Computes ferroalloy mass (kg) needed to raise a steel melt by target_pickup_pct.

        Uses initial steel mass, neglecting alloy addition to final bath mass.

        Args:
            steel_mass: Liquid steel ladle mass (e.g. Q_(100, 'metric_ton')).
            target_pickup_pct: Desired elemental chemistry increase in wt% (e.g. 0.50 for +0.50% Mn).
            recovery_rate: Furnace/ladle metallurgical recovery fraction in (0, 1].
            element: Target alloying element (defaults to primary_element).
        """
        pickup = float(target_pickup_pct)
        rec = float(recovery_rate)
        if not math.isfinite(pickup) or not (0.0 < pickup <= 100.0):
            raise ValueError(f"target_pickup_pct must be in (0, 100], got {target_pickup_pct}")
        if not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError(f"recovery_rate must be in (0, 1], got {recovery_rate}")

        steel_kg = self._validate_mass(steel_mass).to(ureg.kilogram)
        frac = self._element_fraction(element)
        req_alloy_ratio = (pickup / 100.0) / (frac * rec)
        scale = _scale_factor(steel_kg.magnitude, req_alloy_ratio)
        return ureg.Quantity(steel_kg.magnitude * scale, ureg.kilogram)

    # --- Commercial Pricing Parity ($/mt Gross <-> $/kg & $/lb Contained, Pro-Rata) ---
    def gross_mt_price_to_contained_kg_price(self, price_per_mt_gross: float) -> float:
        """Converts $/mt gross alloy quote to $/kg contained primary metal (e.g. $/kg V)."""
        p = float(price_per_mt_gross)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_gross must be finite, got {price_per_mt_gross}")
        return (p / (1000.0 * self.primary_mass_fraction))

    def contained_kg_price_to_gross_mt_price(self, price_per_kg_contained: float) -> float:
        """Converts $/kg contained primary metal quote to $/mt gross alloy."""
        p = float(price_per_kg_contained)
        if not math.isfinite(p):
            raise ValueError(f"price_per_kg_contained must be finite, got {price_per_kg_contained}")
        return (p * 1000.0 * self.primary_mass_fraction)

    def gross_mt_price_to_contained_lb_price(self, price_per_mt_gross: float) -> float:
        """Converts $/mt gross alloy quote to $/lb contained primary metal (e.g. $/lb Mo, $/lb Cr)."""
        p = float(price_per_mt_gross)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_gross must be finite, got {price_per_mt_gross}")
        return (p / (_LB_PER_MT * self.primary_mass_fraction))

    def contained_lb_price_to_gross_mt_price(self, price_per_lb_contained: float) -> float:
        """Converts $/lb contained primary metal quote to $/mt gross alloy."""
        p = float(price_per_lb_contained)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb_contained must be finite, got {price_per_lb_contained}")
        return (p * _LB_PER_MT * self.primary_mass_fraction)

    @staticmethod
    def contained_lb_price_to_contained_kg_price(price_per_lb: float) -> float:
        """Converts $/lb contained metal (US) to $/kg contained metal (Europe/Asia)."""
        p = float(price_per_lb)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb must be finite, got {price_per_lb}")
        return (p * _LB_PER_KG)

    @staticmethod
    def contained_kg_price_to_contained_lb_price(price_per_kg: float) -> float:
        """Converts $/kg contained metal (Europe/Asia) to $/lb contained metal (US)."""
        p = float(price_per_kg)
        if not math.isfinite(p):
            raise ValueError(f"price_per_kg must be finite, got {price_per_kg}")
        return (p / _LB_PER_KG)

    def pro_rata_price(
        self,
        benchmark_price_per_mt: float,
        basis_grade_pct: float | None = None,
    ) -> float:
        """Computes pro-rata $/mt settlement price against a contract basis grade (e.g. FeSi 75%)."""
        p = float(benchmark_price_per_mt)
        basis = self.basis_grade_pct if basis_grade_pct is None else float(basis_grade_pct)
        if basis is None:
            raise ValueError("Provide basis_grade_pct or set it on the FerroalloyAssay.")
        if not math.isfinite(p) or not math.isfinite(basis) or not (0.0 < basis <= 100.0):
            raise ValueError("benchmark_price_per_mt must be finite and basis_grade_pct in (0, 100].")
        return (p * (self.primary_element_pct / basis))

    def blend_with(
        self,
        other: FerroalloyAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> FerroalloyAssay:
        """Blends two ferroalloy lots of the same primary element by gross mass."""
        if not isinstance(other, FerroalloyAssay):
            raise TypeError(f"Expected FerroalloyAssay, got {type(other).__name__}")
        if self.primary_element != other.primary_element:
            raise ValueError(
                f"Cannot blend ferroalloys with different primary elements: "
                f"{self.primary_element!r} vs {other.primary_element!r}"
            )
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1

        b_sec: dict[str, float] | None = None
        if self.secondary_elements is not None or other.secondary_elements is not None:
            d1 = self.secondary_elements or {}
            d2 = other.secondary_elements or {}
            b_sec = {
                k: (d1.get(k, 0.0) * w1) + (d2.get(k, 0.0) * w2)
                for k in sorted(set(d1) | set(d2))
            }

        b_carbon = (
            (self.carbon_pct * w1) + (other.carbon_pct * w2)
            if (self.carbon_pct is not None and other.carbon_pct is not None)
            else None
        )

        b_rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None and other.bulk_density_kg_per_m3 is not None:
            vol_per_kg = (w1 / self.bulk_density_kg_per_m3) + (w2 / other.bulk_density_kg_per_m3)
            b_rho = (1.0 / vol_per_kg)

        return FerroalloyAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            primary_element=self.primary_element,
            primary_element_pct=(self.primary_element_pct * w1) + (other.primary_element_pct * w2),
            basis_grade_pct=self.basis_grade_pct if self.basis_grade_pct == other.basis_grade_pct else None,
            carbon_pct=b_carbon,
            secondary_elements=b_sec,
            bulk_density_kg_per_m3=b_rho,
        )

    # --- Global Ferroalloy Benchmark Presets ---
    @classmethod
    def hc_ferrochrome(cls) -> FerroalloyAssay:
        """High-Carbon / Charge Ferrochrome (52% Cr basis, ~7% C — South African benchmark)."""
        return cls(
            name="HC / Charge Ferrochrome (52% Cr)",
            primary_element="Cr",
            primary_element_pct=52.0,
            basis_grade_pct=52.0,
            carbon_pct=7.0,
            bulk_density_kg_per_m3=3600.0,
        )

    @classmethod
    def lc_ferrochrome(cls) -> FerroalloyAssay:
        """Low-Carbon Ferrochrome (65% Cr, 0.10% max C — stainless steel refining)."""
        return cls(
            name="LC Ferrochrome (65% Cr)",
            primary_element="Cr",
            primary_element_pct=65.0,
            basis_grade_pct=65.0,
            carbon_pct=0.10,
            bulk_density_kg_per_m3=3800.0,
        )

    @classmethod
    def hc_ferromanganese(cls) -> FerroalloyAssay:
        """High-Carbon Ferromanganese (78% Mn basis, ~7.5% C — Europe/US benchmark)."""
        return cls(
            name="HC Ferromanganese (78% Mn)",
            primary_element="Mn",
            primary_element_pct=78.0,
            basis_grade_pct=78.0,
            carbon_pct=7.5,
            bulk_density_kg_per_m3=3800.0,
        )

    @classmethod
    def silicomanganese(cls) -> FerroalloyAssay:
        """Silicomanganese 65/16 (65% Mn basis, 16% Si, 2.0% C)."""
        return cls(
            name="Silicomanganese (65% Mn / 16% Si)",
            primary_element="Mn",
            primary_element_pct=65.0,
            basis_grade_pct=65.0,
            carbon_pct=2.0,
            secondary_elements={"Si": 16.0},
            bulk_density_kg_per_m3=3300.0,
        )

    @classmethod
    def ferrosilicon_75(cls) -> FerroalloyAssay:
        """Ferrosilicon 75 (75% Si basis — Europe/Asia/US steel deoxidiser benchmark)."""
        return cls(
            name="Ferrosilicon (75% Si)",
            primary_element="Si",
            primary_element_pct=75.0,
            basis_grade_pct=75.0,
            carbon_pct=0.10,
            bulk_density_kg_per_m3=1850.0,
        )

    @classmethod
    def ferromolybdenum_60(cls) -> FerroalloyAssay:
        """Ferromolybdenum (60% Mo — quoted in $/lb Mo or $/kg Mo)."""
        return cls(
            name="Ferromolybdenum (60% Mo)",
            primary_element="Mo",
            primary_element_pct=60.0,
            basis_grade_pct=60.0,
            carbon_pct=0.10,
            bulk_density_kg_per_m3=4200.0,
        )

    @classmethod
    def ferrovanadium_80(cls) -> FerroalloyAssay:
        """Ferrovanadium 80 (80% V — European $/kg V & US $/lb V benchmark)."""
        return cls(
            name="Ferrovanadium (80% V)",
            primary_element="V",
            primary_element_pct=80.0,
            basis_grade_pct=80.0,
            carbon_pct=0.20,
            bulk_density_kg_per_m3=3500.0,
        )

    @classmethod
    def ferroniobium_65(cls) -> FerroalloyAssay:
        """Ferroniobium (65% Nb — HSLA structural steel micro-alloy)."""
        return cls(
            name="Ferroniobium (65% Nb)",
            primary_element="Nb",
            primary_element_pct=65.0,
            basis_grade_pct=65.0,
        )



# ---------------------------------------------------------------------------
# Non-Ferrous, Superalloy & Industrial Metal Alloy Specification
# (Brass, Bronze, Aluminium 6061/7075, Monel 400, Inconel 625, Ti-6Al-4V)
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class MetalAlloyAssay:
    """Immutable metallurgical and valuation specification for general metal alloys.

    Bridges solid volume (m^3, cm^3), gross alloy mass (mt, kg, lb), contained
    element mass across arbitrary alloying elements, remelt recovery yields,
    nominal piece/billet/ingot counts, and multi-element melt valuation.
    """

    name: str
    composition_pct: Mapping[str, float]
    solid_density_kg_per_m3: float
    base_metal: str | None = None
    remelt_yield_pct: float = 100.0
    nominal_piece_mass_kg: float | None = None

    def __post_init__(self) -> None:
        if not self.composition_pct:
            raise ValueError("composition_pct must be a non-empty mapping of element -> wt%.")

        cleaned: dict[str, float] = {}
        for elem, pct in self.composition_pct.items():
            if not isinstance(elem, str):
                raise ValueError("Element symbols in composition_pct must be strings.")
            sym = elem.strip()
            val = float(pct)
            if not sym:
                raise ValueError("Element symbol in composition_pct cannot be empty.")
            if not math.isfinite(val) or not (0.0 <= val <= 100.0):
                raise ValueError(f"Weight percentage for {sym!r} must be in [0, 100], got {pct}")
            if sym in cleaned:
                raise ValueError(f"Duplicate element symbol {sym!r} in composition_pct.")
            cleaned[sym] = val

        total_pct = sum(cleaned.values())
        if total_pct <= 0.0:
            raise ValueError("Total composition_pct must be strictly positive.")
        if total_pct > 100.0 + 1e-9:
            raise ValueError(f"Total composition_pct ({total_pct:.4f}%) cannot exceed 100%.")

        if not math.isfinite(self.solid_density_kg_per_m3) or not (1000.0 <= self.solid_density_kg_per_m3 <= 22000.0):
            raise ValueError(
                f"solid_density_kg_per_m3 must be in [1000, 22000] kg/m^3, got {self.solid_density_kg_per_m3}"
            )
        if not math.isfinite(self.remelt_yield_pct) or not (0.0 < self.remelt_yield_pct <= 100.0):
            raise ValueError(f"remelt_yield_pct must be in (0, 100], got {self.remelt_yield_pct}")

        if self.base_metal is not None and not isinstance(self.base_metal, str):
            raise ValueError("base_metal must be an element symbol string or None.")
        resolved_base = self.base_metal.strip() if self.base_metal is not None else max(cleaned, key=cleaned.get)
        if resolved_base not in cleaned:
            raise ValueError(f"base_metal {resolved_base!r} is not present in composition_pct.")
        if cleaned[resolved_base] <= 0.0:
            raise ValueError(f"base_metal {resolved_base!r} must have a positive weight percentage.")
        if self.nominal_piece_mass_kg is not None:
            if not math.isfinite(self.nominal_piece_mass_kg) or self.nominal_piece_mass_kg <= 0.0:
                raise ValueError(
                    f"nominal_piece_mass_kg must be strictly positive, got {self.nominal_piece_mass_kg}"
                )

        object.__setattr__(self, "composition_pct", MappingProxyType(cleaned))
        object.__setattr__(self, "base_metal", resolved_base)

    def __hash__(self) -> int:
        return hash((
            self.name,
            tuple(sorted(self.composition_pct.items())),
            self.solid_density_kg_per_m3,
            self.base_metal,
            self.remelt_yield_pct,
            self.nominal_piece_mass_kg,
        ))

    # --- Alloy Composition & Physical Properties ---
    @property
    def base_metal_pct(self) -> float:
        """Weight percentage of the primary/base metal in the alloy."""
        return self.composition_pct[self.base_metal]  # type: ignore[index]

    @property
    def balance_other_pct(self) -> float:
        """Unlisted/impurity/trace balance (wt%) to 100%."""
        return (max(0.0, 100.0 - sum(self.composition_pct.values())))

    @property
    def remelt_yield_fraction(self) -> float:
        """Remelt recovery yield as a decimal fraction in (0, 1]."""
        return self.remelt_yield_pct / 100.0

    def element_fraction(self, element: str) -> float:
        """Returns mass fraction of `element` in the gross alloy."""
        key = str(element).strip()
        if key not in self.composition_pct:
            raise KeyError(f"Element {key!r} is not defined in {self.name!r}.")
        frac = self.composition_pct[key] / 100.0
        if frac <= 0.0:
            raise ValueError(f"{self.name!r} contains 0% {key!r}.")
        return frac

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging solid [volume] and gross [mass]."""
        return create_metal_alloy_context(self, name=context_name)

    # --- Physical Volume <-> Mass <-> Contained Element <-> Remelt & Piece Conversions ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def solid_volume_to_mass(self, solid_volume: pint.Quantity) -> pint.Quantity:
        """Converts solid metal volume (m^3, cm^3, L) to gross alloy mass (metric_ton)."""
        if not isinstance(solid_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(solid_volume).__name__}")
        m3_q = solid_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.solid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_solid_volume(self, alloy_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross alloy mass (mt, kg, lb) to solid metal volume (m^3)."""
        mt_q = self._validate_mass(alloy_mass)
        mt_per_m3 = self.solid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def alloy_to_contained_element(self, alloy_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Converts gross alloy mass (mt, kg, lb) to contained `element` mass (metric_ton)."""
        mt_q = self._validate_mass(alloy_mass)
        frac = self.element_fraction(element)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def contained_element_to_alloy(self, element_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Converts target contained `element` mass to required gross alloy mass (metric_ton)."""
        mt_q = self._validate_mass(element_mass)
        frac = self.element_fraction(element)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def gross_to_recovered_melt_mass(self, feed_mass: pint.Quantity) -> pint.Quantity:
        """Converts charge/scrap mass (mt, kg) to recoverable remelted ingot/billet mass (metric_ton)."""
        mt_q = self._validate_mass(feed_mass)
        scale = _scale_factor(mt_q.magnitude, self.remelt_yield_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def recovered_melt_to_gross_mass(self, target_melt_mass: pint.Quantity) -> pint.Quantity:
        """Converts target recovered melt mass (mt, kg) to required charge/scrap mass (metric_ton)."""
        mt_q = self._validate_mass(target_melt_mass)
        scale = _scale_factor(mt_q.magnitude, self.remelt_yield_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def mass_to_piece_count(self, alloy_mass: pint.Quantity) -> float:
        """Converts gross alloy mass to nominal ingot/billet/piece count."""
        if self.nominal_piece_mass_kg is None:
            raise ValueError(f"{self.name!r} does not define nominal_piece_mass_kg.")
        kg_val = float(self._validate_mass(alloy_mass).to(ureg.kilogram).magnitude)
        return (kg_val / self.nominal_piece_mass_kg)

    def piece_count_to_mass(self, piece_count: float) -> pint.Quantity:
        """Converts nominal piece/ingot/billet count to gross alloy mass (kilogram)."""
        if self.nominal_piece_mass_kg is None:
            raise ValueError(f"{self.name!r} does not define nominal_piece_mass_kg.")
        n = float(piece_count)
        if not math.isfinite(n) or n < 0.0:
            raise ValueError(f"piece_count must be finite and non-negative, got {piece_count}")
        return ureg.Quantity(n * self.nominal_piece_mass_kg, ureg.kilogram)

    # --- Commercial Melt Valuation & Unit Pricing Parity ---
    def gross_melt_value_per_mt(
        self,
        element_prices_usd_per_mt: Mapping[str, float],
        *,
        apply_remelt_yield: bool = False,
        treatment_deduction_usd_per_mt: float = 0.0,
    ) -> float:
        """Computes gross alloy value ($/mt) from per-element $/mt metal quotes.

        Unquoted elements and the unlisted balance contribute zero value.
        Remelt yield is a uniform bulk recovery estimate.

        Args:
            element_prices_usd_per_mt: Mapping of element symbol -> $/mt contained metal price.
            apply_remelt_yield: If True, scales contained value by remelt_yield_fraction.
            treatment_deduction_usd_per_mt: Fixed tolling/refining/conversion deduction ($/mt gross).
        """
        deduction = float(treatment_deduction_usd_per_mt)
        if not math.isfinite(deduction):
            raise ValueError(f"treatment_deduction_usd_per_mt must be finite, got {treatment_deduction_usd_per_mt}")

        yield_mult = self.remelt_yield_fraction if apply_remelt_yield else 1.0
        gross_val = 0.0
        for elem, pct in self.composition_pct.items():
            if elem in element_prices_usd_per_mt and pct > 0.0:
                price = float(element_prices_usd_per_mt[elem])
                if not math.isfinite(price):
                    raise ValueError(f"Price for {elem!r} must be finite, got {element_prices_usd_per_mt[elem]}")
                gross_val += (pct / 100.0) * yield_mult * price

        return (gross_val - deduction)

    @staticmethod
    def lb_price_to_mt_price(price_per_lb: float) -> float:
        """Converts $/lb quote to $/mt equivalent."""
        p = float(price_per_lb)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb must be finite, got {price_per_lb}")
        return (p * _LB_PER_MT)

    @staticmethod
    def mt_price_to_lb_price(price_per_mt: float) -> float:
        """Converts $/mt quote to $/lb equivalent."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / _LB_PER_MT)

    def blend_with(
        self,
        other: MetalAlloyAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> MetalAlloyAssay:
        """Blends by gross mass; unlisted elements are zero, solid volumes additive."""
        if not isinstance(other, MetalAlloyAssay):
            raise TypeError(f"Expected MetalAlloyAssay, got {type(other).__name__}")
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1

        all_elems = sorted(set(self.composition_pct) | set(other.composition_pct))
        blended_comp = {
            elem: (self.composition_pct.get(elem, 0.0) * w1) + (other.composition_pct.get(elem, 0.0) * w2)
            for elem in all_elems
        }

        vol_per_kg = (w1 / self.solid_density_kg_per_m3) + (w2 / other.solid_density_kg_per_m3)
        b_base = self.base_metal if self.base_metal == other.base_metal else None

        return MetalAlloyAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            composition_pct=blended_comp,
            solid_density_kg_per_m3=(1.0 / vol_per_kg),
            base_metal=b_base,
            remelt_yield_pct=(self.remelt_yield_pct * w1) + (other.remelt_yield_pct * w2),
            nominal_piece_mass_kg=None,
        )

    # --- Global Metal Alloy Benchmark Presets ---
    @classmethod
    def cartridge_brass_c260(cls) -> MetalAlloyAssay:
        """Cartridge Brass (UNS C26000 — 70% Cu, 30% Zn)."""
        return cls(
            name="Cartridge Brass (C26000, 70Cu/30Zn)",
            composition_pct={"Cu": 70.0, "Zn": 30.0},
            solid_density_kg_per_m3=8530.0,
            base_metal="Cu",
        )

    @classmethod
    def aluminium_6061(cls) -> MetalAlloyAssay:
        """Aluminium Alloy 6061 (structural/extrusion grade — ~97.9% Al, 1.0% Mg, 0.6% Si)."""
        return cls(
            name="Aluminium Alloy 6061",
            composition_pct={"Al": 97.9, "Mg": 1.0, "Si": 0.6},
            solid_density_kg_per_m3=2700.0,
            base_metal="Al",
            nominal_piece_mass_kg=22.68,  # ~50 lb standard ingot
        )

    @classmethod
    def aluminium_7075(cls) -> MetalAlloyAssay:
        """Aluminium Alloy 7075 (high-strength aerospace grade — ~90.0% Al, 5.6% Zn, 2.5% Mg, 1.6% Cu)."""
        return cls(
            name="Aluminium Alloy 7075",
            composition_pct={"Al": 90.0, "Zn": 5.6, "Mg": 2.5, "Cu": 1.6},
            solid_density_kg_per_m3=2810.0,
            base_metal="Al",
        )

    @classmethod
    def monel_400(cls) -> MetalAlloyAssay:
        """Monel 400 Nickel-Copper Alloy (UNS N04400 — ~66.5% Ni, 31.5% Cu)."""
        return cls(
            name="Monel 400 (N04400, 66.5Ni/31.5Cu)",
            composition_pct={"Ni": 66.5, "Cu": 31.5},
            solid_density_kg_per_m3=8800.0,
            base_metal="Ni",
        )

    @classmethod
    def inconel_625(cls) -> MetalAlloyAssay:
        """Inconel 625 Nickel-Chromium Superalloy (UNS N06625 — ~61.0% Ni, 21.5% Cr, 9.0% Mo)."""
        return cls(
            name="Inconel 625 (N06625, Ni-Cr-Mo)",
            composition_pct={"Ni": 61.0, "Cr": 21.5, "Mo": 9.0},
            solid_density_kg_per_m3=8440.0,
            base_metal="Ni",
        )

    @classmethod
    def titanium_6al_4v(cls) -> MetalAlloyAssay:
        """Titanium Ti-6Al-4V (Grade 5 aerospace/medical alloy — ~90.0% Ti, 6.0% Al, 4.0% V)."""
        return cls(
            name="Titanium Ti-6Al-4V (Grade 5)",
            composition_pct={"Ti": 90.0, "Al": 6.0, "V": 4.0},
            solid_density_kg_per_m3=4430.0,
            base_metal="Ti",
        )


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
    key = basis.strip().lower()
    if key in {"nar", "lhv", "ncv"}:
        gj_per_mt = assay.gj_per_mt_nar
    elif key in {"gar", "hhv", "gcv"}:
        gj_per_mt = assay.gj_per_mt_gar
    else:
        raise ValueError("Coal context basis must be 'nar'/'lhv'/'ncv' or 'gar'/'hhv'/'gcv'.")
    ctx = pint.Context(name or _default_ctx_name(assay.name, key))

    def mass_to_energy(_u: pint.UnitRegistry, mass: pint.Quantity) -> pint.Quantity:
        tonnes = mass.to(ureg.metric_ton)
        factor = _scale_factor(tonnes.magnitude, gj_per_mt)
        return ureg.Quantity(tonnes.magnitude * factor, ureg.gigajoule)

    def energy_to_mass(_u: pint.UnitRegistry, energy: pint.Quantity) -> pint.Quantity:
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
) -> pint.Quantity:
    """Converts a liquid hydrocarbon quantity (volume, mass, flow rate, or price) to target_unit."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")

    resolved_target = parse_unit(target_unit)
    if quantity.dimensionality == resolved_target.dimensionality:
        return quantity.to(resolved_target)

    ctx = create_hydrocarbon_context(
        assay=assay,
        api_gravity=api_gravity,
        specific_gravity=specific_gravity,
        bbl_per_mt=bbl_per_mt,
    )
    with ureg.context(ctx):
        return _convert_contextual_quantity(quantity, resolved_target, [ctx])


def convert_gas_energy(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: GasAssay | None = None,
    btu_per_scf: float | None = None,
    basis: str = "gcv",
) -> pint.Quantity:
    """Converts a natural gas quantity between reference volumes, energy, flow rates, or prices."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")

    resolved_target = parse_unit(target_unit)
    if quantity.dimensionality == resolved_target.dimensionality:
        return quantity.to(resolved_target)

    ctx = create_gas_heating_context(assay=assay, btu_per_scf=btu_per_scf, basis=basis)
    with ureg.context(ctx):
        return _convert_contextual_quantity(quantity, resolved_target, [ctx])


def convert_hydrogen(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: HydrogenAssay,
    basis: str = "lhv",
    mode: str = "auto",
    cracking_efficiency: float | None = None,
) -> pint.Quantity:
    """Converts hydrogen or carrier quantities across volume, mass, Nm^3, and energy (MWh, GJ, MMBtu)."""
    if not isinstance(assay, HydrogenAssay):
        raise TypeError(f"Expected HydrogenAssay, got {type(assay).__name__}")
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")

    resolved_target = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        ctx = assay.as_context(basis=basis)
        with ureg.context(ctx):
            return _convert_contextual_quantity(quantity, resolved_target, [ctx])

    dispatch = {
        "carrier_to_h2": lambda q: assay.carrier_to_h2_mass(q, cracking_efficiency=cracking_efficiency),
        "h2_to_carrier": lambda q: assay.h2_to_carrier_mass(q, cracking_efficiency=cracking_efficiency),
        "h2_to_nm3": lambda q: HydrogenAssay.h2_mass_to_nm3(q),
        "nm3_to_h2": lambda q: HydrogenAssay.nm3_to_h2_mass(q),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported hydrogen mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(resolved_target)


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
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")

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
STANDARD_ATMOSPHERE_PSI: Final[float] = 14.6959487755134
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
    sea_pa = float(sea_level_pa)
    if not math.isfinite(sea_pa) or sea_pa <= 0.0:
        raise ValueError("sea_level_pa must be finite and strictly positive.")
    height = elevation if isinstance(elevation, pint.Quantity) else ureg.Quantity(elevation, unit or elevation_unit)
    h_m = float(height.to(ureg.meter).magnitude)
    if not math.isfinite(h_m) or not (-500.0 <= h_m <= 11000.0):
        raise ValueError("Elevation in metres must be finite and in [-500, 11000].")
    air_specific_gas_constant = UNIVERSAL_GAS_CONSTANT_J_MOL_K / (DRY_AIR_MOLAR_MASS_KG_KMOL / 1000.0)
    exponent = STANDARD_GRAVITY_M_S2 / (air_specific_gas_constant * BAROMETRIC_LAPSE_RATE_K_M)
    pressure = sea_pa * (1.0 - BAROMETRIC_LAPSE_RATE_K_M * h_m / STANDARD_SEA_LEVEL_TEMP_K) ** exponent
    return ureg.Quantity(pressure, ureg.pascal)


@authorized_conversion
def convert_pressure(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    p_atm: pint.Quantity | float | Decimal | None = None,
    atmospheric_pressure: pint.Quantity | float | Decimal | None = None,
    neutral_is_absolute: bool = True,
    allow_negative_absolute: bool = False,
) -> pint.Quantity:
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
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    resolved_target = parse_unit(target_unit)
    pressure_dim = ureg.pascal.dimensionality
    if quantity.dimensionality != pressure_dim:
        raise pint.DimensionalityError(quantity.units, ureg.pascal, quantity.dimensionality, pressure_dim)
    if resolved_target.dimensionality != pressure_dim:
        raise pint.DimensionalityError(quantity.units, resolved_target, pressure_dim, resolved_target.dimensionality)

    src_datum = _classify_pressure_datum(quantity.units)
    dst_datum = _classify_pressure_datum(resolved_target)
    if not neutral_is_absolute and "neutral" in (src_datum, dst_datum):
        raise ValueError("Unqualified pressure units require neutral_is_absolute=True or explicit datum labels.")
    src_datum = "absolute" if src_datum == "neutral" else src_datum
    dst_datum = "absolute" if dst_datum == "neutral" else dst_datum

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
        raise ValueError("Absolute pressure cannot be negative.")
    if src_datum == dst_datum:
        result = quantity.to(resolved_target)
    else:
        shifted = absolute_magnitude if dst_datum == "absolute" else absolute_magnitude - atm_mag
        result = ureg.Quantity(shifted, ureg.pascal).to(resolved_target)
    if is_decimal and not isinstance(result.magnitude, Decimal):
        return ureg.Quantity(Decimal(str(result.magnitude)), result.units)
    return result


def convert_pressure_datum(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    atmospheric_pressure: pint.Quantity | float | Decimal | None = None,
    p_atm: pint.Quantity | float | Decimal | None = None,
    neutral_is_absolute: bool = True,
    allow_negative_absolute: bool = False,
) -> pint.Quantity:
    """Backward-compatible alias for convert_pressure(), including its options."""
    return convert_pressure(
        quantity, target_unit, p_atm=p_atm, atmospheric_pressure=atmospheric_pressure,
        neutral_is_absolute=neutral_is_absolute, allow_negative_absolute=allow_negative_absolute,
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
    if not isinstance(p1, pint.Quantity) or not isinstance(p2, pint.Quantity):
        raise TypeError("p1 and p2 must be Pint quantities.")
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
    iso = _numerator_currency_iso(ureg.Quantity(1, currency_unit))
    digits = CURRENCY_MINOR_DIGITS.get(iso, 2)
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
) -> pint.Quantity:
    """Converts a currency amount or cross-currency commodity quote to target_currency.

    Preserves physical denominators when converting commodity price quotes
    (e.g. converting "$75.50/bbl" directly into "EUR / oil_barrel" or "GBP_penny / therm_it"
    into "USD / therm_it").

    Parameters
    ----------
    quantity : pint.Quantity
        Source financial amount (e.g. Q_(1000, "EUR")) or quote (e.g. Q_(75.50, "USD / oil_barrel")).
    target_currency : str
        Destination ISO currency code or registered minor unit (e.g. "GBP", "JPY", "USD").
    rates : dict[str, Decimal | float | int]
        Spot rates per 1 unit of base_currency.
    base_currency : str
        Anchor currency for the rates dictionary (default: "USD").
    quantize_settlement : bool
        If True, applies Banker's Rounding (ROUND_HALF_EVEN) using CURRENCY_MINOR_DIGITS
        for the target currency. Float input is converted via its decimal string
        before FX multiplication so settlement rounding avoids float arithmetic.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")

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
        raise ValueError(f"Missing FX rate for source currency {src_iso!r} relative to base {base!r}.")
    if target_iso not in clean_rates:
        raise ValueError(f"Missing FX rate for target currency {target_iso!r} relative to base {base!r}.")

    cross_rate_dec = clean_rates[target_iso] / clean_rates[src_iso]

    conversion_input = quantity
    if quantize_settlement and not isinstance(quantity.magnitude, Decimal):
        conversion_input = ureg.Quantity(_as_decimal(quantity.magnitude), quantity.units)
    converted = _replace_numerator_currency(conversion_input, src_iso, target_unit_obj, cross_rate_dec)

    if quantize_settlement:
        return _quantize_currency_settlement(converted, target_unit_obj)

    return converted


def triangulate_fx_rate(rates: Mapping[tuple[str, str] | str, Any], base: str, quote: str) -> Decimal:
    """Triangulate an FX cross rate from a sparse rate map.
    
    Rates mean 1 BASE = rate QUOTE. Direct, inverse, and multi-hop paths are supported.
    """
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

    raise ValueError(f"No FX triangulation path found from {base_c} to {quote_c}.")


def convert_fx(
    amount: pint.Quantity | float | Decimal,
    target_currency: str,
    rates: Mapping[tuple[str, str] | str, Any],
    *,
    source_currency: str | None = None,
) -> pint.Quantity:
    """Converts currency amounts or commodity quotes using sparse FX rates.

    Quantity inputs infer their currency dimension and retain physical unit
    factors. Registered minor units are supported on both sides. An explicit
    source_currency must agree with a quantity's actual currency dimension.
    """
    tgt_unit = parse_unit(target_currency)
    target_dims = tgt_unit.dimensionality
    if len(target_dims) != 1:
        raise ValueError("target_currency must be an isolated currency unit.")
    tgt = _numerator_currency_iso(ureg.Quantity(1, tgt_unit))

    if isinstance(amount, pint.Quantity):
        src = _numerator_currency_iso(amount)
        if source_currency is not None and _normalize_currency(source_currency) != src:
            raise ValueError("source_currency does not match the quantity's currency dimension.")
        rate = triangulate_fx_rate(rates, src, tgt)
        return _replace_numerator_currency(amount, src, tgt_unit, rate)

    if source_currency is None:
        raise ValueError("source_currency is required when amount is not a pint.Quantity.")
    magnitude = amount if isinstance(amount, Decimal) else float(amount)
    if not _as_decimal(magnitude).is_finite():
        raise ValueError(f"FX amount must be finite, got {amount}")
    src = _normalize_currency(source_currency)
    rate = triangulate_fx_rate(rates, src, tgt)
    return _replace_numerator_currency(ureg.Quantity(magnitude, src), src, tgt_unit, rate)


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
    factors = []
    for name, exponent in unit._units.items():
        candidate = ureg.Unit(name)
        if candidate.dimensionality == ureg.kilogram.dimensionality:
            tokens = _unit_tokens(candidate)
            if tokens & {"wmt", "wet_metric_ton", "wet_metric_tonne", "wet_tonne"}:
                basis = "wet"
            elif tokens & {"dmt", "dry_metric_ton", "dry_metric_tonne", "dry_tonne"}:
                basis = "dry"
            elif tokens & {"dmtu", "dltu", "dry_metric_tonne_unit", "dry_long_ton_unit", "mtu", "stu"}:
                basis = "contained"
            else:
                basis = None
            factors.append((candidate, exponent, basis))
    if len(factors) != 1 or factors[0][1] not in (-1, 1):
        raise ValueError("Ore conversions require one standalone mass factor with exponent +1 or -1.")
    return factors[0]


@authorized_conversion
def convert_wet_to_dry_mass(wet_mass: pint.Quantity, moisture_pct: float | Decimal) -> pint.Quantity:
    """Converts WMT to DMT: dry mass = wet mass * (1 - moisture/100)."""
    if not isinstance(wet_mass, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
    fraction = 1 - _validate_pct_in_range(moisture_pct, low_inclusive=0, high_exclusive=True) / 100
    mass = wet_mass.to(ureg.wmt)
    return ureg.Quantity(_ore_multiply(mass.magnitude, fraction), ureg.dmt)


@authorized_conversion
def convert_dry_to_wet_mass(dry_mass: pint.Quantity, moisture_pct: float | Decimal) -> pint.Quantity:
    """Converts DMT to WMT: wet mass = dry mass / (1 - moisture/100)."""
    if not isinstance(dry_mass, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(dry_mass).__name__}")
    fraction = 1 - _validate_pct_in_range(moisture_pct, low_inclusive=0, high_exclusive=True) / 100
    mass = dry_mass.to(ureg.dmt)
    magnitude = mass.magnitude
    _ore_multiply(magnitude, Decimal(1))  # Validate finiteness before division.
    scale = fraction if isinstance(magnitude, Decimal) else float(fraction)
    return ureg.Quantity(magnitude / scale, ureg.wmt)


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
    if not isinstance(price_quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(price_quantity).__name__}")
    moisture = _validate_pct_in_range(moisture_pct, low_inclusive=0, high_exclusive=True)
    _numerator_currency_iso(price_quantity)
    source_mass, exponent, source_basis = _ore_mass_factor(price_quantity.units)
    target = parse_unit(target_mass_unit)
    if exponent != -1 or target.dimensionality != ureg.kilogram.dimensionality:
        raise ValueError("Price must have a mass denominator and target_mass_unit must be a mass unit.")
    _, _, target_basis = _ore_mass_factor(target)
    if "contained" in (source_basis, target_basis):
        raise ValueError("Moisture prices require wet/dry ore mass units, not contained-metal units.")
    destination = price_quantity.units * source_mass / target
    return convert_bulk_ore(price_quantity, destination, moisture_fraction=moisture / 100)


@authorized_conversion
def extract_contained_metal(
    dry_ore_mass: pint.Quantity,
    grade_pct: float | Decimal,
    *,
    target_unit: str | pint.Unit = "mt",
) -> pint.Quantity:
    """Returns contained metal from dry ore; 1 DMTU represents 10 kg metal."""
    if not isinstance(dry_ore_mass, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(dry_ore_mass).__name__}")
    grade = _validate_pct_in_range(grade_pct, low_inclusive=0, high_exclusive=False) / 100
    dry = dry_ore_mass.to(ureg.dmt)
    metal = ureg.Quantity(_ore_multiply(dry.magnitude, grade), ureg.metric_ton)
    return metal.to(parse_unit(target_unit))


@authorized_conversion
def extract_dry_mass_from_contained(
    contained_metal: pint.Quantity,
    grade_pct: float | Decimal,
    *,
    target_unit: str | pint.Unit = "dmt",
) -> pint.Quantity:
    """Returns dry ore required to supply a contained metal mass at grade_pct."""
    if not isinstance(contained_metal, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(contained_metal).__name__}")
    grade = _validate_pct_in_range(grade_pct, low_inclusive=0, high_exclusive=False) / 100
    if grade == 0:
        raise ValueError("Cannot compute dry ore mass with grade_pct=0.")
    metal = contained_metal.to(ureg.metric_ton)
    _ore_multiply(metal.magnitude, Decimal(1))
    scale = grade if isinstance(metal.magnitude, Decimal) else float(grade)
    dry = ureg.Quantity(metal.magnitude / scale, ureg.dmt)
    return dry.to(parse_unit(target_unit))


@validated_assay
@dataclass(frozen=True, slots=True)
class BulkOreAssay:
    """Lightweight bulk ore moisture and component grades, expressed as fractions.

    Grades are independent assay measurements and need not sum to one.
    """

    name: str = "bulk_ore"
    moisture_fraction: float | Decimal | None = None
    grades: Mapping[str, float | Decimal] = field(default_factory=dict)

    def __post_init__(self) -> None:
        moisture = _as_decimal(self.moisture_fraction) if self.moisture_fraction is not None else None
        if moisture is not None and (not moisture.is_finite() or not (0 <= moisture < 1)):
            raise ValueError("moisture_fraction must be finite and in [0, 1).")
        cleaned = {}
        for component, grade in self.grades.items():
            value = _as_decimal(grade)
            if not value.is_finite() or not (0 <= value <= 1):
                raise ValueError(f"Grade for {component!r} must be finite and in [0, 1].")
            cleaned[str(component)] = grade if isinstance(grade, Decimal) else float(grade)
        object.__setattr__(self, "moisture_fraction", self.moisture_fraction if self.moisture_fraction is None or isinstance(self.moisture_fraction, Decimal) else float(self.moisture_fraction))
        object.__setattr__(self, "grades", MappingProxyType(cleaned))

    def __hash__(self) -> int:
        return hash((self.name, self.moisture_fraction, tuple(sorted(self.grades.items()))))

    @property
    def dry_fraction(self) -> float | Decimal:
        if self.moisture_fraction is None:
            raise ValueError("Explicit ore moisture is required for wet/dry conversion.")
        return 1 - self.moisture_fraction

    def wmt_to_dmt(self, wet_mass: float | Decimal) -> float | Decimal:
        return _ore_multiply(wet_mass, _as_decimal(self.dry_fraction))

    def dmt_to_wmt(self, dry_mass: float | Decimal) -> float | Decimal:
        _ore_multiply(dry_mass, Decimal(1))
        fraction = _as_decimal(self.dry_fraction)
        return dry_mass / (fraction if isinstance(dry_mass, Decimal) else float(fraction))

    def contained_metal_dmt(self, dry_mass: float | Decimal, component: str) -> float | Decimal:
        return _ore_multiply(dry_mass, _as_decimal(self.grades[component]))

    def contained_dmtu(self, dry_mass_dmt: float | Decimal, component: str = "Fe") -> float | Decimal:
        return self.contained_metal_dmt(dry_mass_dmt, component) * 100

    def to_ore_assay(self, component: str = "Fe") -> OreAssay:
        """Builds an OreAssay for one component, subject to its stricter limits.

        OreAssay requires a positive grade and moisture below 50%; independent
        component assays are not automatically treated as additive impurities.
        """
        return OreAssay(
            name=self.name, primary_element=component,
            grade_pct=float(_as_decimal(self.grades[component]) * 100),
            moisture_pct=float((1 - _as_decimal(self.dry_fraction)) * 100),
        )


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
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if isinstance(assay, OreAssay):
        if component != assay.primary_element:
            raise ValueError(f"OreAssay defines {assay.primary_element!r}, not {component!r}.")
        resolved = BulkOreAssay(assay.name, _as_decimal(assay.moisture_pct) / 100, {component: _as_decimal(assay.grade_pct) / 100})
    elif assay is None:
        resolved = BulkOreAssay(
            moisture_fraction=moisture_fraction,
            grades={component: grade} if grade is not None else {},
        )
    elif isinstance(assay, BulkOreAssay):
        resolved = assay
    else:
        raise TypeError("assay must be a BulkOreAssay or OreAssay.")

    target = parse_unit(target_unit)
    _, source_exp, source_basis = _ore_mass_factor(quantity.units)
    _, target_exp, target_basis = _ore_mass_factor(target)
    if source_exp != target_exp:
        raise ValueError("Source and target mass factors must have matching exponents.")
    scaled = quantity.to(target)
    source_basis = source_basis or target_basis
    target_basis = target_basis or source_basis
    magnitude = scaled.magnitude
    _ore_multiply(magnitude, Decimal(1))
    if source_basis != target_basis:
        dry_fraction = _as_decimal(resolved.dry_fraction) if "wet" in (source_basis, target_basis) else Decimal(1)
        fraction = None
        if "contained" in (source_basis, target_basis):
            if component not in resolved.grades:
                raise ValueError(f"Missing dry-basis grade for {component!r}.")
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
        elif source_basis == "contained":
            magnitude = apply_factor(magnitude, fraction, divide=source_exp == 1)
        if target_basis == "wet":
            magnitude = apply_factor(magnitude, dry_fraction, divide=source_exp == 1)
        elif target_basis == "contained":
            magnitude = apply_factor(magnitude, fraction, divide=source_exp == -1)
    return ureg.Quantity(magnitude, target)


# ---------------------------------------------------------------------------
# 8. Unified Cross-Commodity & Cross-Currency Quote Arbitrage Router
# ---------------------------------------------------------------------------


def _split_ratio_units(value: pint.Quantity | pint.Unit) -> tuple[pint.Unit, pint.Unit]:
    """Extracts the actual currency numerator and its physical denominator.

    Minor currencies retain their own unit. A plain currency amount has a
    dimensionless denominator. Compound currency factors are rejected clearly.
    """
    unit = value.units if isinstance(value, pint.Quantity) else value
    iso = _numerator_currency_iso(ureg.Quantity(1, unit))
    for name, exponent in unit._units.items():
        candidate = ureg.Unit(name)
        if exponent == 1 and candidate.dimensionality == {f"[currency_{iso}]": 1}:
            return candidate, candidate / unit
    raise ValueError("Quote must contain a standalone numerator currency unit.")


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
    try:
        return quantity.to(target)
    except pint.DimensionalityError as original_error:
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
    from .heat_content import HeatContentAssay

    contexts = []
    if product is not None and assay is None:
        raise ValueError("product selection requires an assay.")
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
) -> pint.Quantity:
    """Converts quote currency and physical denominator, then optionally settles.

    rates are major-currency fixings per base_currency, as in convert_currency().
    They are required only when ISO currencies differ. Both numerator currencies
    may use registered minor units. Cross-dimensional physical denominators
    require an assay context. Same-dimensional grade, moisture, and remelt
    adjustments still require their explicit conversion methods. product selects
    the density of meal/oil or another multi-product assay when supported.

    energy_basis selects LHV/HHV for biofuel, biomass, petcoke and hydrogen;
    gas_basis selects gas/NGL heating value. active_only=None preserves chemical
    and petrochemical defaults, while True/False overrides them explicitly.
    """
    if not isinstance(quote_quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quote_quantity).__name__}")
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
            raise ValueError(f"FX rates are required to convert {source_iso} to {target_iso}.")
        working = convert_currency(working, target_num, rates=dict(rates), base_currency=base_currency)

    contexts = _build_quote_contexts(
        assay=assay, crude_assay=crude_assay, gas_assay=gas_assay,
        gas_basis=gas_basis, energy_basis=energy_basis, active_only=active_only, product=product,
    )
    cross_dimension = source_den.dimensionality != target_den.dimensionality
    if cross_dimension and not contexts:
        raise ValueError("Changing quote denominator dimensionality requires an assay, crude_assay, or gas_assay.")
    if contexts:
        with ureg.context(*contexts):
            converted = _convert_price_quantity(working, target, contexts=contexts) if cross_dimension else working.to(target)
    else:
        converted = working.to(target)
    if quantize_settlement:
        return _quantize_currency_settlement(converted, target_num)
    return converted


# ---------------------------------------------------------------------------
# 9. Master Conversion Router
# ---------------------------------------------------------------------------


def _is_energy_unit(unit_like: str | pint.Unit) -> bool:
    """Identify energy independently of its market alias or quoting scale."""
    return parse_unit(unit_like).dimensionality == ureg.joule.dimensionality


def _is_thermal_endpoint(unit: pint.Unit) -> bool:
    """Recognise energy amounts, rates, calorific densities and quote legs."""
    if any(dimension.startswith("[currency_") for dimension in unit.dimensionality):
        _, unit = _split_ratio_units(unit)
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
    **heat_overrides: Any,
) -> pint.Quantity:
    """Resolve a composable thermal spec without replacing the commodity assay."""
    from .heat_content import HeatContentAssay, convert_heat_content, heat_content_from_assay

    if heat_content is not None and not isinstance(heat_content, HeatContentAssay):
        raise TypeError("heat_content must be a HeatContentAssay.")
    analytical_basis = heat_overrides.pop("analytical_basis", None)
    # An explicit spec is authoritative, including its density, molar mass and
    # reference state. Commodity adapters are used only when no spec is given.
    spec = heat_content_from_assay(heat_content if heat_content is not None else assay, **heat_overrides)
    return convert_heat_content(
        quantity, target_unit, assay=spec, basis=basis, analytical_basis=analytical_basis,
    )


def convert(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: Any = None,
    heat_content: HeatContentAssay | None = None,
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
    quantize_settlement: bool = False,
    active_only: bool | None = None,
    product: str | None = None,
    **heat_overrides: Any,
) -> pint.Quantity:
    """Dispatches pressure, FX, ore bases and physical assay conversions.

    rates accepts sparse BASE/QUOTE pairs (or tuple pairs), or currency fixings
    per base_currency. Currency-only targets preserve a quote's denominator.
    target_currency overrides the target quote's currency, not its physical leg.
    Settlement rounding happens after both FX and physical conversion.

    basis=None uses GCV for gas, HHV for NGL and LHV for energy assays.
    heat_content selects a shared HeatContentAssay for thermal conversions,
    including its density, alongside an existing commodity assay specification.
    It is engaged only for energy endpoints, rates, prices or calorific densities;
    other commodity conversions retain their original assay. Heat-content keyword
    overrides (gcv/ncv, density, molar_mass, analytical_basis, etc.) are validated
    by the shared thermal layer. Non-fuel assays require explicit calorific data.
    observed_temp corrects volume endpoints using the supplied assay's VCF:
    numeric temperatures are degrees F for crude, degrees C for other supported
    assays; temperature quantities carry their units. Same-dimension ore grade,
    moisture and contained-metal changes use the explicit ore labels.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    magnitude = _as_decimal(quantity.magnitude)
    if not magnitude.is_finite():
        raise ValueError("Conversion input must be finite.")
    physical_amount = quantity.dimensionality in (
        ureg.kilogram.dimensionality, (ureg.meter ** 3).dimensionality, ureg.mole.dimensionality,
    ) or (
        len(quantity.dimensionality) == 1
        and "gas_volume" in next(iter(quantity.dimensionality))
        and next(iter(quantity.dimensionality.values())) == 1
    )
    if physical_amount and magnitude < 0:
        raise ValueError("Physical mass, volume and molar amounts must be non-negative.")
    target = parse_unit(target_unit)
    pressure_dim = ureg.pascal.dimensionality
    if quantity.dimensionality == pressure_dim and target.dimensionality == pressure_dim:
        return convert_pressure(quantity, target, p_atm=p_atm, atmospheric_pressure=atmospheric_pressure)

    resolved_assay = assay
    crude_parameters = any(value is not None for value in (api_gravity, specific_gravity, bbl_per_mt))
    if resolved_assay is None and crude_parameters:
        resolved_assay = CrudeAssay(api_gravity=api_gravity, specific_gravity=specific_gravity, bbl_per_mt=bbl_per_mt)
    if resolved_assay is None and btu_per_scf is not None:
        resolved_assay = GasAssay(name="custom_gas", gcv_btu_per_scf=btu_per_scf)

    # Detect currencies structurally. Non-currency physical amounts skip FX.
    def currency_iso(unit: pint.Unit) -> str | None:
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
    if quantize_settlement and source_iso is not None and not isinstance(working.magnitude, Decimal):
        working = ureg.Quantity(_as_decimal(working.magnitude), working.units)
    if source_iso is not None and source_iso != target_iso:
        if rates is None:
            raise ValueError(f"FX rates are required to convert {source_iso} to {target_iso}.")
        is_fixing_map = all(isinstance(key, str) and "/" not in key for key in rates)
        if is_fixing_map:
            working = convert_currency(working, target_num, rates=dict(rates), base_currency=base_currency)
        else:
            working = convert_fx(working, target_num, rates, source_currency=source_currency)

    # Temperature correction changes volume basis before the standard assay leg.
    vcf = None
    volume_dim = (ureg.meter ** 3).dimensionality
    def is_volume_endpoint(dimensionality: Any) -> bool:
        remainder = dimensionality / volume_dim
        return not remainder or set(remainder) == {"[time]"}

    source_volume = is_volume_endpoint(source_den.dimensionality if source_den is not None else quantity.dimensionality)
    target_volume = is_volume_endpoint(target_den.dimensionality if target_den is not None else target.dimensionality)
    if observed_temp is not None:
        method = getattr(resolved_assay, "volume_correction_factor", None)
        if not callable(method):
            raise ValueError("observed_temp requires an assay with a volume-correction method.")
        vcf = _as_decimal(method(observed_temp))
        if not vcf.is_finite() or vcf <= 0:
            raise ValueError("Volume correction factor must be finite and positive.")
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
        raise ValueError("Heat-content overrides require an energy conversion.")
    bulk_route = isinstance(resolved_assay, (BulkOreAssay, OreAssay)) or moisture_fraction is not None or grade is not None
    if thermal_endpoint and thermal_requested:
        # Preserve existing fuel defaults; explicit shared thermal specs use LHV.
        thermal_basis = basis or (
            "gcv" if heat_content is None and isinstance(resolved_assay, GasAssay)
            else "hhv" if heat_content is None and isinstance(resolved_assay, NGLAssay)
            else "lhv"
        )
        converted = _convert_with_heat_content(
            working, target, assay=resolved_assay, heat_content=heat_content,
            basis=thermal_basis, **heat_overrides,
        )
    elif bulk_route:
        ore_assay = resolved_assay if isinstance(resolved_assay, (BulkOreAssay, OreAssay)) else None
        ore_component = resolved_assay.primary_element if isinstance(resolved_assay, OreAssay) and component == "Fe" else component
        converted = convert_bulk_ore(working, target, assay=ore_assay, moisture_fraction=moisture_fraction, grade=grade, component=ore_component)
    elif source_iso is not None:
        converted = convert_commodity_quote(
            working, target, assay=resolved_assay,
            gas_basis=basis or "gcv", energy_basis=basis or "lhv",
            active_only=active_only, product=product,
        )
    elif resolved_assay is not None:
        contexts = _build_quote_contexts(
            assay=resolved_assay, gas_basis=basis or "gcv",
            energy_basis=basis or "lhv", active_only=active_only, product=product,
        )
        with ureg.context(*contexts):
            converted = _convert_contextual_quantity(working, target, contexts)
    else:
        converted = working.to(target)

    if vcf is not None and target_volume:
        factor = vcf if source_iso is not None else 1 / vcf
        converted = ureg.Quantity(_ore_multiply(converted.magnitude, factor), converted.units)
    if quantize_settlement and target_num is not None:
        return _quantize_currency_settlement(converted, target_num)
    return converted


# ---------------------------------------------------------------------------
# 10. Module Public Exports
# ---------------------------------------------------------------------------

__all__: Final[tuple[str, ...]] = (
    # =========================================================================
    # Physical & Thermodynamic Reference Constants
    # =========================================================================
    "WATER_DENSITY_60F_KG_M3",
    "WATER_DENSITY_60F_LB_GAL",
    "BARREL_TO_M3",
    "MT_PER_BBL_WATER",
    "STANDARD_ATMOSPHERE_PA",
    "STANDARD_ATMOSPHERE_PSI",
    "STANDARD_ATMOSPHERE_BAR",

    # =========================================================================
    # Assay & Specification Models
    # =========================================================================
    "CrudeAssay",
    "DistillationCut",
    "RefineryCrudeAssay",
    "GasAssay",
    "CoalAssay",
    "CleanFuelAssay",
    "BiofuelAssay",
    "NGLAssay",
    "OreAssay",
    "BatteryMaterialAssay",
    "UraniumAssay",
    "BiomassAssay",
    "PetrochemAssay",
    "BullionAssay",
    "FertilizerAssay",
    "BitumenAssay",
    "PetcokeAssay",
    "ChemicalAssay",
    "HydrogenAssay",
    "FerroalloyAssay",
    "MetalAlloyAssay",
    "BulkOreAssay",

    # =========================================================================
    # Pint Transformation Context Builders
    # =========================================================================
    "create_hydrocarbon_context",
    "create_gas_heating_context",
    "create_coal_heating_context",
    "create_clean_fuel_context",
    "create_biofuel_context",
    "create_ngl_context",
    "create_ore_context",
    "create_battery_material_context",
    "create_uranium_context",
    "create_biomass_context",
    "create_petrochem_context",
    "create_bullion_context",
    "create_fertilizer_context",
    "create_bitumen_context",
    "create_petcoke_context",
    "create_chemical_context",
    "create_hydrogen_context",
    "create_ferroalloy_context",
    "create_metal_alloy_context",
    "create_fx_context",

    # =========================================================================
    # High-Level Domain Conversion Functions
    # =========================================================================
    "convert_hydrocarbon",
    "convert_gas_energy",
    "calculate_z_factor",
    "convert_actual_to_standard_gas_volume",
    "convert_standard_to_actual_gas_volume",
    "convert_pressure",
    "convert_pressure_datum",
    "atmospheric_pressure_at_elevation",
    "differential_pressure",
    "convert_currency",
    "convert_fx",
    "triangulate_fx_rate",
    "convert_wet_to_dry_mass",
    "convert_dry_to_wet_mass",
    "convert_ore_moisture_price",
    "extract_contained_metal",
    "extract_dry_mass_from_contained",
    "convert_bulk_ore",
    "convert_biofuel",
    "convert_ngl",
    "convert_ore",
    "convert_battery_material",
    "convert_uranium",
    "convert_biomass",
    "convert_petrochem",
    "convert_bullion",
    "convert_fertilizer",
    "convert_bitumen",
    "convert_petcoke",
    "convert_chemical",
    "convert_ferroalloy",
    "convert_metal_alloy",
    "convert_hydrogen",
    "convert_commodity",
    "convert_commodity_quote",
    "convert",
)
