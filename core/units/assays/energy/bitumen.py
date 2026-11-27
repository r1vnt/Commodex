"""Commodex assays energy bitumen implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final
import inspect
import math
import pint
import warnings

from core.units.conversions import _ASTM_D4311_GROUP0_ALPHA_PER_C
from core.units.validation import validated_assay


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
        from core.units.conversions import _LITRES_PER_BBL
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_BBL / 1000.0)))

    @property
    def gal_per_mt(self) -> float:
        """Standard US gallons per metric tonne at reference_temp_c."""
        from core.units.conversions import _LITRES_PER_US_GAL
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging standard bitumen [volume] and [mass]."""
        from core.units.conversions import create_bitumen_context
        return create_bitumen_context(self, name=context_name)

    # --- ASTM D4311 Thermal Volume Correction & Custody Transfer ---
    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float | None = None,
    ) -> float:
        """Computes ASTM D4311 thermal Volume Correction Factor (VCF) from hot observed temp to reference."""
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(std_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(std_volume).__name__}")
        m3_q = std_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts bitumen mass (mt, short_ton) to standard-temperature volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        std_m3 = self.mass_to_volume(mass)
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        scale = _scale_factor(std_m3.magnitude, vcf)
        return ureg.Quantity(std_m3.magnitude / scale, ureg.meter ** 3)

    # --- Commercial Pricing Parity ($/st <-> $/mt, $/gal <-> $/mt) ---
    @staticmethod
    def st_price_to_mt_price(price_per_st: float) -> float:
        """Converts US $/short_ton asphalt quote to international $/metric_ton."""
        from core.units.conversions import _MT_PER_SHORT_TON
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / _MT_PER_SHORT_TON)

    @staticmethod
    def mt_price_to_st_price(price_per_mt: float) -> float:
        """Converts international $/metric_ton quote to US $/short_ton."""
        from core.units.conversions import _MT_PER_SHORT_TON
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BitumenAssay.pen_60_70")

    @classmethod
    def pen_80_100(cls) -> BitumenAssay:
        """Softer Paving Bitumen 80/100 (cooler-climate road grade)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BitumenAssay.pen_80_100")

    @classmethod
    def pg_64_22(cls) -> BitumenAssay:
        """US Superpave Performance Grade PG 64-22 (primary North American paving binder)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BitumenAssay.pg_64_22")

    @classmethod
    def oxidised_roofing(cls) -> BitumenAssay:
        """Oxidised / Blown Industrial Roofing Bitumen (95/25 grade)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BitumenAssay.oxidised_roofing")


__all__ = ('BitumenAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
