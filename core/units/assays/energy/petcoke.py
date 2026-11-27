"""Commodex assays energy petcoke implementation and public interfaces."""
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

from core.units.validation import validated_assay


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
        from core.units.analytical_basis import dry_matter_fraction
        return float(dry_matter_fraction(moisture_pct=self.moisture_pct))

    @property
    def fixed_carbon_db_pct(self) -> float:
        """Proximate-analysis Fixed Carbon on a dry basis: 100 - VM_db - Ash_db."""
        return (max(0.0, 100.0 - self.volatile_matter_db_pct - self.ash_db_pct))

    @property
    def lhv_gj_per_mt_ar(self) -> float:
        """As-received Net Calorific Value (LHV / NAR) in GJ/mt."""
        from core.units.conversions import _KCAL_KG_TO_GJ_MT
        return (self.lhv_kcal_per_kg_ar * _KCAL_KG_TO_GJ_MT)

    @property
    def hhv_gj_per_mt_ar(self) -> float:
        """As-received Gross Calorific Value (HHV / GAR) in GJ/mt."""
        from core.units.conversions import _KCAL_KG_TO_GJ_MT
        return (self.hhv_kcal_per_kg_ar * _KCAL_KG_TO_GJ_MT)

    @property
    def lhv_mmbtu_per_mt_ar(self) -> float:
        """As-received LHV in MMBtu_IT per metric tonne."""
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.lhv_gj_per_mt_ar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mmbtu_per_mt_ar(self) -> float:
        """As-received HHV in MMBtu_IT per metric tonne."""
        from core.units.conversions import MJ_PER_BTU_IT
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
        from core.units.conversions import create_petcoke_context
        return create_petcoke_context(self, basis=basis, name=context_name)

    # --- Physical Mass, Stowage Volume & Energy Conversions ---
    def wmt_to_dmt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts as-received wet petcoke mass (WMT) to dry metric tonnes (DMT)."""
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.registry import ureg
        if not isinstance(wet_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
        mt_q = wet_mass.to(ureg.metric_ton)
        return wet_to_dry_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.dmt, allow_negative=True)

    def dmt_to_wmt(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts dry petcoke mass (DMT) to as-received wet metric tonnes (WMT)."""
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.registry import ureg
        if not isinstance(dry_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(dry_mass).__name__}")
        mt_q = dry_mass.to(ureg.metric_ton)
        return dry_to_wet_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.wmt, allow_negative=True)

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture_pct: float | Decimal,
        *,
        source_moisture_pct: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity:
        """Converts petcoke mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=self.moisture_pct if source_moisture_pct is None else source_moisture_pct,
            target_moisture_pct=target_moisture_pct,
            target_unit=target_unit or ureg.metric_ton,
        )

    def wmt_to_bulk_volume(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts wet petcoke mass to vessel hold / stockpile bulk volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(wet_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
        mt_q = wet_mass.to(ureg.metric_ton)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def bulk_volume_to_wmt(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts vessel hold / stockpile bulk volume (m^3) to wet petcoke mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(bulk_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bulk_volume).__name__}")
        m3_q = bulk_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_energy(self, wet_mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts as-received petcoke mass to thermal energy (GJ)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(wet_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
        mt_q = wet_mass.to(ureg.metric_ton)
        gj_per_mt = self.lhv_gj_per_mt_ar if basis.strip().lower() in {"lhv", "nar"} else self.hhv_gj_per_mt_ar
        scale = _scale_factor(mt_q.magnitude, gj_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.gigajoule)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts target thermal energy (GJ, MWh, MMBtu) to required petcoke mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "PetcokeAssay.usgc_fuel_grade_65s")

    @classmethod
    def usgc_fuel_grade_45s(cls) -> PetcokeAssay:
        """US Gulf Coast Mid-Sulphur Fuel-Grade Sponge Coke (4.5% S, 7,550 kcal/kg NAR, HGI 50)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "PetcokeAssay.usgc_fuel_grade_45s")

    @classmethod
    def calcined_anode_grade(cls) -> PetcokeAssay:
        """Calcined Petroleum Coke (CPC) Anode Grade for aluminium smelting (>99% FC, 0.2% moisture)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "PetcokeAssay.calcined_anode_grade")

    @classmethod
    def shot_coke(cls) -> PetcokeAssay:
        """Hard USGC Fuel-Grade Shot Coke (6.8% S, low HGI 35)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "PetcokeAssay.shot_coke")


__all__ = ('PetcokeAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
