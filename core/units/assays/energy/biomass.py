"""Commodex assays energy biomass implementation and public interfaces."""
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
        from core.units.analytical_basis import dry_matter_fraction
        from core.units.analytical_basis import rebase_moisture_intensive
        from core.units.conversions import _ISO18125_H_GJ_PER_MT_PCT
        from core.units.conversions import _ISO18125_M_GJ_PER_MT_PCT
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

        dry_frac = float(dry_matter_fraction(moisture_pct=m))
        # ISO 18125 gross-to-net latent heat deduction (GJ/mt as-received)
        latent_deduction = (_ISO18125_H_GJ_PER_MT_PCT * h_db * dry_frac) + (_ISO18125_M_GJ_PER_MT_PCT * m)

        if lhv_gj_per_mt_ar is not None:
            res_lhv_ar = float(lhv_gj_per_mt_ar)
            res_hhv_ar = res_lhv_ar + latent_deduction
            res_hhv_db = float(rebase_moisture_intensive(res_hhv_ar, source_moisture_pct=m, target_moisture_pct=0.0))
        elif hhv_gj_per_mt_ar is not None:
            res_hhv_ar = float(hhv_gj_per_mt_ar)
            res_lhv_ar = res_hhv_ar - latent_deduction
            res_hhv_db = float(rebase_moisture_intensive(res_hhv_ar, source_moisture_pct=m, target_moisture_pct=0.0))
        else:
            assert hhv_gj_per_mt_db is not None
            res_hhv_db = float(hhv_gj_per_mt_db)
            res_hhv_ar = float(rebase_moisture_intensive(res_hhv_db, source_moisture_pct=0.0, target_moisture_pct=m))
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
        from core.units.analytical_basis import dry_matter_fraction
        return float(dry_matter_fraction(moisture_pct=self.moisture_pct))

    @property
    def lhv_mwh_per_mt_ar(self) -> float:
        return (self.lhv_gj_per_mt_ar / 3.6)

    @property
    def hhv_mwh_per_mt_ar(self) -> float:
        return (self.hhv_gj_per_mt_ar / 3.6)

    @property
    def lhv_mmbtu_per_mt_ar(self) -> float:
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.lhv_gj_per_mt_ar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mmbtu_per_mt_ar(self) -> float:
        from core.units.conversions import MJ_PER_BTU_IT
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
        from core.units.conversions import _ISO18125_H_GJ_PER_MT_PCT
        return (self.hhv_gj_per_mt_db - (_ISO18125_H_GJ_PER_MT_PCT * self.hydrogen_db_pct))

    @property
    def lhv_kcal_per_kg_ar(self) -> float:
        """Net As Received calorific value in kcal_IT/kg (Asian PKS / Japanese FIT basis)."""
        from core.units.conversions import _KCAL_KG_TO_GJ_MT
        return (self.lhv_gj_per_mt_ar / _KCAL_KG_TO_GJ_MT)

    @property
    def hhv_kcal_per_kg_ar(self) -> float:
        """Gross As Received calorific value in kcal_IT/kg."""
        from core.units.conversions import _KCAL_KG_TO_GJ_MT
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
        from core.units.conversions import create_biomass_context
        return create_biomass_context(self, basis=basis, name=context_name)

    # --- Physical Conversions (WMT <-> ODT/DMT <-> Bulk m^3 <-> Energy) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def wmt_to_odt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to Oven-Dry Tonnes (ODT / DMT)."""
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(wet_mass)
        return wet_to_dry_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.metric_ton, allow_negative=True)

    def odt_to_wmt(self, oven_dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Oven-Dry Tonnes (ODT / DMT) to Wet Metric Tonnes (WMT)."""
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(oven_dry_mass)
        return dry_to_wet_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.metric_ton, allow_negative=True)

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture_pct: float | Decimal,
        *,
        source_moisture_pct: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity:
        """Converts biomass mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=self.moisture_pct if source_moisture_pct is None else source_moisture_pct,
            target_moisture_pct=target_moisture_pct,
            target_unit=target_unit or ureg.metric_ton,
        )

    def wmt_to_bulk_volume(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to bulk vessel/silo stowage volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(wet_mass)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def bulk_volume_to_wmt(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts bulk stowage volume (e.g. m^3) to Wet Metric Tonnes (WMT)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(bulk_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bulk_volume).__name__}")
        m3_q = bulk_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_energy(self, wet_mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts as-received mass (WMT) to total energy (GJ)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(wet_mass)
        gj_per_mt = self._gj_per_mt(basis)
        scale = _scale_factor(mt_q.magnitude, gj_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.GJ)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. GJ, MWh, MMBtu) to required as-received mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.analytical_basis import rebase_moisture_intensive
        p = float(price_per_wmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_wmt must be finite, got {price_per_wmt}")
        return float(rebase_moisture_intensive(p, source_moisture_pct=self.moisture_pct, target_moisture_pct=0.0))

    def odt_price_to_wmt_price(self, price_per_odt: float) -> float:
        """Converts a $/ODT (Oven-Dry / Bone-Dry Tonne) quote to $/WMT."""
        from core.units.analytical_basis import rebase_moisture_intensive
        p = float(price_per_odt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_odt must be finite, got {price_per_odt}")
        return float(rebase_moisture_intensive(p, source_moisture_pct=0.0, target_moisture_pct=self.moisture_pct))

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
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiomassAssay.industrial_pellets_i2")

    @classmethod
    def enplus_a1_pellets(cls) -> BiomassAssay:
        """ENplus A1 Premium Wood Pellets (17.2 GJ/mt NAR, 7% moisture, <=0.7% ash DB)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiomassAssay.enplus_a1_pellets")

    @classmethod
    def torrefied_pellets(cls) -> BiomassAssay:
        """Torrefied / Black Biomass Pellets (21.0 GJ/mt NAR, 4% moisture, coal-replacement grade)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiomassAssay.torrefied_pellets")

    @classmethod
    def forestry_wood_chips(cls) -> BiomassAssay:
        """Seasoned Forestry Wood Chips (35% moisture, 19.2 GJ/mt HHV dry basis)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiomassAssay.forestry_wood_chips")

    @classmethod
    def palm_kernel_shells(cls) -> BiomassAssay:
        """Southeast Asian Palm Kernel Shells (PKS — ~3,800 kcal/kg NAR, Japanese FIT biomass)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiomassAssay.palm_kernel_shells")


__all__ = ('BiomassAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
