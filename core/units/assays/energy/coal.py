"""Commodex assays energy coal implementation and public interfaces."""
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
        from core.units.analytical_basis import rebase_moisture_intensive
        from core.units.conversions import MJ_PER_BTU_IT
        from core.units.conversions import _ISO1928_H_KCAL_PER_PCT
        from core.units.conversions import _ISO1928_M_KCAL_PER_PCT
        from core.units.conversions import _KCAL_KG_TO_BTU_LB
        from core.units.conversions import _KCAL_KG_TO_GJ_MT
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
            resolved_gar = float(
                rebase_moisture_intensive(
                    resolved_gad, source_moisture_pct=im, target_moisture_pct=tm
                )
            )
            resolved_nar = resolved_gar - latent_deduction

        if not math.isfinite(resolved_nar) or not math.isfinite(resolved_gar) or resolved_nar <= 0 or resolved_gar <= 0:
            raise ValueError(
                f"Resolved calorific values must be finite and positive (NAR={resolved_nar}, GAR={resolved_gar})."
            )

        if resolved_gad is None and im is not None:
            resolved_gad = float(
                rebase_moisture_intensive(
                    resolved_gar, source_moisture_pct=tm, target_moisture_pct=im
                )
            )

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
        from core.units.analytical_basis import rebase_moisture_intensive
        return float(
            rebase_moisture_intensive(
                self.gar_kcal_per_kg, source_moisture_pct=self.total_moisture_pct, target_moisture_pct=0.0
            )
        )

    @property
    def nar_btu_per_lb(self) -> float:
        """Net Calorific Value (NAR) in Btu/lb."""
        from core.units.conversions import _KCAL_KG_TO_BTU_LB
        return (self.nar_kcal_per_kg * _KCAL_KG_TO_BTU_LB)

    @property
    def gj_per_mt_gar(self) -> float:
        """Gross Calorific Value (GAR) in GJ/mt."""
        from core.units.conversions import _KCAL_KG_TO_GJ_MT
        return (self.gar_kcal_per_kg * _KCAL_KG_TO_GJ_MT)

    @property
    def mwh_per_mt_gar(self) -> float:
        """Gross Calorific Value (GAR) in MWh/mt."""
        return (self.gj_per_mt_gar / 3.6)

    @property
    def mmbtu_per_mt_gar(self) -> float:
        """Gross Calorific Value (GAR) in MMBtu_IT/mt."""
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.gj_per_mt_gar / (MJ_PER_BTU_IT * 1000.0))

    @property
    def daf_kcal_per_kg(self) -> float | None:
        """Gross Calorific Value on a Dry Ash-Free (DAF) basis, or None if ash_pct is unset."""
        if self.ash_pct is None:
            return None
        return (self.gar_kcal_per_kg * (100.0 / (100.0 - self.total_moisture_pct - self.ash_pct)))

    def wmt_to_dmt(self, wet_metric_tonnes: float | Decimal | pint.Quantity) -> float | pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to Dry Metric Tonnes (DMT) using total moisture."""
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.registry import ureg
        if isinstance(wet_metric_tonnes, pint.Quantity):
            return wet_to_dry_mass(
                wet_metric_tonnes, moisture_pct=self.total_moisture_pct, target_unit=ureg.dmt
            )
        return float(
            wet_to_dry_mass(
                float(wet_metric_tonnes), moisture_pct=self.total_moisture_pct, allow_negative=True
            )
        )

    def dmt_to_wmt(self, dry_metric_tonnes: float | Decimal | pint.Quantity) -> float | pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT) using total moisture."""
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.registry import ureg
        if isinstance(dry_metric_tonnes, pint.Quantity):
            return dry_to_wet_mass(
                dry_metric_tonnes, moisture_pct=self.total_moisture_pct, target_unit=ureg.wmt
            )
        return float(
            dry_to_wet_mass(
                float(dry_metric_tonnes), moisture_pct=self.total_moisture_pct, allow_negative=True
            )
        )

    def mass_at_moisture(
        self,
        mass: pint.Quantity | float | Decimal,
        target_moisture_pct: float | Decimal,
        *,
        source_moisture_pct: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity | Decimal | float:
        """Converts coal mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        return convert_moisture_mass(
            mass,
            source_moisture_pct=self.total_moisture_pct if source_moisture_pct is None else source_moisture_pct,
            target_moisture_pct=target_moisture_pct,
            target_unit=target_unit,  # type: ignore[arg-type]
        )

    def as_context(self, *, basis: str = "nar", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [mass] and [energy] for this coal."""
        from core.units.conversions import create_coal_heating_context
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "CoalAssay.newcastle_6000")

    @classmethod
    def api2_cif_ara(cls) -> CoalAssay:
        """Argus/McCloskey API 2 CIF ARA 6,000 kcal/kg NAR European benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "CoalAssay.api2_cif_ara")

    @classmethod
    def api4_fob_rbct(cls) -> CoalAssay:
        """Argus/McCloskey API 4 FOB Richards Bay 6,000 kcal/kg NAR benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "CoalAssay.api4_fob_rbct")

    @classmethod
    def indonesian_4200(cls) -> CoalAssay:
        """ICI 4 Indonesian sub-bituminous 4,200 kcal/kg GAR thermal coal."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "CoalAssay.indonesian_4200")

    @classmethod
    def prb_8800(cls) -> CoalAssay:
        """US Powder River Basin 8,800 Btu/lb GAR sub-bituminous benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "CoalAssay.prb_8800")

    @classmethod
    def capp_12500(cls) -> CoalAssay:
        """US Central Appalachian 12,500 Btu/lb GAR bituminous benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "CoalAssay.capp_12500")

    @classmethod
    def plv_hard_coking(cls) -> CoalAssay:
        """Platts Premium Low-Vol (PLV) Hard Coking Coal FOB Australia benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "coal", "CoalAssay.plv_hard_coking")


__all__ = ('CoalAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
