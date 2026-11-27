"""Commodex assays energy refining implementation and public interfaces."""
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


from core.units.validation import validated_assay


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
        from core.units.assays.energy.petroleum import CrudeAssay
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
        from core.units.conversions import _viscosity_blending_number
        from core.units.conversions import _viscosity_from_vbn
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


__all__ = ('DistillationCut', 'RefineryCrudeAssay')
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
    from core.units.assays.energy.petroleum import CrudeAssay
