"""Commodex assays energy nuclear implementation and public interfaces."""
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

from core.units.conversions import _NATURAL_U235_PCT
from core.units.validation import validated_assay


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
        from core.units.conversions import _dirac_value_function
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
        from core.units.conversions import create_uranium_context
        return create_uranium_context(self, name=context_name)

    # --- Front-End Stoichiometric Conversions (U3O8 <-> kgU <-> UF6 <-> UO2) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
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
        from core.units.conversions import _U_FRAC_IN_U3O8
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _U_FRAC_IN_U3O8
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kg_q = self._validate_mass(natural_u_mass)
        loss_mult = (1.0 - self.conversion_loss_pct / 100.0) if include_conversion_loss else 1.0
        factor = _U_FRAC_IN_U3O8 * loss_mult
        scale = _scale_factor(kg_q.magnitude, factor)
        u3o8_kg = ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)
        return u3o8_kg.to(ureg.pound)

    def kgu_to_uf6_mass(self, u_mass: pint.Quantity) -> pint.Quantity:
        """Converts elemental uranium mass (kgU) to total UF6 compound mass (kg UF6)."""
        from core.units.conversions import _U_FRAC_IN_UF6
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kg_q = self._validate_mass(u_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UF6)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def uf6_mass_to_kgu(self, uf6_mass: pint.Quantity) -> pint.Quantity:
        """Converts total UF6 compound mass (kg UF6) to elemental uranium mass (kgU)."""
        from core.units.conversions import _U_FRAC_IN_UF6
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kg_q = self._validate_mass(uf6_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UF6)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def kgu_to_uo2_mass(self, u_mass: pint.Quantity) -> pint.Quantity:
        """Converts elemental uranium mass (kgU) to fabricated UO2 pellet mass (kg UO2)."""
        from core.units.conversions import _U_FRAC_IN_UO2
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kg_q = self._validate_mass(u_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UO2)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def uo2_mass_to_kgu(self, uo2_mass: pint.Quantity) -> pint.Quantity:
        """Converts fabricated UO2 pellet mass (kg UO2) to elemental uranium mass (kgU)."""
        from core.units.conversions import _U_FRAC_IN_UO2
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kg_q = self._validate_mass(uo2_mass)
        scale = _scale_factor(kg_q.magnitude, _U_FRAC_IN_UO2)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    # --- Enrichment Conversions (EUP <-> Feed UF6 <-> U3O8 & SWU) ---
    def eup_to_feed_kgu(self, eup_u_mass: pint.Quantity) -> pint.Quantity:
        """Converts Enriched Uranium Product mass (kgU) to required natural feed (kgU as UF6)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kg_q = self._validate_mass(eup_u_mass)
        scale = _scale_factor(kg_q.magnitude, self.feed_to_product_ratio)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def feed_kgu_to_eup(self, feed_u_mass: pint.Quantity) -> pint.Quantity:
        """Converts natural feed (kgU as UF6) to producible Enriched Uranium Product (kgU)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.registry import ureg
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
        from core.units.conversions import _dirac_value_function
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
            from core.units.conversions import _dirac_value_function
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "UraniumAssay.natural_uranium")

    @classmethod
    def leu_4_5pct(cls, *, tails_u235_pct: float = 0.25) -> UraniumAssay:
        """Standard Light-Water Reactor (LWR) Low-Enriched Uranium (4.5% U-235)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "UraniumAssay.leu_4_5pct", overrides={'tails_u235_pct': tails_u235_pct})

    @classmethod
    def haleu_19_75pct(cls, *, tails_u235_pct: float = 0.20) -> UraniumAssay:
        """High-Assay Low-Enriched Uranium (HALEU, 19.75% U-235) for SMRs."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "UraniumAssay.haleu_19_75pct", overrides={'tails_u235_pct': tails_u235_pct})

    @classmethod
    def depleted_tails(cls, *, tails_u235_pct: float = 0.25) -> UraniumAssay:
        """Depleted Uranium (DU) enrichment tails."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "UraniumAssay.depleted_tails", overrides={'name': f'Depleted Uranium ({tails_u235_pct:.2f}% U-235)', 'product_u235_pct': tails_u235_pct, 'tails_u235_pct': tails_u235_pct})


__all__ = ('UraniumAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
