"""Commodex assays metals alloys implementation and public interfaces."""
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
        from core.units.conversions import create_metal_alloy_context
        return create_metal_alloy_context(self, name=context_name)

    # --- Physical Volume <-> Mass <-> Contained Element <-> Remelt & Piece Conversions ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def solid_volume_to_mass(self, solid_volume: pint.Quantity) -> pint.Quantity:
        """Converts solid metal volume (m^3, cm^3, L) to gross alloy mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(solid_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(solid_volume).__name__}")
        m3_q = solid_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.solid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_solid_volume(self, alloy_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross alloy mass (mt, kg, lb) to solid metal volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(alloy_mass)
        mt_per_m3 = self.solid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def alloy_to_contained_element(self, alloy_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Converts gross alloy mass (mt, kg, lb) to contained `element` mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(alloy_mass)
        frac = self.element_fraction(element)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def contained_element_to_alloy(self, element_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Converts target contained `element` mass to required gross alloy mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(element_mass)
        frac = self.element_fraction(element)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def gross_to_recovered_melt_mass(self, feed_mass: pint.Quantity) -> pint.Quantity:
        """Converts charge/scrap mass (mt, kg) to recoverable remelted ingot/billet mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(feed_mass)
        scale = _scale_factor(mt_q.magnitude, self.remelt_yield_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def recovered_melt_to_gross_mass(self, target_melt_mass: pint.Quantity) -> pint.Quantity:
        """Converts target recovered melt mass (mt, kg) to required charge/scrap mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(target_melt_mass)
        scale = _scale_factor(mt_q.magnitude, self.remelt_yield_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def mass_to_piece_count(self, alloy_mass: pint.Quantity) -> float:
        """Converts gross alloy mass to nominal ingot/billet/piece count."""
        from core.units.registry import ureg
        if self.nominal_piece_mass_kg is None:
            raise ValueError(f"{self.name!r} does not define nominal_piece_mass_kg.")
        kg_val = float(self._validate_mass(alloy_mass).to(ureg.kilogram).magnitude)
        return (kg_val / self.nominal_piece_mass_kg)

    def piece_count_to_mass(self, piece_count: float) -> pint.Quantity:
        """Converts nominal piece/ingot/billet count to gross alloy mass (kilogram)."""
        from core.units.registry import ureg
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
        from core.units.conversions import _LB_PER_MT
        p = float(price_per_lb)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb must be finite, got {price_per_lb}")
        return (p * _LB_PER_MT)

    @staticmethod
    def mt_price_to_lb_price(price_per_mt: float) -> float:
        """Converts $/mt quote to $/lb equivalent."""
        from core.units.conversions import _LB_PER_MT
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "MetalAlloyAssay.cartridge_brass_c260")

    @classmethod
    def aluminium_6061(cls) -> MetalAlloyAssay:
        """Aluminium Alloy 6061 (structural/extrusion grade — ~97.9% Al, 1.0% Mg, 0.6% Si)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "MetalAlloyAssay.aluminium_6061")

    @classmethod
    def aluminium_7075(cls) -> MetalAlloyAssay:
        """Aluminium Alloy 7075 (high-strength aerospace grade — ~90.0% Al, 5.6% Zn, 2.5% Mg, 1.6% Cu)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "MetalAlloyAssay.aluminium_7075")

    @classmethod
    def monel_400(cls) -> MetalAlloyAssay:
        """Monel 400 Nickel-Copper Alloy (UNS N04400 — ~66.5% Ni, 31.5% Cu)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "MetalAlloyAssay.monel_400")

    @classmethod
    def inconel_625(cls) -> MetalAlloyAssay:
        """Inconel 625 Nickel-Chromium Superalloy (UNS N06625 — ~61.0% Ni, 21.5% Cr, 9.0% Mo)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "MetalAlloyAssay.inconel_625")

    @classmethod
    def titanium_6al_4v(cls) -> MetalAlloyAssay:
        """Titanium Ti-6Al-4V (Grade 5 aerospace/medical alloy — ~90.0% Ti, 6.0% Al, 4.0% V)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "MetalAlloyAssay.titanium_6al_4v")


__all__ = ('MetalAlloyAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
