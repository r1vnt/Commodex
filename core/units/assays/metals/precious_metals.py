"""Commodex assays metals precious_metals implementation and public interfaces."""
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
        from core.units.conversions import create_bullion_context
        return create_bullion_context(self, basis=basis, name=context_name)

    # --- Physical Mass Conversions (Gross <-> Fine <-> Payable <-> Bars) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.troy_ounce)

    def gross_to_fine(self, gross_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross bullion/doré mass (e.g. kg, g, ozt) to fine troy ounces."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        ozt_q = self._validate_mass(gross_mass)
        scale = _scale_factor(ozt_q.magnitude, self.fineness)
        return ureg.Quantity(ozt_q.magnitude * scale, ureg.troy_ounce)

    def fine_to_gross(self, fine_mass: pint.Quantity) -> pint.Quantity:
        """Converts fine metal mass (e.g. ozt, kg, g) to required gross troy ounces."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        ozt_q = self._validate_mass(fine_mass)
        scale = _scale_factor(ozt_q.magnitude, self.fineness)
        return ureg.Quantity(ozt_q.magnitude / scale, ureg.troy_ounce)

    def gross_to_payable_fine(self, gross_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross bullion/doré mass to refinery payable fine troy ounces."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        ozt_q = self._validate_mass(gross_mass)
        factor = self.fineness * self.payable_factor
        scale = _scale_factor(ozt_q.magnitude, factor)
        return ureg.Quantity(ozt_q.magnitude * scale, ureg.troy_ounce)

    def payable_fine_to_gross(self, payable_mass: pint.Quantity) -> pint.Quantity:
        """Converts target payable fine metal mass to required gross troy ounces."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        ozt_q = self._validate_mass(payable_mass)
        factor = self.fineness * self.payable_factor
        scale = _scale_factor(ozt_q.magnitude, factor)
        return ureg.Quantity(ozt_q.magnitude / scale, ureg.troy_ounce)

    def gross_to_byproduct_fine(self, gross_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Computes contained fine troy ounces of a doré byproduct metal (e.g. 'Ag')."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.registry import ureg
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
        from core.units.conversions import _GRAMS_PER_TROY_OUNCE
        p = float(price_per_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_ozt must be finite, got {price_per_ozt}")
        return (p / _GRAMS_PER_TROY_OUNCE)

    @staticmethod
    def gram_price_to_ozt_price(price_per_gram: float) -> float:
        """Converts $/g quote to $/ozt equivalent."""
        from core.units.conversions import _GRAMS_PER_TROY_OUNCE
        p = float(price_per_gram)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gram must be finite, got {price_per_gram}")
        return (p * _GRAMS_PER_TROY_OUNCE)

    @staticmethod
    def ozt_price_to_kg_price(price_per_ozt: float) -> float:
        """Converts $/ozt quote to $/kg kilobar equivalent."""
        from core.units.conversions import _TROY_OUNCES_PER_KG
        p = float(price_per_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_ozt must be finite, got {price_per_ozt}")
        return (p * _TROY_OUNCES_PER_KG)

    @staticmethod
    def kg_price_to_ozt_price(price_per_kg: float) -> float:
        """Converts $/kg quote to $/ozt equivalent."""
        from core.units.conversions import _TROY_OUNCES_PER_KG
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
        from core.units.conversions import _GRAMS_PER_TOLA
        from core.units.conversions import _GRAMS_PER_TROY_OUNCE
        p = float(price_per_ozt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_ozt must be finite, got {price_per_ozt}")
        return (p * (_GRAMS_PER_TOLA / _GRAMS_PER_TROY_OUNCE))

    @staticmethod
    def tola_price_to_ozt_price(price_per_tola: float) -> float:
        """Converts India/Dubai $/tola quote to $/ozt equivalent."""
        from core.units.conversions import _GRAMS_PER_TOLA
        from core.units.conversions import _GRAMS_PER_TROY_OUNCE
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BullionAssay.lbma_gold_good_delivery")

    @classmethod
    def comex_gold_100oz(cls) -> BullionAssay:
        """COMEX Good Delivery Gold Bar (995.0 min fineness, 100 ozt nominal)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BullionAssay.comex_gold_100oz")

    @classmethod
    def gold_9999_kilobar(cls) -> BullionAssay:
        """Four-Nines Gold Kilobar (999.9 fineness, 1 kg / 32.1507 ozt — SGE / Dubai / Swiss)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BullionAssay.gold_9999_kilobar")

    @classmethod
    def lbma_silver_good_delivery(cls) -> BullionAssay:
        """LBMA Good Delivery Silver Bar (999.0 min fineness, 1,000 ozt nominal)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BullionAssay.lbma_silver_good_delivery")

    @classmethod
    def lppm_platinum_good_delivery(cls) -> BullionAssay:
        """LPPM Good Delivery Platinum Plate/Ingot (999.5 min fineness)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BullionAssay.lppm_platinum_good_delivery")

    @classmethod
    def lppm_palladium_good_delivery(cls) -> BullionAssay:
        """LPPM Good Delivery Palladium Plate/Ingot (999.5 min fineness)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BullionAssay.lppm_palladium_good_delivery")

    @classmethod
    def gold_dore(cls) -> BullionAssay:
        """Typical Mine-Site Gold Doré Bar (85% Au, 10% Ag, 99.9% refinery payable)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BullionAssay.gold_dore")


__all__ = ('BullionAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
