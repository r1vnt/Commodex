"""Commodex assays metals ferroalloys implementation and public interfaces."""
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
        from core.units.conversions import create_ferroalloy_context
        return create_ferroalloy_context(self, name=context_name)

    # --- Physical Mass, Volume & Ladle Dosing Conversions ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(metal_mass)
        frac = self._element_fraction(element)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def bulk_volume_to_mass(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts bulk stowage volume (m^3) to gross ferroalloy mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.conversions import _LB_PER_MT
        p = float(price_per_mt_gross)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_gross must be finite, got {price_per_mt_gross}")
        return (p / (_LB_PER_MT * self.primary_mass_fraction))

    def contained_lb_price_to_gross_mt_price(self, price_per_lb_contained: float) -> float:
        """Converts $/lb contained primary metal quote to $/mt gross alloy."""
        from core.units.conversions import _LB_PER_MT
        p = float(price_per_lb_contained)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb_contained must be finite, got {price_per_lb_contained}")
        return (p * _LB_PER_MT * self.primary_mass_fraction)

    @staticmethod
    def contained_lb_price_to_contained_kg_price(price_per_lb: float) -> float:
        """Converts $/lb contained metal (US) to $/kg contained metal (Europe/Asia)."""
        from core.units.conversions import _LB_PER_KG
        p = float(price_per_lb)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb must be finite, got {price_per_lb}")
        return (p * _LB_PER_KG)

    @staticmethod
    def contained_kg_price_to_contained_lb_price(price_per_kg: float) -> float:
        """Converts $/kg contained metal (Europe/Asia) to $/lb contained metal (US)."""
        from core.units.conversions import _LB_PER_KG
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.hc_ferrochrome")

    @classmethod
    def lc_ferrochrome(cls) -> FerroalloyAssay:
        """Low-Carbon Ferrochrome (65% Cr, 0.10% max C — stainless steel refining)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.lc_ferrochrome")

    @classmethod
    def hc_ferromanganese(cls) -> FerroalloyAssay:
        """High-Carbon Ferromanganese (78% Mn basis, ~7.5% C — Europe/US benchmark)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.hc_ferromanganese")

    @classmethod
    def silicomanganese(cls) -> FerroalloyAssay:
        """Silicomanganese 65/16 (65% Mn basis, 16% Si, 2.0% C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.silicomanganese")

    @classmethod
    def ferrosilicon_75(cls) -> FerroalloyAssay:
        """Ferrosilicon 75 (75% Si basis — Europe/Asia/US steel deoxidiser benchmark)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.ferrosilicon_75")

    @classmethod
    def ferromolybdenum_60(cls) -> FerroalloyAssay:
        """Ferromolybdenum (60% Mo — quoted in $/lb Mo or $/kg Mo)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.ferromolybdenum_60")

    @classmethod
    def ferrovanadium_80(cls) -> FerroalloyAssay:
        """Ferrovanadium 80 (80% V — European $/kg V & US $/lb V benchmark)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.ferrovanadium_80")

    @classmethod
    def ferroniobium_65(cls) -> FerroalloyAssay:
        """Ferroniobium (65% Nb — HSLA structural steel micro-alloy)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "FerroalloyAssay.ferroniobium_65")


__all__ = ('FerroalloyAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
