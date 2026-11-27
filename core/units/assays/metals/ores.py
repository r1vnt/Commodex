"""Commodex assays metals ores implementation and public interfaces."""
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
        from core.units.analytical_basis import dry_matter_fraction
        return float(dry_matter_fraction(moisture_pct=self.moisture_pct))

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
        from core.units.conversions import create_ore_context
        return create_ore_context(self, basis=basis, name=context_name)

    # --- Physical Mass Conversions (WMT <-> DMT <-> Contained / Payable Metal) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def wmt_to_dmt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to Dry Metric Tonnes (DMT)."""
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(wet_mass)
        return wet_to_dry_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.dmt, allow_negative=True)

    def dmt_to_wmt(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT)."""
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(dry_mass)
        return dry_to_wet_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.wmt, allow_negative=True)

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture_pct: float | Decimal,
        *,
        source_moisture_pct: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity:
        """Converts ore mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=self.moisture_pct if source_moisture_pct is None else source_moisture_pct,
            target_moisture_pct=target_moisture_pct,
            target_unit=target_unit or ureg.metric_ton,
        )

    def dmt_to_contained_metal(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) of ore/concentrate to contained metal mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(dry_mass)
        grade_frac = self.grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, grade_frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.contained_metal_tonne)

    def contained_metal_to_dmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained metal mass (e.g. mt, lb) to required Dry Metric Tonnes (DMT)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(metal_mass)
        grade_frac = self.grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, grade_frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.dmt)

    def wmt_to_contained_metal(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to contained metal mass (metric_ton)."""
        return self.dmt_to_contained_metal(self.wmt_to_dmt(wet_mass))

    def contained_metal_to_wmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained metal mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.contained_metal_to_dmt(metal_mass))

    def dmt_to_payable_metal(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to smelter payable metal mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(dry_mass)
        payable_frac = self.payable_grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, payable_frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.payable_metal_tonne)

    def payable_metal_to_dmt(self, payable_mass: pint.Quantity) -> pint.Quantity:
        """Converts target smelter payable metal mass to required Dry Metric Tonnes (DMT)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.payable_grade_pct <= 0.0:
            raise ValueError(f"{self.name!r} has zero payable grade.")
        mt_q = self._validate_mass(payable_mass)
        payable_frac = self.payable_grade_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, payable_frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.dmt)

    def wmt_to_payable_metal(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to smelter payable metal mass (metric_ton)."""
        return self.dmt_to_payable_metal(self.wmt_to_dmt(wet_mass))

    def payable_metal_to_wmt(self, payable_mass: pint.Quantity) -> pint.Quantity:
        """Converts target payable metal mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.payable_metal_to_dmt(payable_mass))

    def wmt_to_volume(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to bulk vessel stowage volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.exceptions import MissingAssayError
        from core.units.registry import ureg
        if self.bulk_density_wmt_per_m3 is None:
            raise MissingAssayError(f"{self.name!r} does not define bulk_density_wmt_per_m3.")
        mt_q = self._validate_mass(wet_mass)
        scale = _scale_factor(mt_q.magnitude, self.bulk_density_wmt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def volume_to_wmt(self, bulk_volume: pint.Quantity) -> pint.Quantity:
        """Converts bulk vessel stowage volume (e.g. m^3) to Wet Metric Tonnes (WMT)."""
        from core.units.conversions import _scale_factor
        from core.units.exceptions import MissingAssayError
        from core.units.registry import ureg
        if not isinstance(bulk_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bulk_volume).__name__}")
        if self.bulk_density_wmt_per_m3 is None:
            raise MissingAssayError(f"{self.name!r} does not define bulk_density_wmt_per_m3.")
        m3_q = bulk_volume.to(ureg.meter ** 3)
        scale = _scale_factor(m3_q.magnitude, self.bulk_density_wmt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def dmt_to_impurity_mass(self, dry_mass: pint.Quantity, element: str) -> pint.Quantity:
        """Computes contained mass (metric_ton) of a specific impurity element from DMT."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
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
        from core.units.analytical_basis import rebase_moisture_intensive
        p = float(price_per_dmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt must be finite, got {price_per_dmt}")
        return float(rebase_moisture_intensive(p, source_moisture_pct=0.0, target_moisture_pct=self.moisture_pct))

    def wmt_price_to_dmt_price(self, price_per_wmt: float) -> float:
        """Converts a $/wmt quote to $/dmt (dry basis)."""
        from core.units.analytical_basis import rebase_moisture_intensive
        p = float(price_per_wmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_wmt must be finite, got {price_per_wmt}")
        return float(rebase_moisture_intensive(p, source_moisture_pct=self.moisture_pct, target_moisture_pct=0.0))

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
        from core.units.exceptions import IncompatibleCommodityError
        if not isinstance(other, OreAssay):
            raise TypeError(f"Expected OreAssay, got {type(other).__name__}")
        if self.primary_element != other.primary_element:
            raise IncompatibleCommodityError(
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
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "OreAssay.iron_ore_62fe")

    @classmethod
    def iron_ore_65fe(cls) -> OreAssay:
        """Fastmarkets 65% Fe Brazil fines (Carajás) benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "OreAssay.iron_ore_65fe")

    @classmethod
    def iron_ore_58fe(cls) -> OreAssay:
        """58% Fe low-grade iron ore fines benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "OreAssay.iron_ore_58fe")

    @classmethod
    def copper_concentrate(cls) -> OreAssay:
        """Standard clean Copper Concentrate (28% Cu, 96.5% payable / 1.0% min deduction)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "OreAssay.copper_concentrate")

    @classmethod
    def zinc_concentrate(cls) -> OreAssay:
        """Standard Zinc Concentrate (52% Zn, 85% payable / 8.0% min deduction)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "OreAssay.zinc_concentrate")

    @classmethod
    def bauxite_guinea(cls) -> OreAssay:
        """Guinea CBG metallurgical-grade Bauxite (45% Al2O3)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "OreAssay.bauxite_guinea")

    @classmethod
    def nickel_concentrate(cls) -> OreAssay:
        """Standard Nickel Sulphide Concentrate (12% Ni, 75% smelter payable)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "OreAssay.nickel_concentrate")


from core.units.validation import validated_assay


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
        from core.units.conversions import _as_decimal
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
        from core.units.analytical_basis import dry_matter_fraction
        from core.units.exceptions import MissingAssayError
        if self.moisture_fraction is None:
            raise MissingAssayError("Explicit ore moisture is required for wet/dry conversion.")
        dec_dry = dry_matter_fraction(self.moisture_fraction)
        return dec_dry if isinstance(self.moisture_fraction, Decimal) else float(dec_dry)

    def wmt_to_dmt(self, wet_mass: float | Decimal) -> float | Decimal:
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.exceptions import MissingAssayError
        if self.moisture_fraction is None:
            raise MissingAssayError("Explicit ore moisture is required for wet/dry conversion.")
        return wet_to_dry_mass(wet_mass, self.moisture_fraction, allow_negative=True)

    def dmt_to_wmt(self, dry_mass: float | Decimal) -> float | Decimal:
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.exceptions import MissingAssayError
        if self.moisture_fraction is None:
            raise MissingAssayError("Explicit ore moisture is required for wet/dry conversion.")
        return dry_to_wet_mass(dry_mass, self.moisture_fraction, allow_negative=True)

    def mass_at_moisture(
        self,
        mass: pint.Quantity | float | Decimal,
        target_moisture_fraction: float | Decimal,
        *,
        source_moisture_fraction: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity | float | Decimal:
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.exceptions import MissingAssayError
        src_m = self.moisture_fraction if source_moisture_fraction is None else source_moisture_fraction
        if src_m is None:
            raise MissingAssayError("Explicit ore moisture is required for moisture conversion.")
        return convert_moisture_mass(
            mass,
            source_moisture=src_m,
            target_moisture=target_moisture_fraction,
            target_unit=target_unit,  # type: ignore[arg-type]
        )

    def contained_metal_dmt(self, dry_mass: float | Decimal, component: str) -> float | Decimal:
        from core.units.conversions import _as_decimal
        from core.units.conversions import _ore_multiply
        return _ore_multiply(dry_mass, _as_decimal(self.grades[component]))

    def contained_dmtu(self, dry_mass_dmt: float | Decimal, component: str = "Fe") -> float | Decimal:
        return self.contained_metal_dmt(dry_mass_dmt, component) * 100

    def to_ore_assay(self, component: str = "Fe") -> OreAssay:
        """Builds an OreAssay for one component, subject to its stricter limits.

        OreAssay requires a positive grade and moisture below 50%; independent
        component assays are not automatically treated as additive impurities.
        """
        from core.units.conversions import _as_decimal
        return OreAssay(
            name=self.name, primary_element=component,
            grade_pct=float(_as_decimal(self.grades[component]) * 100),
            moisture_pct=float((1 - _as_decimal(self.dry_fraction)) * 100),
        )


__all__ = ('BulkOreAssay', 'OreAssay')
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
