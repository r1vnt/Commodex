"""Commodex assays agriculture rubber implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
from dataclasses import dataclass, replace
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final
from typing import Any, Final, Literal, Self
import math
import pint

from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 10. Natural Rubber Quality Specification (ISO 2000 TSR)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class RubberQualitySpec(AgriculturalQualityBase):
    """Physical ISO 2000 quality specification for Technically Specified / sheet rubber."""

    dirt_pct: Decimal
    ash_pct: Decimal
    volatile_matter_pct: Decimal
    mooney_viscosity_ml1_4_100c: Decimal
    nitrogen_pct: Decimal | None = None
    plasticity_retention_index: Decimal | None = None
    dry_rubber_content_pct: Decimal = Decimal("100.0")

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_pct
        dirt = _agricultural_quality_pct(self.dirt_pct, "dirt_pct", high=Decimal(10))
        ash = _agricultural_quality_pct(self.ash_pct, "ash_pct", high=Decimal(10))
        vm = _agricultural_quality_pct(self.volatile_matter_pct, "volatile_matter_pct", high=Decimal(20))
        mooney = _agricultural_quality_pct(self.mooney_viscosity_ml1_4_100c, "mooney_viscosity_ml1_4_100c", low=Decimal(10), high=Decimal(150))
        nit = _agricultural_quality_pct(self.nitrogen_pct, "nitrogen_pct", high=Decimal(5)) if self.nitrogen_pct is not None else None
        pri = _agricultural_quality_pct(self.plasticity_retention_index, "plasticity_retention_index", high=Decimal(150)) if self.plasticity_retention_index is not None else None
        drc = _agricultural_quality_pct(self.dry_rubber_content_pct, "dry_rubber_content_pct", allow_low=False)

        object.__setattr__(self, "dirt_pct", dirt)
        object.__setattr__(self, "ash_pct", ash)
        object.__setattr__(self, "volatile_matter_pct", vm)
        object.__setattr__(self, "mooney_viscosity_ml1_4_100c", mooney)
        object.__setattr__(self, "nitrogen_pct", nit)
        object.__setattr__(self, "plasticity_retention_index", pri)
        object.__setattr__(self, "dry_rubber_content_pct", drc)

    def evaluate_iso2000_grade(self, grade: Literal["TSR5", "TSR10", "TSR20"] = "TSR20") -> QualityAuditResult:
        """Screen impurity levels against physical ISO 2000 limits for TSR5, TSR10, or TSR20."""
        from core.units.analytical_basis import QualityThreshold
        from core.units.analytical_basis import _agricultural_quality_TSR_ISO2000_LIMITS
        from core.units.analytical_basis import _agricultural_quality_clean_str
        key = _agricultural_quality_clean_str(grade, "grade").upper()
        if key not in _agricultural_quality_TSR_ISO2000_LIMITS:
            raise ValueError(f"Unsupported ISO 2000 grade {grade!r}; expected 'TSR5', 'TSR10', or 'TSR20'.")
        limits = _agricultural_quality_TSR_ISO2000_LIMITS[key]
        rules: dict[str, QualityThreshold] = {
            "dirt_pct": QualityThreshold(max_value=limits["dirt_pct"]),
            "ash_pct": QualityThreshold(max_value=limits["ash_pct"]),
            "volatile_matter_pct": QualityThreshold(max_value=limits["volatile_matter_pct"]),
        }
        if self.nitrogen_pct is not None:
            rules["nitrogen_pct"] = QualityThreshold(max_value=limits["nitrogen_pct"])
        return self.evaluate_thresholds(rules)

    def net_pure_elastomer_mass(self, gross_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        """Compute dry, dirt-free, ash-free, volatile-free elastomer mass."""
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(gross_mass, "gross_mass")
        dry_kg = kg * (self.dry_rubber_content_pct / Decimal(100))
        non_elastomer_frac = (self.dirt_pct + self.ash_pct + self.volatile_matter_pct) / Decimal(100)
        return _agricultural_quality_kg_out(dry_kg * max(Decimal(0), Decimal(1) - non_elastomer_frac), target_unit or gross_mass.units, gross_mass)

    def blend_with(self, other: RubberQualitySpec, *, self_mass_fraction: Decimal | float) -> RubberQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, RubberQualitySpec):
            raise TypeError(f"Expected RubberQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        # Polymer melt/bale log-viscosity mixing rule for Mooney viscosity
        log_m = float(w1) * math.log(float(self.mooney_viscosity_ml1_4_100c)) + float(w2) * math.log(float(other.mooney_viscosity_ml1_4_100c))
        b_mooney = Decimal(str(math.exp(log_m)))
        b_nit = (w1 * self.nitrogen_pct + w2 * other.nitrogen_pct) if (self.nitrogen_pct is not None and other.nitrogen_pct is not None) else None
        b_pri = (w1 * self.plasticity_retention_index + w2 * other.plasticity_retention_index) if (self.plasticity_retention_index is not None and other.plasticity_retention_index is not None) else None
        return RubberQualitySpec(
            dirt_pct=w1 * self.dirt_pct + w2 * other.dirt_pct,
            ash_pct=w1 * self.ash_pct + w2 * other.ash_pct,
            volatile_matter_pct=w1 * self.volatile_matter_pct + w2 * other.volatile_matter_pct,
            mooney_viscosity_ml1_4_100c=b_mooney,
            nitrogen_pct=b_nit,
            plasticity_retention_index=b_pri,
            dry_rubber_content_pct=w1 * self.dry_rubber_content_pct + w2 * other.dry_rubber_content_pct,
        )

    def to_material_spec(self, name: str = "Natural Rubber") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        vm_f = self.volatile_matter_pct / Decimal(100)
        ash_f = self.ash_pct / Decimal(100)
        dirt_f = self.dirt_pct / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                moisture=vm_f,
                purity=self.dry_rubber_content_pct / Decimal(100),
            ),
            composition=MaterialComposition(
                basis=CompositionBasis.AS_RECEIVED,
                moisture_fraction=vm_f,
                ash_fraction=ash_f,
                impurities={"dirt": dirt_f},
            ),
        )


from core.units.agricultural_specs import _STANDARD_TSR_BALE_KG
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class RubberAssay:
    """Immutable specification for natural rubber, latex, cup lump, and TSR bales.

    Bridges wet coagulum/latex mass <-> dry rubber mass via Dry Rubber Content
    (DRC %), bale counts <-> mass, optional bulk/tank density, and USc/kg <-> $/mt
    pricing parity.
    """

    name: str
    grade: str  # e.g. "TSR20", "RSS3", "LATEX", "CUPLUMP"
    dry_rubber_content_pct: float = 100.0
    bale_mass_kg: float = _STANDARD_TSR_BALE_KG
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.grade, str):
            raise ValueError("grade must be a non-empty string.")
        grade_clean = self.grade.strip().upper()
        if not grade_clean:
            raise ValueError("grade must be a non-empty string.")

        drc = _validate_finite_float(self.dry_rubber_content_pct, "dry_rubber_content_pct")
        if not (0.0 < drc <= 100.0):
            raise ValueError(f"dry_rubber_content_pct must be in (0, 100], got {self.dry_rubber_content_pct}")

        bale_kg = _validate_finite_float(self.bale_mass_kg, "bale_mass_kg")
        if bale_kg <= 0.0:
            raise ValueError(f"bale_mass_kg must be strictly positive, got {self.bale_mass_kg}")

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        object.__setattr__(self, "grade", grade_clean)
        object.__setattr__(self, "dry_rubber_content_pct", drc)
        object.__setattr__(self, "bale_mass_kg", bale_kg)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)

    def _effective_drc_fraction(self, drc_pct: float | Decimal | None = None) -> Decimal:
        value = self.dry_rubber_content_pct if drc_pct is None else drc_pct
        percentage = Decimal(str(value))
        if not percentage.is_finite() or not (0 < percentage <= 100):
            raise ValueError(f"drc_pct must be finite and in (0, 100], got {value}")
        return percentage / 100

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        magnitude = Decimal(str(quantity.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Rubber mass must be finite and non-negative.")
        return quantity

    @property
    def drc_fraction(self) -> float:
        """Dry rubber content as a decimal fraction in (0, 1]."""
        return self.dry_rubber_content_pct / 100.0

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_rubber_context
        return create_rubber_context(self, name=context_name)

    # --- Wet / DRC Mass Conversions ---
    def wet_to_dry_mass(
        self,
        wet_mass: pint.Quantity,
        *,
        drc_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts wet latex, cup lump, or coagulum mass to dry rubber mass (metric_ton)."""
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(wet_mass)
        drc_frac = self._effective_drc_fraction(drc_pct)
        return wet_to_dry_mass(mt_q, Decimal(1) - drc_frac, target_unit=ureg.metric_ton)

    def dry_to_wet_mass(
        self,
        dry_mass: pint.Quantity,
        *,
        drc_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts dry rubber mass to required wet latex/cup lump mass (metric_ton)."""
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(dry_mass)
        drc_frac = self._effective_drc_fraction(drc_pct)
        return dry_to_wet_mass(mt_q, Decimal(1) - drc_frac, target_unit=ureg.metric_ton)

    # --- Bale & Bulk Volume Conversions ---
    def mass_to_bales(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal bale count, preserving Decimal mass input."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        scale = _scale_factor(kilograms.magnitude, self.bale_mass_kg)
        return kilograms.magnitude / scale

    def bales_to_mass(self, bale_count: float | Decimal) -> pint.Quantity:
        """Converts nominal bale count to metric tonnes; fractional counts are allowed."""
        from core.units.registry import ureg
        count = Decimal(str(bale_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bale_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bale_mass_kg))
        magnitude = kilograms if isinstance(bale_count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram).to(ureg.metric_ton)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.agricultural_specs import _KG_PER_MT
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        magnitude = Decimal(str(m3_q.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Rubber bulk volume must be finite and non-negative.")
        mt_per_m3 = self.bulk_density_kg_per_m3 / _KG_PER_MT
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.agricultural_specs import _KG_PER_MT
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.bulk_density_kg_per_m3 / _KG_PER_MT
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    # --- Pricing Parity Helpers ---
    @staticmethod
    def _validated_price(price: float | Decimal) -> Decimal:
        value = Decimal(str(price))
        if not value.is_finite():
            raise ValueError("Rubber price must be finite.")
        return value

    @staticmethod
    def cents_per_kg_to_usd_per_mt(cents_per_kg: float | Decimal) -> float | Decimal:
        """Converts US cents/kg to USD/mt at IEEE 754 precision."""
        from core.units.precision import preserve_ieee754_scalar
        value = RubberAssay._validated_price(cents_per_kg) * 10
        return preserve_ieee754_scalar(value, cents_per_kg)

    @staticmethod
    def usd_per_mt_to_cents_per_kg(usd_per_mt: float | Decimal) -> float | Decimal:
        """Converts USD/mt to US cents/kg at IEEE 754 precision."""
        from core.units.precision import preserve_ieee754_scalar
        value = RubberAssay._validated_price(usd_per_mt) / 10
        return preserve_ieee754_scalar(value, usd_per_mt)

    def wet_price_to_dry_price_per_mt(
        self,
        wet_price_per_mt: float | Decimal,
        *,
        drc_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        """Returns dry-rubber USD/mt parity from a wet USD/mt quote at IEEE 754 precision."""
        from core.units.precision import preserve_ieee754_scalar
        value = self._validated_price(wet_price_per_mt) / self._effective_drc_fraction(drc_pct)
        return preserve_ieee754_scalar(value, wet_price_per_mt)

    # --- Standard Natural Rubber Presets ---
    @classmethod
    def tsr20_standard(cls) -> RubberAssay:
        """Technically Specified Rubber 20 (33.33 kg bale, 100% dry rubber basis)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "RubberAssay.tsr20_standard")

    @classmethod
    def rss3_bale_standard(cls, *, bale_mass_kg: float = 111.11) -> RubberAssay:
        """Ribbed Smoked Sheet 3 (RSS3) standard export bale."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "RubberAssay.rss3_bale_standard", overrides={'bale_mass_kg': bale_mass_kg})

    @classmethod
    def field_latex_60drc(cls) -> RubberAssay:
        """Concentrated field/centrifuged natural rubber latex (~60% DRC)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "RubberAssay.field_latex_60drc")

    @classmethod
    def cup_lump_55drc(cls) -> RubberAssay:
        """Farm-gate cup lump coagulum (~55% DRC)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "RubberAssay.cup_lump_55drc")


__all__ = ('RubberAssay', 'RubberQualitySpec')
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
    from core.units.assays.base import MaterialSpec
    from core.units.analytical_basis import QualityAuditResult
    from core.units.registry import Unit
