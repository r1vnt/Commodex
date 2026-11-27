"""Commodex assays agriculture fibres implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, fields
from dataclasses import dataclass, replace
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final
from typing import Any, Final, Literal, Self
from typing import ClassVar
import math
import pint

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class SilkAndNaturalFibreAssay(_agricultural_families_AgriculturalFamilyAssay):
    primary_recovery_fraction: Decimal | None = None
    degumming_recovery_fraction: Decimal | None = None
    _commodities = frozenset({"silk", "jute", "hemp", "sisal", "flax"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        for name in ("primary_recovery_fraction", "degumming_recovery_fraction"):
            object.__setattr__(self, name, _agricultural_families_yield(getattr(self, name), name))

    def recovered_fibre_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.primary_recovery_fraction, "primary fibre recovery")

    def fibre_feedstock_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.primary_recovery_fraction, "primary fibre recovery")

    def degummed_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.degumming_recovery_fraction, "degumming recovery")


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 9. Cotton Lint Quality Specification (HVI Classing)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CottonQualitySpec(AgriculturalQualityBase):
    """Physical HVI quality specification for ginned cotton lint."""

    staple_length_mm: Decimal
    micronaire: Decimal
    strength_gpt: Decimal  # Bundle tenacity in g/tex
    colour_grade: str = "41"  # e.g. USDA HVI Strict Low Middling '41'
    leaf_grade: int | None = None
    trash_area_pct: Decimal = Decimal("0")
    moisture_regain_pct: Decimal = Decimal("8.5")  # Water mass as % of bone-dry fibre

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_clean_str
        from core.units.analytical_basis import _agricultural_quality_pct
        sl = _agricultural_quality_pct(self.staple_length_mm, "staple_length_mm", low=Decimal(15), high=Decimal(50))
        mic = _agricultural_quality_pct(self.micronaire, "micronaire", low=Decimal("1.5"), high=Decimal("8.0"))
        str_gpt = _agricultural_quality_pct(self.strength_gpt, "strength_gpt", low=Decimal(10), high=Decimal(55))
        col = _agricultural_quality_clean_str(self.colour_grade, "colour_grade").upper()
        if self.leaf_grade is not None:
            if isinstance(self.leaf_grade, bool) or not isinstance(self.leaf_grade, int) or not (1 <= self.leaf_grade <= 8):
                raise ValueError("leaf_grade must be an integer in [1, 8].")
        trash = _agricultural_quality_pct(self.trash_area_pct, "trash_area_pct", high=Decimal(25))
        regain = _agricultural_quality_pct(self.moisture_regain_pct, "moisture_regain_pct", high=Decimal(25))

        object.__setattr__(self, "staple_length_mm", sl)
        object.__setattr__(self, "micronaire", mic)
        object.__setattr__(self, "strength_gpt", str_gpt)
        object.__setattr__(self, "colour_grade", col)
        object.__setattr__(self, "trash_area_pct", trash)
        object.__setattr__(self, "moisture_regain_pct", regain)

    @classmethod
    def from_staple_32nds(
        cls,
        *,
        staple_32nds_inch: Decimal | float,
        micronaire: Decimal | float,
        strength_gpt: Decimal | float,
        **kwargs: Any,
    ) -> CottonQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        s32 = _agricultural_quality_nonneg(staple_32nds_inch, "staple_32nds_inch", strict_positive=True)
        mm = s32 * Decimal("25.4") / Decimal(32)
        return cls(staple_length_mm=mm, micronaire=micronaire, strength_gpt=strength_gpt, **kwargs)  # type: ignore[arg-type]

    @property
    def staple_32nds_inch(self) -> Decimal:
        return self.staple_length_mm * Decimal(32) / Decimal("25.4")

    @property
    def staple_inches(self) -> Decimal:
        return self.staple_length_mm / Decimal("25.4")

    @property
    def wet_basis_moisture_pct(self) -> Decimal:
        """Convert textile moisture regain (% of dry fibre) to wet-basis moisture %."""
        return (self.moisture_regain_pct / (Decimal(100) + self.moisture_regain_pct)) * Decimal(100)

    def conditioned_invoice_mass(
        self,
        measured_mass: pint.Quantity,
        *,
        standard_regain_pct: Decimal | float = Decimal("8.5"),
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Normalise cotton bale mass to standard textile moisture regain (ISO/ASTM 8.5%)."""
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        from core.units.analytical_basis import _agricultural_quality_pct
        kg = _agricultural_quality_mass_kg(measured_mass, "measured_mass")
        std_r = _agricultural_quality_pct(standard_regain_pct, "standard_regain_pct", high=Decimal(25))
        factor = (Decimal(100) + std_r) / (Decimal(100) + self.moisture_regain_pct)
        return _agricultural_quality_kg_out(kg * factor, target_unit or measured_mass.units, measured_mass)

    def blend_with(self, other: CottonQualitySpec, *, self_mass_fraction: Decimal | float) -> CottonQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, CottonQualitySpec):
            raise TypeError(f"Expected CottonQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        # Spinning laydown physics: micronaire blends harmonically by fibre mass
        b_mic = Decimal(1) / ((w1 / self.micronaire) + (w2 / other.micronaire))
        b_col = self.colour_grade if self.colour_grade == other.colour_grade else f"{self.colour_grade}/{other.colour_grade}"
        return CottonQualitySpec(
            staple_length_mm=w1 * self.staple_length_mm + w2 * other.staple_length_mm,
            micronaire=b_mic,
            strength_gpt=w1 * self.strength_gpt + w2 * other.strength_gpt,
            colour_grade=b_col,
            leaf_grade=self.leaf_grade if self.leaf_grade == other.leaf_grade else None,
            trash_area_pct=w1 * self.trash_area_pct + w2 * other.trash_area_pct,
            moisture_regain_pct=w1 * self.moisture_regain_pct + w2 * other.moisture_regain_pct,
        )

    def to_material_spec(self, name: str = "Cotton Lint") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        moist_f = self.wet_basis_moisture_pct / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=moist_f),
            composition=MaterialComposition(
                basis=CompositionBasis.AS_RECEIVED,
                moisture_fraction=moist_f,
                impurities={"trash": self.trash_area_pct / Decimal(100)},
            ),
        )


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 13. Wool Quality Specification (IWTO Core Test)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class WoolQualitySpec(AgriculturalQualityBase):
    """Physical IWTO core-test specification for greasy and scoured wool."""

    micron_diameter_um: Decimal
    vegetable_matter_pct: Decimal  # % on greasy basis
    clean_yield_pct: Decimal
    staple_length_mm: Decimal | None = None
    staple_strength_n_per_ktex: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        mic = _agricultural_quality_pct(self.micron_diameter_um, "micron_diameter_um", low=Decimal(10), high=Decimal(60))
        vm = _agricultural_quality_pct(self.vegetable_matter_pct, "vegetable_matter_pct", high=Decimal(40))
        clean = _agricultural_quality_pct(self.clean_yield_pct, "clean_yield_pct", allow_low=False)
        if clean + vm > Decimal(100):
            raise ValueError(f"clean_yield_pct ({clean}%) + vegetable_matter_pct ({vm}%) cannot exceed 100%.")
        sl = _agricultural_quality_nonneg(self.staple_length_mm, "staple_length_mm", strict_positive=True) if self.staple_length_mm is not None else None
        ss = _agricultural_quality_nonneg(self.staple_strength_n_per_ktex, "staple_strength_n_per_ktex", strict_positive=True) if self.staple_strength_n_per_ktex is not None else None

        object.__setattr__(self, "micron_diameter_um", mic)
        object.__setattr__(self, "vegetable_matter_pct", vm)
        object.__setattr__(self, "clean_yield_pct", clean)
        object.__setattr__(self, "staple_length_mm", sl)
        object.__setattr__(self, "staple_strength_n_per_ktex", ss)

    @property
    def grease_suint_moisture_loss_pct(self) -> Decimal:
        return Decimal(100) - self.clean_yield_pct - self.vegetable_matter_pct

    def greasy_to_clean_mass(self, greasy_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(greasy_mass, "greasy_mass")
        return _agricultural_quality_kg_out(kg * self.clean_yield_pct / Decimal(100), target_unit or greasy_mass.units, greasy_mass)

    def clean_to_greasy_mass(self, clean_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(clean_mass, "clean_mass")
        return _agricultural_quality_kg_out(kg / (self.clean_yield_pct / Decimal(100)), target_unit or clean_mass.units, clean_mass)

    def vegetable_matter_mass(self, greasy_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(greasy_mass, "greasy_mass")
        return _agricultural_quality_kg_out(kg * self.vegetable_matter_pct / Decimal(100), target_unit or greasy_mass.units, greasy_mass)

    def blend_with(self, other: WoolQualitySpec, *, self_mass_fraction: Decimal | float) -> WoolQualitySpec:
        """Blend greasy wool lots conserving clean wool mass, VM mass, and IWTO fibre fineness."""
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, WoolQualitySpec):
            raise TypeError(f"Expected WoolQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        c1 = w1 * self.clean_yield_pct
        c2 = w2 * other.clean_yield_pct
        total_c = c1 + c2
        cw1, cw2 = c1 / total_c, c2 / total_c

        # IWTO Core Test blending physics: d_blend = sqrt(1 / (cw1 / d1^2 + cw2 / d2^2))
        inv_d2 = (cw1 / (self.micron_diameter_um ** 2)) + (cw2 / (other.micron_diameter_um ** 2))
        b_mic = (Decimal(1) / inv_d2).sqrt()
        b_sl = (cw1 * self.staple_length_mm + cw2 * other.staple_length_mm) if (self.staple_length_mm is not None and other.staple_length_mm is not None) else None
        b_ss = (cw1 * self.staple_strength_n_per_ktex + cw2 * other.staple_strength_n_per_ktex) if (self.staple_strength_n_per_ktex is not None and other.staple_strength_n_per_ktex is not None) else None

        return WoolQualitySpec(
            micron_diameter_um=b_mic,
            vegetable_matter_pct=w1 * self.vegetable_matter_pct + w2 * other.vegetable_matter_pct,
            clean_yield_pct=total_c,
            staple_length_mm=b_sl,
            staple_strength_n_per_ktex=b_ss,
        )

    def to_material_spec(self, name: str = "Greasy Wool") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        from core.units.processing import ProcessingSpecification
        clean_f = self.clean_yield_pct / Decimal(100)
        vm_f = self.vegetable_matter_pct / Decimal(100)
        loss_f = Decimal(1) - clean_f
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(purity=clean_f),
            composition=MaterialComposition(
                components={"clean_wool": clean_f},
                basis=CompositionBasis.AS_RECEIVED,
                impurities={"vegetable_matter": vm_f},
            ),
            processing=ProcessingSpecification(
                input_material="greasy_wool",
                output_yields={"scoured_clean_wool": clean_f},
                process_losses={"scouring_vm_and_grease_loss": loss_f} if loss_f > Decimal(0) else {},
            ),
        )


from core.units.agricultural_specs import _STANDARD_GREASY_WOOL_BALE_KG
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class WoolAssay:
    """Immutable specification for greasy wool, scoured clean wool, and tops.

    Bridges greasy fleece mass <-> clean scoured mass via clean yield (%),
    scoured-to-top combing yield, bale counts, bulk stowage density, and clean
    auction vs. greasy farm-gate pricing parity.
    """

    name: str
    wool_type: str  # e.g. "MERINO_FINE", "MERINO_MEDIUM", "CROSSBRED", "SCOURED"
    clean_yield_pct: float = 65.0          # Greasy -> clean scoured wool yield (%)
    combing_top_yield_pct: float = 90.0    # Clean scoured -> combed wool top yield (%)
    micron: float | None = None            # Mean fibre diameter in micrometres (e.g. 18.5, 21.0, 28.0)
    bale_mass_kg: float = _STANDARD_GREASY_WOOL_BALE_KG
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.wool_type, str):
            raise ValueError("wool_type must be a non-empty string.")
        wtype = self.wool_type.strip().upper()
        if not wtype:
            raise ValueError("wool_type must be a non-empty string.")

        clean_y = _validate_finite_float(self.clean_yield_pct, "clean_yield_pct")
        if not (0.0 < clean_y <= 100.0):
            raise ValueError(f"clean_yield_pct must be in (0, 100], got {self.clean_yield_pct}")

        top_y = _validate_finite_float(self.combing_top_yield_pct, "combing_top_yield_pct")
        if not (0.0 < top_y <= 100.0):
            raise ValueError(f"combing_top_yield_pct must be in (0, 100], got {self.combing_top_yield_pct}")

        mic: float | None = None
        if self.micron is not None:
            mic = _validate_finite_float(self.micron, "micron")
            if not (10.0 <= mic <= 60.0):
                raise ValueError(f"micron must be in [10.0, 60.0], got {self.micron}")

        bale_kg = _validate_finite_float(self.bale_mass_kg, "bale_mass_kg")
        if bale_kg <= 0.0:
            raise ValueError(f"bale_mass_kg must be strictly positive, got {self.bale_mass_kg}")

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        object.__setattr__(self, "wool_type", wtype)
        object.__setattr__(self, "clean_yield_pct", clean_y)
        object.__setattr__(self, "combing_top_yield_pct", top_y)
        object.__setattr__(self, "micron", mic)
        object.__setattr__(self, "bale_mass_kg", bale_kg)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)

    def _effective_clean_fraction(self, clean_yield_pct: float | Decimal | None = None) -> Decimal:
        value = self.clean_yield_pct if clean_yield_pct is None else clean_yield_pct
        percentage = Decimal(str(value))
        if not percentage.is_finite() or not (0 < percentage <= 100):
            raise ValueError("clean_yield_pct must be finite and in (0, 100].")
        return percentage / 100

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.kilogram)
        value = Decimal(str(quantity.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Wool mass must be finite and non-negative.")
        return quantity

    @property
    def clean_fraction(self) -> float:
        return self.clean_yield_pct / 100.0

    @property
    def greasy_to_top_fraction(self) -> float:
        """Overall greasy-to-combed-top yield fraction."""
        return (self.clean_yield_pct / 100.0) * (self.combing_top_yield_pct / 100.0)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_wool_context
        return create_wool_context(self, name=context_name)

    # --- Greasy <-> Clean <-> Top Mass Conversions ---
    def greasy_to_clean_mass(
        self,
        greasy_mass: pint.Quantity,
        *,
        clean_yield_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts greasy fleece mass to clean scoured wool mass (kilogram)."""
        from core.units.registry import ureg
        if not isinstance(greasy_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(greasy_mass).__name__}")
        kg_q = self._validate_mass(greasy_mass)
        frac = self._effective_clean_fraction(clean_yield_pct)
        scale = frac if isinstance(kg_q.magnitude, Decimal) else float(frac)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def clean_to_greasy_mass(
        self,
        clean_mass: pint.Quantity,
        *,
        clean_yield_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts clean scoured wool mass back to required greasy fleece mass (kilogram)."""
        from core.units.registry import ureg
        if not isinstance(clean_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(clean_mass).__name__}")
        kg_q = self._validate_mass(clean_mass)
        frac = self._effective_clean_fraction(clean_yield_pct)
        scale = frac if isinstance(kg_q.magnitude, Decimal) else float(frac)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def clean_to_top_mass(self, clean_mass: pint.Quantity) -> pint.Quantity:
        """Converts clean scoured wool mass to combed wool top mass (kilogram)."""
        from core.units.registry import ureg
        if not isinstance(clean_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(clean_mass).__name__}")
        kg_q = self._validate_mass(clean_mass)
        fraction = Decimal(str(self.combing_top_yield_pct)) / 100
        scale = fraction if isinstance(kg_q.magnitude, Decimal) else float(fraction)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def greasy_to_top_mass(
        self,
        greasy_mass: pint.Quantity,
        *,
        clean_yield_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts greasy wool mass directly to combed wool top mass (kilogram)."""
        clean_q = self.greasy_to_clean_mass(greasy_mass, clean_yield_pct=clean_yield_pct)
        return self.clean_to_top_mass(clean_q)

    def top_to_clean_mass(self, top_mass: pint.Quantity) -> pint.Quantity:
        """Returns required clean scoured input from combed top output."""
        from core.units.registry import ureg
        quantity = self._validate_mass(top_mass)
        fraction = Decimal(str(self.combing_top_yield_pct)) / 100
        scale = fraction if isinstance(quantity.magnitude, Decimal) else float(fraction)
        return ureg.Quantity(quantity.magnitude / scale, ureg.kilogram)

    def top_to_greasy_mass(
        self, top_mass: pint.Quantity, *, clean_yield_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        return self.clean_to_greasy_mass(self.top_to_clean_mass(top_mass), clean_yield_pct=clean_yield_pct)

    # --- Bale & Bulk Volume Conversions ---
    def mass_to_bales(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal bale count, preserving Decimal mass input."""
        from core.units.conversions import _scale_factor
        kilograms = self._validate_mass(mass)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bale_mass_kg)

    def bales_to_mass(self, bale_count: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        count = Decimal(str(bale_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bale_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bale_mass_kg))
        magnitude = kilograms if isinstance(bale_count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        value = Decimal(str(m3_q.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Wool volume must be finite and non-negative.")
        scale = _scale_factor(m3_q.magnitude, self.bulk_density_kg_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.kilogram)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        kg_q = self._validate_mass(mass)
        scale = _scale_factor(kg_q.magnitude, self.bulk_density_kg_per_m3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.meter ** 3)

    # --- Pricing Parity Helpers ---
    @staticmethod
    def _validated_price(price: float | Decimal) -> Decimal:
        value = Decimal(str(price))
        if not value.is_finite():
            raise ValueError("Wool price must be finite.")
        return value

    @staticmethod
    def _price_result(value: Decimal, original: float | Decimal) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        return preserve_ieee754_scalar(value, original)

    def clean_price_to_greasy_price(
        self, clean_price_per_kg: float | Decimal, *, clean_yield_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        """Returns same-currency greasy/kg parity from clean/kg, before processing costs."""
        value = self._validated_price(clean_price_per_kg) * self._effective_clean_fraction(clean_yield_pct)
        return self._price_result(value, clean_price_per_kg)

    def greasy_price_to_clean_price(
        self, greasy_price_per_kg: float | Decimal, *, clean_yield_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        value = self._validated_price(greasy_price_per_kg) / self._effective_clean_fraction(clean_yield_pct)
        return self._price_result(value, greasy_price_per_kg)

    def clean_cents_per_kg_to_greasy_price(
        self, clean_cents_per_kg: float | Decimal, *, clean_yield_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        """Returns greasy currency-units/kg from clean cents/kg, without FX conversion."""
        value = self._validated_price(clean_cents_per_kg) / 100 * self._effective_clean_fraction(clean_yield_pct)
        return self._price_result(value, clean_cents_per_kg)

    def greasy_price_to_clean_cents_per_kg(
        self, greasy_price_per_kg: float | Decimal, *, clean_yield_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        value = self._validated_price(greasy_price_per_kg) / self._effective_clean_fraction(clean_yield_pct) * 100
        return self._price_result(value, greasy_price_per_kg)

    def greasy_price_per_kg_to_bale_price(self, greasy_price_per_kg: float | Decimal) -> float | Decimal:
        return self._price_result(self._validated_price(greasy_price_per_kg) * Decimal(str(self.bale_mass_kg)), greasy_price_per_kg)

    def bale_price_to_greasy_price_per_kg(self, bale_price: float | Decimal) -> float | Decimal:
        return self._price_result(self._validated_price(bale_price) / Decimal(str(self.bale_mass_kg)), bale_price)

    @staticmethod
    def cents_per_lb_to_price_per_kg(cents_per_lb: float | Decimal) -> float | Decimal:
        """Converts cents/lb to currency-units/kg in the same currency."""
        value = WoolAssay._validated_price(cents_per_lb) / 100 / Decimal("0.45359237")
        return WoolAssay._price_result(value, cents_per_lb)

    @staticmethod
    def price_per_kg_to_cents_per_lb(price_per_kg: float | Decimal) -> float | Decimal:
        value = WoolAssay._validated_price(price_per_kg) * Decimal("0.45359237") * 100
        return WoolAssay._price_result(value, price_per_kg)

    # --- Standard Wool Presets ---
    @classmethod
    def australian_merino_19micron(cls) -> WoolAssay:
        """Australian fine Merino fleece (~19.0 micron, 66% clean yield, 175 kg bale)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "WoolAssay.australian_merino_19micron")

    @classmethod
    def medium_merino_21micron(cls) -> WoolAssay:
        """AWEX 21-micron Merino benchmark (~65% clean yield, 175 kg bale)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "WoolAssay.medium_merino_21micron")

    @classmethod
    def crossbred_28micron(cls) -> WoolAssay:
        """Crossbred apparel/carpet wool (~28.0 micron, 72% clean yield, 180 kg bale)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "WoolAssay.crossbred_28micron")


from core.units.agricultural_specs import _US_COTTON_BALE_KG
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class CottonSeedLintAssay:
    """Immutable specification for cotton ginning turnout, lint bales, and cottonseed.

    Bridges raw harvested seed cotton with ginned cotton lint mass, whole
    cottonseed mass, ginning trash loss, international bale weights (480 lb US,
    170 kg Indian, 227 kg AU/BR), bulk seed density, and ginning netback parity.
    """

    name: str
    origin: str  # e.g. "US_UPLAND", "INDIA_SHANKAR6", "BRAZIL", "AUSTRALIA"
    lint_turnout_pct: float = 40.0
    cottonseed_turnout_pct: float = 54.0
    bale_mass_kg: float = _US_COTTON_BALE_KG
    seed_bulk_density_kg_per_m3: float = 420.0

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.origin, str):
            raise ValueError("origin must be a non-empty string.")
        orig = self.origin.strip().upper()
        if not orig:
            raise ValueError("origin must be a non-empty string.")

        lint_y = _validate_finite_float(self.lint_turnout_pct, "lint_turnout_pct")
        if not (0.0 < lint_y <= 100.0):
            raise ValueError(f"lint_turnout_pct must be in (0, 100], got {self.lint_turnout_pct}")

        seed_y = _validate_finite_float(self.cottonseed_turnout_pct, "cottonseed_turnout_pct")
        if not (0.0 <= seed_y <= 100.0) or (lint_y + seed_y > 100.0 + 1e-9):
            raise ValueError(
                f"Invalid ginning turnout: lint_turnout_pct={lint_y}, cottonseed_turnout_pct={seed_y}."
            )

        bale_kg = _validate_finite_float(self.bale_mass_kg, "bale_mass_kg")
        if bale_kg <= 0.0:
            raise ValueError(f"bale_mass_kg must be strictly positive, got {self.bale_mass_kg}")

        seed_rho = _validate_finite_float(self.seed_bulk_density_kg_per_m3, "seed_bulk_density_kg_per_m3")
        if not (200.0 <= seed_rho <= 1000.0):
            raise ValueError(
                f"seed_bulk_density_kg_per_m3 must be in [200, 1000], got {self.seed_bulk_density_kg_per_m3}"
            )

        object.__setattr__(self, "origin", orig)
        object.__setattr__(self, "lint_turnout_pct", lint_y)
        object.__setattr__(self, "cottonseed_turnout_pct", seed_y)
        object.__setattr__(self, "bale_mass_kg", bale_kg)
        object.__setattr__(self, "seed_bulk_density_kg_per_m3", seed_rho)

    @property
    def gin_trash_loss_pct(self) -> float:
        """Unrecovered trash, motes, and moisture loss percentage during ginning."""
        from core.units.precision import to_ieee754_float
        return to_ieee754_float(max(Decimal(0), Decimal(100) - Decimal(str(self.lint_turnout_pct)) - Decimal(str(self.cottonseed_turnout_pct))))

    @property
    def seed_to_lint_ratio(self) -> float:
        """Mass ratio of recovered cottonseed per unit of ginned cotton lint."""
        from core.units.precision import to_ieee754_float
        return to_ieee754_float(Decimal(str(self.cottonseed_turnout_pct)) / Decimal(str(self.lint_turnout_pct)))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_cotton_seed_lint_context
        return create_cotton_seed_lint_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        number = Decimal(str(quantity.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Cotton/lint/seed mass must be finite and non-negative.")
        return quantity

    @staticmethod
    def _yield_scale(magnitude: Any, percentage: float) -> float | Decimal:
        factor = Decimal(str(percentage)) / 100
        return factor if isinstance(magnitude, Decimal) else float(factor)

    # --- Ginning Turnout Conversions (Seed Cotton <-> Lint & Cottonseed) ---
    def seed_cotton_to_lint_mass(self, seed_cotton_mass: pint.Quantity) -> pint.Quantity:
        """Converts raw harvested seed cotton mass to ginned cotton lint mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(seed_cotton_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_cotton_mass).__name__}")
        mt_q = self._validate_mass(seed_cotton_mass)
        scale = self._yield_scale(mt_q.magnitude, self.lint_turnout_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def lint_to_seed_cotton_mass(self, lint_mass: pint.Quantity) -> pint.Quantity:
        """Converts ginned cotton lint mass back to required raw seed cotton mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(lint_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(lint_mass).__name__}")
        mt_q = self._validate_mass(lint_mass)
        scale = self._yield_scale(mt_q.magnitude, self.lint_turnout_pct)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def seed_cotton_to_cottonseed_mass(self, seed_cotton_mass: pint.Quantity) -> pint.Quantity:
        """Converts raw seed cotton mass to recovered whole cottonseed mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(seed_cotton_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_cotton_mass).__name__}")
        mt_q = self._validate_mass(seed_cotton_mass)
        scale = self._yield_scale(mt_q.magnitude, self.cottonseed_turnout_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def cottonseed_to_seed_cotton_mass(self, cottonseed_mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if self.cottonseed_turnout_pct == 0:
            raise ValueError("Cannot infer seed-cotton requirements from zero cottonseed turnout.")
        quantity = self._validate_mass(cottonseed_mass)
        scale = self._yield_scale(quantity.magnitude, self.cottonseed_turnout_pct)
        return ureg.Quantity(quantity.magnitude / scale, ureg.metric_ton)

    def gin_seed_cotton(self, seed_cotton_mass: pint.Quantity) -> dict[str, pint.Quantity]:
        """Returns lint, cottonseed and remaining trash without intermediate rounding."""
        from core.units.registry import ureg
        quantity = self._validate_mass(seed_cotton_mass)
        lint = self.seed_cotton_to_lint_mass(quantity)
        seed = self.seed_cotton_to_cottonseed_mass(quantity)
        residual = quantity.magnitude - lint.magnitude - seed.magnitude
        if residual < 0:
            residual = Decimal(0) if isinstance(quantity.magnitude, Decimal) else 0.0
        return {"lint": lint, "cottonseed": seed, "trash": ureg.Quantity(residual, ureg.metric_ton)}

    # --- Lint Bales & Bulk Whole Cottonseed ---
    def lint_mass_to_bales(self, lint_mass: pint.Quantity) -> float | Decimal:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kilograms = self._validate_mass(lint_mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bale_mass_kg)

    def bales_to_lint_mass(self, bale_count: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        count = Decimal(str(bale_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bale_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bale_mass_kg))
        value = kilograms if isinstance(bale_count, Decimal) else float(kilograms)
        return ureg.Quantity(value, ureg.kilogram).to(ureg.metric_ton)

    def bales_to_seed_cotton_mass(self, bale_count: float | Decimal) -> pint.Quantity:
        """Returns required raw seed cotton to supply a nominal lint-bale order."""
        return self.lint_to_seed_cotton_mass(self.bales_to_lint_mass(bale_count))

    def cottonseed_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts bulk whole cottonseed volume (m^3, ft^3) to mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        number = Decimal(str(m3_q.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Cottonseed volume must be finite and non-negative.")
        mt_per_m3 = Decimal(str(self.seed_bulk_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(m3_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def cottonseed_mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts bulk whole cottonseed mass (mt, short_ton) to volume (m^3)."""
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = Decimal(str(self.seed_bulk_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    # --- Pricing & Ginning Netback ---
    @staticmethod
    def _validated_price(price: float | Decimal) -> Decimal:
        value = Decimal(str(price))
        if not value.is_finite():
            raise ValueError("Cotton/seed price must be finite.")
        return value

    @staticmethod
    def _price_result(value: Decimal, *originals: Any) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        return preserve_ieee754_scalar(value, *originals)

    @staticmethod
    def lint_cents_per_lb_to_usd_per_mt(cents_per_lb: float | Decimal) -> float | Decimal:
        value = CottonSeedLintAssay._validated_price(cents_per_lb) * 10 / Decimal("0.45359237")
        return CottonSeedLintAssay._price_result(value, cents_per_lb)

    @staticmethod
    def lint_usd_per_mt_to_cents_per_lb(usd_per_mt: float | Decimal) -> float | Decimal:
        value = CottonSeedLintAssay._validated_price(usd_per_mt) * Decimal("0.45359237") / 10
        return CottonSeedLintAssay._price_result(value, usd_per_mt)

    def lint_cents_per_lb_to_usd_per_bale(self, cents_per_lb: float | Decimal) -> float | Decimal:
        pounds = Decimal(str(self.bale_mass_kg)) / Decimal("0.45359237")
        value = self._validated_price(cents_per_lb) / 100 * pounds
        return self._price_result(value, cents_per_lb)

    def lint_usd_per_bale_to_cents_per_lb(self, usd_per_bale: float | Decimal) -> float | Decimal:
        pounds = Decimal(str(self.bale_mass_kg)) / Decimal("0.45359237")
        value = self._validated_price(usd_per_bale) / pounds * 100
        return self._price_result(value, usd_per_bale)

    def gross_ginning_value_usd_per_mt_seed_cotton(
        self, lint_usd_per_mt: float | Decimal, cottonseed_usd_per_mt: float | Decimal,
    ) -> float | Decimal:
        """Gross lint/seed revenue before ginning, transport and other costs."""
        value = (Decimal(str(self.lint_turnout_pct)) / 100 * self._validated_price(lint_usd_per_mt)
                 + Decimal(str(self.cottonseed_turnout_pct)) / 100 * self._validated_price(cottonseed_usd_per_mt))
        return self._price_result(value, lint_usd_per_mt, cottonseed_usd_per_mt)

    # --- Standard Cotton Origin Presets ---
    @classmethod
    def us_upland_480lb(cls) -> CottonSeedLintAssay:
        """US Upland cotton (40% lint / 54% seed turnout, 480 lb bale)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "CottonSeedLintAssay.us_upland_480lb")

    @classmethod
    def indian_shankar6_170kg(cls) -> CottonSeedLintAssay:
        """Indian Shankar-6 cotton (~35% lint / 61% seed turnout, 170 kg bale)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "CottonSeedLintAssay.indian_shankar6_170kg")

    @classmethod
    def brazil_australian_227kg(cls) -> CottonSeedLintAssay:
        """Brazilian / Australian high-turnout cotton (42% lint / 52% seed, 227 kg bale)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "CottonSeedLintAssay.brazil_australian_227kg")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (SilkAndNaturalFibreAssay,))
__all__ = ('CottonQualitySpec', 'CottonSeedLintAssay', 'SilkAndNaturalFibreAssay', 'WoolAssay', 'WoolQualitySpec', 'convert_silk_and_natural_fibre', 'create_silk_and_natural_fibre_context')
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
    from core.units.registry import Unit
    from core.units.assays.agriculture.fibres import convert_silk_and_natural_fibre
    from core.units.assays.agriculture.fibres import create_silk_and_natural_fibre_context
