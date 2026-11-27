"""Commodex assays agriculture softs implementation and public interfaces."""
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
# 6. Coffee Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CoffeeQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for green Arabica / Robusta coffee."""

    screen_size: Decimal  # Sieve size in 1/64ths of an inch (e.g. 15, 16, 17, 18)
    moisture_pct: Decimal
    primary_defects_count: Decimal = Decimal("0")
    secondary_defects_count: Decimal = Decimal("0")
    processing_loss_pct: Decimal = Decimal("18.0")  # Roasting / cleaning shrink %
    water_activity: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        scr = _agricultural_quality_pct(self.screen_size, "screen_size", low=Decimal(8), high=Decimal(24))
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        prim = _agricultural_quality_nonneg(self.primary_defects_count, "primary_defects_count")
        sec = _agricultural_quality_nonneg(self.secondary_defects_count, "secondary_defects_count")
        loss = _agricultural_quality_pct(self.processing_loss_pct, "processing_loss_pct", allow_high=False)
        aw = _agricultural_quality_pct(self.water_activity, "water_activity", low=Decimal(0), high=Decimal(1)) if self.water_activity is not None else None

        object.__setattr__(self, "screen_size", scr)
        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "primary_defects_count", prim)
        object.__setattr__(self, "secondary_defects_count", sec)
        object.__setattr__(self, "processing_loss_pct", loss)
        object.__setattr__(self, "water_activity", aw)

    @property
    def screen_aperture_mm(self) -> Decimal:
        """Physical sieve aperture in millimetres (``screen_size * 25.4 / 64``)."""
        return self.screen_size * Decimal("25.4") / Decimal(64)

    @property
    def processing_yield_pct(self) -> Decimal:
        return Decimal(100) - self.processing_loss_pct

    def full_defect_equivalents(self, *, secondary_per_full_defect: Decimal | float = Decimal("5")) -> Decimal:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        ratio = _agricultural_quality_nonneg(secondary_per_full_defect, "secondary_per_full_defect", strict_positive=True)
        return self.primary_defects_count + (self.secondary_defects_count / ratio)

    def green_to_roasted_mass(self, green_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(green_mass, "green_mass")
        return _agricultural_quality_kg_out(kg * self.processing_yield_pct / Decimal(100), target_unit or green_mass.units, green_mass)

    def roasted_to_green_mass(self, roasted_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(roasted_mass, "roasted_mass")
        return _agricultural_quality_kg_out(kg / (self.processing_yield_pct / Decimal(100)), target_unit or roasted_mass.units, roasted_mass)

    def blend_with(self, other: CoffeeQualitySpec, *, self_mass_fraction: Decimal | float) -> CoffeeQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, CoffeeQualitySpec):
            raise TypeError(f"Expected CoffeeQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_aw = (w1 * self.water_activity + w2 * other.water_activity) if (self.water_activity is not None and other.water_activity is not None) else None
        return CoffeeQualitySpec(
            screen_size=w1 * self.screen_size + w2 * other.screen_size,
            moisture_pct=w1 * self.moisture_pct + w2 * other.moisture_pct,
            primary_defects_count=w1 * self.primary_defects_count + w2 * other.primary_defects_count,
            secondary_defects_count=w1 * self.secondary_defects_count + w2 * other.secondary_defects_count,
            processing_loss_pct=w1 * self.processing_loss_pct + w2 * other.processing_loss_pct,
            water_activity=b_aw,
        )

    def to_material_spec(self, name: str = "Green Coffee") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import MaterialProperties
        from core.units.processing import ProcessingSpecification
        yld = self.processing_yield_pct / Decimal(100)
        loss = self.processing_loss_pct / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=self.moisture_pct / Decimal(100)),
            processing=ProcessingSpecification(
                input_material="green_coffee",
                output_yields={"roasted_coffee": yld},
                process_losses={"roast_and_cleaning_loss": loss} if loss > Decimal(0) else {},
            ),
        )


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 7. Cocoa Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CocoaQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for fermented cocoa beans."""

    bean_count_per_100g: Decimal
    fermentation_pct: Decimal
    fat_content_pct: Decimal  # Fat % of nib on dry/as-is nib basis
    moisture_pct: Decimal = Decimal("7.5")
    defective_beans_pct: Decimal = Decimal("0")  # Mouldy + slaty + insect/germinated %
    shell_pct: Decimal = Decimal("12.0")
    free_fatty_acid_pct: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        cnt = _agricultural_quality_nonneg(self.bean_count_per_100g, "bean_count_per_100g", strict_positive=True)
        ferm = _agricultural_quality_pct(self.fermentation_pct, "fermentation_pct")
        fat = _agricultural_quality_pct(self.fat_content_pct, "fat_content_pct", low=Decimal(20), high=Decimal(70))
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        defects = _agricultural_quality_pct(self.defective_beans_pct, "defective_beans_pct")
        shell = _agricultural_quality_pct(self.shell_pct, "shell_pct", high=Decimal(30))
        ffa = _agricultural_quality_pct(self.free_fatty_acid_pct, "free_fatty_acid_pct", high=Decimal(15)) if self.free_fatty_acid_pct is not None else None

        object.__setattr__(self, "bean_count_per_100g", cnt)
        object.__setattr__(self, "fermentation_pct", ferm)
        object.__setattr__(self, "fat_content_pct", fat)
        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "defective_beans_pct", defects)
        object.__setattr__(self, "shell_pct", shell)
        object.__setattr__(self, "free_fatty_acid_pct", ffa)

    @property
    def mean_bean_mass_grams(self) -> Decimal:
        return Decimal(100) / self.bean_count_per_100g

    @property
    def nib_yield_pct(self) -> Decimal:
        """Clean cocoa nib yield percentage after deshelling/winnowing."""
        return max(Decimal(0), Decimal(100) - self.shell_pct - self.moisture_pct)

    def bean_to_nib_mass(self, bean_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(bean_mass, "bean_mass")
        return _agricultural_quality_kg_out(kg * self.nib_yield_pct / Decimal(100), target_unit or bean_mass.units, bean_mass)

    def extractable_butter_mass(
        self,
        bean_mass: pint.Quantity,
        *,
        press_fat_recovery_pct: Decimal | float = Decimal("88.0"),
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        from core.units.analytical_basis import _agricultural_quality_pct
        rec = _agricultural_quality_pct(press_fat_recovery_pct, "press_fat_recovery_pct") / Decimal(100)
        nib_kg = _agricultural_quality_mass_kg(self.bean_to_nib_mass(bean_mass, target_unit="kilogram"), "nib_mass")
        butter_kg = nib_kg * (self.fat_content_pct / Decimal(100)) * rec
        return _agricultural_quality_kg_out(butter_kg, target_unit or bean_mass.units, bean_mass)

    def blend_with(self, other: CocoaQualitySpec, *, self_mass_fraction: Decimal | float) -> CocoaQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, CocoaQualitySpec):
            raise TypeError(f"Expected CocoaQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        # Bean count per 100 g has mass in the denominator, so mass-weighted average strictly conserves total bean count!
        b_ffa = (w1 * self.free_fatty_acid_pct + w2 * other.free_fatty_acid_pct) if (self.free_fatty_acid_pct is not None and other.free_fatty_acid_pct is not None) else None
        return CocoaQualitySpec(
            bean_count_per_100g=w1 * self.bean_count_per_100g + w2 * other.bean_count_per_100g,
            fermentation_pct=w1 * self.fermentation_pct + w2 * other.fermentation_pct,
            fat_content_pct=w1 * self.fat_content_pct + w2 * other.fat_content_pct,
            moisture_pct=w1 * self.moisture_pct + w2 * other.moisture_pct,
            defective_beans_pct=w1 * self.defective_beans_pct + w2 * other.defective_beans_pct,
            shell_pct=w1 * self.shell_pct + w2 * other.shell_pct,
            free_fatty_acid_pct=b_ffa,
        )

    def to_material_spec(self, name: str = "Cocoa Beans") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import MaterialProperties
        from core.units.processing import ProcessingSpecification
        nib_f = self.nib_yield_pct / Decimal(100)
        loss_f = Decimal(1) - nib_f
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=self.moisture_pct / Decimal(100)),
            processing=ProcessingSpecification(
                input_material="cocoa_beans",
                output_yields={"cocoa_liquor_nib": nib_f},
                process_losses={"shell_and_roast_moisture": loss_f} if loss_f > Decimal(0) else {},
            ),
        )


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 8. Sugar Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SugarQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for raw (VHP / 96° Pol) and refined sugar."""

    polarisation_deg_z: Decimal
    icumsa_colour_iu: Decimal
    brix_pct: Decimal | None = None
    moisture_pct: Decimal = Decimal("0.05")
    ash_pct: Decimal = Decimal("0.04")
    refining_yield_pct: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        pol = _agricultural_quality_pct(self.polarisation_deg_z, "polarisation_deg_z", allow_low=False)
        col = _agricultural_quality_nonneg(self.icumsa_colour_iu, "icumsa_colour_iu")
        brix = _agricultural_quality_pct(self.brix_pct, "brix_pct", allow_low=False) if self.brix_pct is not None else None
        if brix is not None and pol > brix + Decimal("1e-6"):
            raise ValueError(f"polarisation_deg_z ({pol}) cannot exceed total dissolved solids brix_pct ({brix}).")
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        ash = _agricultural_quality_pct(self.ash_pct, "ash_pct", allow_high=False)
        ref_y = _agricultural_quality_pct(self.refining_yield_pct, "refining_yield_pct", allow_low=False) if self.refining_yield_pct is not None else None

        object.__setattr__(self, "polarisation_deg_z", pol)
        object.__setattr__(self, "icumsa_colour_iu", col)
        object.__setattr__(self, "brix_pct", brix)
        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "ash_pct", ash)
        object.__setattr__(self, "refining_yield_pct", ref_y)

    @property
    def apparent_purity_pct(self) -> Decimal | None:
        """Apparent sucrose purity ``100 * Pol / Brix`` when ``brix_pct`` is supplied."""
        if self.brix_pct is None:
            return None
        return (self.polarisation_deg_z / self.brix_pct) * Decimal(100)

    def commercial_rendement_pct(self, *, ash_molasses_factor: Decimal | float = Decimal("4.5")) -> Decimal:
        """Theoretical refinery white-sugar rendement ``Pol - k_ash * Ash``."""
        from core.units.analytical_basis import _agricultural_quality_nonneg
        k = _agricultural_quality_nonneg(ash_molasses_factor, "ash_molasses_factor")
        return max(Decimal(0), self.polarisation_deg_z - k * self.ash_pct)

    @property
    def effective_refining_yield_pct(self) -> Decimal:
        return self.refining_yield_pct if self.refining_yield_pct is not None else self.commercial_rendement_pct()

    def normalize_to_reference_pol(
        self,
        mass: pint.Quantity,
        *,
        reference_pol_deg_z: Decimal | float = Decimal("96.0"),
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Convert raw sugar mass to physical sucrose-equivalent mass at ``reference_pol_deg_z``."""
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        from core.units.analytical_basis import _agricultural_quality_pct
        kg = _agricultural_quality_mass_kg(mass, "mass")
        ref_pol = _agricultural_quality_pct(reference_pol_deg_z, "reference_pol_deg_z", allow_low=False)
        return _agricultural_quality_kg_out(kg * (self.polarisation_deg_z / ref_pol), target_unit or mass.units, mass)

    def raw_to_refined_mass(self, raw_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(raw_mass, "raw_mass")
        return _agricultural_quality_kg_out(kg * self.effective_refining_yield_pct / Decimal(100), target_unit or raw_mass.units, raw_mass)

    def blend_with(self, other: SugarQualitySpec, *, self_mass_fraction: Decimal | float) -> SugarQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, SugarQualitySpec):
            raise TypeError(f"Expected SugarQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_brix = (w1 * self.brix_pct + w2 * other.brix_pct) if (self.brix_pct is not None and other.brix_pct is not None) else None
        b_ref = (w1 * self.refining_yield_pct + w2 * other.refining_yield_pct) if (self.refining_yield_pct is not None and other.refining_yield_pct is not None) else None
        return SugarQualitySpec(
            polarisation_deg_z=w1 * self.polarisation_deg_z + w2 * other.polarisation_deg_z,
            icumsa_colour_iu=w1 * self.icumsa_colour_iu + w2 * other.icumsa_colour_iu,
            brix_pct=b_brix,
            moisture_pct=w1 * self.moisture_pct + w2 * other.moisture_pct,
            ash_pct=w1 * self.ash_pct + w2 * other.ash_pct,
            refining_yield_pct=b_ref,
        )

    def to_material_spec(self, name: str = "Sugar") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        from core.units.processing import ProcessingSpecification
        pol_f = self.polarisation_deg_z / Decimal(100)
        moist_f = self.moisture_pct / Decimal(100)
        ash_f = self.ash_pct / Decimal(100)
        ref_f = self.effective_refining_yield_pct / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(purity=pol_f, moisture=moist_f),
            composition=(
                MaterialComposition(
                    components={"sucrose": pol_f},
                    basis=CompositionBasis.AS_RECEIVED,
                    moisture_fraction=moist_f,
                    ash_fraction=ash_f,
                )
                if (pol_f + moist_f + ash_f) <= Decimal(1)
                else None
            ),
            processing=ProcessingSpecification(
                input_material="raw_sugar",
                output_yields={"refined_white_sugar": ref_f},
                process_losses={"molasses_and_refining_loss": max(Decimal(0), Decimal(1) - ref_f)},
            ),
        )


from core.units.validation import validated_assay


# ---------------------------------------------------------------------------
# 1. Softs Specification (Sugar, Coffee, Cocoa, Cotton)
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class SoftsAssay:
    """Immutable physical, yield, and commercial pricing specification for softs.

    Supports sugar (raw-to-white refining yield & polarisation), coffee (60 kg
    bags & green-to-roasted yield), cocoa (62.5 kg bags & bean-to-butter/powder
    press ratios), and cotton (480 lb / 170 kg bales).
    """

    name: str
    commodity: str
    bag_or_bale_kg: float | None = None
    processing_yield_pct: float = 100.0
    secondary_yield_pct: float = 0.0
    reference_polarisation_deg: float | None = None
    actual_polarisation_deg: float | None = None
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        comm = self.commodity.strip().lower()
        if comm not in {"sugar", "coffee", "cocoa", "cotton"}:
            raise ValueError(f"commodity must be 'sugar', 'coffee', 'cocoa', or 'cotton', got {self.commodity!r}")

        yld = _validate_finite_float(self.processing_yield_pct, "processing_yield_pct")
        if not (0.0 < yld <= 100.0):
            raise ValueError(f"processing_yield_pct must be in (0, 100], got {self.processing_yield_pct}")

        sec_yld = _validate_finite_float(self.secondary_yield_pct, "secondary_yield_pct")
        if not (0.0 <= sec_yld <= 100.0):
            raise ValueError(f"secondary_yield_pct must be in [0, 100], got {self.secondary_yield_pct}")

        if self.bag_or_bale_kg is not None:
            pkg = _validate_finite_float(self.bag_or_bale_kg, "bag_or_bale_kg")
            if pkg <= 0.0:
                raise ValueError(f"bag_or_bale_kg must be strictly positive, got {self.bag_or_bale_kg}")
            object.__setattr__(self, "bag_or_bale_kg", pkg)

        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")
            object.__setattr__(self, "bulk_density_kg_per_m3", rho)

        object.__setattr__(self, "commodity", comm)
        object.__setattr__(self, "processing_yield_pct", yld)
        object.__setattr__(self, "secondary_yield_pct", sec_yld)

    @property
    def processing_yield_fraction(self) -> float:
        return self.processing_yield_pct / 100.0

    @property
    def polarisation_factor(self) -> float:
        """Linear sugar polarisation adjustment factor (actual_pol / reference_pol)."""
        from core.units.agricultural_specs import _validate_finite_float
        if self.reference_polarisation_deg is None or self.actual_polarisation_deg is None:
            return 1.0
        ref = _validate_finite_float(self.reference_polarisation_deg, "reference_polarisation_deg")
        act = _validate_finite_float(self.actual_polarisation_deg, "actual_polarisation_deg")
        if ref <= 0.0 or act <= 0.0:
            raise ValueError("Polarisation degrees must be strictly positive.")
        return act / ref

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_softs_context
        return create_softs_context(self, name=context_name)

    def raw_to_processed_mass(self, raw_mass: pint.Quantity) -> pint.Quantity:
        """Converts raw/unprocessed mass (e.g. raw sugar, green coffee, cocoa beans) to output mass."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(raw_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(raw_mass).__name__}")
        mt_q = raw_mass.to(ureg.metric_ton)
        eff_frac = self.processing_yield_fraction * self.polarisation_factor
        scale = _scale_factor(mt_q.magnitude, eff_frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def processed_to_raw_mass(self, processed_mass: pint.Quantity) -> pint.Quantity:
        """Converts processed output mass (e.g. refined white sugar, roasted coffee) to raw input mass."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(processed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(processed_mass).__name__}")
        mt_q = processed_mass.to(ureg.metric_ton)
        eff_frac = self.processing_yield_fraction * self.polarisation_factor
        scale = _scale_factor(mt_q.magnitude, eff_frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def mass_to_packages(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal package count, preserving Decimal input."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.bag_or_bale_kg is None:
            raise ValueError(f"{self.name!r} does not define bag_or_bale_kg.")
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        kilograms = mass.to(ureg.kilogram)
        number = Decimal(str(kilograms.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Package mass must be finite and non-negative.")
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bag_or_bale_kg)

    def packages_to_mass(self, count: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        if self.bag_or_bale_kg is None:
            raise ValueError(f"{self.name!r} does not define bag_or_bale_kg.")
        number = Decimal(str(count))
        if not number.is_finite() or number < 0:
            raise ValueError("Package count must be finite and non-negative.")
        kilograms = number * Decimal(str(self.bag_or_bale_kg))
        value = kilograms if isinstance(count, Decimal) else float(kilograms)
        return ureg.Quantity(value, ureg.kilogram).to(ureg.metric_ton)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    @staticmethod
    def cents_per_lb_to_usd_per_mt(cents_per_lb: float | Decimal) -> float | Decimal:
        """Converts ICE quotes in US cents/lb (e.g. No. 11 Sugar, Cotton) to $/metric_ton."""
        from core.units.agricultural_specs import _validate_finite_float
        from core.units.precision import preserve_ieee754_scalar
        _validate_finite_float(cents_per_lb, "cents_per_lb")
        val = (Decimal(str(cents_per_lb)) / Decimal(100)) * (Decimal(1000) / Decimal("0.45359237"))
        return preserve_ieee754_scalar(val, cents_per_lb)

    @staticmethod
    def usd_per_mt_to_cents_per_lb(usd_per_mt: float | Decimal) -> float | Decimal:
        """Converts $/metric_ton quote to US cents/lb."""
        from core.units.agricultural_specs import _validate_finite_float
        from core.units.precision import preserve_ieee754_scalar
        _validate_finite_float(usd_per_mt, "usd_per_mt")
        val = Decimal(str(usd_per_mt)) * Decimal("0.45359237") / Decimal(10)
        return preserve_ieee754_scalar(val, usd_per_mt)

    def raw_mt_price_to_refined_mt_price(self, raw_price_per_mt: float | Decimal) -> float | Decimal:
        """Converts raw/unprocessed $/mt quote to implied processed $/mt parity."""
        from core.units.agricultural_specs import _validate_finite_float
        from core.units.precision import preserve_ieee754_scalar
        _validate_finite_float(raw_price_per_mt, "raw_price_per_mt")
        eff_frac = (Decimal(str(self.processing_yield_pct)) / Decimal(100)) * Decimal(str(self.polarisation_factor))
        return preserve_ieee754_scalar(Decimal(str(raw_price_per_mt)) / eff_frac, raw_price_per_mt)

    @classmethod
    def raw_sugar_no11(cls) -> SoftsAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SoftsAssay.raw_sugar_no11")

    @classmethod
    def green_coffee_60kg_bag(cls) -> SoftsAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SoftsAssay.green_coffee_60kg_bag")

    @classmethod
    def cocoa_beans_62_5kg_bag(cls) -> SoftsAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SoftsAssay.cocoa_beans_62_5kg_bag")

    @classmethod
    def us_cotton_480lb_bale(cls) -> SoftsAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SoftsAssay.us_cotton_480lb_bale")

    @classmethod
    def indian_cotton_170kg_bale(cls) -> SoftsAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SoftsAssay.indian_cotton_170kg_bale")


from core.units.agricultural_specs import _HONEY_EXPORT_DRUM_KG
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class HoneyApicultureAssay:
    """Immutable specification for bulk honey, moisture/Brix solids, drums, and beeswax.

    Bridges tank volume <-> mass using a supplied density for the measured
    product state. Density is not inferred from moisture or temperature.
    Normalises mass to a contractual moisture basis, converts net package counts,
    and estimates beeswax cappings relative to extracted honey mass.
    Presets are examples, not a certification of a regulatory honey grade.
    """

    name: str
    floral_source: str  # e.g. "POLYFLORAL", "CLOVER", "ACACIA", "MANUKA"
    moisture_pct: float = 17.5
    base_moisture_pct: float = 18.0
    density_kg_per_m3: float = 1420.0
    drum_mass_kg: float = _HONEY_EXPORT_DRUM_KG
    beeswax_ratio_pct_of_honey: float = 1.8  # Typical cappings beeswax mass as % of extracted honey

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.floral_source, str):
            raise ValueError("floral_source must be a non-empty string.")
        floral = self.floral_source.strip().upper()
        if not floral:
            raise ValueError("floral_source must be a non-empty string.")

        moist = _validate_finite_float(self.moisture_pct, "moisture_pct")
        if not (10.0 <= moist <= 30.0):
            raise ValueError(f"moisture_pct must be in [10.0, 30.0], got {self.moisture_pct}")

        base_m = _validate_finite_float(self.base_moisture_pct, "base_moisture_pct")
        if not (10.0 <= base_m <= 30.0):
            raise ValueError(f"base_moisture_pct must be in [10.0, 30.0], got {self.base_moisture_pct}")

        rho = _validate_finite_float(self.density_kg_per_m3, "density_kg_per_m3")
        if not (1300.0 <= rho <= 1500.0):
            raise ValueError(f"density_kg_per_m3 must be in [1300, 1500], got {self.density_kg_per_m3}")

        drum_kg = _validate_finite_float(self.drum_mass_kg, "drum_mass_kg")
        if drum_kg <= 0.0:
            raise ValueError(f"drum_mass_kg must be strictly positive, got {self.drum_mass_kg}")

        wax_pct = _validate_finite_float(self.beeswax_ratio_pct_of_honey, "beeswax_ratio_pct_of_honey")
        if not (0.0 <= wax_pct <= 25.0):
            raise ValueError(
                f"beeswax_ratio_pct_of_honey must be in [0.0, 25.0], got {self.beeswax_ratio_pct_of_honey}"
            )

        object.__setattr__(self, "floral_source", floral)
        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "base_moisture_pct", base_m)
        object.__setattr__(self, "density_kg_per_m3", rho)
        object.__setattr__(self, "drum_mass_kg", drum_kg)
        object.__setattr__(self, "beeswax_ratio_pct_of_honey", wax_pct)

    @staticmethod
    def _validate_moisture(value: float | Decimal) -> Decimal:
        percentage = Decimal(str(value))
        if not percentage.is_finite() or not (10 <= percentage <= 30):
            raise ValueError("Honey moisture percentage must be finite and in [10, 30].")
        return percentage

    @staticmethod
    def moisture_to_brix(moisture_pct: float | Decimal) -> float | Decimal:
        """Approximate soluble-solids percentage; not a refractometer correction."""
        value = 100 - HoneyApicultureAssay._validate_moisture(moisture_pct)
        return value if isinstance(moisture_pct, Decimal) else float(value)

    @staticmethod
    def brix_to_moisture(brix_pct: float | Decimal) -> float | Decimal:
        """Inverse of the approximate 100 - moisture relationship for honey."""
        value = HoneyApicultureAssay._validate_moisture(100 - Decimal(str(brix_pct)))
        return value if isinstance(brix_pct, Decimal) else float(value)

    @property
    def brix_solids_pct(self) -> float:
        return float(self.moisture_to_brix(self.moisture_pct))

    def _solids_factor(self, actual_moisture_pct: float | Decimal | None = None) -> Decimal:
        from core.units.analytical_basis import moisture_level_ratio
        actual = self._validate_moisture(self.moisture_pct if actual_moisture_pct is None else actual_moisture_pct)
        return moisture_level_ratio(
            source_moisture_pct=actual,
            target_moisture_pct=Decimal(str(self.base_moisture_pct)),
        )

    @property
    def solids_normalisation_factor(self) -> float:
        return float(self._solids_factor())

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        value = Decimal(str(quantity.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Honey mass must be finite and non-negative.")
        return quantity

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_honey_apiculture_context
        return create_honey_apiculture_context(self, name=context_name)

    # --- Tank Volume <-> Mass ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        quantity = volume.to(ureg.meter ** 3)
        value = Decimal(str(quantity.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Honey volume must be finite and non-negative.")
        density = Decimal(str(self.density_kg_per_m3)) / 1000
        factor = density if isinstance(quantity.magnitude, Decimal) else float(density)
        return ureg.Quantity(quantity.magnitude * factor, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        quantity = self._validate_mass(mass)
        density = Decimal(str(self.density_kg_per_m3)) / 1000
        factor = density if isinstance(quantity.magnitude, Decimal) else float(density)
        return ureg.Quantity(quantity.magnitude / factor, ureg.meter ** 3)

    # --- Solids Normalisation & Beeswax By-Product ---
    def normalise_to_base_moisture(
        self, honey_mass: pint.Quantity, *, actual_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Returns constant-solids equivalent mass; does not alter stored density."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        quantity = self._validate_mass(honey_mass)
        actual = self._validate_moisture(self.moisture_pct if actual_moisture_pct is None else actual_moisture_pct)
        return convert_moisture_mass(
            quantity,
            source_moisture_pct=actual,
            target_moisture_pct=Decimal(str(self.base_moisture_pct)),
            target_unit=ureg.metric_ton,
        )

    def reverse_moisture_normalisation(
        self, base_mass: pint.Quantity, *, actual_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        quantity = self._validate_mass(base_mass)
        actual = self._validate_moisture(self.moisture_pct if actual_moisture_pct is None else actual_moisture_pct)
        return convert_moisture_mass(
            quantity,
            source_moisture_pct=Decimal(str(self.base_moisture_pct)),
            target_moisture_pct=actual,
            target_unit=ureg.metric_ton,
        )

    def honey_to_beeswax_mass(self, honey_mass: pint.Quantity) -> pint.Quantity:
        """Estimates associated cappings wax; wax is not deducted from honey mass."""
        from core.units.registry import ureg
        kilograms = self._validate_mass(honey_mass).to(ureg.kilogram)
        ratio = Decimal(str(self.beeswax_ratio_pct_of_honey)) / 100
        scale = ratio if isinstance(kilograms.magnitude, Decimal) else float(ratio)
        return ureg.Quantity(kilograms.magnitude * scale, ureg.kilogram)

    # --- Net Honey Mass & Nominal Packaging Counts ---
    def _mass_to_count(self, mass: pint.Quantity, package_kg: float) -> float | Decimal:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, package_kg)

    @staticmethod
    def _count_to_mass(count: float | Decimal, package_kg: float) -> pint.Quantity:
        from core.units.registry import ureg
        number = Decimal(str(count))
        if not number.is_finite() or number < 0:
            raise ValueError("Package count must be finite and non-negative.")
        kilograms = number * Decimal(str(package_kg))
        value = kilograms if isinstance(count, Decimal) else float(kilograms)
        return ureg.Quantity(value, ureg.kilogram).to(ureg.metric_ton)

    def mass_to_drums(self, mass: pint.Quantity) -> float | Decimal:
        return self._mass_to_count(mass, self.drum_mass_kg)

    def drums_to_mass(self, drum_count: float | Decimal) -> pint.Quantity:
        return self._count_to_mass(drum_count, self.drum_mass_kg)

    def mass_to_ibc_totes(self, mass: pint.Quantity) -> float | Decimal:
        from core.units.agricultural_specs import _HONEY_IBC_TOTE_KG
        return self._mass_to_count(mass, _HONEY_IBC_TOTE_KG)

    def ibc_totes_to_mass(self, tote_count: float | Decimal) -> pint.Quantity:
        from core.units.agricultural_specs import _HONEY_IBC_TOTE_KG
        return self._count_to_mass(tote_count, _HONEY_IBC_TOTE_KG)

    def mass_to_pails(self, mass: pint.Quantity) -> float | Decimal:
        from core.units.agricultural_specs import _HONEY_PAIL_60LB_KG
        return self._mass_to_count(mass, _HONEY_PAIL_60LB_KG)

    def pails_to_mass(self, pail_count: float | Decimal) -> pint.Quantity:
        from core.units.agricultural_specs import _HONEY_PAIL_60LB_KG
        return self._count_to_mass(pail_count, _HONEY_PAIL_60LB_KG)

    # --- Pricing Parity ---
    @staticmethod
    def _price_parity(price: float | Decimal, *, per_pound_to_tonne: bool) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        value = Decimal(str(price))
        if not value.is_finite():
            raise ValueError("Honey price must be finite.")
        pounds_per_tonne = Decimal(1000) / Decimal("0.45359237")
        output = value * pounds_per_tonne if per_pound_to_tonne else value / pounds_per_tonne
        return preserve_ieee754_scalar(output, price)

    @staticmethod
    def usd_per_lb_to_usd_per_mt(usd_per_lb: float | Decimal) -> float | Decimal:
        return HoneyApicultureAssay._price_parity(usd_per_lb, per_pound_to_tonne=True)

    @staticmethod
    def usd_per_mt_to_usd_per_lb(usd_per_mt: float | Decimal) -> float | Decimal:
        return HoneyApicultureAssay._price_parity(usd_per_mt, per_pound_to_tonne=False)

    # --- Standard Honey Presets ---
    @classmethod
    def bulk_polyfloral_export(cls) -> HoneyApicultureAssay:
        """Standard bulk export polyfloral honey (17.5% moisture, 1,420 kg/m^3, 300 kg drum)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "HoneyApicultureAssay.bulk_polyfloral_export")

    @classmethod
    def us_grade_a_clover(cls) -> HoneyApicultureAssay:
        """US Grade A clover honey (18.0% moisture, 1,417 kg/m^3, 60 lb pail option)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "HoneyApicultureAssay.us_grade_a_clover")


__all__ = ('CocoaQualitySpec', 'CoffeeQualitySpec', 'HoneyApicultureAssay', 'SoftsAssay', 'SugarQualitySpec')
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
