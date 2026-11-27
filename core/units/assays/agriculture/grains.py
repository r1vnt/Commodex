"""Commodex assays agriculture grains implementation and public interfaces."""
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
from core.units.validation import _safety_QuantitySnapshot


# ---------------------------------------------------------------------------
# 1. Wheat Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class WheatQualitySpec(_safety_QuantitySnapshot, AgriculturalQualityBase):
    """Physical quality specification for milling and feed wheat.

    Attributes cover protein (with explicit reference moisture basis, e.g. 12.0% US FGIS,
    13.5% CGC Canada, or 0.0% dry basis), Hagberg-Perten falling number, wet/dry gluten,
    gluten index, moisture, and hectolitre / Winchester bushel test weight.
    """

    moisture_pct: Decimal
    protein_pct: Decimal
    protein_moisture_basis_pct: Decimal = Decimal("12.0")
    falling_number_s: Decimal | None = None
    wet_gluten_pct: Decimal | None = None
    dry_gluten_pct: Decimal | None = None
    gluten_index: Decimal | None = None
    test_weight: Quantity | None = None  # Accepts kg/hL, kg/m^3, or lb/bu via from_test_weight_lb_bu

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_dec
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        prot = _agricultural_quality_pct(self.protein_pct, "protein_pct", high=Decimal(35))
        p_basis = _agricultural_quality_pct(self.protein_moisture_basis_pct, "protein_moisture_basis_pct", allow_high=False)
        fn = _agricultural_quality_nonneg(self.falling_number_s, "falling_number_s", strict_positive=True) if self.falling_number_s is not None else None
        if fn is not None and fn < Decimal(55):
            raise ValueError(f"Hagberg-Perten falling_number_s must be >= 55 s, got {fn}.")
        wg = _agricultural_quality_pct(self.wet_gluten_pct, "wet_gluten_pct", high=Decimal(60)) if self.wet_gluten_pct is not None else None
        dg = _agricultural_quality_pct(self.dry_gluten_pct, "dry_gluten_pct", high=Decimal(30)) if self.dry_gluten_pct is not None else None
        if wg is not None and dg is not None and dg > wg:
            raise ValueError(f"dry_gluten_pct ({dg}) cannot exceed wet_gluten_pct ({wg}).")
        gi = _agricultural_quality_pct(self.gluten_index, "gluten_index") if self.gluten_index is not None else None

        tw: Quantity | None = None
        if self.test_weight is not None:
            if not isinstance(self.test_weight, pint.Quantity):
                raise TypeError("test_weight must be a pint.Quantity with mass/volume dimensionality (e.g. 'kg/hL').")
            validate_registry_quantity(self.test_weight)
            kg_m3 = _agricultural_quality_nonneg(self.test_weight.to("kg / meter ** 3").magnitude, "test_weight", strict_positive=True)
            if not (Decimal(400) <= kg_m3 <= Decimal(1000)):
                raise ValueError(f"Wheat test_weight ({kg_m3} kg/m^3) is outside physical range [400, 1000] kg/m^3.")
            tw = ureg.Quantity(_agricultural_quality_dec(self.test_weight.magnitude, "test_weight"), self.test_weight.units)

        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "protein_pct", prot)
        object.__setattr__(self, "protein_moisture_basis_pct", p_basis)
        object.__setattr__(self, "falling_number_s", fn)
        object.__setattr__(self, "wet_gluten_pct", wg)
        object.__setattr__(self, "dry_gluten_pct", dg)
        object.__setattr__(self, "gluten_index", gi)
        object.__setattr__(self, "test_weight", tw)

    @classmethod
    def from_winchester_test_weight(
        cls,
        *,
        moisture_pct: Decimal | float,
        protein_pct: Decimal | float,
        test_weight_lb_per_bu: Decimal | float,
        **kwargs: Any,
    ) -> WheatQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_KG_PER_HL_PER_LB_WINCHESTER_BU
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.registry import ureg
        lb_bu = _agricultural_quality_nonneg(test_weight_lb_per_bu, "test_weight_lb_per_bu", strict_positive=True)
        kg_hl = lb_bu * _agricultural_quality_KG_PER_HL_PER_LB_WINCHESTER_BU
        return cls(
            moisture_pct=moisture_pct,  # type: ignore[arg-type]
            protein_pct=protein_pct,  # type: ignore[arg-type]
            test_weight=ureg.Quantity(kg_hl, "kg / hectoliter"),
            **kwargs,
        )

    @property
    def test_weight_kg_per_hl(self) -> Decimal | None:
        from core.units.analytical_basis import _agricultural_quality_dec
        if self.test_weight is None:
            return None
        return _agricultural_quality_dec(self.test_weight.to("kg / hectoliter").magnitude, "test_weight")

    @property
    def test_weight_lb_per_winchester_bu(self) -> Decimal | None:
        from core.units.analytical_basis import _agricultural_quality_KG_PER_HL_PER_LB_WINCHESTER_BU
        kg_hl = self.test_weight_kg_per_hl
        return None if kg_hl is None else kg_hl / _agricultural_quality_KG_PER_HL_PER_LB_WINCHESTER_BU

    @property
    def liquefaction_number(self) -> Decimal | None:
        """Perten Liquefaction Number ``LN = 6000 / (FN - 50)``, additive by mass for blending."""
        if self.falling_number_s is None:
            return None
        return Decimal(6000) / (self.falling_number_s - Decimal(50))

    def protein_at_moisture_basis(self, target_moisture_basis_pct: Decimal | float) -> Decimal:
        """Project protein percentage onto ``target_moisture_basis_pct`` (0 for dry basis)."""
        from core.units.analytical_basis import rebase_constituent_pct
        return rebase_constituent_pct(
            self.protein_pct,
            source_moisture_pct=self.protein_moisture_basis_pct,
            target_moisture_pct=target_moisture_basis_pct,
        )

    @property
    def protein_dry_basis_pct(self) -> Decimal:
        return self.protein_at_moisture_basis(Decimal(0))

    @property
    def protein_as_is_pct(self) -> Decimal:
        return self.protein_at_moisture_basis(self.moisture_pct)

    def normalized_mass(
        self,
        mass: pint.Quantity,
        *,
        base_moisture_pct: Decimal | float = Decimal("13.5"),
        allow_moisture_gain: bool = False,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        from core.units.analytical_basis import normalize_mass_to_moisture
        return normalize_mass_to_moisture(
            mass,
            actual_moisture_pct=self.moisture_pct,
            base_moisture_pct=base_moisture_pct,
            allow_moisture_gain=allow_moisture_gain,
            target_unit=target_unit,
        )

    def contained_protein_mass(self, mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(mass, "mass")
        prot_kg = kg * (self.protein_as_is_pct / Decimal(100))
        return _agricultural_quality_kg_out(prot_kg, target_unit or mass.units, mass)

    def compare_to_benchmark(self, benchmark: WheatQualitySpec) -> dict[str, Decimal]:
        from core.units.analytical_basis import AgriculturalQualityBase
        if not isinstance(benchmark, WheatQualitySpec):
            raise TypeError(f"Expected WheatQualitySpec, got {type(benchmark).__name__}.")
        deltas = AgriculturalQualityBase.compare_to_benchmark(self, benchmark)
        aligned_prot = self.protein_at_moisture_basis(benchmark.protein_moisture_basis_pct)
        deltas["protein_pct"] = aligned_prot - benchmark.protein_pct
        deltas.pop("protein_moisture_basis_pct", None)
        if self.test_weight_kg_per_hl is not None and benchmark.test_weight_kg_per_hl is not None:
            deltas["test_weight_kg_per_hl"] = self.test_weight_kg_per_hl - benchmark.test_weight_kg_per_hl
        return deltas

    def blend_with(self, other: WheatQualitySpec, *, self_mass_fraction: Decimal | float) -> WheatQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        from core.units.analytical_basis import _agricultural_quality_dec
        from core.units.analytical_basis import rebase_constituent_pct
        from core.units.registry import ureg
        if not isinstance(other, WheatQualitySpec):
            raise TypeError(f"Expected WheatQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_moist = w1 * self.moisture_pct + w2 * other.moisture_pct

        # Conserve true physical protein mass across wet lots, then express on self's basis
        prot_as_is_blend = w1 * self.protein_as_is_pct + w2 * other.protein_as_is_pct
        b_prot = rebase_constituent_pct(
            prot_as_is_blend,
            source_moisture_pct=b_moist,
            target_moisture_pct=self.protein_moisture_basis_pct,
        )

        # Falling number blends linearly in Perten Liquefaction Number space (dry-flour basis)
        b_fn: Decimal | None = None
        if self.liquefaction_number is not None and other.liquefaction_number is not None:
            d1 = w1 * (Decimal(100) - self.moisture_pct)
            d2 = w2 * (Decimal(100) - other.moisture_pct)
            ln_blend = (d1 * self.liquefaction_number + d2 * other.liquefaction_number) / (d1 + d2)
            b_fn = (Decimal(6000) / ln_blend) + Decimal(50)

        b_wg = (w1 * self.wet_gluten_pct + w2 * other.wet_gluten_pct) if (self.wet_gluten_pct is not None and other.wet_gluten_pct is not None) else None
        b_dg = (w1 * self.dry_gluten_pct + w2 * other.dry_gluten_pct) if (self.dry_gluten_pct is not None and other.dry_gluten_pct is not None) else None
        b_gi = (w1 * self.gluten_index + w2 * other.gluten_index) if (self.gluten_index is not None and other.gluten_index is not None) else None

        b_tw: Quantity | None = None
        if self.test_weight is not None and other.test_weight is not None:
            rho1 = _agricultural_quality_dec(self.test_weight.to("kg / hectoliter").magnitude, "tw1")
            rho2 = _agricultural_quality_dec(other.test_weight.to("kg / hectoliter").magnitude, "tw2")
            rho_blend = Decimal(1) / ((w1 / rho1) + (w2 / rho2))
            b_tw = ureg.Quantity(rho_blend, "kg / hectoliter").to(self.test_weight.units)

        return WheatQualitySpec(
            moisture_pct=b_moist,
            protein_pct=b_prot,
            protein_moisture_basis_pct=self.protein_moisture_basis_pct,
            falling_number_s=b_fn,
            wet_gluten_pct=b_wg,
            dry_gluten_pct=b_dg,
            gluten_index=b_gi,
            test_weight=b_tw,
        )

    def to_material_spec(self, name: str = "Wheat") -> MaterialSpec:
        from core.units.analytical_basis import rebase_constituent_pct
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        moist_f = self.moisture_pct / Decimal(100)
        comps = {"protein": self.protein_dry_basis_pct / Decimal(100)}
        if self.dry_gluten_pct is not None:
            comps["dry_gluten"] = rebase_constituent_pct(
                self.dry_gluten_pct,
                source_moisture_pct=self.protein_moisture_basis_pct,
                target_moisture_pct=Decimal(0),
            ) / Decimal(100)
        bulk_rho = self.test_weight.to("kg / meter ** 3") if self.test_weight is not None else None
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(bulk_density=bulk_rho, moisture=moist_f),
            composition=MaterialComposition(components=comps, basis=CompositionBasis.DRY, moisture_fraction=moist_f),
        )


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 2. Corn / Maize Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class CornQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for yellow/white corn (maize)."""

    moisture_pct: Decimal
    broken_kernels_pct: Decimal = Decimal("0")
    foreign_material_pct: Decimal = Decimal("0")
    damaged_kernels_pct: Decimal = Decimal("0")
    heat_damaged_pct: Decimal = Decimal("0")
    aflatoxin_ppb: Decimal | None = None
    test_weight_lb_per_bu: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        bk = _agricultural_quality_pct(self.broken_kernels_pct, "broken_kernels_pct")
        fm = _agricultural_quality_pct(self.foreign_material_pct, "foreign_material_pct")
        dmg = _agricultural_quality_pct(self.damaged_kernels_pct, "damaged_kernels_pct")
        hdmg = _agricultural_quality_pct(self.heat_damaged_pct, "heat_damaged_pct")
        if hdmg > dmg:
            raise ValueError(f"heat_damaged_pct ({hdmg}) cannot exceed total damaged_kernels_pct ({dmg}).")
        if bk + fm + moist > Decimal(100):
            raise ValueError("broken_kernels_pct + foreign_material_pct + moisture_pct cannot exceed 100%.")
        afla = _agricultural_quality_nonneg(self.aflatoxin_ppb, "aflatoxin_ppb") if self.aflatoxin_ppb is not None else None
        tw = _agricultural_quality_nonneg(self.test_weight_lb_per_bu, "test_weight_lb_per_bu", strict_positive=True) if self.test_weight_lb_per_bu is not None else None

        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "broken_kernels_pct", bk)
        object.__setattr__(self, "foreign_material_pct", fm)
        object.__setattr__(self, "damaged_kernels_pct", dmg)
        object.__setattr__(self, "heat_damaged_pct", hdmg)
        object.__setattr__(self, "aflatoxin_ppb", afla)
        object.__setattr__(self, "test_weight_lb_per_bu", tw)

    @property
    def bcfm_pct(self) -> Decimal:
        """Combined Broken Corn and Foreign Material (BCFM %)."""
        return self.broken_kernels_pct + self.foreign_material_pct

    def meets_aflatoxin_limit(self, max_ppb: Decimal | float = Decimal("20.0")) -> bool:
        from core.units.analytical_basis import _agricultural_quality_nonneg
        if self.aflatoxin_ppb is None:
            raise ValueError("aflatoxin_ppb is not assayed on this CornQualitySpec.")
        return self.aflatoxin_ppb <= _agricultural_quality_nonneg(max_ppb, "max_ppb")

    def normalized_clean_mass(
        self,
        mass: pint.Quantity,
        *,
        base_moisture_pct: Decimal | float = Decimal("15.5"),
        deductible_fm_threshold_pct: Decimal | float = Decimal("0.0"),
        allow_moisture_gain: bool = False,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Deduct excess foreign-material dockage above threshold, then normalise to ``base_moisture_pct``."""
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        from core.units.analytical_basis import _agricultural_quality_pct
        from core.units.analytical_basis import normalize_mass_to_moisture
        from core.units.registry import ureg
        kg = _agricultural_quality_mass_kg(mass, "mass")
        fm_thresh = _agricultural_quality_pct(deductible_fm_threshold_pct, "deductible_fm_threshold_pct")
        excess_fm = max(Decimal(0), self.foreign_material_pct - fm_thresh) / Decimal(100)
        clean_kg = kg * (Decimal(1) - excess_fm)
        return normalize_mass_to_moisture(
            ureg.Quantity(clean_kg, ureg.kilogram),
            actual_moisture_pct=self.moisture_pct,
            base_moisture_pct=base_moisture_pct,
            allow_moisture_gain=allow_moisture_gain,
            target_unit=target_unit or mass.units,
        )

    def blend_with(self, other: CornQualitySpec, *, self_mass_fraction: Decimal | float) -> CornQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, CornQualitySpec):
            raise TypeError(f"Expected CornQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_afla = (w1 * self.aflatoxin_ppb + w2 * other.aflatoxin_ppb) if (self.aflatoxin_ppb is not None and other.aflatoxin_ppb is not None) else None
        b_tw = (Decimal(1) / ((w1 / self.test_weight_lb_per_bu) + (w2 / other.test_weight_lb_per_bu))) if (self.test_weight_lb_per_bu is not None and other.test_weight_lb_per_bu is not None) else None
        return CornQualitySpec(
            moisture_pct=w1 * self.moisture_pct + w2 * other.moisture_pct,
            broken_kernels_pct=w1 * self.broken_kernels_pct + w2 * other.broken_kernels_pct,
            foreign_material_pct=w1 * self.foreign_material_pct + w2 * other.foreign_material_pct,
            damaged_kernels_pct=w1 * self.damaged_kernels_pct + w2 * other.damaged_kernels_pct,
            heat_damaged_pct=w1 * self.heat_damaged_pct + w2 * other.heat_damaged_pct,
            aflatoxin_ppb=b_afla,
            test_weight_lb_per_bu=b_tw,
        )

    def to_material_spec(self, name: str = "Corn") -> MaterialSpec:
        from core.units.analytical_basis import _agricultural_quality_KG_PER_HL_PER_LB_WINCHESTER_BU
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        from core.units.registry import ureg
        moist_f = self.moisture_pct / Decimal(100)
        bulk_rho = (
            ureg.Quantity(self.test_weight_lb_per_bu * _agricultural_quality_KG_PER_HL_PER_LB_WINCHESTER_BU * Decimal(10), "kg / meter ** 3")
            if self.test_weight_lb_per_bu is not None
            else None
        )
        imps = {
            "broken_kernels": self.broken_kernels_pct / Decimal(100),
            "foreign_material": self.foreign_material_pct / Decimal(100),
        }
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(bulk_density=bulk_rho, moisture=moist_f),
            composition=MaterialComposition(basis=CompositionBasis.AS_RECEIVED, moisture_fraction=moist_f, impurities=imps),
        )


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 5. Barley Quality Specification (Malting vs Feed)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class BarleyQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for two-row/six-row malting and feed barley."""

    moisture_pct: Decimal
    protein_pct: Decimal
    germination_energy_pct: Decimal
    protein_dry_basis: bool = True
    plump_kernels_pct: Decimal | None = None
    skinnings_pct: Decimal = Decimal("0")
    fine_extract_db_pct: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_pct
        if not isinstance(self.protein_dry_basis, bool):
            raise TypeError("protein_dry_basis must be a bool.")
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        prot = _agricultural_quality_pct(self.protein_pct, "protein_pct", high=Decimal(30))
        germ = _agricultural_quality_pct(self.germination_energy_pct, "germination_energy_pct")
        plump = _agricultural_quality_pct(self.plump_kernels_pct, "plump_kernels_pct") if self.plump_kernels_pct is not None else None
        skin = _agricultural_quality_pct(self.skinnings_pct, "skinnings_pct")
        ext = _agricultural_quality_pct(self.fine_extract_db_pct, "fine_extract_db_pct") if self.fine_extract_db_pct is not None else None

        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "protein_pct", prot)
        object.__setattr__(self, "germination_energy_pct", germ)
        object.__setattr__(self, "plump_kernels_pct", plump)
        object.__setattr__(self, "skinnings_pct", skin)
        object.__setattr__(self, "fine_extract_db_pct", ext)

    @property
    def protein_db_pct(self) -> Decimal:
        from core.units.analytical_basis import rebase_constituent_pct
        if self.protein_dry_basis:
            return self.protein_pct
        return rebase_constituent_pct(self.protein_pct, source_moisture_pct=self.moisture_pct, target_moisture_pct=Decimal(0))

    @property
    def protein_as_is_pct(self) -> Decimal:
        from core.units.analytical_basis import rebase_constituent_pct
        if not self.protein_dry_basis:
            return self.protein_pct
        return rebase_constituent_pct(self.protein_pct, source_moisture_pct=Decimal(0), target_moisture_pct=self.moisture_pct)

    def meets_malting_thresholds(
        self,
        *,
        min_germination_pct: Decimal | float = Decimal("95.0"),
        min_protein_db_pct: Decimal | float = Decimal("9.5"),
        max_protein_db_pct: Decimal | float = Decimal("11.5"),
        max_moisture_pct: Decimal | float = Decimal("14.5"),
        min_plump_kernels_pct: Decimal | float | None = None,
    ) -> QualityAuditResult:
        """Evaluate physical malting suitability bounds (returns signed deviations, no pricing)."""
        from core.units.analytical_basis import QualityThreshold
        from core.units.analytical_basis import _agricultural_quality_pct
        rules: dict[str, QualityThreshold] = {
            "germination_energy_pct": QualityThreshold(min_value=_agricultural_quality_pct(min_germination_pct, "min_germination_pct")),
            "protein_db_pct": QualityThreshold(
                min_value=_agricultural_quality_pct(min_protein_db_pct, "min_protein_db_pct"),
                max_value=_agricultural_quality_pct(max_protein_db_pct, "max_protein_db_pct"),
            ),
            "moisture_pct": QualityThreshold(max_value=_agricultural_quality_pct(max_moisture_pct, "max_moisture_pct")),
        }
        if min_plump_kernels_pct is not None:
            rules["plump_kernels_pct"] = QualityThreshold(min_value=_agricultural_quality_pct(min_plump_kernels_pct, "min_plump_kernels_pct"))
        return self.evaluate_thresholds(rules)

    def extractable_malt_solids_mass(self, barley_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_dec
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import wet_to_dry_mass
        if self.fine_extract_db_pct is None:
            raise ValueError("fine_extract_db_pct is required to compute extractable malt solids.")
        dry_q = wet_to_dry_mass(barley_mass, moisture_pct=self.moisture_pct, target_unit="kilogram")
        dry_kg = _agricultural_quality_dec(dry_q.magnitude, "dry_barley_mass")
        return _agricultural_quality_kg_out(dry_kg * self.fine_extract_db_pct / Decimal(100), target_unit or barley_mass.units, barley_mass)

    def blend_with(self, other: BarleyQualitySpec, *, self_mass_fraction: Decimal | float) -> BarleyQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        from core.units.analytical_basis import rebase_constituent_pct
        if not isinstance(other, BarleyQualitySpec):
            raise TypeError(f"Expected BarleyQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_moist = w1 * self.moisture_pct + w2 * other.moisture_pct
        prot_as_is = w1 * self.protein_as_is_pct + w2 * other.protein_as_is_pct
        b_prot = (
            rebase_constituent_pct(prot_as_is, source_moisture_pct=b_moist, target_moisture_pct=Decimal(0))
            if self.protein_dry_basis
            else prot_as_is
        )
        b_plump = (w1 * self.plump_kernels_pct + w2 * other.plump_kernels_pct) if (self.plump_kernels_pct is not None and other.plump_kernels_pct is not None) else None
        d1 = w1 * (Decimal(100) - self.moisture_pct)
        d2 = w2 * (Decimal(100) - other.moisture_pct)
        b_ext = ((d1 * self.fine_extract_db_pct + d2 * other.fine_extract_db_pct) / (d1 + d2)) if (self.fine_extract_db_pct is not None and other.fine_extract_db_pct is not None) else None
        return BarleyQualitySpec(
            moisture_pct=b_moist,
            protein_pct=b_prot,
            germination_energy_pct=w1 * self.germination_energy_pct + w2 * other.germination_energy_pct,
            protein_dry_basis=self.protein_dry_basis,
            plump_kernels_pct=b_plump,
            skinnings_pct=w1 * self.skinnings_pct + w2 * other.skinnings_pct,
            fine_extract_db_pct=b_ext,
        )

    def to_material_spec(self, name: str = "Barley") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        moist_f = self.moisture_pct / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=moist_f),
            composition=MaterialComposition(
                components={"protein": self.protein_db_pct / Decimal(100)},
                basis=CompositionBasis.DRY,
                moisture_fraction=moist_f,
            ),
        )


from core.units.validation import validated_assay


# ---------------------------------------------------------------------------
# 4. Grain & Oilseed Bushel Weights and Moisture Shrinkage Specification
# (Corn 56 lb/bu, Wheat & Soybeans 60 lb/bu, Barley 48 lb/bu, Oats 32 lb/bu)
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class GrainAssay:
    """Immutable specification for grain bushel weights and moisture deductions.

    Bridges mass (mt, kg, lb) <-> bushels at crop-specific bushel weights,
    applies moisture shrinkage / base-moisture normalisation, and computes
    drying charges or price adjustments for excess moisture.
    """

    name: str
    grain: str
    bushel_weight_lb: float
    base_moisture_pct: float = 14.0
    actual_moisture_pct: float | None = None
    handling_loss_pct: float = 0.0
    allow_moisture_gain: bool = False
    drying_cost_per_mt_per_pct: float = 0.0
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.grain, str) or not self.grain.strip():
            raise ValueError("grain must be a non-empty string.")
        object.__setattr__(self, "grain", self.grain.strip().lower())
        bu_lb = _validate_finite_float(self.bushel_weight_lb, "bushel_weight_lb")
        if bu_lb <= 0.0:
            raise ValueError(f"bushel_weight_lb must be strictly positive, got {self.bushel_weight_lb}")

        base_m = _validate_finite_float(self.base_moisture_pct, "base_moisture_pct")
        if not (0.0 <= base_m < 100.0):
            raise ValueError(f"base_moisture_pct must be in [0, 100), got {self.base_moisture_pct}")

        act_m: float | None = None
        if self.actual_moisture_pct is not None:
            act_m = _validate_finite_float(self.actual_moisture_pct, "actual_moisture_pct")
            if not (0.0 <= act_m < 100.0):
                raise ValueError(f"actual_moisture_pct must be in [0, 100), got {self.actual_moisture_pct}")

        loss_pct = _validate_finite_float(self.handling_loss_pct, "handling_loss_pct")
        if not (0.0 <= loss_pct < 50.0):
            raise ValueError(f"handling_loss_pct must be in [0, 50), got {self.handling_loss_pct}")

        dry_cost = _validate_finite_float(self.drying_cost_per_mt_per_pct, "drying_cost_per_mt_per_pct")
        if dry_cost < 0.0:
            raise ValueError(
                f"drying_cost_per_mt_per_pct must be non-negative, got {self.drying_cost_per_mt_per_pct}"
            )

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        object.__setattr__(self, "bushel_weight_lb", bu_lb)
        object.__setattr__(self, "base_moisture_pct", base_m)
        object.__setattr__(self, "actual_moisture_pct", act_m)
        object.__setattr__(self, "handling_loss_pct", loss_pct)
        object.__setattr__(self, "drying_cost_per_mt_per_pct", dry_cost)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)

    def _effective_actual_moisture(self, actual_moisture_pct: float | Decimal | None = None) -> float:
        from core.units.agricultural_specs import _validate_finite_float
        if actual_moisture_pct is not None:
            m = _validate_finite_float(actual_moisture_pct, "actual_moisture_pct")
            if not (0.0 <= m < 100.0):
                raise ValueError(f"actual_moisture_pct must be in [0, 100), got {actual_moisture_pct}")
            return m
        if self.actual_moisture_pct is not None:
            return self.actual_moisture_pct
        return self.base_moisture_pct

    def _moisture_adjustment_decimal(self, actual_moisture_pct: float | Decimal | None = None) -> Decimal:
        from core.units.analytical_basis import moisture_level_ratio
        moisture = self._effective_actual_moisture(actual_moisture_pct)
        actual = Decimal(str(actual_moisture_pct if actual_moisture_pct is not None else moisture))
        return moisture_level_ratio(
            source_moisture_pct=actual,
            target_moisture_pct=Decimal(str(self.base_moisture_pct)),
            allow_moisture_gain=self.allow_moisture_gain,
            handling_loss_fraction=Decimal(str(self.handling_loss_pct)) / 100,
        )

    def moisture_adjustment_factor(self, actual_moisture_pct: float | Decimal | None = None) -> float:
        """Returns delivered-to-base mass multiplier under this contract's rules.

        Handling loss is subtracted when a moisture adjustment applies. Grain
        exactly at base moisture is unchanged. Below-base adjustments are
        disabled unless allow_moisture_gain=True. Contract presets are examples;
        handling losses and drying rates must reflect the actual contract.
        """
        return float(self._moisture_adjustment_decimal(actual_moisture_pct))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_grain_context
        return create_grain_context(self, name=context_name)

    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        result = mass.to(ureg.metric_ton)
        magnitude = Decimal(str(result.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Grain mass must be finite and non-negative.")
        return result

    def mass_to_bushels(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded crop-weight bushel count, preserving Decimal input.

        This explicit count is distinct from Pint's generic volumetric bushel.
        """
        from core.units.registry import ureg
        pounds = self._validate_mass(mass).to(ureg.pound)
        weight = Decimal(str(self.bushel_weight_lb)) if isinstance(pounds.magnitude, Decimal) else self.bushel_weight_lb
        return pounds.magnitude / weight

    def bushels_to_mass(self, bushels: float | Decimal) -> pint.Quantity:
        """Returns mass from crop-weight bushel count, preserving Decimal input."""
        from core.units.registry import ureg
        count = Decimal(str(bushels))
        if not count.is_finite() or count < 0:
            raise ValueError("bushels must be finite and non-negative.")
        if isinstance(bushels, Decimal):
            pounds = bushels * Decimal(str(self.bushel_weight_lb))
        else:
            pounds = float(bushels) * self.bushel_weight_lb
        return ureg.Quantity(pounds, ureg.pound).to(ureg.metric_ton)

    def apply_moisture_deduction(
        self,
        mass: pint.Quantity,
        *,
        actual_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Adjusts grain mass from actual moisture to `base_moisture_pct` (metric_ton)."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        moisture = self._effective_actual_moisture(actual_moisture_pct)
        actual = Decimal(str(actual_moisture_pct if actual_moisture_pct is not None else moisture))
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=actual,
            target_moisture_pct=Decimal(str(self.base_moisture_pct)),
            target_unit=ureg.metric_ton,
            allow_moisture_gain=self.allow_moisture_gain,
            handling_loss_fraction=Decimal(str(self.handling_loss_pct)) / 100,
        )

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture_pct: float | Decimal,
        *,
        source_moisture_pct: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity:
        """Physical two-moisture mass conversion ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        src_m = self._effective_actual_moisture(source_moisture_pct)
        src_dec = Decimal(str(source_moisture_pct if source_moisture_pct is not None else src_m))
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=src_dec,
            target_moisture_pct=target_moisture_pct,
            target_unit=target_unit or ureg.metric_ton,
        )

    def reverse_moisture_deduction(
        self,
        base_moisture_mass: pint.Quantity,
        *,
        actual_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts base-moisture-equivalent mass back to pre-shrinkage mass at actual moisture."""
        from core.units.registry import ureg
        if not isinstance(base_moisture_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(base_moisture_mass).__name__}")
        mt_q = self._validate_mass(base_moisture_mass)
        factor = self._moisture_adjustment_decimal(actual_moisture_pct)
        if factor <= 0.0:
            raise ValueError("Moisture adjustment factor is zero; cannot reverse deduction.")
        scale = factor if isinstance(mt_q.magnitude, Decimal) else float(factor)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def drying_charge_usd(
        self,
        mass: pint.Quantity,
        *,
        actual_moisture_pct: float | Decimal | None = None,
        charge_per_mt_per_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        """Returns drying fee on delivered tonnes, rounded to four USD decimals."""
        from core.units.precision import preserve_ieee754_scalar
        tonnes = self._validate_mass(mass)
        actual = self._effective_actual_moisture(actual_moisture_pct)
        actual_dec = Decimal(str(actual_moisture_pct if actual_moisture_pct is not None else actual))
        excess = max(Decimal(0), actual_dec - Decimal(str(self.base_moisture_pct)))
        rate = Decimal(str(self.drying_cost_per_mt_per_pct if charge_per_mt_per_pct is None else charge_per_mt_per_pct))
        if not rate.is_finite() or rate < 0:
            raise ValueError("Drying charge rate must be finite and non-negative.")
        charge = Decimal(str(tonnes.magnitude)) * excess * rate
        return preserve_ieee754_scalar(charge, tonnes.magnitude)

    def adjust_price_for_moisture(
        self,
        price_per_mt: float | Decimal,
        *,
        actual_moisture_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        """Returns base-moisture price scaled per delivered tonne at IEEE 754 precision."""
        from core.units.precision import preserve_ieee754_scalar
        price = Decimal(str(price_per_mt))
        if not price.is_finite():
            raise ValueError("price_per_mt must be finite.")
        adjusted = price * self._moisture_adjustment_decimal(actual_moisture_pct)
        return preserve_ieee754_scalar(adjusted, price_per_mt)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        if not Decimal(str(m3_q.magnitude)).is_finite() or m3_q.magnitude < 0:
            raise ValueError("Grain bulk volume must be finite and non-negative.")
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
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    @classmethod
    def corn_56lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "GrainAssay.corn_56lb", overrides={'actual_moisture_pct': actual_moisture_pct})

    @classmethod
    def wheat_60lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "GrainAssay.wheat_60lb", overrides={'actual_moisture_pct': actual_moisture_pct})

    @classmethod
    def soybeans_60lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "GrainAssay.soybeans_60lb", overrides={'actual_moisture_pct': actual_moisture_pct})

    @classmethod
    def barley_48lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "GrainAssay.barley_48lb", overrides={'actual_moisture_pct': actual_moisture_pct})

    @classmethod
    def oats_32lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "GrainAssay.oats_32lb", overrides={'actual_moisture_pct': actual_moisture_pct})


__all__ = ('BarleyQualitySpec', 'CornQualitySpec', 'GrainAssay', 'WheatQualitySpec')
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
    from core.units.registry import Quantity
    from core.units.registry import Unit
