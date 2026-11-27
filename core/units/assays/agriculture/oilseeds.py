"""Commodex assays agriculture oilseeds implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, fields
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final, Literal, Self
from typing import ClassVar
import math
import pint

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class OilseedAssay(_agricultural_families_AgriculturalFamilyAssay):
    oil_recovery_fraction: Decimal | None = None
    dry_meal_yield_fraction: Decimal | None = None
    _commodities = frozenset({"sunflower", "sesame", "flaxseed", "mustard", "groundnuts"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        for name in ("oil_recovery_fraction", "dry_meal_yield_fraction"):
            object.__setattr__(self, name, _agricultural_families_yield(getattr(self, name), name))
        if self.oil_recovery_fraction is not None and self.dry_meal_yield_fraction is not None:
            if self.composition is None:
                raise ValueError("Oil and meal mass balance requires dry-basis oil grade.")
            if self.composition.fraction("oil") * self.oil_recovery_fraction + self.dry_meal_yield_fraction > 1:
                raise ValueError("Oil plus dry meal cannot exceed dry seed feed.")

    def oil_mass(self, quantity, *, mass_basis="received"):
        from core.units.analytical_basis import _agricultural_families_forward
        contained = self.contained_mass(quantity, component="oil", mass_basis=mass_basis)
        return _agricultural_families_forward(contained, self.oil_recovery_fraction, "oil recovery")

    def dry_meal_mass(self, quantity, *, mass_basis="received"):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(self._dry_feed(quantity, mass_basis), self.dry_meal_yield_fraction, "dry meal yield")

    def meal_mass(self, quantity, *, meal_moisture_fraction, mass_basis="received"):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        fraction = _converters_petroleum_fraction(meal_moisture_fraction, "meal moisture", allow_one=False)
        return scale_mass(self.dry_meal_mass(quantity, mass_basis=mass_basis), 1 / (1 - fraction))


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 3. Soybean Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SoybeanQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for crush and export soybeans."""

    moisture_pct: Decimal
    protein_pct: Decimal
    oil_content_pct: Decimal
    basis_moisture_pct: Decimal | None = None  # None => reported on as-received moisture_pct
    damaged_beans_pct: Decimal = Decimal("0")
    heat_damaged_pct: Decimal = Decimal("0")
    foreign_material_pct: Decimal = Decimal("0")
    splits_pct: Decimal = Decimal("0")

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_pct
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        b_moist = _agricultural_quality_pct(self.basis_moisture_pct, "basis_moisture_pct", allow_high=False) if self.basis_moisture_pct is not None else moist
        prot = _agricultural_quality_pct(self.protein_pct, "protein_pct", high=Decimal(60))
        oil = _agricultural_quality_pct(self.oil_content_pct, "oil_content_pct", high=Decimal(40))
        if prot + oil + b_moist > Decimal(100):
            raise ValueError("protein_pct + oil_content_pct + basis_moisture_pct cannot exceed 100%.")
        dmg = _agricultural_quality_pct(self.damaged_beans_pct, "damaged_beans_pct")
        hdmg = _agricultural_quality_pct(self.heat_damaged_pct, "heat_damaged_pct")
        if hdmg > dmg:
            raise ValueError("heat_damaged_pct cannot exceed total damaged_beans_pct.")
        fm = _agricultural_quality_pct(self.foreign_material_pct, "foreign_material_pct")
        splits = _agricultural_quality_pct(self.splits_pct, "splits_pct")

        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "protein_pct", prot)
        object.__setattr__(self, "oil_content_pct", oil)
        object.__setattr__(self, "basis_moisture_pct", b_moist)
        object.__setattr__(self, "damaged_beans_pct", dmg)
        object.__setattr__(self, "heat_damaged_pct", hdmg)
        object.__setattr__(self, "foreign_material_pct", fm)
        object.__setattr__(self, "splits_pct", splits)

    def protein_at_moisture(self, target_moisture_pct: Decimal | float) -> Decimal:
        from core.units.analytical_basis import rebase_constituent_pct
        assert self.basis_moisture_pct is not None
        return rebase_constituent_pct(
            self.protein_pct,
            source_moisture_pct=self.basis_moisture_pct,
            target_moisture_pct=target_moisture_pct,
        )

    def oil_at_moisture(self, target_moisture_pct: Decimal | float) -> Decimal:
        from core.units.analytical_basis import rebase_constituent_pct
        assert self.basis_moisture_pct is not None
        return rebase_constituent_pct(
            self.oil_content_pct,
            source_moisture_pct=self.basis_moisture_pct,
            target_moisture_pct=target_moisture_pct,
        )

    def profat_pct(self, *, target_moisture_pct: Decimal | float | None = None) -> Decimal:
        """Combined Protein + Oil ('ProFat') percentage on ``target_moisture_pct`` (default reported basis)."""
        if target_moisture_pct is None:
            return self.protein_pct + self.oil_content_pct
        return self.protein_at_moisture(target_moisture_pct) + self.oil_at_moisture(target_moisture_pct)

    def contained_oil_mass(self, mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(mass, "mass")
        oil_as_is = self.oil_at_moisture(self.moisture_pct) / Decimal(100)
        return _agricultural_quality_kg_out(kg * oil_as_is, target_unit or mass.units, mass)

    def contained_protein_mass(self, mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(mass, "mass")
        prot_as_is = self.protein_at_moisture(self.moisture_pct) / Decimal(100)
        return _agricultural_quality_kg_out(kg * prot_as_is, target_unit or mass.units, mass)

    def normalized_mass(
        self,
        mass: pint.Quantity,
        *,
        base_moisture_pct: Decimal | float = Decimal("14.0"),
        base_foreign_material_pct: Decimal | float = Decimal("1.0"),
        allow_moisture_gain: bool = False,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        from core.units.analytical_basis import _agricultural_quality_pct
        from core.units.analytical_basis import normalize_mass_to_moisture
        from core.units.registry import ureg
        kg = _agricultural_quality_mass_kg(mass, "mass")
        fm_base = _agricultural_quality_pct(base_foreign_material_pct, "base_foreign_material_pct")
        excess_fm = max(Decimal(0), self.foreign_material_pct - fm_base) / Decimal(100)
        clean_kg = kg * (Decimal(1) - excess_fm)
        return normalize_mass_to_moisture(
            ureg.Quantity(clean_kg, ureg.kilogram),
            actual_moisture_pct=self.moisture_pct,
            base_moisture_pct=base_moisture_pct,
            allow_moisture_gain=allow_moisture_gain,
            target_unit=target_unit or mass.units,
        )

    def compare_to_benchmark(self, benchmark: SoybeanQualitySpec) -> dict[str, Decimal]:
        from core.units.analytical_basis import AgriculturalQualityBase
        if not isinstance(benchmark, SoybeanQualitySpec):
            raise TypeError(f"Expected SoybeanQualitySpec, got {type(benchmark).__name__}.")
        assert benchmark.basis_moisture_pct is not None
        deltas = AgriculturalQualityBase.compare_to_benchmark(self, benchmark)
        deltas["protein_pct"] = self.protein_at_moisture(benchmark.basis_moisture_pct) - benchmark.protein_pct
        deltas["oil_content_pct"] = self.oil_at_moisture(benchmark.basis_moisture_pct) - benchmark.oil_content_pct
        deltas["profat_pct"] = deltas["protein_pct"] + deltas["oil_content_pct"]
        deltas.pop("basis_moisture_pct", None)
        return deltas

    def blend_with(self, other: SoybeanQualitySpec, *, self_mass_fraction: Decimal | float) -> SoybeanQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        from core.units.analytical_basis import rebase_constituent_pct
        if not isinstance(other, SoybeanQualitySpec):
            raise TypeError(f"Expected SoybeanQualitySpec, got {type(other).__name__}.")
        assert self.basis_moisture_pct is not None
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_moist = w1 * self.moisture_pct + w2 * other.moisture_pct
        prot_as_is = w1 * self.protein_at_moisture(self.moisture_pct) + w2 * other.protein_at_moisture(other.moisture_pct)
        oil_as_is = w1 * self.oil_at_moisture(self.moisture_pct) + w2 * other.oil_at_moisture(other.moisture_pct)
        return SoybeanQualitySpec(
            moisture_pct=b_moist,
            protein_pct=rebase_constituent_pct(prot_as_is, source_moisture_pct=b_moist, target_moisture_pct=self.basis_moisture_pct),
            oil_content_pct=rebase_constituent_pct(oil_as_is, source_moisture_pct=b_moist, target_moisture_pct=self.basis_moisture_pct),
            basis_moisture_pct=self.basis_moisture_pct,
            damaged_beans_pct=w1 * self.damaged_beans_pct + w2 * other.damaged_beans_pct,
            heat_damaged_pct=w1 * self.heat_damaged_pct + w2 * other.heat_damaged_pct,
            foreign_material_pct=w1 * self.foreign_material_pct + w2 * other.foreign_material_pct,
            splits_pct=w1 * self.splits_pct + w2 * other.splits_pct,
        )

    def to_material_spec(self, name: str = "Soybeans") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        moist_f = self.moisture_pct / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=moist_f),
            composition=MaterialComposition(
                components={
                    "protein": self.protein_at_moisture(Decimal(0)) / Decimal(100),
                    "oil": self.oil_at_moisture(Decimal(0)) / Decimal(100),
                },
                basis=CompositionBasis.DRY,
                moisture_fraction=moist_f,
            ),
        )


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (OilseedAssay,))
__all__ = ('OilseedAssay', 'SoybeanQualitySpec', 'convert_oilseed', 'create_oilseed_context')
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
    from core.units.assays.agriculture.oilseeds import convert_oilseed
    from core.units.assays.agriculture.oilseeds import create_oilseed_context
