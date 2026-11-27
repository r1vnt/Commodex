"""Commodex assays agriculture dairy implementation and public interfaces."""
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
class DairyProductAssay(_agricultural_families_AgriculturalFamilyAssay):
    _commodities = frozenset({"smp", "wmp", "butter", "cheese", "whey", "lactose"})

    def milk_solids_mass(self, quantity):
        return self.wet_to_dry_mass(quantity)

    def milk_fat_mass(self, quantity, *, mass_basis="received"):
        return self.contained_mass(quantity, component="fat", mass_basis=mass_basis)

    def protein_mass(self, quantity, *, mass_basis="received"):
        return self.contained_mass(quantity, component="protein", mass_basis=mass_basis)

    def fat_equivalent_mass(self, quantity, *, target_fat_fraction, fat_recovery_fraction, mass_basis="received"):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        target = _converters_petroleum_fraction(target_fat_fraction, "target product fat")
        recovery = _converters_petroleum_fraction(fat_recovery_fraction, "fat recovery")
        if target == 0:
            raise ValueError("Target fat fraction must be positive.")
        return scale_mass(self.milk_fat_mass(quantity, mass_basis=mass_basis), recovery / target)


from core.units.analytical_basis import AgriculturalQualityBase
from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# 11. Milk & Dairy Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MilkQualitySpec(AgriculturalQualityBase):
    """Physical composition and standardisation specification for raw milk, skim, and cream."""

    butterfat_pct: Decimal
    protein_pct: Decimal
    lactose_pct: Decimal
    minerals_ash_pct: Decimal = Decimal("0.72")
    density_kg_per_m3: Decimal = Decimal("1030.0")

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_pct
        fat = _agricultural_quality_pct(self.butterfat_pct, "butterfat_pct", high=Decimal(60))
        prot = _agricultural_quality_pct(self.protein_pct, "protein_pct", high=Decimal(40))
        lac = _agricultural_quality_pct(self.lactose_pct, "lactose_pct", high=Decimal(40))
        ash = _agricultural_quality_pct(self.minerals_ash_pct, "minerals_ash_pct", high=Decimal(10))
        if fat + prot + lac + ash >= Decimal(100):
            raise ValueError("Total milk solids (butterfat + protein + lactose + minerals) must be < 100%.")
        rho = _agricultural_quality_pct(self.density_kg_per_m3, "density_kg_per_m3", low=Decimal(900), high=Decimal(1200))

        object.__setattr__(self, "butterfat_pct", fat)
        object.__setattr__(self, "protein_pct", prot)
        object.__setattr__(self, "lactose_pct", lac)
        object.__setattr__(self, "minerals_ash_pct", ash)
        object.__setattr__(self, "density_kg_per_m3", rho)

    @property
    def solids_not_fat_pct(self) -> Decimal:
        """Solids-Not-Fat (SNF % = protein + lactose + minerals)."""
        return self.protein_pct + self.lactose_pct + self.minerals_ash_pct

    @property
    def milk_solids_pct(self) -> Decimal:
        """Commercial dairy Milk Solids percentage (kgMS basis = butterfat % + protein %)."""
        return self.butterfat_pct + self.protein_pct

    @property
    def total_solids_pct(self) -> Decimal:
        """Total milk solids percentage (butterfat + SNF)."""
        return self.butterfat_pct + self.solids_not_fat_pct

    @property
    def moisture_pct(self) -> Decimal:
        return Decimal(100) - self.total_solids_pct

    @authorized_conversion
    def _milk_mass_kg(self, milk_qty: pint.Quantity) -> Decimal:
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        if not isinstance(milk_qty, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(milk_qty).__name__}.")
        validate_registry_quantity(milk_qty)
        if milk_qty.dimensionality == (ureg.meter ** 3).dimensionality:
            m3 = _agricultural_quality_nonneg(milk_qty.to("meter ** 3").magnitude, "milk_qty")
            return m3 * self.density_kg_per_m3
        return _agricultural_quality_mass_kg(milk_qty, "milk_qty")

    def component_masses(self, milk_qty: pint.Quantity, *, target_unit: str | Unit = "kilogram") -> dict[str, pint.Quantity]:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        kg = self._milk_mass_kg(milk_qty)
        return {
            "butterfat": _agricultural_quality_kg_out(kg * self.butterfat_pct / Decimal(100), target_unit, milk_qty),
            "protein": _agricultural_quality_kg_out(kg * self.protein_pct / Decimal(100), target_unit, milk_qty),
            "lactose": _agricultural_quality_kg_out(kg * self.lactose_pct / Decimal(100), target_unit, milk_qty),
            "solids_not_fat": _agricultural_quality_kg_out(kg * self.solids_not_fat_pct / Decimal(100), target_unit, milk_qty),
            "milk_solids": _agricultural_quality_kg_out(kg * self.milk_solids_pct / Decimal(100), target_unit, milk_qty),
            "total_solids": _agricultural_quality_kg_out(kg * self.total_solids_pct / Decimal(100), target_unit, milk_qty),
        }

    def energy_corrected_milk(self, milk_qty: pint.Quantity, *, target_unit: str | Unit = "kilogram") -> pint.Quantity:
        """IDF standard Energy-Corrected Milk (ECM, standardised to 4.0% fat, 3.3% protein)."""
        from core.units.analytical_basis import _agricultural_quality_kg_out
        kg = self._milk_mass_kg(milk_qty)
        factor = (
            Decimal("0.383") * self.butterfat_pct
            + Decimal("0.242") * self.protein_pct
            + Decimal("0.7832")
        ) / Decimal("3.140")
        return _agricultural_quality_kg_out(kg * factor, target_unit, milk_qty)

    def fat_corrected_milk_4pct(self, milk_qty: pint.Quantity, *, target_unit: str | Unit = "kilogram") -> pint.Quantity:
        """Gaines 4% Fat-Corrected Milk: ``0.4 * M_milk + 15 * M_fat``."""
        from core.units.analytical_basis import _agricultural_quality_kg_out
        kg = self._milk_mass_kg(milk_qty)
        fcm_kg = Decimal("0.4") * kg + Decimal("15") * (kg * self.butterfat_pct / Decimal(100))
        return _agricultural_quality_kg_out(fcm_kg, target_unit, milk_qty)

    def pearson_square_fat_weights(
        self,
        other: MilkQualitySpec,
        *,
        target_butterfat_pct: Decimal | float,
    ) -> tuple[Decimal, Decimal]:
        """Return ``(self_mass_fraction, other_mass_fraction)`` to standardise butterfat to ``target_butterfat_pct``."""
        from core.units.analytical_basis import _agricultural_quality_pct
        if not isinstance(other, MilkQualitySpec):
            raise TypeError(f"Expected MilkQualitySpec, got {type(other).__name__}.")
        tgt = _agricultural_quality_pct(target_butterfat_pct, "target_butterfat_pct", high=Decimal(60))
        f1, f2 = self.butterfat_pct, other.butterfat_pct
        if f1 == f2:
            if tgt == f1:
                return Decimal("0.5"), Decimal("0.5")
            raise ValueError("Cannot shift butterfat percentage using two streams of identical fat content.")
        if not (min(f1, f2) <= tgt <= max(f1, f2)):
            raise ValueError(f"target_butterfat_pct ({tgt}%) must lie between {min(f1, f2)}% and {max(f1, f2)}%.")
        w1 = (tgt - f2) / (f1 - f2)
        return w1, Decimal(1) - w1

    def blend_with(self, other: MilkQualitySpec, *, self_mass_fraction: Decimal | float) -> MilkQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, MilkQualitySpec):
            raise TypeError(f"Expected MilkQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_rho = Decimal(1) / ((w1 / self.density_kg_per_m3) + (w2 / other.density_kg_per_m3))
        return MilkQualitySpec(
            butterfat_pct=w1 * self.butterfat_pct + w2 * other.butterfat_pct,
            protein_pct=w1 * self.protein_pct + w2 * other.protein_pct,
            lactose_pct=w1 * self.lactose_pct + w2 * other.lactose_pct,
            minerals_ash_pct=w1 * self.minerals_ash_pct + w2 * other.minerals_ash_pct,
            density_kg_per_m3=b_rho,
        )

    def to_material_spec(self, name: str = "Milk") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        from core.units.registry import ureg
        moist_f = self.moisture_pct / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                density=ureg.Quantity(self.density_kg_per_m3, "kg / meter ** 3"),
                moisture=moist_f,
            ),
            composition=MaterialComposition(
                components={
                    "butterfat": self.butterfat_pct / Decimal(100),
                    "protein": self.protein_pct / Decimal(100),
                    "lactose": self.lactose_pct / Decimal(100),
                },
                basis=CompositionBasis.AS_RECEIVED,
                moisture_fraction=moist_f,
                ash_fraction=self.minerals_ash_pct / Decimal(100),
            ),
        )


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (DairyProductAssay,))
__all__ = ('DairyProductAssay', 'MilkQualitySpec', 'convert_dairy_product', 'create_dairy_product_context')
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
    from core.units.assays.agriculture.dairy import convert_dairy_product
    from core.units.assays.agriculture.dairy import create_dairy_product_context
