"""Commodex assays agriculture livestock implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, replace
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final
from typing import ClassVar
import math
import pint

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class MeatProductAssay(_agricultural_families_AgriculturalFamilyAssay):
    saleable_yield_fraction: Decimal | None = None
    _commodities = frozenset({"beef", "pork", "lamb", "offal", "processed-meat"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "saleable_yield_fraction", _agricultural_families_yield(self.saleable_yield_fraction, "saleable meat yield"))

    def saleable_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.saleable_yield_fraction, "saleable meat yield")

    def raw_meat_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.saleable_yield_fraction, "saleable meat yield")


from core.units.validation import validated_assay


# ---------------------------------------------------------------------------
# 3. Livestock & Dairy Specification (Carcass Yield, $/cwt, Milk Solids)
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class LivestockDairyAssay:
    """Immutable yield and solids specification for livestock carcasses and milk.

    Supports liveweight <-> dressed carcass weight (dressing %), $/cwt <-> $/mt,
    and milk volume/mass <-> milk solids (kgMS = fat % + protein %).
    """

    name: str
    category: str
    dressing_pct: float | None = None
    fat_pct: float | None = None
    protein_pct: float | None = None
    milk_density_kg_per_m3: float = 1030.0

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        cat = self.category.strip().lower()
        if cat not in {"cattle", "hogs", "milk"}:
            raise ValueError(f"category must be 'cattle', 'hogs', or 'milk', got {self.category!r}")

        if self.dressing_pct is not None:
            dr = _validate_finite_float(self.dressing_pct, "dressing_pct")
            if not (0.0 < dr <= 100.0):
                raise ValueError(f"dressing_pct must be in (0, 100], got {self.dressing_pct}")
            object.__setattr__(self, "dressing_pct", dr)

        if self.fat_pct is not None:
            fat = _validate_finite_float(self.fat_pct, "fat_pct")
            if not (0.0 <= fat <= 100.0):
                raise ValueError(f"fat_pct must be in [0, 100], got {self.fat_pct}")
            object.__setattr__(self, "fat_pct", fat)

        if self.protein_pct is not None:
            prot = _validate_finite_float(self.protein_pct, "protein_pct")
            if not (0.0 <= prot <= 100.0):
                raise ValueError(f"protein_pct must be in [0, 100], got {self.protein_pct}")
            object.__setattr__(self, "protein_pct", prot)

        rho = _validate_finite_float(self.milk_density_kg_per_m3, "milk_density_kg_per_m3")
        if not (900.0 <= rho <= 1200.0):
            raise ValueError(f"milk_density_kg_per_m3 must be in [900, 1200], got {self.milk_density_kg_per_m3}")

        object.__setattr__(self, "category", cat)
        object.__setattr__(self, "milk_density_kg_per_m3", rho)

    @property
    def milk_solids_pct(self) -> float:
        from core.units.precision import to_ieee754_float
        return to_ieee754_float(Decimal(str(self.fat_pct or 0.0)) + Decimal(str(self.protein_pct or 0.0)))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_livestock_dairy_context
        return create_livestock_dairy_context(self, name=context_name)

    def live_to_carcass_mass(self, live_mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(live_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(live_mass).__name__}")
        if self.dressing_pct is None:
            raise ValueError(f"{self.name!r} does not define dressing_pct.")
        kg_q = live_mass.to(ureg.kilogram)
        frac = self.dressing_pct / 100.0
        scale = _scale_factor(kg_q.magnitude, frac)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def carcass_to_live_mass(self, carcass_mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(carcass_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(carcass_mass).__name__}")
        if self.dressing_pct is None:
            raise ValueError(f"{self.name!r} does not define dressing_pct.")
        kg_q = carcass_mass.to(ureg.kilogram)
        frac = self.dressing_pct / 100.0
        scale = _scale_factor(kg_q.magnitude, frac)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def milk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        kg_per_m3 = self.milk_density_kg_per_m3
        scale = _scale_factor(m3_q.magnitude, kg_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.kilogram)

    def milk_mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        kg_q = mass.to(ureg.kilogram)
        kg_per_m3 = self.milk_density_kg_per_m3
        scale = _scale_factor(kg_q.magnitude, kg_per_m3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.meter ** 3)

    def milk_to_solids_kg(self, milk_qty: pint.Quantity) -> pint.Quantity:
        """Converts milk mass or volume to Milk Solids (kgMS = fat + protein)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(milk_qty, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(milk_qty).__name__}")
        if self.milk_solids_pct <= 0.0:
            raise ValueError(f"{self.name!r} does not define positive milk solids (fat_pct / protein_pct).")
        if milk_qty.dimensionality == (ureg.meter ** 3).dimensionality:
            kg_q = self.milk_volume_to_mass(milk_qty)
        else:
            kg_q = milk_qty.to(ureg.kilogram)
        frac = self.milk_solids_pct / 100.0
        scale = _scale_factor(kg_q.magnitude, frac)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    @staticmethod
    def usd_per_cwt_to_usd_per_mt(usd_per_cwt: float | Decimal, *, metric_cwt: bool = False) -> float | Decimal:
        """Converts $/cwt (US 100 lb cwt by default, or 100 kg metric cwt) to $/metric_ton."""
        from core.units.agricultural_specs import _validate_finite_float
        from core.units.precision import preserve_ieee754_scalar
        _validate_finite_float(usd_per_cwt, "usd_per_cwt")
        p_dec = Decimal(str(usd_per_cwt))
        if metric_cwt:
            return preserve_ieee754_scalar(p_dec * Decimal(10), usd_per_cwt)
        return preserve_ieee754_scalar(p_dec * (Decimal(1000) / Decimal("0.45359237")) / Decimal(100), usd_per_cwt)

    @staticmethod
    def usd_per_mt_to_usd_per_cwt(usd_per_mt: float | Decimal, *, metric_cwt: bool = False) -> float | Decimal:
        """Converts $/metric_ton to $/cwt."""
        from core.units.agricultural_specs import _validate_finite_float
        from core.units.precision import preserve_ieee754_scalar
        _validate_finite_float(usd_per_mt, "usd_per_mt")
        p_dec = Decimal(str(usd_per_mt))
        if metric_cwt:
            return preserve_ieee754_scalar(p_dec / Decimal(10), usd_per_mt)
        return preserve_ieee754_scalar(p_dec * Decimal(100) / (Decimal(1000) / Decimal("0.45359237")), usd_per_mt)

    @classmethod
    def live_cattle_standard(cls) -> LivestockDairyAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "LivestockDairyAssay.live_cattle_standard")

    @classmethod
    def live_hogs_standard(cls) -> LivestockDairyAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "LivestockDairyAssay.live_hogs_standard")

    @classmethod
    def whole_milk_standard(cls) -> LivestockDairyAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "LivestockDairyAssay.whole_milk_standard")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (MeatProductAssay,))
__all__ = ('LivestockDairyAssay', 'MeatProductAssay', 'convert_meat_product', 'create_meat_product_context')
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
    from core.units.assays.agriculture.livestock import convert_meat_product
    from core.units.assays.agriculture.livestock import create_meat_product_context
