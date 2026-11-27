"""Commodex assays agriculture rice implementation and public interfaces."""
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
# 4. Rice Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class RiceQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for rough/paddy and milled rice.

    ``head_rice_basis`` explicitly distinguishes USDA rough-rice convention
    (``'of_paddy'``, e.g. 55/70 milling yield = 55% whole kernels of paddy,
    70% total milled of paddy) from milled-export convention (``'of_milled'``,
    e.g. Thai 5% Broken = 95% head rice of milled rice).
    """

    milling_recovery_pct: Decimal = Decimal("70.0")
    head_rice_yield_pct: Decimal | None = None
    broken_rice_pct: Decimal | None = None  # Share of milled rice that is broken
    head_rice_basis: Literal["of_paddy", "of_milled"] = "of_paddy"
    moisture_pct: Decimal = Decimal("14.0")
    chalky_kernels_pct: Decimal = Decimal("0")

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_pct
        mill = _agricultural_quality_pct(self.milling_recovery_pct, "milling_recovery_pct", allow_low=False)
        moist = _agricultural_quality_pct(self.moisture_pct, "moisture_pct", allow_high=False)
        chalk = _agricultural_quality_pct(self.chalky_kernels_pct, "chalky_kernels_pct")
        if self.head_rice_basis not in {"of_paddy", "of_milled"}:
            raise ValueError("head_rice_basis must be 'of_paddy' or 'of_milled'.")

        hr_in = _agricultural_quality_pct(self.head_rice_yield_pct, "head_rice_yield_pct") if self.head_rice_yield_pct is not None else None
        br_in = _agricultural_quality_pct(self.broken_rice_pct, "broken_rice_pct") if self.broken_rice_pct is not None else None
        if hr_in is None and br_in is None:
            raise ValueError("Supply at least one of head_rice_yield_pct or broken_rice_pct.")

        if hr_in is not None:
            if self.head_rice_basis == "of_paddy":
                if hr_in > mill:
                    raise ValueError(f"head_rice_yield_pct of paddy ({hr_in}) cannot exceed milling_recovery_pct ({mill}).")
                derived_broken = (Decimal(1) - (hr_in / mill)) * Decimal(100)
            else:
                derived_broken = Decimal(100) - hr_in
            if br_in is not None and abs(br_in - derived_broken) > Decimal("1e-6"):
                raise ValueError(
                    f"broken_rice_pct ({br_in}%) conflicts with head_rice_yield_pct={hr_in}% ({self.head_rice_basis})."
                )
            br_final = derived_broken if br_in is None else br_in
            hr_final = hr_in
        else:
            assert br_in is not None
            br_final = br_in
            hr_of_milled = Decimal(100) - br_final
            hr_final = (mill * hr_of_milled / Decimal(100)) if self.head_rice_basis == "of_paddy" else hr_of_milled

        object.__setattr__(self, "milling_recovery_pct", mill)
        object.__setattr__(self, "head_rice_yield_pct", hr_final)
        object.__setattr__(self, "broken_rice_pct", br_final)
        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "chalky_kernels_pct", chalk)

    @property
    def head_rice_pct_of_paddy(self) -> Decimal:
        assert self.head_rice_yield_pct is not None
        if self.head_rice_basis == "of_paddy":
            return self.head_rice_yield_pct
        return self.milling_recovery_pct * self.head_rice_yield_pct / Decimal(100)

    @property
    def head_rice_pct_of_milled(self) -> Decimal:
        assert self.broken_rice_pct is not None
        return Decimal(100) - self.broken_rice_pct

    def milling_breakdown(self, paddy_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> dict[str, pint.Quantity]:
        """Mass-conserving paddy -> head rice + broken rice + hull/bran loss breakdown."""
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(paddy_mass, "paddy_mass")
        dst = target_unit or paddy_mass.units
        milled_kg = kg * self.milling_recovery_pct / Decimal(100)
        head_kg = kg * self.head_rice_pct_of_paddy / Decimal(100)
        broken_kg = milled_kg - head_kg
        loss_kg = kg - milled_kg
        return {
            "milled_total": _agricultural_quality_kg_out(milled_kg, dst, paddy_mass),
            "head_rice": _agricultural_quality_kg_out(head_kg, dst, paddy_mass),
            "broken_rice": _agricultural_quality_kg_out(broken_kg, dst, paddy_mass),
            "hull_and_bran_loss": _agricultural_quality_kg_out(loss_kg, dst, paddy_mass),
        }

    def required_paddy_for_head_rice(self, head_rice_mass: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg_hr = _agricultural_quality_mass_kg(head_rice_mass, "head_rice_mass")
        hr_frac = self.head_rice_pct_of_paddy / Decimal(100)
        if hr_frac <= Decimal(0):
            raise ValueError("head_rice_pct_of_paddy is zero; required paddy mass is undefined.")
        return _agricultural_quality_kg_out(kg_hr / hr_frac, target_unit or head_rice_mass.units, head_rice_mass)

    def blend_with(self, other: RiceQualitySpec, *, self_mass_fraction: Decimal | float) -> RiceQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, RiceQualitySpec):
            raise TypeError(f"Expected RiceQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        b_mill = w1 * self.milling_recovery_pct + w2 * other.milling_recovery_pct
        b_hr_paddy = w1 * self.head_rice_pct_of_paddy + w2 * other.head_rice_pct_of_paddy
        b_hr = b_hr_paddy if self.head_rice_basis == "of_paddy" else (b_hr_paddy / b_mill * Decimal(100))
        return RiceQualitySpec(
            milling_recovery_pct=b_mill,
            head_rice_yield_pct=b_hr,
            head_rice_basis=self.head_rice_basis,
            moisture_pct=w1 * self.moisture_pct + w2 * other.moisture_pct,
            chalky_kernels_pct=w1 * self.chalky_kernels_pct + w2 * other.chalky_kernels_pct,
        )

    def to_material_spec(self, name: str = "Rice") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import MaterialProperties
        from core.units.processing import ProcessingSpecification
        assert self.broken_rice_pct is not None
        mill_f = self.milling_recovery_pct / Decimal(100)
        hr_f = self.head_rice_pct_of_paddy / Decimal(100)
        br_f = mill_f - hr_f
        loss_f = Decimal(1) - mill_f
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=self.moisture_pct / Decimal(100)),
            processing=ProcessingSpecification(
                input_material="paddy_rice",
                output_yields={"head_rice": hr_f, "broken_rice": br_f},
                process_losses={"hull_and_bran": loss_f} if loss_f > Decimal(0) else {},
            ),
        )


from core.units.agricultural_specs import _RICE_BUSHEL_LB
from core.units.agricultural_specs import _RICE_EXPORT_BAG_KG
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class RiceAssay:
    """Immutable milling, bag, and pricing specification for rough and milled rice.

    Bridges rough/paddy rice with brown and milled white rice yields, head-rice
    vs. broken fractions, 50 kg export bags, US rough-rice bushel/cwt parity,
    and optional bulk density.
    """

    name: str
    variety: str  # e.g. "LONG_GRAIN", "JASMINE", "BASMATI", "MEDIUM_GRAIN"
    milling_yield_pct: float = 70.0       # Paddy -> total milled rice recovery (%)
    head_rice_pct: float = 80.0         # Share of milled rice that is whole/head rice (%)
    brown_rice_yield_pct: float = 78.0  # Paddy -> dehulled brown/cargo rice recovery (%)
    bushel_weight_lb: float = _RICE_BUSHEL_LB
    bag_mass_kg: float = _RICE_EXPORT_BAG_KG
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.variety, str):
            raise ValueError("variety must be a non-empty string.")
        var = self.variety.strip().upper()
        if not var:
            raise ValueError("variety must be a non-empty string.")

        mill_y = _validate_finite_float(self.milling_yield_pct, "milling_yield_pct")
        if not (0.0 < mill_y <= 100.0):
            raise ValueError(f"milling_yield_pct must be in (0, 100], got {self.milling_yield_pct}")

        head_y = _validate_finite_float(self.head_rice_pct, "head_rice_pct")
        if not (0.0 <= head_y <= 100.0):
            raise ValueError(f"head_rice_pct must be in [0, 100], got {self.head_rice_pct}")

        brown_y = _validate_finite_float(self.brown_rice_yield_pct, "brown_rice_yield_pct")
        if not (0.0 < brown_y <= 100.0) or brown_y < mill_y:
            raise ValueError(
                f"brown_rice_yield_pct ({brown_y}) must be in (0, 100] and >= milling_yield_pct ({mill_y})."
            )

        bu_lb = _validate_finite_float(self.bushel_weight_lb, "bushel_weight_lb")
        if bu_lb <= 0.0:
            raise ValueError(f"bushel_weight_lb must be strictly positive, got {self.bushel_weight_lb}")

        bag_kg = _validate_finite_float(self.bag_mass_kg, "bag_mass_kg")
        if bag_kg <= 0.0:
            raise ValueError(f"bag_mass_kg must be strictly positive, got {self.bag_mass_kg}")

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        object.__setattr__(self, "variety", var)
        object.__setattr__(self, "milling_yield_pct", mill_y)
        object.__setattr__(self, "head_rice_pct", head_y)
        object.__setattr__(self, "brown_rice_yield_pct", brown_y)
        object.__setattr__(self, "bushel_weight_lb", bu_lb)
        object.__setattr__(self, "bag_mass_kg", bag_kg)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)

    @property
    def milling_fraction(self) -> float:
        return self.milling_yield_pct / 100.0

    @property
    def head_rice_fraction_of_paddy(self) -> float:
        """Whole head-rice output as a fraction of original paddy mass."""
        return (self.milling_yield_pct / 100.0) * (self.head_rice_pct / 100.0)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_rice_context
        return create_rice_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        number = Decimal(str(quantity.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Rice mass must be finite and non-negative.")
        return quantity

    @staticmethod
    def _percentage_scale(magnitude: Any, percentage: float) -> float | Decimal:
        fraction = Decimal(str(percentage)) / 100
        return fraction if isinstance(magnitude, Decimal) else float(fraction)

    # --- Paddy <-> Brown <-> Milled Conversions ---
    def paddy_to_brown_mass(self, paddy_mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(paddy_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(paddy_mass).__name__}")
        mt_q = self._validate_mass(paddy_mass)
        scale = self._percentage_scale(mt_q.magnitude, self.brown_rice_yield_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def paddy_to_milled_mass(self, paddy_mass: pint.Quantity) -> pint.Quantity:
        """Converts rough/paddy rice mass to total milled rice mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(paddy_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(paddy_mass).__name__}")
        mt_q = self._validate_mass(paddy_mass)
        scale = self._percentage_scale(mt_q.magnitude, self.milling_yield_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def milled_to_paddy_mass(self, milled_mass: pint.Quantity) -> pint.Quantity:
        """Converts milled rice mass back to required rough/paddy rice mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(milled_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(milled_mass).__name__}")
        mt_q = self._validate_mass(milled_mass)
        scale = self._percentage_scale(mt_q.magnitude, self.milling_yield_pct)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def brown_to_paddy_mass(self, brown_mass: pint.Quantity) -> pint.Quantity:
        """Returns required paddy mass from dehulled brown/cargo rice mass."""
        from core.units.registry import ureg
        quantity = self._validate_mass(brown_mass)
        scale = self._percentage_scale(quantity.magnitude, self.brown_rice_yield_pct)
        return ureg.Quantity(quantity.magnitude / scale, ureg.metric_ton)

    def brown_to_milled_mass(self, brown_mass: pint.Quantity) -> pint.Quantity:
        """Applies the relative milling recovery from brown rice, not paddy yield twice."""
        return self.paddy_to_milled_mass(self.brown_to_paddy_mass(brown_mass))

    def milled_to_brown_mass(self, milled_mass: pint.Quantity) -> pint.Quantity:
        return self.paddy_to_brown_mass(self.milled_to_paddy_mass(milled_mass))

    def milled_to_head_and_broken(self, milled_mass: pint.Quantity) -> dict[str, pint.Quantity]:
        """Splits total milled mass; head_rice_pct is a share of milled output."""
        from core.units.registry import ureg
        quantity = self._validate_mass(milled_mass)
        scale = self._percentage_scale(quantity.magnitude, self.head_rice_pct)
        head = quantity.magnitude * scale
        return {
            "head_rice": ureg.Quantity(head, ureg.metric_ton),
            "broken_rice": ureg.Quantity(quantity.magnitude - head, ureg.metric_ton),
        }

    # --- US Bushels, Bags & Bulk Volume ---
    def mass_to_bushels(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded rough-rice bushels by the assay's specified weight."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        pounds = self._validate_mass(mass).to(ureg.pound)
        return pounds.magnitude / _scale_factor(pounds.magnitude, self.bushel_weight_lb)

    def bushels_to_mass(self, bushels: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        count = Decimal(str(bushels))
        if not count.is_finite() or count < 0:
            raise ValueError("bushels must be finite and non-negative.")
        pounds = count * Decimal(str(self.bushel_weight_lb))
        magnitude = pounds if isinstance(bushels, Decimal) else float(pounds)
        return ureg.Quantity(magnitude, ureg.pound).to(ureg.metric_ton)

    def mass_to_bags(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal export-bag count, preserving Decimal input."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bag_mass_kg)

    def bags_to_mass(self, bag_count: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        count = Decimal(str(bag_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bag_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bag_mass_kg))
        magnitude = kilograms if isinstance(bag_count, Decimal) else float(kilograms)
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
            raise ValueError("Rice volume must be finite and non-negative.")
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
            raise ValueError("Rice price must be finite.")
        return value

    @staticmethod
    def _price_result(value: Decimal, original: float | Decimal) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        return preserve_ieee754_scalar(value, original)

    @staticmethod
    def usd_per_cwt_to_usd_per_mt(usd_per_cwt: float | Decimal) -> float | Decimal:
        """Converts US 100-lb hundredweight price to USD/metric tonne."""
        from core.units.registry import ureg
        pounds_per_tonne = ureg.Quantity(Decimal("1"), ureg.metric_ton).to(ureg.pound).magnitude
        price = RiceAssay._validated_price(usd_per_cwt) * pounds_per_tonne / 100
        return RiceAssay._price_result(price, usd_per_cwt)

    @staticmethod
    def usd_per_mt_to_usd_per_cwt(usd_per_mt: float | Decimal) -> float | Decimal:
        """Converts USD/metric tonne to US 100-lb hundredweight price."""
        from core.units.registry import ureg
        pounds_per_tonne = ureg.Quantity(Decimal("1"), ureg.metric_ton).to(ureg.pound).magnitude
        price = RiceAssay._validated_price(usd_per_mt) * 100 / pounds_per_tonne
        return RiceAssay._price_result(price, usd_per_mt)

    def paddy_price_to_milled_parity_per_mt(self, paddy_usd_per_mt: float | Decimal) -> float | Decimal:
        """Returns milled price covering paddy cost only; excludes milling costs/credits."""
        price = self._validated_price(paddy_usd_per_mt) / (Decimal(str(self.milling_yield_pct)) / 100)
        return self._price_result(price, paddy_usd_per_mt)

    # --- Standard Rice Presets ---
    @classmethod
    def thai_white_5pct_broken(cls) -> RiceAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "RiceAssay.thai_white_5pct_broken")

    @classmethod
    def us_cme_rough_rice(cls) -> RiceAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "RiceAssay.us_cme_rough_rice")


__all__ = ('RiceAssay', 'RiceQualitySpec')
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
