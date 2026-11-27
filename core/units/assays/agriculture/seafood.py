"""Commodex assays agriculture seafood implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final
import math
import pint

from core.units.validation import validated_assay


# ---------------------------------------------------------------------------
# 5. Seafood Yields, Glaze & Marine Reduction Specification
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class SeafoodAssay:
    """Specifications for finfish, crustaceans, and marine reduction bulk.

    Dressed and edible yields are independent percentages of whole input mass;
    edible yield is not applied sequentially to dressed output. Glaze percentage
    is ice as a share of gross glazed mass. Presets are illustrative recovery
    assumptions; use the actual supplier/contract specification for settlement.
    """

    name: str
    species_group: str
    dressed_yield_pct: float = 100.0
    edible_yield_pct: float = 100.0
    oil_yield_pct: float = 0.0
    glaze_pct: float = 0.0
    count_per_lb: float | None = None
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.species_group, str):
            raise ValueError("species_group must be a string.")
        grp = self.species_group.strip().lower()
        if grp not in {"finfish", "lobster", "crab", "shrimp", "reduction"}:
            raise ValueError(
                f"species_group must be 'finfish', 'lobster', 'crab', 'shrimp', or 'reduction', got {self.species_group!r}"
            )
        dressed = _validate_finite_float(self.dressed_yield_pct, "dressed_yield_pct")
        edible = _validate_finite_float(self.edible_yield_pct, "edible_yield_pct")
        oil = _validate_finite_float(self.oil_yield_pct, "oil_yield_pct")
        glaze = _validate_finite_float(self.glaze_pct, "glaze_pct")

        if not (0.0 < dressed <= 100.0):
            raise ValueError(f"dressed_yield_pct must be in (0, 100], got {self.dressed_yield_pct}")
        if not (0.0 < edible <= 100.0):
            raise ValueError(f"edible_yield_pct must be in (0, 100], got {self.edible_yield_pct}")
        if not (0.0 <= oil <= 100.0):
            raise ValueError(f"oil_yield_pct must be in [0, 100], got {self.oil_yield_pct}")
        if not (0.0 <= glaze < 100.0):
            raise ValueError(f"glaze_pct must be in [0, 100), got {self.glaze_pct}")

        if grp == "reduction" and dressed + oil > 100.0 + 1e-9:
            raise ValueError("Reduction meal and oil yields cannot sum to more than 100%.")

        cnt: float | None = None
        if self.count_per_lb is not None:
            cnt = _validate_finite_float(self.count_per_lb, "count_per_lb")
            if cnt <= 0.0:
                raise ValueError(f"count_per_lb must be strictly positive, got {self.count_per_lb}")

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        object.__setattr__(self, "species_group", grp)
        object.__setattr__(self, "dressed_yield_pct", dressed)
        object.__setattr__(self, "edible_yield_pct", edible)
        object.__setattr__(self, "oil_yield_pct", oil)
        object.__setattr__(self, "glaze_pct", glaze)
        object.__setattr__(self, "count_per_lb", cnt)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)

    @property
    def glaze_retention_fraction(self) -> float:
        """Net deglazed mass divided by gross glazed mass."""
        return 1.0 - (self.glaze_pct / 100.0)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Bridges bulk volume and gross product mass if density is supplied."""
        from core.units.agricultural_specs import create_seafood_context
        return create_seafood_context(self, name=context_name)

    def _validate_mass(self, mass: pint.Quantity, unit: pint.Unit) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(unit)
        magnitude = Decimal(str(quantity.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Seafood mass must be finite and non-negative.")
        return quantity

    @staticmethod
    def _yield_scale(magnitude: Any, percentage: float) -> float | Decimal:
        fraction = Decimal(str(percentage)) / 100
        return fraction if isinstance(magnitude, Decimal) else float(fraction)

    def _glaze_scale(self, magnitude: Any) -> float | Decimal:
        fraction = 1 - Decimal(str(self.glaze_pct)) / 100
        return fraction if isinstance(magnitude, Decimal) else float(fraction)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        cubic_metres = volume.to(ureg.meter ** 3)
        magnitude = Decimal(str(cubic_metres.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Seafood bulk volume must be finite and non-negative.")
        scale = _scale_factor(cubic_metres.magnitude, self.bulk_density_kg_per_m3)
        return ureg.Quantity(cubic_metres.magnitude * scale, ureg.kilogram)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        kilograms = self._validate_mass(mass, ureg.kilogram)
        scale = _scale_factor(kilograms.magnitude, self.bulk_density_kg_per_m3)
        return ureg.Quantity(kilograms.magnitude / scale, ureg.meter ** 3)

    def whole_to_dressed(self, whole_mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        q = self._validate_mass(whole_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.dressed_yield_pct)
        return ureg.Quantity(q.magnitude * scale, ureg.kilogram)

    def dressed_to_whole(self, dressed_mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(dressed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(dressed_mass).__name__}")
        q = self._validate_mass(dressed_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.dressed_yield_pct)
        return ureg.Quantity(q.magnitude / scale, ureg.kilogram)

    def whole_to_edible(self, whole_mass: pint.Quantity) -> pint.Quantity:
        """Converts whole live/round mass to fillet, picked meat, or P&D shrimp mass."""
        from core.units.registry import ureg
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        q = self._validate_mass(whole_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.edible_yield_pct)
        return ureg.Quantity(q.magnitude * scale, ureg.kilogram)

    def edible_to_whole(self, edible_mass: pint.Quantity) -> pint.Quantity:
        """Converts edible/fillet/picked mass back to whole live/round mass."""
        from core.units.registry import ureg
        if not isinstance(edible_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(edible_mass).__name__}")
        q = self._validate_mass(edible_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.edible_yield_pct)
        return ureg.Quantity(q.magnitude / scale, ureg.kilogram)

    def net_to_glazed_mass(self, net_mass: pint.Quantity) -> pint.Quantity:
        """Converts net deglazed/drained mass to gross glazed frozen mass."""
        from core.units.registry import ureg
        if not isinstance(net_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(net_mass).__name__}")
        q = self._validate_mass(net_mass, ureg.kilogram)
        scale = self._glaze_scale(q.magnitude)
        return ureg.Quantity(q.magnitude / scale, ureg.kilogram)

    def glazed_to_net_mass(self, glazed_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross glazed frozen mass to net deglazed/drained mass."""
        from core.units.registry import ureg
        if not isinstance(glazed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(glazed_mass).__name__}")
        q = self._validate_mass(glazed_mass, ureg.kilogram)
        scale = self._glaze_scale(q.magnitude)
        return ureg.Quantity(q.magnitude * scale, ureg.kilogram)

    def whole_to_reduction_yields(self, whole_mass: pint.Quantity) -> dict[str, pint.Quantity]:
        """Returns meal (dressed_yield_pct) and marine oil (oil_yield_pct) in metric tonnes."""
        from core.units.exceptions import IncompatibleCommodityError
        from core.units.registry import ureg
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        if self.species_group != "reduction":
            raise IncompatibleCommodityError("whole_to_reduction_yields requires species_group='reduction'.")
        mt_q = self._validate_mass(whole_mass, ureg.metric_ton)
        meal_scale = self._yield_scale(mt_q.magnitude, self.dressed_yield_pct)
        oil_scale = self._yield_scale(mt_q.magnitude, self.oil_yield_pct)
        return {
            "meal": ureg.Quantity(mt_q.magnitude * meal_scale, ureg.metric_ton),
            "oil": ureg.Quantity(mt_q.magnitude * oil_scale, ureg.metric_ton),
        }

    def piece_mass_grams(self) -> float:
        """Nominal single-piece weight in grams from `count_per_lb`."""
        from core.units.precision import to_ieee754_float
        if self.count_per_lb is None:
            raise ValueError(f"{self.name!r} does not define count_per_lb.")
        return to_ieee754_float(453.59237 / self.count_per_lb)

    def whole_price_to_edible_price(self, whole_price: float | Decimal) -> float | Decimal:
        """Converts $/unit whole price to implied $/unit edible/fillet/meat parity."""
        from core.units.precision import preserve_ieee754_scalar
        px = Decimal(str(whole_price))
        if not px.is_finite():
            raise ValueError("whole_price must be finite.")
        result = px / (Decimal(str(self.edible_yield_pct)) / 100)
        return preserve_ieee754_scalar(result, whole_price)

    @classmethod
    def atlantic_salmon(cls, *, glaze_pct: float = 0.0) -> SeafoodAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SeafoodAssay.atlantic_salmon", overrides={'glaze_pct': glaze_pct})

    @classmethod
    def american_lobster(cls, *, glaze_pct: float = 0.0) -> SeafoodAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SeafoodAssay.american_lobster", overrides={'glaze_pct': glaze_pct})

    @classmethod
    def snow_crab(cls, *, glaze_pct: float = 8.0) -> SeafoodAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SeafoodAssay.snow_crab", overrides={'glaze_pct': glaze_pct})

    @classmethod
    def vannamei_shrimp(cls, *, count_per_lb: float = 21.0, glaze_pct: float = 10.0) -> SeafoodAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SeafoodAssay.vannamei_shrimp", overrides={'name': f'Vannamei Shrimp ({count_per_lb:g} ct/lb, HLSO 65% / P&D 50%)', 'glaze_pct': glaze_pct, 'count_per_lb': count_per_lb})

    @classmethod
    def anchoveta_reduction(cls) -> SeafoodAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "SeafoodAssay.anchoveta_reduction")


from core.units.agricultural_specs import _BULK_FISHMEAL_KG_PER_M3
from core.units.agricultural_specs import _CRUDE_MARINE_OIL_KG_PER_M3
from core.units.agricultural_specs import _FISHMEAL_BAG_KG
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class FishmealMarineOilAssay:
    """Immutable reduction, protein, and stowage specification for fishmeal and marine oil.

    Bridges whole pelagic fish catch with fishmeal and crude marine oil yields,
    protein content tiers, bag counts, and bulk stowage/tank volume <-> mass.
    """

    name: str
    species: str  # e.g. "ANCHOVETA", "MENHADEN", "CAPELIN", "SARDINE"
    meal_yield_pct: float = 22.5
    oil_yield_pct: float = 4.5
    protein_pct: float = 65.0
    oil_density_kg_per_m3: float = _CRUDE_MARINE_OIL_KG_PER_M3
    meal_bulk_density_kg_per_m3: float = _BULK_FISHMEAL_KG_PER_M3
    bag_mass_kg: float = _FISHMEAL_BAG_KG

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.species, str):
            raise ValueError("species must be a non-empty string.")
        sp = self.species.strip().upper()
        if not sp:
            raise ValueError("species must be a non-empty string.")

        meal_y = _validate_finite_float(self.meal_yield_pct, "meal_yield_pct")
        if not (0.0 < meal_y <= 100.0):
            raise ValueError(f"meal_yield_pct must be in (0, 100], got {self.meal_yield_pct}")

        oil_y = _validate_finite_float(self.oil_yield_pct, "oil_yield_pct")
        if not (0.0 <= oil_y <= 100.0) or (meal_y + oil_y > 100.0 + 1e-9):
            raise ValueError(
                f"Invalid reduction yields: meal_yield_pct={meal_y}, oil_yield_pct={oil_y}."
            )

        prot = _validate_finite_float(self.protein_pct, "protein_pct")
        if not (40.0 <= prot <= 80.0):
            raise ValueError(f"protein_pct must be in [40.0, 80.0], got {self.protein_pct}")

        oil_rho = _validate_finite_float(self.oil_density_kg_per_m3, "oil_density_kg_per_m3")
        if not (750.0 <= oil_rho <= 1100.0):
            raise ValueError(f"oil_density_kg_per_m3 must be in [750, 1100], got {self.oil_density_kg_per_m3}")

        meal_rho = _validate_finite_float(self.meal_bulk_density_kg_per_m3, "meal_bulk_density_kg_per_m3")
        if not (300.0 <= meal_rho <= 1200.0):
            raise ValueError(
                f"meal_bulk_density_kg_per_m3 must be in [300, 1200], got {self.meal_bulk_density_kg_per_m3}"
            )

        bag_kg = _validate_finite_float(self.bag_mass_kg, "bag_mass_kg")
        if bag_kg <= 0.0:
            raise ValueError(f"bag_mass_kg must be strictly positive, got {self.bag_mass_kg}")

        object.__setattr__(self, "species", sp)
        object.__setattr__(self, "meal_yield_pct", meal_y)
        object.__setattr__(self, "oil_yield_pct", oil_y)
        object.__setattr__(self, "protein_pct", prot)
        object.__setattr__(self, "oil_density_kg_per_m3", oil_rho)
        object.__setattr__(self, "meal_bulk_density_kg_per_m3", meal_rho)
        object.__setattr__(self, "bag_mass_kg", bag_kg)

    @property
    def combined_yield_pct(self) -> float:
        """Total meal plus oil mass recovery percentage from whole raw catch."""
        from core.units.precision import to_ieee754_float
        return to_ieee754_float(Decimal(str(self.meal_yield_pct)) + Decimal(str(self.oil_yield_pct)))

    def as_context(self, *, context_name: str | None = None, product: str = "meal") -> pint.Context:
        from core.units.agricultural_specs import create_fishmeal_marine_oil_context
        return create_fishmeal_marine_oil_context(self, name=context_name, product=product)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        number = Decimal(str(quantity.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Whole catch and product mass must be finite and non-negative.")
        return quantity

    @staticmethod
    def _ratio_scale(magnitude: Any, value: float, denominator: int) -> float | Decimal:
        fraction = Decimal(str(value)) / denominator
        return fraction if isinstance(magnitude, Decimal) else float(fraction)

    # --- Whole Catch <-> Meal & Marine Oil Conversions ---
    def whole_to_meal_mass(self, whole_mass: pint.Quantity) -> pint.Quantity:
        """Converts whole raw fish catch to fishmeal mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        mt_q = self._validate_mass(whole_mass)
        scale = self._ratio_scale(mt_q.magnitude, self.meal_yield_pct, 100)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def meal_to_whole_mass(self, meal_mass: pint.Quantity) -> pint.Quantity:
        """Converts fishmeal mass back to required whole raw fish catch (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(meal_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(meal_mass).__name__}")
        mt_q = self._validate_mass(meal_mass)
        scale = self._ratio_scale(mt_q.magnitude, self.meal_yield_pct, 100)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def whole_to_oil_mass(self, whole_mass: pint.Quantity) -> pint.Quantity:
        """Converts whole raw fish catch to crude marine oil mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        mt_q = self._validate_mass(whole_mass)
        scale = self._ratio_scale(mt_q.magnitude, self.oil_yield_pct, 100)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def oil_to_whole_mass(self, oil_mass: pint.Quantity) -> pint.Quantity:
        """Converts crude marine oil mass back to implied whole raw fish catch (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(oil_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(oil_mass).__name__}")
        if self.oil_yield_pct <= 0.0:
            raise ValueError(f"{self.name!r} has zero oil_yield_pct; cannot convert oil to whole catch.")
        mt_q = self._validate_mass(oil_mass)
        scale = self._ratio_scale(mt_q.magnitude, self.oil_yield_pct, 100)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def whole_to_meal_and_oil(self, whole_mass: pint.Quantity) -> dict[str, pint.Quantity]:
        """Returns both fishmeal and marine oil output quantities in metric tonnes."""
        return {
            "meal": self.whole_to_meal_mass(whole_mass),
            "oil": self.whole_to_oil_mass(whole_mass),
        }

    # --- Volume <-> Mass & Bag Count Conversions ---
    def oil_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        number = Decimal(str(m3_q.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Meal/oil volume must be finite and non-negative.")
        mt_per_m3 = Decimal(str(self.oil_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(m3_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def oil_mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = Decimal(str(self.oil_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def meal_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        number = Decimal(str(m3_q.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Meal/oil volume must be finite and non-negative.")
        mt_per_m3 = Decimal(str(self.meal_bulk_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(m3_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def meal_mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = Decimal(str(self.meal_bulk_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def meal_mass_to_bags(self, meal_mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal bag count, preserving Decimal input."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kilograms = self._validate_mass(meal_mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bag_mass_kg)

    def bags_to_meal_mass(self, bag_count: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        count = Decimal(str(bag_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bag_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bag_mass_kg))
        magnitude = kilograms if isinstance(bag_count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram).to(ureg.metric_ton)

    def gross_product_value_per_mt_whole(
        self, meal_price_per_mt: float | Decimal, oil_price_per_mt: float | Decimal,
    ) -> float | Decimal:
        """Gross meal/oil revenue per whole-catch tonne, before processing costs."""
        from core.units.precision import preserve_ieee754_scalar
        meal_price = Decimal(str(meal_price_per_mt))
        oil_price = Decimal(str(oil_price_per_mt))
        if not meal_price.is_finite() or not oil_price.is_finite():
            raise ValueError("Meal and oil prices must be finite.")
        value = Decimal(str(self.meal_yield_pct)) / 100 * meal_price + Decimal(str(self.oil_yield_pct)) / 100 * oil_price
        return preserve_ieee754_scalar(value, meal_price_per_mt, oil_price_per_mt)

    # --- Standard Reduction Presets ---
    @classmethod
    def peruvian_anchoveta_super_prime(cls) -> FishmealMarineOilAssay:
        """Peruvian anchoveta super-prime fishmeal (~68% protein, ~23% meal, ~5% oil)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "FishmealMarineOilAssay.peruvian_anchoveta_super_prime")

    @classmethod
    def menhaden_standard(cls) -> FishmealMarineOilAssay:
        """Atlantic/Gulf menhaden reduction benchmark (~62% protein, ~21.5% meal, ~6.5% oil)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "FishmealMarineOilAssay.menhaden_standard")


__all__ = ('FishmealMarineOilAssay', 'SeafoodAssay')
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
