"""Commodex assays agriculture poultry implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final
import math
import pint

from core.units.agricultural_specs import _DEFAULT_LARGE_EGG_GRAMS
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class PoultryEggAssay:
    """Immutable yield and packaging specification for poultry meat and shell/breaker eggs.

    Supports live-to-RTC (ready-to-cook) eviscerated carcass and boneless meat
    yields for broilers/turkeys, as well as egg/dozen/30-dozen case counts to
    shell, liquid, and dried mass equivalents. Both egg recovery percentages
    are relative to shell-egg input, not sequential processing yields.
    Bulk density describes the specified product state; the egg preset uses
    liquid whole egg, not a case or shell-egg packing density.
    """

    name: str
    category: str  # "broiler", "turkey", or "eggs"
    rtc_dressing_pct: float = 74.0          # Live -> Ready-to-Cook (RTC) eviscerated carcass (%)
    boneless_meat_pct_of_live: float = 52.0  # Live -> boneless/skinless meat yield (%)
    egg_mass_grams: float = _DEFAULT_LARGE_EGG_GRAMS
    liquid_egg_yield_pct: float = 88.0      # Shell egg -> liquid whole egg recovery (%)
    dried_egg_yield_pct: float = 24.0       # Shell egg -> spray-dried whole egg powder yield (%)
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.category, str):
            raise ValueError("category must be 'broiler', 'turkey', or 'eggs'.")
        cat = self.category.strip().lower()
        if cat not in {"broiler", "turkey", "eggs"}:
            raise ValueError(f"category must be 'broiler', 'turkey', or 'eggs', got {self.category!r}")

        rtc = _validate_finite_float(self.rtc_dressing_pct, "rtc_dressing_pct")
        if not (0.0 < rtc <= 100.0):
            raise ValueError(f"rtc_dressing_pct must be in (0, 100], got {self.rtc_dressing_pct}")

        boneless = _validate_finite_float(self.boneless_meat_pct_of_live, "boneless_meat_pct_of_live")
        if not (0.0 < boneless <= 100.0) or boneless > rtc:
            raise ValueError(
                f"boneless_meat_pct_of_live ({boneless}) must be in (0, 100] and <= rtc_dressing_pct ({rtc})."
            )

        egg_g = _validate_finite_float(self.egg_mass_grams, "egg_mass_grams")
        if not (20.0 <= egg_g <= 150.0):
            raise ValueError(f"egg_mass_grams must be in [20, 150], got {self.egg_mass_grams}")

        liq_y = _validate_finite_float(self.liquid_egg_yield_pct, "liquid_egg_yield_pct")
        if not (0.0 < liq_y <= 100.0):
            raise ValueError(f"liquid_egg_yield_pct must be in (0, 100], got {self.liquid_egg_yield_pct}")

        dry_y = _validate_finite_float(self.dried_egg_yield_pct, "dried_egg_yield_pct")
        if not (0.0 < dry_y <= 100.0) or dry_y > liq_y:
            raise ValueError(
                f"dried_egg_yield_pct ({dry_y}) must be in (0, 100] and <= liquid_egg_yield_pct ({liq_y})."
            )

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        object.__setattr__(self, "category", cat)
        object.__setattr__(self, "rtc_dressing_pct", rtc)
        object.__setattr__(self, "boneless_meat_pct_of_live", boneless)
        object.__setattr__(self, "egg_mass_grams", egg_g)
        object.__setattr__(self, "liquid_egg_yield_pct", liq_y)
        object.__setattr__(self, "dried_egg_yield_pct", dry_y)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)

    @property
    def dozen_mass_kg(self) -> float:
        """Net shell-egg mass per dozen in kilograms."""
        from core.units.precision import to_ieee754_float
        return to_ieee754_float((Decimal(str(self.egg_mass_grams)) * Decimal("12")) / Decimal("1000"))

    @property
    def case_mass_kg(self) -> float:
        """Net shell-egg mass per 30-dozen wholesale case in kilograms."""
        from core.units.precision import to_ieee754_float
        return to_ieee754_float((Decimal(str(self.egg_mass_grams)) * Decimal("360")) / Decimal("1000"))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_poultry_egg_context
        return create_poultry_egg_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.kilogram)
        value = Decimal(str(quantity.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Poultry/egg mass must be finite and non-negative.")
        return quantity

    def _require_eggs(self) -> None:
        from core.units.exceptions import IncompatibleCommodityError
        if self.category != "eggs":
            raise IncompatibleCommodityError("Egg count and processing methods require category='eggs'.")

    def _require_poultry(self) -> None:
        from core.units.exceptions import IncompatibleCommodityError
        if self.category not in {"broiler", "turkey"}:
            raise IncompatibleCommodityError("Poultry yield methods require a broiler or turkey assay.")

    def _yield_mass(self, mass: pint.Quantity, percentage: float, *, reverse: bool = False) -> pint.Quantity:
        from core.units.registry import ureg
        quantity = self._validate_mass(mass)
        fraction = Decimal(str(percentage)) / 100
        scale = fraction if isinstance(quantity.magnitude, Decimal) else float(fraction)
        value = quantity.magnitude / scale if reverse else quantity.magnitude * scale
        return ureg.Quantity(value, ureg.kilogram)

    # --- Poultry Carcass Conversions ---
    def live_to_rtc_mass(self, live_mass: pint.Quantity) -> pint.Quantity:
        self._require_poultry()
        return self._yield_mass(live_mass, self.rtc_dressing_pct)

    def rtc_to_live_mass(self, rtc_mass: pint.Quantity) -> pint.Quantity:
        self._require_poultry()
        return self._yield_mass(rtc_mass, self.rtc_dressing_pct, reverse=True)

    def live_to_boneless_mass(self, live_mass: pint.Quantity) -> pint.Quantity:
        self._require_poultry()
        return self._yield_mass(live_mass, self.boneless_meat_pct_of_live)

    def boneless_to_live_mass(self, boneless_mass: pint.Quantity) -> pint.Quantity:
        self._require_poultry()
        return self._yield_mass(boneless_mass, self.boneless_meat_pct_of_live, reverse=True)

    # --- Individual Egg / Dozen / Case / Liquid / Dried Conversions ---
    def eggs_to_shell_mass(self, eggs: float | Decimal) -> pint.Quantity:
        """Returns nominal shell-egg mass in kg; fractional counts are permitted."""
        from core.units.registry import ureg
        self._require_eggs()
        count = Decimal(str(eggs))
        if not count.is_finite() or count < 0:
            raise ValueError("Egg count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.egg_mass_grams)) / 1000
        value = kilograms if isinstance(eggs, Decimal) else float(kilograms)
        return ureg.Quantity(value, ureg.kilogram)

    def shell_mass_to_eggs(self, shell_mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal individual-egg count."""
        self._require_eggs()
        mass = self._validate_mass(shell_mass)
        per_egg = Decimal(str(self.egg_mass_grams)) / 1000
        scale = per_egg if isinstance(mass.magnitude, Decimal) else float(per_egg)
        return mass.magnitude / scale

    def dozens_to_shell_mass(self, dozens: float | Decimal) -> pint.Quantity:
        return self.eggs_to_shell_mass(dozens * 12)

    def cases_to_shell_mass(self, cases: float | Decimal) -> pint.Quantity:
        """Returns shell-egg mass for nominal 30-dozen/360-egg cases."""
        return self.eggs_to_shell_mass(cases * 360)

    def shell_mass_to_dozens(self, shell_mass: pint.Quantity) -> float | Decimal:
        return self.shell_mass_to_eggs(shell_mass) / 12

    def shell_mass_to_cases(self, shell_mass: pint.Quantity) -> float | Decimal:
        return self.shell_mass_to_eggs(shell_mass) / 360

    def shell_to_liquid_egg_mass(self, shell_mass: pint.Quantity) -> pint.Quantity:
        self._require_eggs()
        return self._yield_mass(shell_mass, self.liquid_egg_yield_pct)

    def liquid_to_shell_egg_mass(self, liquid_mass: pint.Quantity) -> pint.Quantity:
        self._require_eggs()
        return self._yield_mass(liquid_mass, self.liquid_egg_yield_pct, reverse=True)

    def shell_to_dried_egg_mass(self, shell_mass: pint.Quantity) -> pint.Quantity:
        self._require_eggs()
        return self._yield_mass(shell_mass, self.dried_egg_yield_pct)

    def dried_to_shell_egg_mass(self, dried_mass: pint.Quantity) -> pint.Quantity:
        self._require_eggs()
        return self._yield_mass(dried_mass, self.dried_egg_yield_pct, reverse=True)

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
            raise ValueError("Bulk volume must be finite and non-negative.")
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
            raise ValueError("Poultry/egg price must be finite.")
        return value

    @staticmethod
    def _price_result(value: Decimal, original: float | Decimal) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        return preserve_ieee754_scalar(value, original)

    def live_price_to_rtc_price(self, live_price: float | Decimal) -> float | Decimal:
        """Returns RTC price parity on the same currency/mass unit, before processing costs."""
        self._require_poultry()
        price = self._validated_price(live_price) / (Decimal(str(self.rtc_dressing_pct)) / 100)
        return self._price_result(price, live_price)

    def rtc_price_to_live_price(self, rtc_price: float | Decimal) -> float | Decimal:
        self._require_poultry()
        price = self._validated_price(rtc_price) * Decimal(str(self.rtc_dressing_pct)) / 100
        return self._price_result(price, rtc_price)

    def usd_per_dozen_to_liquid_usd_per_lb(self, usd_per_dozen: float | Decimal) -> float | Decimal:
        """Returns liquid whole-egg USD/lb parity, before breaking/processing costs."""
        from core.units.registry import ureg
        self._require_eggs()
        liquid = self.shell_to_liquid_egg_mass(self.dozens_to_shell_mass(Decimal("1")))
        pounds = liquid.to(ureg.pound).magnitude
        return self._price_result(self._validated_price(usd_per_dozen) / pounds, usd_per_dozen)

    def liquid_usd_per_lb_to_usd_per_dozen(self, usd_per_lb: float | Decimal) -> float | Decimal:
        from core.units.registry import ureg
        self._require_eggs()
        liquid = self.shell_to_liquid_egg_mass(self.dozens_to_shell_mass(Decimal("1")))
        pounds = liquid.to(ureg.pound).magnitude
        return self._price_result(self._validated_price(usd_per_lb) * pounds, usd_per_lb)

    # --- Standard Poultry & Egg Presets ---
    @classmethod
    def broiler_chicken_standard(cls) -> PoultryEggAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "PoultryEggAssay.broiler_chicken_standard")

    @classmethod
    def turkey_standard(cls) -> PoultryEggAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "PoultryEggAssay.turkey_standard")

    @classmethod
    def us_large_shell_eggs(cls) -> PoultryEggAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "PoultryEggAssay.us_large_shell_eggs")


__all__ = ('PoultryEggAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
