"""Commodex assays agriculture feedstocks implementation and public interfaces."""
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
class StarchAssay(_agricultural_families_AgriculturalFamilyAssay):
    _commodities = frozenset({"corn", "wheat", "potato", "tapioca"})

    def starch_mass(self, quantity, *, mass_basis="received"):
        return self.contained_mass(quantity, component="starch", mass_basis=mass_basis)

    def recovered_starch_mass(self, quantity, *, mass_basis="received"):
        return self.recovered_mass(quantity, component="starch", mass_basis=mass_basis)


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class AnimalFeedAssay(_agricultural_families_AgriculturalFamilyAssay):
    _commodities = frozenset({"compound-feed", "protein-meal", "bran", "silage"})

    def nutrient_mass(self, quantity, *, nutrient, mass_basis="received"):
        return self.contained_mass(quantity, component=nutrient, mass_basis=mass_basis)

    def blend(self, other, own_mass, other_mass, *, name):
        """Wet mass and dry nutrient balance; blend density must be measured."""
        from core.units.assays.base import mass
        from core.units.composition import ElementComposition
        if not isinstance(other, AnimalFeedAssay) or self.composition is None or other.composition is None:
            raise TypeError("Blend two AnimalFeedAssay instances with declared compositions.")
        if self.composition.mass_fractions.keys() != other.composition.mass_fractions.keys():
            raise ValueError("Declare every nutrient on both feeds, including measured zero.")
        a, b = mass(own_mass), mass(other_mass)
        if a + b == 0:
            raise ValueError("Blend mass must be positive.")
        dry_a = mass(self.wet_to_dry_mass(own_mass))
        dry_b = mass(other.wet_to_dry_mass(other_mass))
        dry_total = dry_a + dry_b
        if dry_total == 0:
            raise ValueError("Blend must contain dry matter.")
        components = {k: (self.composition.fraction(k) * dry_a + other.composition.fraction(k) * dry_b) / dry_total for k in self.composition.mass_fractions}
        # Retain the weighted unassayed balance, correcting only the Decimal
        # division residual rather than normalizing partial analyses to 100%.
        total_fraction = (sum(self.composition.mass_fractions.values(), Decimal(0)) * dry_a
                          + sum(other.composition.mass_fractions.values(), Decimal(0)) * dry_b) / dry_total
        largest = max(components, key=components.get)
        components[largest] += total_fraction - sum(components.values(), Decimal(0))
        return AnimalFeedAssay(name=name, commodity="compound-feed", composition=ElementComposition(components),
                               moisture_fraction=1 - dry_total / (a + b))


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class AgriculturalResidueAssay(_agricultural_families_AgriculturalFamilyAssay):
    collectable_fraction: Decimal | None = None
    _commodities = frozenset({"straw", "husks", "bagasse", "other"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "collectable_fraction", _agricultural_families_yield(self.collectable_fraction, "collectable residue fraction"))

    def collectable_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.collectable_fraction, "collectable residue fraction")

    def ash_mass(self, quantity, *, mass_basis="received"):
        return self.contained_mass(quantity, component="ash", mass_basis=mass_basis)

    def fuel_energy(self, quantity, *, basis="lhv", target_unit="MWh"):
        from core.units.heat_content import convert_heat_content
        return convert_heat_content(quantity, target_unit, assay=self.as_heat_content(), basis=basis)


from core.units.agricultural_specs import _CPO_DENSITY_KG_PER_M3_50C
from core.units.validation import validated_assay


# ---------------------------------------------------------------------------
# 2. Oilseed Crush Specification (Soybeans, Canola/Rapeseed, Crude Palm Oil)
# ---------------------------------------------------------------------------


@validated_assay
@dataclass(frozen=True, slots=True)
class OilseedCrushAssay:
    """Immutable crush and volumetric specification for oilseeds and vegetable oils.

    Models bushel weight (e.g. 60 lb/bu soy, 50 lb/bu canola), oil and meal
    extraction yields, board/gross crush margins, and CPO bulk/thermal density.
    """

    name: str
    seed: str
    bushel_weight_lb: float
    oil_yield_pct: float
    meal_yield_pct: float
    oil_density_kg_per_m3: float = _CPO_DENSITY_KG_PER_M3_50C
    reference_temp_c: float = 50.0
    bulk_seed_density_kg_per_m3: float = 720.0

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        bu_lb = _validate_finite_float(self.bushel_weight_lb, "bushel_weight_lb")
        oil_y = _validate_finite_float(self.oil_yield_pct, "oil_yield_pct")
        meal_y = _validate_finite_float(self.meal_yield_pct, "meal_yield_pct")
        oil_rho = _validate_finite_float(self.oil_density_kg_per_m3, "oil_density_kg_per_m3")
        ref_t = _validate_finite_float(self.reference_temp_c, "reference_temp_c")
        seed_rho = _validate_finite_float(self.bulk_seed_density_kg_per_m3, "bulk_seed_density_kg_per_m3")

        if bu_lb <= 0.0:
            raise ValueError(f"bushel_weight_lb must be strictly positive, got {self.bushel_weight_lb}")
        if not (0.0 <= oil_y <= 100.0) or not (0.0 <= meal_y <= 100.0) or (oil_y + meal_y > 100.0 + 1e-6):
            raise ValueError(f"Invalid oil/meal yields: oil={oil_y}%, meal={meal_y}%.")
        if not (500.0 <= oil_rho <= 1500.0) or not (300.0 <= seed_rho <= 1500.0):
            raise ValueError("Densities must be physically realistic.")

        object.__setattr__(self, "bushel_weight_lb", bu_lb)
        object.__setattr__(self, "oil_yield_pct", oil_y)
        object.__setattr__(self, "meal_yield_pct", meal_y)
        object.__setattr__(self, "oil_density_kg_per_m3", oil_rho)
        object.__setattr__(self, "reference_temp_c", ref_t)
        object.__setattr__(self, "bulk_seed_density_kg_per_m3", seed_rho)

    @property
    def oil_lb_per_bushel(self) -> float:
        from core.units.precision import to_ieee754_float
        return to_ieee754_float(Decimal(str(self.bushel_weight_lb)) * Decimal(str(self.oil_yield_pct)) / Decimal(100))

    @property
    def meal_lb_per_bushel(self) -> float:
        from core.units.precision import to_ieee754_float
        return to_ieee754_float(Decimal(str(self.bushel_weight_lb)) * Decimal(str(self.meal_yield_pct)) / Decimal(100))

    def as_context(self, *, context_name: str | None = None, temp_c: float | Decimal | pint.Quantity | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_oilseed_crush_context
        return create_oilseed_crush_context(self, name=context_name, temp_c=temp_c)

    def volume_correction_factor(self, observed_temp: float | Decimal | pint.Quantity) -> Decimal:
        """Approximate observed-to-reference oil VCF within [-20, 100] degrees C."""
        from core.units.agricultural_specs import _CPO_EXPANSION_PER_C
        from core.units.exceptions import StandardsApplicabilityError
        from core.units.registry import ureg
        raw = observed_temp.to(ureg.degC).magnitude if isinstance(observed_temp, pint.Quantity) else observed_temp
        value = Decimal(str(raw))
        if not value.is_finite() or not (-20 <= value <= 100):
            raise StandardsApplicabilityError("Observed temperature must be finite and in [-20, 100] degrees C.")
        return 1 - Decimal(str(_CPO_EXPANSION_PER_C)) * (value - Decimal(str(self.reference_temp_c)))

    def seed_mass_to_bushels(self, seed_mass: pint.Quantity) -> float | Decimal:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(seed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_mass).__name__}")
        pounds = seed_mass.to(ureg.pound)
        number = Decimal(str(pounds.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Seed mass must be finite and non-negative.")
        return pounds.magnitude / _scale_factor(pounds.magnitude, self.bushel_weight_lb)

    def bushels_to_seed_mass(self, bushels: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        count = Decimal(str(bushels))
        if not count.is_finite() or count < 0:
            raise ValueError("bushels must be finite and non-negative.")
        pounds = count * Decimal(str(self.bushel_weight_lb))
        value = pounds if isinstance(bushels, Decimal) else float(pounds)
        return ureg.Quantity(value, ureg.pound).to(ureg.metric_ton)

    def seed_mass_to_oil_mass(self, seed_mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(seed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_mass).__name__}")
        mt_q = seed_mass.to(ureg.metric_ton)
        frac = self.oil_yield_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def seed_mass_to_meal_mass(self, seed_mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(seed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_mass).__name__}")
        mt_q = seed_mass.to(ureg.metric_ton)
        frac = self.meal_yield_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def oil_volume_to_mass(self, volume: pint.Quantity, *, temp_c: float | None = None) -> pint.Quantity:
        """Converts vegetable oil volume (m^3, gal) to mass (metric_ton) with optional temp adjustment."""
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        rho = self.oil_density_kg_per_m3
        if temp_c is not None:
            vcf = self.volume_correction_factor(temp_c)
            rho = Decimal(str(rho)) * vcf
        mt_per_m3 = Decimal(str(rho)) / 1000
        scale = mt_per_m3 if isinstance(m3_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def oil_mass_to_volume(self, mass: pint.Quantity, *, temp_c: float | None = None) -> pint.Quantity:
        """Converts vegetable oil mass (mt) to volume (m^3) with optional temp adjustment."""
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        rho = self.oil_density_kg_per_m3
        if temp_c is not None:
            vcf = self.volume_correction_factor(temp_c)
            rho = Decimal(str(rho)) * vcf
        mt_per_m3 = Decimal(str(rho)) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def gross_crush_margin_usd_per_mt(
        self,
        seed_usd_per_mt: float | Decimal,
        oil_usd_per_mt: float | Decimal,
        meal_usd_per_mt: float | Decimal,
    ) -> float:
        """Gross crush margin ($/mt seed): (oil_yield * oil_px + meal_yield * meal_px) - seed_px."""
        from core.units.agricultural_specs import _validate_finite_float
        from core.units.precision import to_ieee754_float
        s_px = _validate_finite_float(seed_usd_per_mt, "seed_usd_per_mt")
        o_px = _validate_finite_float(oil_usd_per_mt, "oil_usd_per_mt")
        m_px = _validate_finite_float(meal_usd_per_mt, "meal_usd_per_mt")
        out_val = (self.oil_yield_pct / 100.0) * o_px + (self.meal_yield_pct / 100.0) * m_px
        return to_ieee754_float(out_val - s_px)

    @classmethod
    def soybeans_cb_60lb(cls) -> OilseedCrushAssay:
        # Standard 60 lb/bu: ~11 lb oil (~18.33%) + ~44 lb meal (~73.33%)
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "OilseedCrushAssay.soybeans_cb_60lb")

    @classmethod
    def canola_50lb(cls) -> OilseedCrushAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "OilseedCrushAssay.canola_50lb")

    @classmethod
    def crude_palm_oil_cpo(cls) -> OilseedCrushAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "OilseedCrushAssay.crude_palm_oil_cpo")


from core.units.agricultural_specs import _ETHANOL_DENSITY_KG_PER_M3
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class EthanolFeedstockAssay:
    """Immutable specification for grain-to-ethanol conversion and DDGS yield.

    Bridges grain feedstock mass and bushels with denatured fuel ethanol volume,
    ethanol mass, Distillers Dried Grains with Solubles (DDGS) co-product mass,
    and gross ethanol crush economics.
    """

    name: str
    feedstock: str
    bushel_weight_lb: float
    ethanol_gal_per_bushel: float
    ddgs_lb_per_bushel: float
    ethanol_density_kg_per_m3: float = _ETHANOL_DENSITY_KG_PER_M3

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.feedstock, str):
            raise ValueError("feedstock must be a non-empty string.")
        feed = self.feedstock.strip().lower()
        if not feed:
            raise ValueError("feedstock must be a non-empty string.")

        bu_lb = _validate_finite_float(self.bushel_weight_lb, "bushel_weight_lb")
        if bu_lb <= 0.0:
            raise ValueError(f"bushel_weight_lb must be strictly positive, got {self.bushel_weight_lb}")

        etoh_gal = _validate_finite_float(self.ethanol_gal_per_bushel, "ethanol_gal_per_bushel")
        if etoh_gal <= 0.0:
            raise ValueError(
                f"ethanol_gal_per_bushel must be strictly positive, got {self.ethanol_gal_per_bushel}"
            )

        ddgs_lb = _validate_finite_float(self.ddgs_lb_per_bushel, "ddgs_lb_per_bushel")
        if ddgs_lb < 0.0 or ddgs_lb > bu_lb:
            raise ValueError(
                f"ddgs_lb_per_bushel must be in [0, {bu_lb}], got {self.ddgs_lb_per_bushel}"
            )

        rho = _validate_finite_float(self.ethanol_density_kg_per_m3, "ethanol_density_kg_per_m3")
        if not (650.0 <= rho <= 950.0):
            raise ValueError(
                f"ethanol_density_kg_per_m3 must be in [650, 950], got {self.ethanol_density_kg_per_m3}"
            )

        object.__setattr__(self, "feedstock", feed)
        object.__setattr__(self, "bushel_weight_lb", bu_lb)
        object.__setattr__(self, "ethanol_gal_per_bushel", etoh_gal)
        object.__setattr__(self, "ddgs_lb_per_bushel", ddgs_lb)
        object.__setattr__(self, "ethanol_density_kg_per_m3", rho)

    @property
    def ethanol_litres_per_mt_feedstock(self) -> float:
        """Litres of denatured ethanol produced per metric tonne of grain feedstock."""
        from core.units.agricultural_specs import _LB_PER_MT
        from core.units.agricultural_specs import _LITRES_PER_US_GALLON
        from core.units.precision import to_ieee754_float
        bushels_per_mt = _LB_PER_MT / self.bushel_weight_lb
        litres_per_gal = _LITRES_PER_US_GALLON
        return to_ieee754_float(bushels_per_mt * self.ethanol_gal_per_bushel * litres_per_gal)

    @property
    def ddgs_yield_pct_by_mass(self) -> float:
        """DDGS output mass as a percentage of grain feedstock mass."""
        from core.units.precision import to_ieee754_float
        return to_ieee754_float((self.ddgs_lb_per_bushel / self.bushel_weight_lb) * 100.0)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_ethanol_feedstock_context
        return create_ethanol_feedstock_context(self, name=context_name)

    @staticmethod
    def _validated_number(value: float | Decimal, label: str, *, nonnegative: bool = True) -> Decimal:
        number = Decimal(str(value))
        if not number.is_finite() or (nonnegative and number < 0):
            raise ValueError(f"{label} must be finite" + (" and non-negative." if nonnegative else "."))
        return number

    @staticmethod
    def _validate_quantity(quantity: pint.Quantity, unit: pint.Unit) -> pint.Quantity:
        if not isinstance(quantity, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
        converted = quantity.to(unit)
        EthanolFeedstockAssay._validated_number(converted.magnitude, "mass or volume")
        return converted

    @staticmethod
    def _number_like(number: Decimal, original: Any) -> float | Decimal:
        return number if isinstance(original, Decimal) else float(number)

    # --- Feedstock Mass <-> Bushels Conversions ---
    def feedstock_mass_to_bushels(self, feedstock_mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded crop-weight bushels, preserving Decimal mass input."""
        from core.units.registry import ureg
        pounds = self._validate_quantity(feedstock_mass, ureg.pound)
        count = Decimal(str(pounds.magnitude)) / Decimal(str(self.bushel_weight_lb))
        return self._number_like(count, pounds.magnitude)

    def bushels_to_feedstock_mass(self, bushels: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        count = self._validated_number(bushels, "bushels")
        pounds = count * Decimal(str(self.bushel_weight_lb))
        return ureg.Quantity(self._number_like(pounds, bushels), ureg.pound).to(ureg.metric_ton)

    # --- Ethanol & DDGS Output Conversions ---
    def feedstock_mass_to_ethanol_volume(self, feedstock_mass: pint.Quantity) -> pint.Quantity:
        """Returns US gallons of ethanol produced from grain feedstock mass."""
        from core.units.registry import ureg
        bushels = self.feedstock_mass_to_bushels(feedstock_mass)
        gallons = Decimal(str(bushels)) * Decimal(str(self.ethanol_gal_per_bushel))
        return ureg.Quantity(self._number_like(gallons, bushels), ureg.gallon)

    def ethanol_volume_to_feedstock_mass(self, ethanol_volume: pint.Quantity) -> pint.Quantity:
        """Returns required feedstock tonnes from ethanol volume, not ethanol mass."""
        from core.units.registry import ureg
        gallons = self._validate_quantity(ethanol_volume, ureg.gallon)
        bushels = Decimal(str(gallons.magnitude)) / Decimal(str(self.ethanol_gal_per_bushel))
        return self.bushels_to_feedstock_mass(self._number_like(bushels, gallons.magnitude))

    def feedstock_mass_to_ddgs_mass(self, feedstock_mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        bushels = self.feedstock_mass_to_bushels(feedstock_mass)
        pounds = Decimal(str(bushels)) * Decimal(str(self.ddgs_lb_per_bushel))
        return ureg.Quantity(self._number_like(pounds, bushels), ureg.pound).to(ureg.metric_ton)

    def ddgs_mass_to_feedstock_mass(self, ddgs_mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if self.ddgs_lb_per_bushel == 0:
            raise ValueError("Cannot infer feedstock requirements from zero DDGS yield.")
        pounds = self._validate_quantity(ddgs_mass, ureg.pound)
        bushels = Decimal(str(pounds.magnitude)) / Decimal(str(self.ddgs_lb_per_bushel))
        return self.bushels_to_feedstock_mass(self._number_like(bushels, pounds.magnitude))

    def ethanol_volume_to_mass(self, ethanol_volume: pint.Quantity) -> pint.Quantity:
        """Returns physical ethanol mass using density, independently of grain yield."""
        from core.units.registry import ureg
        volume = self._validate_quantity(ethanol_volume, ureg.meter ** 3)
        kilograms = Decimal(str(volume.magnitude)) * Decimal(str(self.ethanol_density_kg_per_m3))
        return ureg.Quantity(self._number_like(kilograms, volume.magnitude), ureg.kilogram).to(ureg.metric_ton)

    def ethanol_mass_to_volume(self, ethanol_mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        kilograms = self._validate_quantity(ethanol_mass, ureg.kilogram)
        volume = Decimal(str(kilograms.magnitude)) / Decimal(str(self.ethanol_density_kg_per_m3))
        return ureg.Quantity(self._number_like(volume, kilograms.magnitude), ureg.meter ** 3)

    # --- Crush Margin & Breakeven Parity Helpers ---
    def _gross_margin_decimal(self, grain_price: float | Decimal, ethanol_price: float | Decimal, ddgs_price: float | Decimal) -> Decimal:
        grain = self._validated_number(grain_price, "grain price", nonnegative=False)
        ethanol = self._validated_number(ethanol_price, "ethanol price", nonnegative=False)
        ddgs = self._validated_number(ddgs_price, "DDGS price", nonnegative=False)
        ddgs_short_tons = Decimal(str(self.ddgs_lb_per_bushel)) / 2000
        return Decimal(str(self.ethanol_gal_per_bushel)) * ethanol + ddgs_short_tons * ddgs - grain

    @staticmethod
    def _rounded_economics(value: Decimal, *prices: Any) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        return preserve_ieee754_scalar(value, *prices)

    def gross_crush_margin_usd_per_bushel(
        self, grain_usd_per_bushel: float | Decimal, ethanol_usd_per_gallon: float | Decimal,
        ddgs_usd_per_short_ton: float | Decimal,
    ) -> float | Decimal:
        """Gross ethanol + DDGS revenue less grain cost; excludes processing costs."""
        prices = (grain_usd_per_bushel, ethanol_usd_per_gallon, ddgs_usd_per_short_ton)
        return self._rounded_economics(self._gross_margin_decimal(*prices), *prices)

    def gross_crush_margin_usd_per_mt(
        self, grain_usd_per_bushel: float | Decimal, ethanol_usd_per_gallon: float | Decimal,
        ddgs_usd_per_short_ton: float | Decimal,
    ) -> float | Decimal:
        """Gross margin per feedstock tonne, without rounding the per-bushel leg."""
        from core.units.registry import ureg
        prices = (grain_usd_per_bushel, ethanol_usd_per_gallon, ddgs_usd_per_short_ton)
        bushels = self.feedstock_mass_to_bushels(ureg.Quantity(Decimal("1"), ureg.metric_ton))
        return self._rounded_economics(self._gross_margin_decimal(*prices) * bushels, *prices)

    def breakeven_ethanol_usd_per_gallon(
        self, grain_usd_per_bushel: float | Decimal, ddgs_usd_per_short_ton: float | Decimal,
    ) -> float | Decimal:
        """Ethanol price covering grain cost after DDGS credit; excludes processing."""
        deficit = -self._gross_margin_decimal(grain_usd_per_bushel, 0, ddgs_usd_per_short_ton)
        price = deficit / Decimal(str(self.ethanol_gal_per_bushel))
        return self._rounded_economics(price, grain_usd_per_bushel, ddgs_usd_per_short_ton)

    def breakeven_grain_usd_per_bushel(
        self, ethanol_usd_per_gallon: float | Decimal, ddgs_usd_per_short_ton: float | Decimal,
    ) -> float | Decimal:
        """Maximum grain price supported by gross ethanol + DDGS revenue."""
        revenue = self._gross_margin_decimal(0, ethanol_usd_per_gallon, ddgs_usd_per_short_ton)
        return self._rounded_economics(revenue, ethanol_usd_per_gallon, ddgs_usd_per_short_ton)

    @classmethod
    def us_corn_ethanol(cls) -> EthanolFeedstockAssay:
        """Standard US dry-grind corn ethanol (56 lb/bu, 2.85 gal/bu, 17.0 lb DDGS/bu)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "EthanolFeedstockAssay.us_corn_ethanol")

    @classmethod
    def wheat_ethanol(cls) -> EthanolFeedstockAssay:
        """Standard wheat ethanol benchmark (60 lb/bu, 2.70 gal/bu, 18.0 lb DDGS/bu)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "EthanolFeedstockAssay.wheat_ethanol")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (StarchAssay, AnimalFeedAssay, AgriculturalResidueAssay,))
__all__ = ('AgriculturalResidueAssay', 'AnimalFeedAssay', 'EthanolFeedstockAssay', 'OilseedCrushAssay', 'StarchAssay', 'convert_agricultural_residue', 'convert_animal_feed', 'convert_starch', 'create_agricultural_residue_context', 'create_animal_feed_context', 'create_starch_context')
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
    from core.units.assays.agriculture.feedstocks import convert_agricultural_residue
    from core.units.assays.agriculture.feedstocks import convert_animal_feed
    from core.units.assays.agriculture.feedstocks import convert_starch
    from core.units.assays.agriculture.feedstocks import create_agricultural_residue_context
    from core.units.assays.agriculture.feedstocks import create_animal_feed_context
    from core.units.assays.agriculture.feedstocks import create_starch_context
