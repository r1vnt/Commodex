from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
import math
from types import MappingProxyType
from typing import Any, Final

import pint

from .conversions import _scale_factor, convert as convert_units
from .registry import parse_unit
from .registry import ureg

__all__: Final[tuple[str, ...]] = (
    "SoftsAssay",
    "OilseedCrushAssay",
    "LivestockDairyAssay",
    "GrainAssay",
    "SeafoodAssay",
    "PulpLumberAssay",
    "RubberAssay",
    "EthanolFeedstockAssay",
    "RiceAssay",
    "WoolAssay",
    "PoultryEggAssay",
    "FishmealMarineOilAssay",
    "HoneyApicultureAssay",
    "OliveOilAssay",
    "CottonSeedLintAssay",
    "create_softs_context",
    "create_oilseed_crush_context",
    "create_livestock_dairy_context",
    "create_grain_context",
    "create_seafood_context",
    "create_pulp_lumber_context",
    "create_rubber_context",
    "create_ethanol_feedstock_context",
    "create_rice_context",
    "create_wool_context",
    "create_poultry_egg_context",
    "create_fishmeal_marine_oil_context",
    "create_honey_apiculture_context",
    "create_olive_oil_context",
    "create_cotton_seed_lint_context",
    "convert_softs",
    "convert_oilseed_crush",
    "convert_livestock_dairy",
    "convert_grain",
    "convert_seafood",
    "convert_pulp_lumber",
    "convert_rubber",
    "convert_ethanol_feedstock",
    "convert_rice",
    "convert_wool",
    "convert_poultry_egg",
    "convert_fishmeal_marine_oil",
    "convert_honey_apiculture",
    "convert_olive_oil",
    "convert_cotton_seed_lint",
)

_LB_PER_MT: Final[float] = 2204.6226218487757
_CWT_US_LB: Final[float] = 100.0
_CWT_METRIC_KG: Final[float] = 100.0
_CPO_DENSITY_KG_PER_M3_50C: Final[float] = 890.0  # Typical heated CPO bulk density (~50°C)
_CPO_EXPANSION_PER_C: Final[float] = 0.0008     # Approximate volumetric expansion coefficient (/°C)


def _validate_finite_float(val: float | Decimal, name: str) -> float:
    fval = float(val)
    if not math.isfinite(fval):
        raise ValueError(f"{name} must be finite, got {val!r}")
    return fval


# ---------------------------------------------------------------------------
# 1. Softs Specification (Sugar, Coffee, Cocoa, Cotton)
# ---------------------------------------------------------------------------


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
        if self.reference_polarisation_deg is None or self.actual_polarisation_deg is None:
            return 1.0
        ref = _validate_finite_float(self.reference_polarisation_deg, "reference_polarisation_deg")
        act = _validate_finite_float(self.actual_polarisation_deg, "actual_polarisation_deg")
        if ref <= 0.0 or act <= 0.0:
            raise ValueError("Polarisation degrees must be strictly positive.")
        return act / ref

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_softs_context(self, name=context_name)

    def raw_to_processed_mass(self, raw_mass: pint.Quantity) -> pint.Quantity:
        """Converts raw/unprocessed mass (e.g. raw sugar, green coffee, cocoa beans) to output mass."""
        if not isinstance(raw_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(raw_mass).__name__}")
        mt_q = raw_mass.to(ureg.metric_ton)
        eff_frac = self.processing_yield_fraction * self.polarisation_factor
        scale = _scale_factor(mt_q.magnitude, eff_frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def processed_to_raw_mass(self, processed_mass: pint.Quantity) -> pint.Quantity:
        """Converts processed output mass (e.g. refined white sugar, roasted coffee) to raw input mass."""
        if not isinstance(processed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(processed_mass).__name__}")
        mt_q = processed_mass.to(ureg.metric_ton)
        eff_frac = self.processing_yield_fraction * self.polarisation_factor
        scale = _scale_factor(mt_q.magnitude, eff_frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def mass_to_packages(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal package count, preserving Decimal input."""
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
        if self.bag_or_bale_kg is None:
            raise ValueError(f"{self.name!r} does not define bag_or_bale_kg.")
        number = Decimal(str(count))
        if not number.is_finite() or number < 0:
            raise ValueError("Package count must be finite and non-negative.")
        kilograms = number * Decimal(str(self.bag_or_bale_kg))
        value = kilograms if isinstance(count, Decimal) else float(kilograms)
        return ureg.Quantity(value, ureg.kilogram).to(ureg.metric_ton)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    @staticmethod
    def cents_per_lb_to_usd_per_mt(cents_per_lb: float | Decimal) -> float:
        """Converts ICE quotes in US cents/lb (e.g. No. 11 Sugar, Cotton) to $/metric_ton."""
        c = _validate_finite_float(cents_per_lb, "cents_per_lb")
        return round((c / 100.0) * _LB_PER_MT, 4)

    @staticmethod
    def usd_per_mt_to_cents_per_lb(usd_per_mt: float | Decimal) -> float:
        """Converts $/metric_ton quote to US cents/lb."""
        p = _validate_finite_float(usd_per_mt, "usd_per_mt")
        return round((p / _LB_PER_MT) * 100.0, 4)

    def raw_mt_price_to_refined_mt_price(self, raw_price_per_mt: float | Decimal) -> float:
        """Converts raw/unprocessed $/mt quote to implied processed $/mt parity."""
        p = _validate_finite_float(raw_price_per_mt, "raw_price_per_mt")
        eff_frac = self.processing_yield_fraction * self.polarisation_factor
        return round(p / eff_frac, 4)

    @classmethod
    def raw_sugar_no11(cls) -> SoftsAssay:
        return cls(
            name="Raw Sugar (ICE No. 11, 96° Pol)",
            commodity="sugar",
            processing_yield_pct=92.5,
            reference_polarisation_deg=96.0,
            actual_polarisation_deg=96.0,
            bulk_density_kg_per_m3=850.0,
        )

    @classmethod
    def green_coffee_60kg_bag(cls) -> SoftsAssay:
        return cls(
            name="Green Coffee Beans (60 kg Bag)",
            commodity="coffee",
            bag_or_bale_kg=60.0,
            processing_yield_pct=82.0,  # Green-to-roasted yield
            bulk_density_kg_per_m3=650.0,
        )

    @classmethod
    def cocoa_beans_62_5kg_bag(cls) -> SoftsAssay:
        return cls(
            name="Cocoa Beans (62.5 kg Bag)",
            commodity="cocoa",
            bag_or_bale_kg=62.5,
            processing_yield_pct=45.0,   # Typical cocoa butter extraction
            secondary_yield_pct=40.0,  # Typical cocoa cake/powder yield
            bulk_density_kg_per_m3=620.0,
        )

    @classmethod
    def us_cotton_480lb_bale(cls) -> SoftsAssay:
        return cls(
            name="US Cotton (480 lb Bale)",
            commodity="cotton",
            bag_or_bale_kg=480.0 * 0.45359237,
            processing_yield_pct=100.0,
        )

    @classmethod
    def indian_cotton_170kg_bale(cls) -> SoftsAssay:
        return cls(
            name="Indian Cotton (170 kg Bale)",
            commodity="cotton",
            bag_or_bale_kg=170.0,
            processing_yield_pct=100.0,
        )


# ---------------------------------------------------------------------------
# 2. Oilseed Crush Specification (Soybeans, Canola/Rapeseed, Crude Palm Oil)
# ---------------------------------------------------------------------------


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
        return round(self.bushel_weight_lb * (self.oil_yield_pct / 100.0), 4)

    @property
    def meal_lb_per_bushel(self) -> float:
        return round(self.bushel_weight_lb * (self.meal_yield_pct / 100.0), 4)

    def as_context(self, *, context_name: str | None = None, temp_c: float | Decimal | pint.Quantity | None = None) -> pint.Context:
        return create_oilseed_crush_context(self, name=context_name, temp_c=temp_c)

    def volume_correction_factor(self, observed_temp: float | Decimal | pint.Quantity) -> Decimal:
        """Approximate observed-to-reference oil VCF within [-20, 100] degrees C."""
        raw = observed_temp.to(ureg.degC).magnitude if isinstance(observed_temp, pint.Quantity) else observed_temp
        value = Decimal(str(raw))
        if not value.is_finite() or not (-20 <= value <= 100):
            raise ValueError("Observed temperature must be finite and in [-20, 100] degrees C.")
        return 1 - Decimal(str(_CPO_EXPANSION_PER_C)) * (value - Decimal(str(self.reference_temp_c)))

    def seed_mass_to_bushels(self, seed_mass: pint.Quantity) -> float | Decimal:
        if not isinstance(seed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_mass).__name__}")
        pounds = seed_mass.to(ureg.pound)
        number = Decimal(str(pounds.magnitude))
        if not number.is_finite() or number < 0:
            raise ValueError("Seed mass must be finite and non-negative.")
        return pounds.magnitude / _scale_factor(pounds.magnitude, self.bushel_weight_lb)

    def bushels_to_seed_mass(self, bushels: float | Decimal) -> pint.Quantity:
        count = Decimal(str(bushels))
        if not count.is_finite() or count < 0:
            raise ValueError("bushels must be finite and non-negative.")
        pounds = count * Decimal(str(self.bushel_weight_lb))
        value = pounds if isinstance(bushels, Decimal) else float(pounds)
        return ureg.Quantity(value, ureg.pound).to(ureg.metric_ton)

    def seed_mass_to_oil_mass(self, seed_mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(seed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_mass).__name__}")
        mt_q = seed_mass.to(ureg.metric_ton)
        frac = self.oil_yield_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def seed_mass_to_meal_mass(self, seed_mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(seed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_mass).__name__}")
        mt_q = seed_mass.to(ureg.metric_ton)
        frac = self.meal_yield_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def oil_volume_to_mass(self, volume: pint.Quantity, *, temp_c: float | None = None) -> pint.Quantity:
        """Converts vegetable oil volume (m^3, gal) to mass (metric_ton) with optional temp adjustment."""
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
        s_px = _validate_finite_float(seed_usd_per_mt, "seed_usd_per_mt")
        o_px = _validate_finite_float(oil_usd_per_mt, "oil_usd_per_mt")
        m_px = _validate_finite_float(meal_usd_per_mt, "meal_usd_per_mt")
        out_val = (self.oil_yield_pct / 100.0) * o_px + (self.meal_yield_pct / 100.0) * m_px
        return round(out_val - s_px, 4)

    @classmethod
    def soybeans_cb_60lb(cls) -> OilseedCrushAssay:
        # Standard 60 lb/bu: ~11 lb oil (~18.33%) + ~44 lb meal (~73.33%)
        return cls(
            name="US Soybeans (60 lb/bu, ~18.3% Oil / ~73.3% Meal)",
            seed="soybean",
            bushel_weight_lb=60.0,
            oil_yield_pct=(11.0 / 60.0) * 100.0,
            meal_yield_pct=(44.0 / 60.0) * 100.0,
            oil_density_kg_per_m3=920.0,
            reference_temp_c=25.0,
            bulk_seed_density_kg_per_m3=770.0,
        )

    @classmethod
    def canola_50lb(cls) -> OilseedCrushAssay:
        return cls(
            name="Canola / Rapeseed (50 lb/bu, ~42% Oil / ~54% Meal)",
            seed="canola",
            bushel_weight_lb=50.0,
            oil_yield_pct=42.0,
            meal_yield_pct=54.0,
            oil_density_kg_per_m3=915.0,
            reference_temp_c=25.0,
            bulk_seed_density_kg_per_m3=670.0,
        )

    @classmethod
    def crude_palm_oil_cpo(cls) -> OilseedCrushAssay:
        return cls(
            name="Crude Palm Oil (CPO, 50°C Reference)",
            seed="palm",
            bushel_weight_lb=50.0,
            oil_yield_pct=100.0,
            meal_yield_pct=0.0,
            oil_density_kg_per_m3=_CPO_DENSITY_KG_PER_M3_50C,
            reference_temp_c=50.0,
            bulk_seed_density_kg_per_m3=_CPO_DENSITY_KG_PER_M3_50C,
        )


# ---------------------------------------------------------------------------
# 3. Livestock & Dairy Specification (Carcass Yield, $/cwt, Milk Solids)
# ---------------------------------------------------------------------------


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
        return round((self.fat_pct or 0.0) + (self.protein_pct or 0.0), 4)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_livestock_dairy_context(self, name=context_name)

    def live_to_carcass_mass(self, live_mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(live_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(live_mass).__name__}")
        if self.dressing_pct is None:
            raise ValueError(f"{self.name!r} does not define dressing_pct.")
        kg_q = live_mass.to(ureg.kilogram)
        frac = self.dressing_pct / 100.0
        scale = _scale_factor(kg_q.magnitude, frac)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def carcass_to_live_mass(self, carcass_mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(carcass_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(carcass_mass).__name__}")
        if self.dressing_pct is None:
            raise ValueError(f"{self.name!r} does not define dressing_pct.")
        kg_q = carcass_mass.to(ureg.kilogram)
        frac = self.dressing_pct / 100.0
        scale = _scale_factor(kg_q.magnitude, frac)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def milk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        kg_per_m3 = self.milk_density_kg_per_m3
        scale = _scale_factor(m3_q.magnitude, kg_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.kilogram)

    def milk_mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        kg_q = mass.to(ureg.kilogram)
        kg_per_m3 = self.milk_density_kg_per_m3
        scale = _scale_factor(kg_q.magnitude, kg_per_m3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.meter ** 3)

    def milk_to_solids_kg(self, milk_qty: pint.Quantity) -> pint.Quantity:
        """Converts milk mass or volume to Milk Solids (kgMS = fat + protein)."""
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
    def usd_per_cwt_to_usd_per_mt(usd_per_cwt: float | Decimal, *, metric_cwt: bool = False) -> float:
        """Converts $/cwt (US 100 lb cwt by default, or 100 kg metric cwt) to $/metric_ton."""
        p = _validate_finite_float(usd_per_cwt, "usd_per_cwt")
        if metric_cwt:
            return round(p * (_CWT_METRIC_KG / 1000.0) * 1000.0 / _CWT_METRIC_KG * 10.0, 4)
        return round(p * (_LB_PER_MT / _CWT_US_LB), 4)

    @staticmethod
    def usd_per_mt_to_usd_per_cwt(usd_per_mt: float | Decimal, *, metric_cwt: bool = False) -> float:
        """Converts $/metric_ton to $/cwt."""
        p = _validate_finite_float(usd_per_mt, "usd_per_mt")
        if metric_cwt:
            return round(p / 10.0, 4)
        return round(p / (_LB_PER_MT / _CWT_US_LB), 4)

    @classmethod
    def live_cattle_standard(cls) -> LivestockDairyAssay:
        return cls(name="Live Cattle (62% Dressing)", category="cattle", dressing_pct=62.0)

    @classmethod
    def live_hogs_standard(cls) -> LivestockDairyAssay:
        return cls(name="Lean / Live Hogs (74% Dressing)", category="hogs", dressing_pct=74.0)

    @classmethod
    def whole_milk_standard(cls) -> LivestockDairyAssay:
        return cls(
            name="Whole Milk (4.0% Fat / 3.3% Protein)",
            category="milk",
            fat_pct=4.0,
            protein_pct=3.3,
            milk_density_kg_per_m3=1030.0,
        )


# ---------------------------------------------------------------------------
# 4. Grain & Oilseed Bushel Weights and Moisture Shrinkage Specification
# (Corn 56 lb/bu, Wheat & Soybeans 60 lb/bu, Barley 48 lb/bu, Oats 32 lb/bu)
# ---------------------------------------------------------------------------


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
        if actual_moisture_pct is not None:
            m = _validate_finite_float(actual_moisture_pct, "actual_moisture_pct")
            if not (0.0 <= m < 100.0):
                raise ValueError(f"actual_moisture_pct must be in [0, 100), got {actual_moisture_pct}")
            return m
        if self.actual_moisture_pct is not None:
            return self.actual_moisture_pct
        return self.base_moisture_pct

    def _moisture_adjustment_decimal(self, actual_moisture_pct: float | Decimal | None = None) -> Decimal:
        moisture = self._effective_actual_moisture(actual_moisture_pct)
        actual = Decimal(str(actual_moisture_pct if actual_moisture_pct is not None else moisture))
        base = Decimal(str(self.base_moisture_pct))
        if actual == base or (actual < base and not self.allow_moisture_gain):
            return Decimal(1)
        ratio = (100 - actual) / (100 - base)
        return max(Decimal(0), ratio - Decimal(str(self.handling_loss_pct)) / 100)

    def moisture_adjustment_factor(self, actual_moisture_pct: float | Decimal | None = None) -> float:
        """Returns delivered-to-base mass multiplier under this contract's rules.

        Handling loss is subtracted when a moisture adjustment applies. Grain
        exactly at base moisture is unchanged. Below-base adjustments are
        disabled unless allow_moisture_gain=True. Contract presets are examples;
        handling losses and drying rates must reflect the actual contract.
        """
        return float(self._moisture_adjustment_decimal(actual_moisture_pct))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_grain_context(self, name=context_name)

    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
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
        pounds = self._validate_mass(mass).to(ureg.pound)
        weight = Decimal(str(self.bushel_weight_lb)) if isinstance(pounds.magnitude, Decimal) else self.bushel_weight_lb
        return pounds.magnitude / weight

    def bushels_to_mass(self, bushels: float | Decimal) -> pint.Quantity:
        """Returns mass from crop-weight bushel count, preserving Decimal input."""
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
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        factor = self._moisture_adjustment_decimal(actual_moisture_pct)
        scale = factor if isinstance(mt_q.magnitude, Decimal) else float(factor)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def reverse_moisture_deduction(
        self,
        base_moisture_mass: pint.Quantity,
        *,
        actual_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts base-moisture-equivalent mass back to pre-shrinkage mass at actual moisture."""
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
        tonnes = self._validate_mass(mass)
        actual = self._effective_actual_moisture(actual_moisture_pct)
        actual_dec = Decimal(str(actual_moisture_pct if actual_moisture_pct is not None else actual))
        excess = max(Decimal(0), actual_dec - Decimal(str(self.base_moisture_pct)))
        rate = Decimal(str(self.drying_cost_per_mt_per_pct if charge_per_mt_per_pct is None else charge_per_mt_per_pct))
        if not rate.is_finite() or rate < 0:
            raise ValueError("Drying charge rate must be finite and non-negative.")
        charge = (Decimal(str(tonnes.magnitude)) * excess * rate).quantize(Decimal("0.0001"))
        return charge if isinstance(tonnes.magnitude, Decimal) else float(charge)

    def adjust_price_for_moisture(
        self,
        price_per_mt: float | Decimal,
        *,
        actual_moisture_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        """Returns base-moisture price scaled per delivered tonne (four decimals)."""
        price = Decimal(str(price_per_mt))
        if not price.is_finite():
            raise ValueError("price_per_mt must be finite.")
        adjusted = (price * self._moisture_adjustment_decimal(actual_moisture_pct)).quantize(Decimal("0.0001"))
        return adjusted if isinstance(price_per_mt, Decimal) else float(adjusted)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
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
        return cls(
            name="Yellow Corn / Maize (56 lb/bu, 15.5% Base Moisture)",
            grain="corn",
            bushel_weight_lb=56.0,
            base_moisture_pct=15.5,
            actual_moisture_pct=actual_moisture_pct,
            bulk_density_kg_per_m3=720.0,
        )

    @classmethod
    def wheat_60lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        return cls(
            name="Wheat (60 lb/bu, 13.5% Base Moisture)",
            grain="wheat",
            bushel_weight_lb=60.0,
            base_moisture_pct=13.5,
            actual_moisture_pct=actual_moisture_pct,
            bulk_density_kg_per_m3=770.0,
        )

    @classmethod
    def soybeans_60lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        return cls(
            name="Soybeans (60 lb/bu, 13.0% Base Moisture)",
            grain="soybeans",
            bushel_weight_lb=60.0,
            base_moisture_pct=13.0,
            actual_moisture_pct=actual_moisture_pct,
            bulk_density_kg_per_m3=770.0,
        )

    @classmethod
    def barley_48lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        return cls(
            name="Barley (48 lb/bu, 14.5% Base Moisture)",
            grain="barley",
            bushel_weight_lb=48.0,
            base_moisture_pct=14.5,
            actual_moisture_pct=actual_moisture_pct,
            bulk_density_kg_per_m3=620.0,
        )

    @classmethod
    def oats_32lb(cls, *, actual_moisture_pct: float | None = None) -> GrainAssay:
        return cls(
            name="Oats (32 lb/bu, 13.5% Base Moisture)",
            grain="oats",
            bushel_weight_lb=32.0,
            base_moisture_pct=13.5,
            actual_moisture_pct=actual_moisture_pct,
            bulk_density_kg_per_m3=410.0,
        )


# ---------------------------------------------------------------------------
# 5. Seafood Yields, Glaze & Marine Reduction Specification
# ---------------------------------------------------------------------------


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
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        kilograms = self._validate_mass(mass, ureg.kilogram)
        scale = _scale_factor(kilograms.magnitude, self.bulk_density_kg_per_m3)
        return ureg.Quantity(kilograms.magnitude / scale, ureg.meter ** 3)

    def whole_to_dressed(self, whole_mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        q = self._validate_mass(whole_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.dressed_yield_pct)
        return ureg.Quantity(q.magnitude * scale, ureg.kilogram)

    def dressed_to_whole(self, dressed_mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(dressed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(dressed_mass).__name__}")
        q = self._validate_mass(dressed_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.dressed_yield_pct)
        return ureg.Quantity(q.magnitude / scale, ureg.kilogram)

    def whole_to_edible(self, whole_mass: pint.Quantity) -> pint.Quantity:
        """Converts whole live/round mass to fillet, picked meat, or P&D shrimp mass."""
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        q = self._validate_mass(whole_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.edible_yield_pct)
        return ureg.Quantity(q.magnitude * scale, ureg.kilogram)

    def edible_to_whole(self, edible_mass: pint.Quantity) -> pint.Quantity:
        """Converts edible/fillet/picked mass back to whole live/round mass."""
        if not isinstance(edible_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(edible_mass).__name__}")
        q = self._validate_mass(edible_mass, ureg.kilogram)
        scale = self._yield_scale(q.magnitude, self.edible_yield_pct)
        return ureg.Quantity(q.magnitude / scale, ureg.kilogram)

    def net_to_glazed_mass(self, net_mass: pint.Quantity) -> pint.Quantity:
        """Converts net deglazed/drained mass to gross glazed frozen mass."""
        if not isinstance(net_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(net_mass).__name__}")
        q = self._validate_mass(net_mass, ureg.kilogram)
        scale = self._glaze_scale(q.magnitude)
        return ureg.Quantity(q.magnitude / scale, ureg.kilogram)

    def glazed_to_net_mass(self, glazed_mass: pint.Quantity) -> pint.Quantity:
        """Converts gross glazed frozen mass to net deglazed/drained mass."""
        if not isinstance(glazed_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(glazed_mass).__name__}")
        q = self._validate_mass(glazed_mass, ureg.kilogram)
        scale = self._glaze_scale(q.magnitude)
        return ureg.Quantity(q.magnitude * scale, ureg.kilogram)

    def whole_to_reduction_yields(self, whole_mass: pint.Quantity) -> dict[str, pint.Quantity]:
        """Returns meal (dressed_yield_pct) and marine oil (oil_yield_pct) in metric tonnes."""
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        if self.species_group != "reduction":
            raise ValueError("whole_to_reduction_yields requires species_group='reduction'.")
        mt_q = self._validate_mass(whole_mass, ureg.metric_ton)
        meal_scale = self._yield_scale(mt_q.magnitude, self.dressed_yield_pct)
        oil_scale = self._yield_scale(mt_q.magnitude, self.oil_yield_pct)
        return {
            "meal": ureg.Quantity(mt_q.magnitude * meal_scale, ureg.metric_ton),
            "oil": ureg.Quantity(mt_q.magnitude * oil_scale, ureg.metric_ton),
        }

    def piece_mass_grams(self) -> float:
        """Nominal single-piece weight in grams from `count_per_lb`."""
        if self.count_per_lb is None:
            raise ValueError(f"{self.name!r} does not define count_per_lb.")
        return round(453.59237 / self.count_per_lb, 4)

    def whole_price_to_edible_price(self, whole_price: float | Decimal) -> float | Decimal:
        """Converts $/unit whole price to implied $/unit edible/fillet/meat parity."""
        px = Decimal(str(whole_price))
        if not px.is_finite():
            raise ValueError("whole_price must be finite.")
        result = (px / (Decimal(str(self.edible_yield_pct)) / 100)).quantize(Decimal("0.0001"))
        return result if isinstance(whole_price, Decimal) else float(result)

    @classmethod
    def atlantic_salmon(cls, *, glaze_pct: float = 0.0) -> SeafoodAssay:
        return cls(
            name="Atlantic Salmon (HOG 84% / Skinless Fillet 48%)",
            species_group="finfish",
            dressed_yield_pct=84.0,
            edible_yield_pct=48.0,
            glaze_pct=glaze_pct,
        )

    @classmethod
    def american_lobster(cls, *, glaze_pct: float = 0.0) -> SeafoodAssay:
        return cls(
            name="American Lobster (Cooked 90% / Picked Meat 23%)",
            species_group="lobster",
            dressed_yield_pct=90.0,
            edible_yield_pct=23.0,
            glaze_pct=glaze_pct,
        )

    @classmethod
    def snow_crab(cls, *, glaze_pct: float = 8.0) -> SeafoodAssay:
        return cls(
            name="Snow Crab (Clusters 62% / Picked Meat 21%)",
            species_group="crab",
            dressed_yield_pct=62.0,
            edible_yield_pct=21.0,
            glaze_pct=glaze_pct,
        )

    @classmethod
    def vannamei_shrimp(cls, *, count_per_lb: float = 21.0, glaze_pct: float = 10.0) -> SeafoodAssay:
        return cls(
            name=f"Vannamei Shrimp ({count_per_lb:g} ct/lb, HLSO 65% / P&D 50%)",
            species_group="shrimp",
            dressed_yield_pct=65.0,
            edible_yield_pct=50.0,
            glaze_pct=glaze_pct,
            count_per_lb=count_per_lb,
        )

    @classmethod
    def anchoveta_reduction(cls) -> SeafoodAssay:
        return cls(
            name="Anchoveta Reduction (Meal 22.5% / Crude Oil 4.5%)",
            species_group="reduction",
            dressed_yield_pct=22.5,
            edible_yield_pct=22.5,
            oil_yield_pct=4.5,
            bulk_density_kg_per_m3=640.0,
        )


# ---------------------------------------------------------------------------
# Pulp & Lumber Specification (Chemical/Mechanical Pulp & Sawn Timber)
# (NBSK, BEK, ADMT 90% dry / BDMT 100% dry, Softwood Lumber MBF <-> m^3)
# ---------------------------------------------------------------------------

_ADMT_DRY_FRACTION: Final[float] = 0.90  # Standard ISO/TAPPI air-dry pulp: 90% bone-dry fibre, 10% moisture
_CUBIC_FEET_PER_MBF: Final[float] = 1000.0 / 12.0  # Exact: 83.33333333333333 ft^3 per MBF
_CUBIC_METRES_PER_MBF: Final[float] = 2.359737216  # Exact: 1000 * (0.3048 m)^2 * 0.0254 m


@dataclass(frozen=True, slots=True)
class PulpLumberAssay:
    """Immutable physical, moisture, and volumetric specification for pulp and lumber.

    Bridges pulp moisture bases (ADMT at 90% fibre, BDMT at 100% fibre, or custom
    moisture_pct), sawn timber volume (MBF <-> m^3 with optional dressed recovery),
    bulk stowage density, unit bale/bundle counts, and price parity.
    """

    name: str
    product_type: str  # "pulp" or "lumber"
    moisture_pct: float = 10.0  # 10.0% = standard ADMT (90% bone-dry fibre); 0.0% = BDMT
    dressed_recovery_pct: float = 100.0  # Dressed/finished volume relative to nominal MBF
    bulk_density_kg_per_m3: float | None = None
    unit_package_kg: float | None = None  # e.g. pulp bale (~250 kg) or lumber pack

    def __post_init__(self) -> None:
        if not isinstance(self.product_type, str):
            raise ValueError("product_type must be 'pulp' or 'lumber'.")
        prod = self.product_type.strip().lower()
        if prod not in {"pulp", "lumber"}:
            raise ValueError(f"product_type must be 'pulp' or 'lumber', got {self.product_type!r}")

        moist = _validate_finite_float(self.moisture_pct, "moisture_pct")
        if not (0.0 <= moist < 100.0):
            raise ValueError(f"moisture_pct must be in [0, 100), got {self.moisture_pct}")

        dressed_rec = _validate_finite_float(self.dressed_recovery_pct, "dressed_recovery_pct")
        if not (0.0 < dressed_rec <= 100.0):
            raise ValueError(f"dressed_recovery_pct must be in (0, 100], got {self.dressed_recovery_pct}")

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        pkg_kg: float | None = None
        if self.unit_package_kg is not None:
            pkg_kg = _validate_finite_float(self.unit_package_kg, "unit_package_kg")
            if pkg_kg <= 0.0:
                raise ValueError(f"unit_package_kg must be strictly positive, got {self.unit_package_kg}")

        object.__setattr__(self, "product_type", prod)
        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "dressed_recovery_pct", dressed_rec)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)
        object.__setattr__(self, "unit_package_kg", pkg_kg)

    @property
    def bone_dry_fraction(self) -> float:
        """Bone-dry fibre fraction in [0, 1] (0.90 for standard 10% moisture ADMT)."""
        return 1.0 - (self.moisture_pct / 100.0)

    @property
    def admt_factor(self) -> float:
        """Multiplier converting mass at this assay's moisture to standard ADMT (90% dry)."""
        return self.bone_dry_fraction / _ADMT_DRY_FRACTION

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_pulp_lumber_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        magnitude = Decimal(str(quantity.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Pulp/lumber mass must be finite and non-negative.")
        return quantity

    @staticmethod
    def _validated_number(value: float | Decimal, label: str, *, nonnegative: bool = True) -> Decimal:
        number = Decimal(str(value))
        if not number.is_finite() or (nonnegative and number < 0):
            raise ValueError(f"{label} must be finite" + (" and non-negative." if nonnegative else "."))
        return number

    def _dry_scale(self, magnitude: Any, *, air_dry: bool = False) -> float | Decimal:
        fraction = 1 - Decimal(str(self.moisture_pct)) / 100
        if air_dry:
            fraction /= Decimal("0.9")
        return fraction if isinstance(magnitude, Decimal) else float(fraction)

    # --- Pulp Moisture Conversions (ADMT <-> BDMT <-> Actual Moisture) ---
    def mass_to_bdmt(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts pulp/fibre mass (mt, kg, short_ton) to Bone-Dry Metric Tonnes (BDMT)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        scale = self._dry_scale(mt_q.magnitude)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.bdmt)

    def bdmt_to_mass(self, bdmt_mass: pint.Quantity) -> pint.Quantity:
        """Converts Bone-Dry Metric Tonnes (BDMT) to mass at this assay's moisture_pct."""
        if not isinstance(bdmt_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(bdmt_mass).__name__}")
        if self.bone_dry_fraction <= 0.0:
            raise ValueError("Cannot convert BDMT to mass when bone_dry_fraction is zero.")
        mt_q = self._validate_mass(bdmt_mass)
        scale = self._dry_scale(mt_q.magnitude)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def mass_to_admt(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts mass at this assay's moisture_pct to Air-Dry Metric Tonnes (90% bone-dry)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        scale = self._dry_scale(mt_q.magnitude, air_dry=True)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.admt)

    def admt_to_mass(self, admt_mass: pint.Quantity) -> pint.Quantity:
        """Converts Air-Dry Metric Tonnes (90% bone-dry) to mass at this assay's moisture_pct."""
        if not isinstance(admt_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(admt_mass).__name__}")
        if self.admt_factor <= 0.0:
            raise ValueError("Cannot convert ADMT to mass when admt_factor is zero.")
        mt_q = self._validate_mass(admt_mass)
        scale = self._dry_scale(mt_q.magnitude, air_dry=True)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    # --- Lumber Volume Conversions (MBF <-> m^3) ---
    def mbf_to_cubic_metres(self, mbf: float | Decimal, *, apply_dressed_recovery: bool = False) -> pint.Quantity:
        """Converts nominal MBF count to m^3, optionally applying finished recovery."""
        count = self._validated_number(mbf, "mbf")
        recovery = Decimal(str(self.dressed_recovery_pct)) / 100 if apply_dressed_recovery else Decimal(1)
        volume = count * Decimal(str(_CUBIC_METRES_PER_MBF)) * recovery
        magnitude = volume if isinstance(mbf, Decimal) else float(volume)
        return ureg.Quantity(magnitude, ureg.meter ** 3)

    def cubic_metres_to_mbf(self, volume: pint.Quantity, *, apply_dressed_recovery: bool = False) -> float | Decimal:
        """Returns nominal MBF count, reversing any requested dressed recovery."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        quantity = volume.to(ureg.meter ** 3)
        cubic_metres = self._validated_number(quantity.magnitude, "volume")
        recovery = Decimal(str(self.dressed_recovery_pct)) / 100 if apply_dressed_recovery else Decimal(1)
        count = cubic_metres / (Decimal(str(_CUBIC_METRES_PER_MBF)) * recovery)
        return count if isinstance(quantity.magnitude, Decimal) else float(count)

    # --- Bulk Stowage & Package Conversions ---
    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        self._validated_number(m3_q.magnitude, "volume")
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def mass_to_packages(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded bale/bundle count, preserving Decimal magnitudes."""
        if self.unit_package_kg is None:
            raise ValueError(f"{self.name!r} does not define unit_package_kg.")
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        scale = _scale_factor(kilograms.magnitude, self.unit_package_kg)
        return kilograms.magnitude / scale

    def packages_to_mass(self, count: float | Decimal) -> pint.Quantity:
        if self.unit_package_kg is None:
            raise ValueError(f"{self.name!r} does not define unit_package_kg.")
        number = self._validated_number(count, "count")
        kilograms = number * Decimal(str(self.unit_package_kg))
        magnitude = kilograms if isinstance(count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram).to(ureg.metric_ton)

    @staticmethod
    def _price_parity(price: float | Decimal, factor: Decimal, *, divide: bool) -> float | Decimal:
        number = PulpLumberAssay._validated_number(price, "price", nonnegative=False)
        value = number / factor if divide else number * factor
        rounded = value.quantize(Decimal("0.0001"))
        return rounded if isinstance(price, Decimal) else float(rounded)

    @staticmethod
    def admt_price_to_bdmt_price(price_per_admt: float | Decimal) -> float | Decimal:
        """Converts nominal $/ADMT to $/BDMT, rounded to four decimals."""
        return PulpLumberAssay._price_parity(price_per_admt, Decimal("0.9"), divide=True)

    @staticmethod
    def bdmt_price_to_admt_price(price_per_bdmt: float | Decimal) -> float | Decimal:
        """Converts nominal $/BDMT to $/ADMT, rounded to four decimals."""
        return PulpLumberAssay._price_parity(price_per_bdmt, Decimal("0.9"), divide=False)

    @staticmethod
    def mbf_price_to_m3_price(price_per_mbf: float | Decimal) -> float | Decimal:
        """Converts nominal $/MBF to $/m^3; no dressed-recovery adjustment."""
        return PulpLumberAssay._price_parity(price_per_mbf, Decimal(str(_CUBIC_METRES_PER_MBF)), divide=True)

    @staticmethod
    def m3_price_to_mbf_price(price_per_m3: float | Decimal) -> float | Decimal:
        """Converts $/m^3 to nominal $/MBF; no dressed-recovery adjustment."""
        return PulpLumberAssay._price_parity(price_per_m3, Decimal(str(_CUBIC_METRES_PER_MBF)), divide=False)

    # --- Standard Forestry Presets ---
    @classmethod
    def nbsk_softwood_pulp(cls) -> PulpLumberAssay:
        """Northern Bleached Softwood Kraft (NBSK) standard air-dry pulp (90% dry / 10% moisture)."""
        return cls(
            name="NBSK Softwood Kraft Pulp (ADMT)",
            product_type="pulp",
            moisture_pct=10.0,
            bulk_density_kg_per_m3=750.0,
            unit_package_kg=250.0,
        )

    @classmethod
    def bek_hardwood_pulp(cls) -> PulpLumberAssay:
        """Bleached Eucalyptus Kraft (BEK) standard air-dry pulp (90% dry / 10% moisture)."""
        return cls(
            name="BEK Hardwood Kraft Pulp (ADMT)",
            product_type="pulp",
            moisture_pct=10.0,
            bulk_density_kg_per_m3=780.0,
            unit_package_kg=250.0,
        )

    @classmethod
    def spf_dimensional_lumber(cls, *, dressed_recovery_pct: float = 92.0) -> PulpLumberAssay:
        """Spruce-Pine-Fir (SPF) dimensional lumber benchmark."""
        return cls(
            name="SPF Dimensional Lumber",
            product_type="lumber",
            moisture_pct=19.0,  # Kiln-dried structural lumber max moisture
            dressed_recovery_pct=dressed_recovery_pct,
            bulk_density_kg_per_m3=480.0,
        )


# ---------------------------------------------------------------------------
# Natural Rubber Specification (TSR20, RSS3, Field Latex & Cup Lump)
# ---------------------------------------------------------------------------

_KG_PER_MT: Final[float] = 1000.0
_STANDARD_TSR_BALE_KG: Final[float] = 1000.0 / 30.0  # 33.3333 kg (30 bales per metric tonne)


@dataclass(frozen=True, slots=True)
class RubberAssay:
    """Immutable specification for natural rubber, latex, cup lump, and TSR bales.

    Bridges wet coagulum/latex mass <-> dry rubber mass via Dry Rubber Content
    (DRC %), bale counts <-> mass, optional bulk/tank density, and USc/kg <-> $/mt
    pricing parity.
    """

    name: str
    grade: str  # e.g. "TSR20", "RSS3", "LATEX", "CUPLUMP"
    dry_rubber_content_pct: float = 100.0
    bale_mass_kg: float = _STANDARD_TSR_BALE_KG
    bulk_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        if not isinstance(self.grade, str):
            raise ValueError("grade must be a non-empty string.")
        grade_clean = self.grade.strip().upper()
        if not grade_clean:
            raise ValueError("grade must be a non-empty string.")

        drc = _validate_finite_float(self.dry_rubber_content_pct, "dry_rubber_content_pct")
        if not (0.0 < drc <= 100.0):
            raise ValueError(f"dry_rubber_content_pct must be in (0, 100], got {self.dry_rubber_content_pct}")

        bale_kg = _validate_finite_float(self.bale_mass_kg, "bale_mass_kg")
        if bale_kg <= 0.0:
            raise ValueError(f"bale_mass_kg must be strictly positive, got {self.bale_mass_kg}")

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        object.__setattr__(self, "grade", grade_clean)
        object.__setattr__(self, "dry_rubber_content_pct", drc)
        object.__setattr__(self, "bale_mass_kg", bale_kg)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)

    def _effective_drc_fraction(self, drc_pct: float | Decimal | None = None) -> Decimal:
        value = self.dry_rubber_content_pct if drc_pct is None else drc_pct
        percentage = Decimal(str(value))
        if not percentage.is_finite() or not (0 < percentage <= 100):
            raise ValueError(f"drc_pct must be finite and in (0, 100], got {value}")
        return percentage / 100

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        magnitude = Decimal(str(quantity.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Rubber mass must be finite and non-negative.")
        return quantity

    @property
    def drc_fraction(self) -> float:
        """Dry rubber content as a decimal fraction in (0, 1]."""
        return self.dry_rubber_content_pct / 100.0

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_rubber_context(self, name=context_name)

    # --- Wet / DRC Mass Conversions ---
    def wet_to_dry_mass(
        self,
        wet_mass: pint.Quantity,
        *,
        drc_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts wet latex, cup lump, or coagulum mass to dry rubber mass (metric_ton)."""
        if not isinstance(wet_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(wet_mass).__name__}")
        mt_q = self._validate_mass(wet_mass)
        frac = self._effective_drc_fraction(drc_pct)
        scale = frac if isinstance(mt_q.magnitude, Decimal) else float(frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def dry_to_wet_mass(
        self,
        dry_mass: pint.Quantity,
        *,
        drc_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts dry rubber mass to required wet latex/cup lump mass (metric_ton)."""
        if not isinstance(dry_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(dry_mass).__name__}")
        mt_q = self._validate_mass(dry_mass)
        frac = self._effective_drc_fraction(drc_pct)
        scale = frac if isinstance(mt_q.magnitude, Decimal) else float(frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    # --- Bale & Bulk Volume Conversions ---
    def mass_to_bales(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal bale count, preserving Decimal mass input."""
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        scale = _scale_factor(kilograms.magnitude, self.bale_mass_kg)
        return kilograms.magnitude / scale

    def bales_to_mass(self, bale_count: float | Decimal) -> pint.Quantity:
        """Converts nominal bale count to metric tonnes; fractional counts are allowed."""
        count = Decimal(str(bale_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bale_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bale_mass_kg))
        magnitude = kilograms if isinstance(bale_count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram).to(ureg.metric_ton)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        magnitude = Decimal(str(m3_q.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Rubber bulk volume must be finite and non-negative.")
        mt_per_m3 = self.bulk_density_kg_per_m3 / _KG_PER_MT
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
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
            raise ValueError("Rubber price must be finite.")
        return value

    @staticmethod
    def cents_per_kg_to_usd_per_mt(cents_per_kg: float | Decimal) -> float | Decimal:
        """Converts US cents/kg to USD/mt, rounded to four decimals."""
        value = (RubberAssay._validated_price(cents_per_kg) * 10).quantize(Decimal("0.0001"))
        return value if isinstance(cents_per_kg, Decimal) else float(value)

    @staticmethod
    def usd_per_mt_to_cents_per_kg(usd_per_mt: float | Decimal) -> float | Decimal:
        """Converts USD/mt to US cents/kg, rounded to four decimals."""
        value = (RubberAssay._validated_price(usd_per_mt) / 10).quantize(Decimal("0.0001"))
        return value if isinstance(usd_per_mt, Decimal) else float(value)

    def wet_price_to_dry_price_per_mt(
        self,
        wet_price_per_mt: float | Decimal,
        *,
        drc_pct: float | Decimal | None = None,
    ) -> float | Decimal:
        """Returns dry-rubber USD/mt parity from a wet USD/mt quote (four decimals)."""
        value = (self._validated_price(wet_price_per_mt) / self._effective_drc_fraction(drc_pct)).quantize(Decimal("0.0001"))
        return value if isinstance(wet_price_per_mt, Decimal) else float(value)

    # --- Standard Natural Rubber Presets ---
    @classmethod
    def tsr20_standard(cls) -> RubberAssay:
        """Technically Specified Rubber 20 (33.33 kg bale, 100% dry rubber basis)."""
        return cls(
            name="TSR20 Natural Rubber (33.33 kg Bale)",
            grade="TSR20",
            dry_rubber_content_pct=100.0,
            bale_mass_kg=_STANDARD_TSR_BALE_KG,
            bulk_density_kg_per_m3=920.0,
        )

    @classmethod
    def rss3_bale_standard(cls, *, bale_mass_kg: float = 111.11) -> RubberAssay:
        """Ribbed Smoked Sheet 3 (RSS3) standard export bale."""
        return cls(
            name="RSS3 Ribbed Smoked Sheet",
            grade="RSS3",
            dry_rubber_content_pct=100.0,
            bale_mass_kg=bale_mass_kg,
            bulk_density_kg_per_m3=910.0,
        )

    @classmethod
    def field_latex_60drc(cls) -> RubberAssay:
        """Concentrated field/centrifuged natural rubber latex (~60% DRC)."""
        return cls(
            name="Concentrated Latex (60% DRC)",
            grade="LATEX",
            dry_rubber_content_pct=60.0,
            bale_mass_kg=_STANDARD_TSR_BALE_KG,
            bulk_density_kg_per_m3=950.0,
        )

    @classmethod
    def cup_lump_55drc(cls) -> RubberAssay:
        """Farm-gate cup lump coagulum (~55% DRC)."""
        return cls(
            name="Cup Lump Coagulum (55% DRC)",
            grade="CUPLUMP",
            dry_rubber_content_pct=55.0,
            bale_mass_kg=_STANDARD_TSR_BALE_KG,
            bulk_density_kg_per_m3=780.0,
        )


# ---------------------------------------------------------------------------
# Agricultural Ethanol Feedstock & DDGS Co-Product Specification
# (Corn 56 lb/bu, Wheat 60 lb/bu, Sorghum 56 lb/bu, Fuel Ethanol & DDGS)
# ---------------------------------------------------------------------------

_LITRES_PER_US_GALLON: Final[float] = 3.785411784
_KG_PER_SHORT_TON: Final[float] = 907.18474
_ETHANOL_DENSITY_KG_PER_M3: Final[float] = 789.0  # Denatured fuel ethanol (~20°C)


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
        bushels_per_mt = _LB_PER_MT / self.bushel_weight_lb
        litres_per_gal = _LITRES_PER_US_GALLON
        return round(bushels_per_mt * self.ethanol_gal_per_bushel * litres_per_gal, 4)

    @property
    def ddgs_yield_pct_by_mass(self) -> float:
        """DDGS output mass as a percentage of grain feedstock mass."""
        return round((self.ddgs_lb_per_bushel / self.bushel_weight_lb) * 100.0, 4)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
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
        pounds = self._validate_quantity(feedstock_mass, ureg.pound)
        count = Decimal(str(pounds.magnitude)) / Decimal(str(self.bushel_weight_lb))
        return self._number_like(count, pounds.magnitude)

    def bushels_to_feedstock_mass(self, bushels: float | Decimal) -> pint.Quantity:
        count = self._validated_number(bushels, "bushels")
        pounds = count * Decimal(str(self.bushel_weight_lb))
        return ureg.Quantity(self._number_like(pounds, bushels), ureg.pound).to(ureg.metric_ton)

    # --- Ethanol & DDGS Output Conversions ---
    def feedstock_mass_to_ethanol_volume(self, feedstock_mass: pint.Quantity) -> pint.Quantity:
        """Returns US gallons of ethanol produced from grain feedstock mass."""
        bushels = self.feedstock_mass_to_bushels(feedstock_mass)
        gallons = Decimal(str(bushels)) * Decimal(str(self.ethanol_gal_per_bushel))
        return ureg.Quantity(self._number_like(gallons, bushels), ureg.gallon)

    def ethanol_volume_to_feedstock_mass(self, ethanol_volume: pint.Quantity) -> pint.Quantity:
        """Returns required feedstock tonnes from ethanol volume, not ethanol mass."""
        gallons = self._validate_quantity(ethanol_volume, ureg.gallon)
        bushels = Decimal(str(gallons.magnitude)) / Decimal(str(self.ethanol_gal_per_bushel))
        return self.bushels_to_feedstock_mass(self._number_like(bushels, gallons.magnitude))

    def feedstock_mass_to_ddgs_mass(self, feedstock_mass: pint.Quantity) -> pint.Quantity:
        bushels = self.feedstock_mass_to_bushels(feedstock_mass)
        pounds = Decimal(str(bushels)) * Decimal(str(self.ddgs_lb_per_bushel))
        return ureg.Quantity(self._number_like(pounds, bushels), ureg.pound).to(ureg.metric_ton)

    def ddgs_mass_to_feedstock_mass(self, ddgs_mass: pint.Quantity) -> pint.Quantity:
        if self.ddgs_lb_per_bushel == 0:
            raise ValueError("Cannot infer feedstock requirements from zero DDGS yield.")
        pounds = self._validate_quantity(ddgs_mass, ureg.pound)
        bushels = Decimal(str(pounds.magnitude)) / Decimal(str(self.ddgs_lb_per_bushel))
        return self.bushels_to_feedstock_mass(self._number_like(bushels, pounds.magnitude))

    def ethanol_volume_to_mass(self, ethanol_volume: pint.Quantity) -> pint.Quantity:
        """Returns physical ethanol mass using density, independently of grain yield."""
        volume = self._validate_quantity(ethanol_volume, ureg.meter ** 3)
        kilograms = Decimal(str(volume.magnitude)) * Decimal(str(self.ethanol_density_kg_per_m3))
        return ureg.Quantity(self._number_like(kilograms, volume.magnitude), ureg.kilogram).to(ureg.metric_ton)

    def ethanol_mass_to_volume(self, ethanol_mass: pint.Quantity) -> pint.Quantity:
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
        rounded = value.quantize(Decimal("0.0001"))
        return rounded if any(isinstance(price, Decimal) for price in prices) else float(rounded)

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
        return cls(
            name="US Corn Ethanol (56 lb/bu, 2.85 gal/bu, 17 lb DDGS/bu)",
            feedstock="corn",
            bushel_weight_lb=56.0,
            ethanol_gal_per_bushel=2.85,
            ddgs_lb_per_bushel=17.0,
        )

    @classmethod
    def wheat_ethanol(cls) -> EthanolFeedstockAssay:
        """Standard wheat ethanol benchmark (60 lb/bu, 2.70 gal/bu, 18.0 lb DDGS/bu)."""
        return cls(
            name="Wheat Ethanol (60 lb/bu, 2.70 gal/bu, 18 lb DDGS/bu)",
            feedstock="wheat",
            bushel_weight_lb=60.0,
            ethanol_gal_per_bushel=2.70,
            ddgs_lb_per_bushel=18.0,
        )


# ---------------------------------------------------------------------------
# Rice Milling & Export Specification (Paddy, Brown, Milled White & Bags)
# ---------------------------------------------------------------------------

_RICE_BUSHEL_LB: Final[float] = 45.0        # Standard US rough rice bushel weight (45 lb/bu)
_RICE_EXPORT_BAG_KG: Final[float] = 50.0   # Standard international milled rice export bag (50 kg)


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
        return create_rice_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
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
        if not isinstance(paddy_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(paddy_mass).__name__}")
        mt_q = self._validate_mass(paddy_mass)
        scale = self._percentage_scale(mt_q.magnitude, self.brown_rice_yield_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def paddy_to_milled_mass(self, paddy_mass: pint.Quantity) -> pint.Quantity:
        """Converts rough/paddy rice mass to total milled rice mass (metric_ton)."""
        if not isinstance(paddy_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(paddy_mass).__name__}")
        mt_q = self._validate_mass(paddy_mass)
        scale = self._percentage_scale(mt_q.magnitude, self.milling_yield_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def milled_to_paddy_mass(self, milled_mass: pint.Quantity) -> pint.Quantity:
        """Converts milled rice mass back to required rough/paddy rice mass (metric_ton)."""
        if not isinstance(milled_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(milled_mass).__name__}")
        mt_q = self._validate_mass(milled_mass)
        scale = self._percentage_scale(mt_q.magnitude, self.milling_yield_pct)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def brown_to_paddy_mass(self, brown_mass: pint.Quantity) -> pint.Quantity:
        """Returns required paddy mass from dehulled brown/cargo rice mass."""
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
        pounds = self._validate_mass(mass).to(ureg.pound)
        return pounds.magnitude / _scale_factor(pounds.magnitude, self.bushel_weight_lb)

    def bushels_to_mass(self, bushels: float | Decimal) -> pint.Quantity:
        count = Decimal(str(bushels))
        if not count.is_finite() or count < 0:
            raise ValueError("bushels must be finite and non-negative.")
        pounds = count * Decimal(str(self.bushel_weight_lb))
        magnitude = pounds if isinstance(bushels, Decimal) else float(pounds)
        return ureg.Quantity(magnitude, ureg.pound).to(ureg.metric_ton)

    def mass_to_bags(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal export-bag count, preserving Decimal input."""
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bag_mass_kg)

    def bags_to_mass(self, bag_count: float | Decimal) -> pint.Quantity:
        count = Decimal(str(bag_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bag_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bag_mass_kg))
        magnitude = kilograms if isinstance(bag_count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram).to(ureg.metric_ton)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
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
        rounded = value.quantize(Decimal("0.0001"))
        return rounded if isinstance(original, Decimal) else float(rounded)

    @staticmethod
    def usd_per_cwt_to_usd_per_mt(usd_per_cwt: float | Decimal) -> float | Decimal:
        """Converts US 100-lb hundredweight price to USD/metric tonne."""
        pounds_per_tonne = ureg.Quantity(Decimal("1"), ureg.metric_ton).to(ureg.pound).magnitude
        price = RiceAssay._validated_price(usd_per_cwt) * pounds_per_tonne / 100
        return RiceAssay._price_result(price, usd_per_cwt)

    @staticmethod
    def usd_per_mt_to_usd_per_cwt(usd_per_mt: float | Decimal) -> float | Decimal:
        """Converts USD/metric tonne to US 100-lb hundredweight price."""
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
        return cls(
            name="Thai White Rice 5% Broken (50 kg Bag)",
            variety="LONG_GRAIN",
            milling_yield_pct=70.0,
            head_rice_pct=95.0,
            brown_rice_yield_pct=79.0,
            bag_mass_kg=50.0,
            bulk_density_kg_per_m3=780.0,
        )

    @classmethod
    def us_cme_rough_rice(cls) -> RiceAssay:
        return cls(
            name="US CME Rough Rice (45 lb/bu, 100 lb/cwt)",
            variety="LONG_GRAIN",
            milling_yield_pct=68.0,
            head_rice_pct=82.0,
            brown_rice_yield_pct=77.0,
            bushel_weight_lb=_RICE_BUSHEL_LB,
            bag_mass_kg=50.0,
            bulk_density_kg_per_m3=580.0,
        )


# ---------------------------------------------------------------------------
# Wool & Natural Fibre Specification (Greasy, Scoured Clean, Tops & Bales)
# ---------------------------------------------------------------------------

_STANDARD_GREASY_WOOL_BALE_KG: Final[float] = 175.0  # Typical Australian/NZ greasy wool bale (~170-180 kg)


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
        return create_wool_context(self, name=context_name)

    # --- Greasy <-> Clean <-> Top Mass Conversions ---
    def greasy_to_clean_mass(
        self,
        greasy_mass: pint.Quantity,
        *,
        clean_yield_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts greasy fleece mass to clean scoured wool mass (kilogram)."""
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
        if not isinstance(clean_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(clean_mass).__name__}")
        kg_q = self._validate_mass(clean_mass)
        frac = self._effective_clean_fraction(clean_yield_pct)
        scale = frac if isinstance(kg_q.magnitude, Decimal) else float(frac)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram)

    def clean_to_top_mass(self, clean_mass: pint.Quantity) -> pint.Quantity:
        """Converts clean scoured wool mass to combed wool top mass (kilogram)."""
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
        kilograms = self._validate_mass(mass)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bale_mass_kg)

    def bales_to_mass(self, bale_count: float | Decimal) -> pint.Quantity:
        count = Decimal(str(bale_count))
        if not count.is_finite() or count < 0:
            raise ValueError("bale_count must be finite and non-negative.")
        kilograms = count * Decimal(str(self.bale_mass_kg))
        magnitude = kilograms if isinstance(bale_count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram)

    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
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
        rounded = value.quantize(Decimal("0.0001"))
        return rounded if isinstance(original, Decimal) else float(rounded)

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
        return cls(
            name="Australian Merino 19 Micron (66% Clean Yield)",
            wool_type="MERINO_FINE",
            clean_yield_pct=66.0,
            combing_top_yield_pct=91.0,
            micron=19.0,
            bale_mass_kg=_STANDARD_GREASY_WOOL_BALE_KG,
            bulk_density_kg_per_m3=320.0,
        )

    @classmethod
    def medium_merino_21micron(cls) -> WoolAssay:
        """AWEX 21-micron Merino benchmark (~65% clean yield, 175 kg bale)."""
        return cls(
            name="Merino 21 Micron Benchmark (65% Clean Yield)",
            wool_type="MERINO_MEDIUM",
            clean_yield_pct=65.0,
            combing_top_yield_pct=90.0,
            micron=21.0,
            bale_mass_kg=_STANDARD_GREASY_WOOL_BALE_KG,
            bulk_density_kg_per_m3=330.0,
        )

    @classmethod
    def crossbred_28micron(cls) -> WoolAssay:
        """Crossbred apparel/carpet wool (~28.0 micron, 72% clean yield, 180 kg bale)."""
        return cls(
            name="Crossbred Wool 28 Micron (72% Clean Yield)",
            wool_type="CROSSBRED",
            clean_yield_pct=72.0,
            combing_top_yield_pct=88.0,
            micron=28.0,
            bale_mass_kg=180.0,
            bulk_density_kg_per_m3=340.0,
        )


# ---------------------------------------------------------------------------
# Poultry & Egg Specification (Broilers, Turkeys, Shell, Liquid & Dried Eggs)
# ---------------------------------------------------------------------------

_EGGS_PER_DOZEN: Final[float] = 12.0
_DOZENS_PER_CASE: Final[float] = 30.0
_EGGS_PER_CASE: Final[float] = _EGGS_PER_DOZEN * _DOZENS_PER_CASE  # 360 eggs per standard wholesale case
_DEFAULT_LARGE_EGG_GRAMS: Final[float] = 56.7  # USDA Large egg minimum (~24 oz/dozen = 56.7 g/egg)


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
        return round((self.egg_mass_grams * _EGGS_PER_DOZEN) / 1000.0, 6)

    @property
    def case_mass_kg(self) -> float:
        """Net shell-egg mass per 30-dozen wholesale case in kilograms."""
        return round((self.egg_mass_grams * _EGGS_PER_CASE) / 1000.0, 6)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_poultry_egg_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.kilogram)
        value = Decimal(str(quantity.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Poultry/egg mass must be finite and non-negative.")
        return quantity

    def _require_eggs(self) -> None:
        if self.category != "eggs":
            raise ValueError("Egg count and processing methods require category='eggs'.")

    def _require_poultry(self) -> None:
        if self.category not in {"broiler", "turkey"}:
            raise ValueError("Poultry yield methods require a broiler or turkey assay.")

    def _yield_mass(self, mass: pint.Quantity, percentage: float, *, reverse: bool = False) -> pint.Quantity:
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
        rounded = value.quantize(Decimal("0.0001"))
        return rounded if isinstance(original, Decimal) else float(rounded)

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
        self._require_eggs()
        liquid = self.shell_to_liquid_egg_mass(self.dozens_to_shell_mass(Decimal("1")))
        pounds = liquid.to(ureg.pound).magnitude
        return self._price_result(self._validated_price(usd_per_dozen) / pounds, usd_per_dozen)

    def liquid_usd_per_lb_to_usd_per_dozen(self, usd_per_lb: float | Decimal) -> float | Decimal:
        self._require_eggs()
        liquid = self.shell_to_liquid_egg_mass(self.dozens_to_shell_mass(Decimal("1")))
        pounds = liquid.to(ureg.pound).magnitude
        return self._price_result(self._validated_price(usd_per_lb) * pounds, usd_per_lb)

    # --- Standard Poultry & Egg Presets ---
    @classmethod
    def broiler_chicken_standard(cls) -> PoultryEggAssay:
        return cls(
            name="Broiler Chicken (74% RTC / 52% Boneless)",
            category="broiler",
            rtc_dressing_pct=74.0,
            boneless_meat_pct_of_live=52.0,
        )

    @classmethod
    def turkey_standard(cls) -> PoultryEggAssay:
        return cls(
            name="Commercial Turkey (79% RTC / 58% Boneless)",
            category="turkey",
            rtc_dressing_pct=79.0,
            boneless_meat_pct_of_live=58.0,
        )

    @classmethod
    def us_large_shell_eggs(cls) -> PoultryEggAssay:
        return cls(
            name="US Large Shell Eggs (56.7 g/egg, 88% Liquid / 24% Dried)",
            category="eggs",
            egg_mass_grams=_DEFAULT_LARGE_EGG_GRAMS,
            liquid_egg_yield_pct=88.0,
            dried_egg_yield_pct=24.0,
            bulk_density_kg_per_m3=1035.0,  # Liquid whole egg density
        )


# ---------------------------------------------------------------------------
# Fishmeal & Marine Oil Reduction Specification
# (Anchoveta, Menhaden, Capelin, Sardine, Meal & Crude Marine Oil)
# ---------------------------------------------------------------------------

_FISHMEAL_BAG_KG: Final[float] = 50.0          # Standard export fishmeal bag (50 kg)
_CRUDE_MARINE_OIL_KG_PER_M3: Final[float] = 925.0  # Typical crude fish/marine oil density (~920-930 kg/m^3)
_BULK_FISHMEAL_KG_PER_M3: Final[float] = 640.0     # Typical stowed bulk/bagged fishmeal density (~600-680 kg/m^3)


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
        return round(self.meal_yield_pct + self.oil_yield_pct, 4)

    def as_context(self, *, context_name: str | None = None, product: str = "meal") -> pint.Context:
        return create_fishmeal_marine_oil_context(self, name=context_name, product=product)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
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
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        mt_q = self._validate_mass(whole_mass)
        scale = self._ratio_scale(mt_q.magnitude, self.meal_yield_pct, 100)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def meal_to_whole_mass(self, meal_mass: pint.Quantity) -> pint.Quantity:
        """Converts fishmeal mass back to required whole raw fish catch (metric_ton)."""
        if not isinstance(meal_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(meal_mass).__name__}")
        mt_q = self._validate_mass(meal_mass)
        scale = self._ratio_scale(mt_q.magnitude, self.meal_yield_pct, 100)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def whole_to_oil_mass(self, whole_mass: pint.Quantity) -> pint.Quantity:
        """Converts whole raw fish catch to crude marine oil mass (metric_ton)."""
        if not isinstance(whole_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(whole_mass).__name__}")
        mt_q = self._validate_mass(whole_mass)
        scale = self._ratio_scale(mt_q.magnitude, self.oil_yield_pct, 100)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def oil_to_whole_mass(self, oil_mass: pint.Quantity) -> pint.Quantity:
        """Converts crude marine oil mass back to implied whole raw fish catch (metric_ton)."""
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
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = Decimal(str(self.oil_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def meal_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
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
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = Decimal(str(self.meal_bulk_density_kg_per_m3)) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def meal_mass_to_bags(self, meal_mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded nominal bag count, preserving Decimal input."""
        kilograms = self._validate_mass(meal_mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bag_mass_kg)

    def bags_to_meal_mass(self, bag_count: float | Decimal) -> pint.Quantity:
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
        meal_price = Decimal(str(meal_price_per_mt))
        oil_price = Decimal(str(oil_price_per_mt))
        if not meal_price.is_finite() or not oil_price.is_finite():
            raise ValueError("Meal and oil prices must be finite.")
        value = (Decimal(str(self.meal_yield_pct)) / 100 * meal_price + Decimal(str(self.oil_yield_pct)) / 100 * oil_price).quantize(Decimal("0.0001"))
        return value if isinstance(meal_price_per_mt, Decimal) or isinstance(oil_price_per_mt, Decimal) else float(value)

    # --- Standard Reduction Presets ---
    @classmethod
    def peruvian_anchoveta_super_prime(cls) -> FishmealMarineOilAssay:
        """Peruvian anchoveta super-prime fishmeal (~68% protein, ~23% meal, ~5% oil)."""
        return cls(
            name="Peruvian Anchoveta Super Prime (68% Protein)",
            species="ANCHOVETA",
            meal_yield_pct=23.0,
            oil_yield_pct=5.0,
            protein_pct=68.0,
        )

    @classmethod
    def menhaden_standard(cls) -> FishmealMarineOilAssay:
        """Atlantic/Gulf menhaden reduction benchmark (~62% protein, ~21.5% meal, ~6.5% oil)."""
        return cls(
            name="Menhaden Reduction (62% Protein)",
            species="MENHADEN",
            meal_yield_pct=21.5,
            oil_yield_pct=6.5,
            protein_pct=62.0,
        )


# ---------------------------------------------------------------------------
# Honey & Apiculture Specification (Bulk Honey, Drums, IBCs & Beeswax)
# ---------------------------------------------------------------------------

_HONEY_EXPORT_DRUM_KG: Final[float] = 300.0  # Standard 200 L steel export honey drum (~290-300 kg net)
_HONEY_IBC_TOTE_KG: Final[float] = 1400.0    # Standard 1,000 L IBC tote (~1,400-1,420 kg net)
_HONEY_PAIL_60LB_KG: Final[float] = 27.2155422  # 27.2155 kg (US 5-gallon / 60 lb pail)


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
        actual = self._validate_moisture(self.moisture_pct if actual_moisture_pct is None else actual_moisture_pct)
        return (100 - actual) / (100 - Decimal(str(self.base_moisture_pct)))

    @property
    def solids_normalisation_factor(self) -> float:
        return float(self._solids_factor())

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        value = Decimal(str(quantity.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Honey mass must be finite and non-negative.")
        return quantity

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_honey_apiculture_context(self, name=context_name)

    # --- Tank Volume <-> Mass ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
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
        quantity = self._validate_mass(mass)
        density = Decimal(str(self.density_kg_per_m3)) / 1000
        factor = density if isinstance(quantity.magnitude, Decimal) else float(density)
        return ureg.Quantity(quantity.magnitude / factor, ureg.meter ** 3)

    # --- Solids Normalisation & Beeswax By-Product ---
    def normalise_to_base_moisture(
        self, honey_mass: pint.Quantity, *, actual_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Returns constant-solids equivalent mass; does not alter stored density."""
        quantity = self._validate_mass(honey_mass)
        factor = self._solids_factor(actual_moisture_pct)
        scale = factor if isinstance(quantity.magnitude, Decimal) else float(factor)
        return ureg.Quantity(quantity.magnitude * scale, ureg.metric_ton)

    def reverse_moisture_normalisation(
        self, base_mass: pint.Quantity, *, actual_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        quantity = self._validate_mass(base_mass)
        factor = self._solids_factor(actual_moisture_pct)
        scale = factor if isinstance(quantity.magnitude, Decimal) else float(factor)
        return ureg.Quantity(quantity.magnitude / scale, ureg.metric_ton)

    def honey_to_beeswax_mass(self, honey_mass: pint.Quantity) -> pint.Quantity:
        """Estimates associated cappings wax; wax is not deducted from honey mass."""
        kilograms = self._validate_mass(honey_mass).to(ureg.kilogram)
        ratio = Decimal(str(self.beeswax_ratio_pct_of_honey)) / 100
        scale = ratio if isinstance(kilograms.magnitude, Decimal) else float(ratio)
        return ureg.Quantity(kilograms.magnitude * scale, ureg.kilogram)

    # --- Net Honey Mass & Nominal Packaging Counts ---
    def _mass_to_count(self, mass: pint.Quantity, package_kg: float) -> float | Decimal:
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, package_kg)

    @staticmethod
    def _count_to_mass(count: float | Decimal, package_kg: float) -> pint.Quantity:
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
        return self._mass_to_count(mass, _HONEY_IBC_TOTE_KG)

    def ibc_totes_to_mass(self, tote_count: float | Decimal) -> pint.Quantity:
        return self._count_to_mass(tote_count, _HONEY_IBC_TOTE_KG)

    def mass_to_pails(self, mass: pint.Quantity) -> float | Decimal:
        return self._mass_to_count(mass, _HONEY_PAIL_60LB_KG)

    def pails_to_mass(self, pail_count: float | Decimal) -> pint.Quantity:
        return self._count_to_mass(pail_count, _HONEY_PAIL_60LB_KG)

    # --- Pricing Parity ---
    @staticmethod
    def _price_parity(price: float | Decimal, *, per_pound_to_tonne: bool) -> float | Decimal:
        value = Decimal(str(price))
        if not value.is_finite():
            raise ValueError("Honey price must be finite.")
        pounds_per_tonne = Decimal(1000) / Decimal("0.45359237")
        output = value * pounds_per_tonne if per_pound_to_tonne else value / pounds_per_tonne
        rounded = output.quantize(Decimal("0.0001"))
        return rounded if isinstance(price, Decimal) else float(rounded)

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
        return cls(
            name="Bulk Export Polyfloral Honey (17.5% Moisture, 300 kg Drum)",
            floral_source="POLYFLORAL",
            moisture_pct=17.5,
            base_moisture_pct=18.0,
            density_kg_per_m3=1420.0,
            drum_mass_kg=_HONEY_EXPORT_DRUM_KG,
        )

    @classmethod
    def us_grade_a_clover(cls) -> HoneyApicultureAssay:
        """US Grade A clover honey (18.0% moisture, 1,417 kg/m^3, 60 lb pail option)."""
        return cls(
            name="US Grade A Clover Honey (18.0% Moisture)",
            floral_source="CLOVER",
            moisture_pct=18.0,
            base_moisture_pct=18.6,
            density_kg_per_m3=1417.0,
            drum_mass_kg=290.0,
        )


# ---------------------------------------------------------------------------
# Olive Oil Specification (EVOO, Virgin, Lampante, Refined & Pomace)
# ---------------------------------------------------------------------------

_OLIVE_OIL_DENSITY_KG_PER_M3_20C: Final[float] = 915.0  # Typical olive oil density at 20°C (~913-916 kg/m^3)
_OLIVE_OIL_EXPANSION_PER_C: Final[float] = 0.00072    # Approximate volumetric expansion coefficient (/°C)


@dataclass(frozen=True, slots=True)
class OliveOilAssay:
    """Immutable specification for bulk olive oil grades, thermal density, and pricing.

    Records a declared grade and free acidity, with acidity-only screening for
    EVOO, virgin and refined oil. Grade certification also requires other
    chemical/sensory criteria, not modelled here. Includes explicit refining
    yield, approximate linear thermal correction, and pricing parity.
    Density is specified at reference_temp_c; the thermal model is limited
    to observed temperatures in [-20, 100] degrees C.
    """

    name: str
    grade: str  # "EVOO", "VIRGIN", "LAMPANTE", "REFINED", "POMACE"
    free_fatty_acid_pct: float = 0.4
    refining_yield_pct: float = 100.0
    density_kg_per_m3: float = _OLIVE_OIL_DENSITY_KG_PER_M3_20C
    reference_temp_c: float = 20.0

    def __post_init__(self) -> None:
        if not isinstance(self.grade, str):
            raise ValueError("grade must be a supported olive-oil grade string.")
        grd = self.grade.strip().upper()
        if grd not in {"EVOO", "VIRGIN", "LAMPANTE", "REFINED", "POMACE"}:
            raise ValueError(
                f"grade must be 'EVOO', 'VIRGIN', 'LAMPANTE', 'REFINED', or 'POMACE', got {self.grade!r}"
            )

        ffa = _validate_finite_float(self.free_fatty_acid_pct, "free_fatty_acid_pct")
        if not (0.0 <= ffa <= 30.0):
            raise ValueError(f"free_fatty_acid_pct must be in [0.0, 30.0], got {self.free_fatty_acid_pct}")

        ref_y = _validate_finite_float(self.refining_yield_pct, "refining_yield_pct")
        if not (0.0 < ref_y <= 100.0):
            raise ValueError(f"refining_yield_pct must be in (0.0, 100.0], got {self.refining_yield_pct}")

        rho = _validate_finite_float(self.density_kg_per_m3, "density_kg_per_m3")
        if not (850.0 <= rho <= 980.0):
            raise ValueError(f"density_kg_per_m3 must be in [850, 980], got {self.density_kg_per_m3}")

        ref_t = _validate_finite_float(self.reference_temp_c, "reference_temp_c")
        if not (-20.0 <= ref_t <= 100.0):
            raise ValueError(f"reference_temp_c must be in [-20, 100], got {self.reference_temp_c}")

        object.__setattr__(self, "grade", grd)
        object.__setattr__(self, "free_fatty_acid_pct", ffa)
        object.__setattr__(self, "refining_yield_pct", ref_y)
        object.__setattr__(self, "density_kg_per_m3", rho)
        object.__setattr__(self, "reference_temp_c", ref_t)

    @property
    def applicable_ffa_limit_pct(self) -> float | None:
        """Acidity ceiling for explicitly modelled grades; None for broad categories.

        https://eur-lex.europa.eu/legal-content/EN/TXT/HTML/?uri=CELEX:32022R2104
        Lampante cannot be identified from acidity alone; pomace subgrades have
        different limits. Neither category is classified by this property.
        """
        return {"EVOO": 0.8, "VIRGIN": 2.0, "REFINED": 0.3}.get(self.grade)

    @property
    def meets_grade_acidity_limit(self) -> bool | None:
        """Acidity-only screen, not a full grade compliance determination."""
        limit = self.applicable_ffa_limit_pct
        return None if limit is None else self.free_fatty_acid_pct <= limit

    def volume_correction_factor(
        self, observed_temp: float | Decimal | pint.Quantity,
    ) -> Decimal:
        """Returns linear observed-to-reference VCF within [-20, 100] degrees C."""
        value = observed_temp.to(ureg.degC).magnitude if isinstance(observed_temp, pint.Quantity) else observed_temp
        temperature = Decimal(str(value))
        if not temperature.is_finite() or not (-20 <= temperature <= 100):
            raise ValueError("Observed temperature must be finite and in [-20, 100] degrees C.")
        delta = temperature - Decimal(str(self.reference_temp_c))
        return 1 - Decimal(str(_OLIVE_OIL_EXPANSION_PER_C)) * delta

    def _effective_density_kg_per_m3(
        self, temp_c: float | Decimal | pint.Quantity | None = None,
    ) -> Decimal:
        factor = Decimal(1) if temp_c is None else self.volume_correction_factor(temp_c)
        return Decimal(str(self.density_kg_per_m3)) * factor

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        value = Decimal(str(quantity.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Olive oil mass must be finite and non-negative.")
        return quantity

    def as_context(
        self, *, context_name: str | None = None,
        temp_c: float | Decimal | pint.Quantity | None = None,
    ) -> pint.Context:
        return create_olive_oil_context(self, name=context_name, temp_c=temp_c)

    # --- Tank Volume <-> Mass Conversions ---
    def volume_to_mass(
        self,
        volume: pint.Quantity,
        *,
        temp_c: float | Decimal | pint.Quantity | None = None,
    ) -> pint.Quantity:
        """Converts bulk olive oil volume (L, gal, m^3) to mass (metric_ton)."""
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        value = Decimal(str(m3_q.magnitude))
        if not value.is_finite() or value < 0:
            raise ValueError("Olive oil volume must be finite and non-negative.")
        mt_per_m3 = self._effective_density_kg_per_m3(temp_c) / 1000
        scale = mt_per_m3 if isinstance(m3_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(
        self,
        mass: pint.Quantity,
        *,
        temp_c: float | Decimal | pint.Quantity | None = None,
    ) -> pint.Quantity:
        """Converts bulk olive oil mass (mt, kg, lb) to volume (m^3)."""
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self._effective_density_kg_per_m3(temp_c) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    # --- Refining Yield Conversions ---
    def crude_to_refined_mass(self, crude_mass: pint.Quantity) -> pint.Quantity:
        """Converts crude/lampante olive oil mass to refined oil mass (metric_ton)."""
        if not isinstance(crude_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(crude_mass).__name__}")
        mt_q = self._validate_mass(crude_mass)
        fraction = Decimal(str(self.refining_yield_pct)) / 100
        scale = fraction if isinstance(mt_q.magnitude, Decimal) else float(fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def refined_to_crude_mass(self, refined_mass: pint.Quantity) -> pint.Quantity:
        """Converts refined olive oil mass back to required crude/lampante mass (metric_ton)."""
        if not isinstance(refined_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(refined_mass).__name__}")
        mt_q = self._validate_mass(refined_mass)
        fraction = Decimal(str(self.refining_yield_pct)) / 100
        scale = fraction if isinstance(mt_q.magnitude, Decimal) else float(fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    # --- Pricing Parity Helpers (same currency; FX is handled by the routers) ---
    @staticmethod
    def _price_result(price: float | Decimal, factor: Decimal, *, divide: bool = False) -> float | Decimal:
        value = Decimal(str(price))
        if not value.is_finite():
            raise ValueError("Olive oil price must be finite.")
        value = value / factor if divide else value * factor
        rounded = value.quantize(Decimal("0.0001"))
        return rounded if isinstance(price, Decimal) else float(rounded)

    @staticmethod
    def price_per_100kg_to_price_per_mt(price_per_100kg: float | Decimal) -> float | Decimal:
        return OliveOilAssay._price_result(price_per_100kg, Decimal(10))

    @staticmethod
    def price_per_mt_to_price_per_100kg(price_per_mt: float | Decimal) -> float | Decimal:
        return OliveOilAssay._price_result(price_per_mt, Decimal(10), divide=True)

    def price_per_kg_to_price_per_us_gallon(
        self, price_per_kg: float | Decimal, *, temp_c: float | Decimal | pint.Quantity | None = None,
    ) -> float | Decimal:
        kilograms_per_gallon = self._effective_density_kg_per_m3(temp_c) * Decimal(str(_LITRES_PER_US_GALLON)) / 1000
        return self._price_result(price_per_kg, kilograms_per_gallon)

    def price_per_us_gallon_to_price_per_kg(
        self, price_per_gallon: float | Decimal, *, temp_c: float | Decimal | pint.Quantity | None = None,
    ) -> float | Decimal:
        kilograms_per_gallon = self._effective_density_kg_per_m3(temp_c) * Decimal(str(_LITRES_PER_US_GALLON)) / 1000
        return self._price_result(price_per_gallon, kilograms_per_gallon, divide=True)

    # --- Standard Olive Oil Presets ---
    @classmethod
    def extra_virgin_evoo(cls) -> OliveOilAssay:
        """Extra Virgin Olive Oil (FFA <= 0.8%, ~915 kg/m^3 at 20°C)."""
        return cls(
            name="Extra Virgin Olive Oil (EVOO, <=0.8% FFA)",
            grade="EVOO",
            free_fatty_acid_pct=0.4,
            refining_yield_pct=100.0,
            density_kg_per_m3=_OLIVE_OIL_DENSITY_KG_PER_M3_20C,
            reference_temp_c=20.0,
        )

    @classmethod
    def lampante_crude(cls, *, free_fatty_acid_pct: float = 3.0, refining_yield_pct: float = 95.0) -> OliveOilAssay:
        """Lampante crude olive oil requiring physical/chemical refining."""
        return cls(
            name=f"Lampante Olive Oil ({free_fatty_acid_pct:g}% FFA)",
            grade="LAMPANTE",
            free_fatty_acid_pct=free_fatty_acid_pct,
            refining_yield_pct=refining_yield_pct,
            density_kg_per_m3=916.0,
            reference_temp_c=20.0,
        )

    @classmethod
    def refined_olive_oil(cls) -> OliveOilAssay:
        """Refined olive oil benchmark (low FFA, ~914 kg/m^3 at 20°C)."""
        return cls(
            name="Refined Olive Oil",
            grade="REFINED",
            free_fatty_acid_pct=0.2,
            refining_yield_pct=100.0,
            density_kg_per_m3=914.0,
            reference_temp_c=20.0,
        )


# ---------------------------------------------------------------------------
# Cotton Ginning & Seed-Lint Specification (Seed Cotton, Lint Bales & Seed)
# ---------------------------------------------------------------------------

_US_COTTON_BALE_KG: Final[float] = 217.7243376  # 217.7243376 kg (standard 480 lb US bale)
_INDIAN_COTTON_BALE_KG: Final[float] = 170.0           # Standard 170 kg Indian press bale
_AU_BR_COTTON_BALE_KG: Final[float] = 227.0            # Standard 227 kg (~500 lb) Australian/Brazilian bale


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
        return round(max(0.0, 100.0 - self.lint_turnout_pct - self.cottonseed_turnout_pct), 4)

    @property
    def seed_to_lint_ratio(self) -> float:
        """Mass ratio of recovered cottonseed per unit of ginned cotton lint."""
        return round(self.cottonseed_turnout_pct / self.lint_turnout_pct, 6)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        return create_cotton_seed_lint_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
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
        if not isinstance(seed_cotton_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_cotton_mass).__name__}")
        mt_q = self._validate_mass(seed_cotton_mass)
        scale = self._yield_scale(mt_q.magnitude, self.lint_turnout_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def lint_to_seed_cotton_mass(self, lint_mass: pint.Quantity) -> pint.Quantity:
        """Converts ginned cotton lint mass back to required raw seed cotton mass (metric_ton)."""
        if not isinstance(lint_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(lint_mass).__name__}")
        mt_q = self._validate_mass(lint_mass)
        scale = self._yield_scale(mt_q.magnitude, self.lint_turnout_pct)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def seed_cotton_to_cottonseed_mass(self, seed_cotton_mass: pint.Quantity) -> pint.Quantity:
        """Converts raw seed cotton mass to recovered whole cottonseed mass (metric_ton)."""
        if not isinstance(seed_cotton_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(seed_cotton_mass).__name__}")
        mt_q = self._validate_mass(seed_cotton_mass)
        scale = self._yield_scale(mt_q.magnitude, self.cottonseed_turnout_pct)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def cottonseed_to_seed_cotton_mass(self, cottonseed_mass: pint.Quantity) -> pint.Quantity:
        if self.cottonseed_turnout_pct == 0:
            raise ValueError("Cannot infer seed-cotton requirements from zero cottonseed turnout.")
        quantity = self._validate_mass(cottonseed_mass)
        scale = self._yield_scale(quantity.magnitude, self.cottonseed_turnout_pct)
        return ureg.Quantity(quantity.magnitude / scale, ureg.metric_ton)

    def gin_seed_cotton(self, seed_cotton_mass: pint.Quantity) -> dict[str, pint.Quantity]:
        """Returns lint, cottonseed and remaining trash without intermediate rounding."""
        quantity = self._validate_mass(seed_cotton_mass)
        lint = self.seed_cotton_to_lint_mass(quantity)
        seed = self.seed_cotton_to_cottonseed_mass(quantity)
        residual = quantity.magnitude - lint.magnitude - seed.magnitude
        if residual < 0:
            residual = Decimal(0) if isinstance(quantity.magnitude, Decimal) else 0.0
        return {"lint": lint, "cottonseed": seed, "trash": ureg.Quantity(residual, ureg.metric_ton)}

    # --- Lint Bales & Bulk Whole Cottonseed ---
    def lint_mass_to_bales(self, lint_mass: pint.Quantity) -> float | Decimal:
        kilograms = self._validate_mass(lint_mass).to(ureg.kilogram)
        return kilograms.magnitude / _scale_factor(kilograms.magnitude, self.bale_mass_kg)

    def bales_to_lint_mass(self, bale_count: float | Decimal) -> pint.Quantity:
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
        rounded = value.quantize(Decimal("0.0001"))
        return rounded if any(isinstance(value, Decimal) for value in originals) else float(rounded)

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
        return cls(
            name="US Upland Cotton (40% Lint / 54% Seed, 480 lb Bale)",
            origin="US_UPLAND",
            lint_turnout_pct=40.0,
            cottonseed_turnout_pct=54.0,
            bale_mass_kg=_US_COTTON_BALE_KG,
        )

    @classmethod
    def indian_shankar6_170kg(cls) -> CottonSeedLintAssay:
        """Indian Shankar-6 cotton (~35% lint / 61% seed turnout, 170 kg bale)."""
        return cls(
            name="Indian Shankar-6 Cotton (35% Lint / 61% Seed, 170 kg Bale)",
            origin="INDIA_SHANKAR6",
            lint_turnout_pct=35.0,
            cottonseed_turnout_pct=61.0,
            bale_mass_kg=_INDIAN_COTTON_BALE_KG,
        )

    @classmethod
    def brazil_australian_227kg(cls) -> CottonSeedLintAssay:
        """Brazilian / Australian high-turnout cotton (42% lint / 52% seed, 227 kg bale)."""
        return cls(
            name="Brazil/Australia Cotton (42% Lint / 52% Seed, 227 kg Bale)",
            origin="BRAZIL",
            lint_turnout_pct=42.0,
            cottonseed_turnout_pct=52.0,
            bale_mass_kg=_AU_BR_COTTON_BALE_KG,
        )


# ---------------------------------------------------------------------------
# 16. Pint Context Factories & High-Level Converters
# ---------------------------------------------------------------------------


def create_softs_context(assay: SoftsAssay, *, name: str | None = None) -> pint.Context:
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def create_oilseed_crush_context(
    assay: OilseedCrushAssay, *, name: str | None = None,
    temp_c: float | Decimal | pint.Quantity | None = None,
) -> pint.Context:
    if temp_c is not None:
        assay.volume_correction_factor(temp_c)
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.oil_volume_to_mass(q, temp_c=temp_c))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.oil_mass_to_volume(q, temp_c=temp_c))
    return ctx


def create_livestock_dairy_context(assay: LivestockDairyAssay, *, name: str | None = None) -> pint.Context:
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.category == "milk":
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.milk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.milk_mass_to_volume(q))
    return ctx


def convert_softs(quantity: pint.Quantity, target_unit: str | pint.Unit, *, assay: SoftsAssay, mode: str = "auto",) -> pint.Quantity:
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, SoftsAssay):
        raise TypeError(f"Expected SoftsAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if quantity.dimensionality == tgt.dimensionality:
            return quantity.to(tgt)
        with ureg.context(assay.as_context()):
            return quantity.to(tgt)
    dispatch = {"raw_to_processed": lambda q: assay.raw_to_processed_mass(q), "processed_to_raw": lambda q: assay.processed_to_raw_mass(q),}
    if m not in dispatch:
        raise ValueError(f"Unsupported softs mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def convert_oilseed_crush(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: OilseedCrushAssay,
    mode: str = "auto",
    temp_c: float | None = None,
) -> pint.Quantity:
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, OilseedCrushAssay):
        raise TypeError(f"Expected OilseedCrushAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay, observed_temp=temp_c)
    dispatch = {
        "seed_to_oil": lambda q: assay.seed_mass_to_oil_mass(q),
        "seed_to_meal": lambda q: assay.seed_mass_to_meal_mass(q),
        "oil_volume_to_mass": lambda q: assay.oil_volume_to_mass(q, temp_c=temp_c),
        "oil_mass_to_volume": lambda q: assay.oil_mass_to_volume(q, temp_c=temp_c),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported crush mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def convert_livestock_dairy(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: LivestockDairyAssay,
    mode: str = "auto",
) -> pint.Quantity:
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, LivestockDairyAssay):
        raise TypeError(f"Expected LivestockDairyAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if quantity.dimensionality == tgt.dimensionality:
            return quantity.to(tgt)
        with ureg.context(assay.as_context()):
            return quantity.to(tgt)
    dispatch = {
        "live_to_carcass": lambda q: assay.live_to_carcass_mass(q),
        "carcass_to_live": lambda q: assay.carcass_to_live_mass(q),
        "milk_to_solids": lambda q: assay.milk_to_solids_kg(q),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported livestock/dairy mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)

def create_grain_context(assay: GrainAssay, *, name: str | None = None) -> pint.Context:
    """Bridges actual bulk volume and mass when density is supplied.

    Use explicit bushel-count methods or registered crop-weight bushel units;
    generic Pint bushel remains a volume measure and uses the bulk density.
    """
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def convert_grain(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: GrainAssay,
    mode: str = "auto",
    actual_moisture_pct: float | Decimal | None = None,
) -> pint.Quantity:
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, GrainAssay):
        raise TypeError(f"Expected GrainAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)
    dispatch = {
        "apply_moisture": lambda q: assay.apply_moisture_deduction(q, actual_moisture_pct=actual_moisture_pct),
        "reverse_moisture": lambda q: assay.reverse_moisture_deduction(q, actual_moisture_pct=actual_moisture_pct),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported grain mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def create_seafood_context(assay: SeafoodAssay, *, name: str | None = None) -> pint.Context:
    """Bridges bulk volume and gross mass; yields and glaze use explicit methods."""
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def convert_seafood(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: SeafoodAssay,
    mode: str = "auto",
) -> pint.Quantity:
    """Converts physical units or an explicit seafood processing/glaze mode."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, SeafoodAssay):
        raise TypeError(f"Expected SeafoodAssay, got {type(assay).__name__}")
    target = parse_unit(target_unit)
    key = mode.strip().lower()
    if key == "auto":
        return convert_units(quantity, target, assay=assay)
    dispatch = {
        "whole_to_dressed": assay.whole_to_dressed,
        "dressed_to_whole": assay.dressed_to_whole,
        "whole_to_edible": assay.whole_to_edible,
        "edible_to_whole": assay.edible_to_whole,
        "net_to_glazed": assay.net_to_glazed_mass,
        "glazed_to_net": assay.glazed_to_net_mass,
        "whole_to_meal": lambda q: assay.whole_to_reduction_yields(q)["meal"],
        "whole_to_oil": lambda q: assay.whole_to_reduction_yields(q)["oil"],
    }
    if key not in dispatch:
        raise ValueError(f"Unsupported seafood mode {mode!r}. Expected {sorted(dispatch)} or 'auto'.")
    return dispatch[key](quantity).to(target)



def create_pulp_lumber_context(assay: PulpLumberAssay, *, name: str | None = None) -> pint.Context:
    """Bridges [volume] <-> [mass] for PulpLumberAssay when bulk density is defined."""
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def convert_pulp_lumber(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: PulpLumberAssay,
    mode: str = "auto",
) -> pint.Quantity:
    """Converts pulp or lumber quantities across volume, mass, BDMT, and ADMT."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, PulpLumberAssay):
        raise TypeError(f"Expected PulpLumberAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "to_bdmt": lambda q: assay.mass_to_bdmt(q),
        "from_bdmt": lambda q: assay.bdmt_to_mass(q),
        "to_admt": lambda q: assay.mass_to_admt(q),
        "from_admt": lambda q: assay.admt_to_mass(q),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported pulp/lumber mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def create_rubber_context(assay: RubberAssay, *, name: str | None = None) -> pint.Context:
    """Bridges [volume] <-> [mass] for RubberAssay when bulk density is defined."""
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def convert_rubber(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: RubberAssay,
    mode: str = "auto",
    drc_pct: float | Decimal | None = None,
) -> pint.Quantity:
    """Converts rubber physical units or explicit wet/dry processing modes.

    auto changes units only; drc_pct applies to wet_to_dry and dry_to_wet.
    Bale weights are nominal specification inputs, not a moisture deduction.
    """
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, RubberAssay):
        raise TypeError(f"Expected RubberAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "wet_to_dry": lambda q: assay.wet_to_dry_mass(q, drc_pct=drc_pct),
        "dry_to_wet": lambda q: assay.dry_to_wet_mass(q, drc_pct=drc_pct),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported rubber mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def create_ethanol_feedstock_context(
    assay: EthanolFeedstockAssay,
    *,
    name: str | None = None,
) -> pint.Context:
    """Bridges grain feedstock mass and ethanol output volume using process yield.

    Context mass always means grain feedstock, never ethanol or DDGS product.
    Use explicit converter modes for product mass/density conversions.
    """
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.feedstock_mass_to_ethanol_volume(q))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.ethanol_volume_to_feedstock_mass(q))
    return ctx


def convert_ethanol_feedstock(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: EthanolFeedstockAssay,
    mode: str = "auto",
) -> pint.Quantity:
    """Converts feedstock, ethanol, and DDGS quantities across mass and volume."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, EthanolFeedstockAssay):
        raise TypeError(f"Expected EthanolFeedstockAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "feedstock_to_ethanol": lambda q: assay.feedstock_mass_to_ethanol_volume(q),
        "ethanol_to_feedstock": lambda q: assay.ethanol_volume_to_feedstock_mass(q),
        "feedstock_to_ddgs": lambda q: assay.feedstock_mass_to_ddgs_mass(q),
        "ddgs_to_feedstock": lambda q: assay.ddgs_mass_to_feedstock_mass(q),
        "ethanol_vol_to_mass": lambda q: assay.ethanol_volume_to_mass(q),
        "ethanol_mass_to_vol": lambda q: assay.ethanol_mass_to_volume(q),
    }
    if m not in dispatch:
        raise ValueError(
            f"Unsupported ethanol feedstock mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    return dispatch[m](quantity).to(tgt)


def create_rice_context(assay: RiceAssay, *, name: str | None = None) -> pint.Context:
    """Bridges [volume] <-> [mass] for RiceAssay when bulk density is defined."""
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def convert_rice(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: RiceAssay,
    mode: str = "auto",
) -> pint.Quantity:
    """Converts rice quantities across volume, mass, paddy, brown, and milled basis."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, RiceAssay):
        raise TypeError(f"Expected RiceAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "paddy_to_brown": lambda q: assay.paddy_to_brown_mass(q),
        "brown_to_paddy": lambda q: assay.brown_to_paddy_mass(q),
        "brown_to_milled": lambda q: assay.brown_to_milled_mass(q),
        "milled_to_brown": lambda q: assay.milled_to_brown_mass(q),
        "paddy_to_milled": lambda q: assay.paddy_to_milled_mass(q),
        "milled_to_paddy": lambda q: assay.milled_to_paddy_mass(q),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported rice mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def create_wool_context(assay: WoolAssay, *, name: str | None = None) -> pint.Context:
    """Bridges [volume] <-> [mass] for WoolAssay when bulk density is defined."""
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def convert_wool(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: WoolAssay,
    mode: str = "auto",
    clean_yield_pct: float | Decimal | None = None,
) -> pint.Quantity:
    """Converts wool quantities across volume, greasy mass, clean scoured mass, and tops."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, WoolAssay):
        raise TypeError(f"Expected WoolAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "greasy_to_clean": lambda q: assay.greasy_to_clean_mass(q, clean_yield_pct=clean_yield_pct),
        "clean_to_greasy": lambda q: assay.clean_to_greasy_mass(q, clean_yield_pct=clean_yield_pct),
        "clean_to_top": lambda q: assay.clean_to_top_mass(q),
        "top_to_clean": lambda q: assay.top_to_clean_mass(q),
        "top_to_greasy": lambda q: assay.top_to_greasy_mass(q, clean_yield_pct=clean_yield_pct),
        "greasy_to_top": lambda q: assay.greasy_to_top_mass(q, clean_yield_pct=clean_yield_pct),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported wool mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def create_poultry_egg_context(assay: PoultryEggAssay, *, name: str | None = None) -> pint.Context:
    """Bridges volume/mass for the density's product state (liquid egg in its preset).

    Counts, carcass and egg processing yields require explicit methods/modes.
    """
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    if assay.bulk_density_kg_per_m3 is not None:
        ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.bulk_volume_to_mass(q))
        ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_bulk_volume(q))
    return ctx


def convert_poultry_egg(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: PoultryEggAssay,
    mode: str = "auto",
) -> pint.Quantity:
    """Converts poultry and egg quantities across volume, live, RTC, boneless, shell, liquid, and dried mass."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, PoultryEggAssay):
        raise TypeError(f"Expected PoultryEggAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "live_to_rtc": lambda q: assay.live_to_rtc_mass(q),
        "rtc_to_live": lambda q: assay.rtc_to_live_mass(q),
        "live_to_boneless": lambda q: assay.live_to_boneless_mass(q),
        "boneless_to_live": lambda q: assay.boneless_to_live_mass(q),
        "liquid_to_shell": lambda q: assay.liquid_to_shell_egg_mass(q),
        "dried_to_shell": lambda q: assay.dried_to_shell_egg_mass(q),
        "shell_to_liquid": lambda q: assay.shell_to_liquid_egg_mass(q),
        "shell_to_dried": lambda q: assay.shell_to_dried_egg_mass(q),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported poultry/egg mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def create_fishmeal_marine_oil_context(
    assay: FishmealMarineOilAssay,
    *,
    name: str | None = None,
    product: str = "meal",
) -> pint.Context:
    """Bridges [volume] <-> [mass] for either fishmeal ('meal') or marine oil ('oil')."""
    if not isinstance(product, str):
        raise ValueError("product must be 'meal' or 'oil'.")
    prod = product.strip().lower()
    if prod == "meal":
        to_mass = assay.meal_volume_to_mass
        to_vol = assay.meal_mass_to_volume
    elif prod == "oil":
        to_mass = assay.oil_volume_to_mass
        to_vol = assay.oil_mass_to_volume
    else:
        raise ValueError(f"Unsupported product {product!r}; expected 'meal' or 'oil'.")

    ctx = pint.Context(name or f"{assay.name.strip().lower().replace(' ', '_')}_{prod}")
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: to_vol(q))
    return ctx


def convert_fishmeal_marine_oil(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: FishmealMarineOilAssay,
    mode: str = "auto",
    product: str = "meal",
) -> pint.Quantity:
    """Converts whole catch, fishmeal, and crude marine oil across volume and mass."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, FishmealMarineOilAssay):
        raise TypeError(f"Expected FishmealMarineOilAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay, product=product)

    dispatch = {
        "whole_to_meal": lambda q: assay.whole_to_meal_mass(q),
        "meal_to_whole": lambda q: assay.meal_to_whole_mass(q),
        "whole_to_oil": lambda q: assay.whole_to_oil_mass(q),
        "oil_to_whole": lambda q: assay.oil_to_whole_mass(q),
        "oil_vol_to_mass": lambda q: assay.oil_volume_to_mass(q),
        "oil_mass_to_vol": lambda q: assay.oil_mass_to_volume(q),
        "meal_vol_to_mass": lambda q: assay.meal_volume_to_mass(q),
        "meal_mass_to_vol": lambda q: assay.meal_mass_to_volume(q),
    }
    if m not in dispatch:
        raise ValueError(
            f"Unsupported fishmeal/marine oil mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    return dispatch[m](quantity).to(tgt)


def create_honey_apiculture_context(assay: HoneyApicultureAssay, *, name: str | None = None) -> pint.Context:
    """Bridges [volume] <-> [mass] for HoneyApicultureAssay via honey density."""
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q))
    return ctx


def convert_honey_apiculture(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: HoneyApicultureAssay,
    mode: str = "auto",
    actual_moisture_pct: float | Decimal | None = None,
) -> pint.Quantity:
    """Converts bulk honey and beeswax quantities across volume, mass, and moisture basis."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, HoneyApicultureAssay):
        raise TypeError(f"Expected HoneyApicultureAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "normalise_moisture": lambda q: assay.normalise_to_base_moisture(
            q, actual_moisture_pct=actual_moisture_pct
        ),
        "reverse_moisture": lambda q: assay.reverse_moisture_normalisation(q, actual_moisture_pct=actual_moisture_pct),
        "honey_to_beeswax": lambda q: assay.honey_to_beeswax_mass(q),
    }
    if m not in dispatch:
        raise ValueError(
            f"Unsupported honey/apiculture mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    return dispatch[m](quantity).to(tgt)


def create_olive_oil_context(
    assay: OliveOilAssay,
    *,
    name: str | None = None,
    temp_c: float | Decimal | pint.Quantity | None = None,
) -> pint.Context:
    """Bridges [volume] <-> [mass] for OliveOilAssay with optional thermal VCF."""
    assay._effective_density_kg_per_m3(temp_c)  # Validate temperature before activation.
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.volume_to_mass(q, temp_c=temp_c))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.mass_to_volume(q, temp_c=temp_c))
    return ctx


def convert_olive_oil(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: OliveOilAssay,
    mode: str = "auto",
    temp_c: float | Decimal | pint.Quantity | None = None,
) -> pint.Quantity:
    """Converts olive oil quantities across volume, mass, and refining yield."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, OliveOilAssay):
        raise TypeError(f"Expected OliveOilAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay, observed_temp=temp_c)

    dispatch = {
        "vol_to_mass": lambda q: assay.volume_to_mass(q, temp_c=temp_c),
        "mass_to_vol": lambda q: assay.mass_to_volume(q, temp_c=temp_c),
        "crude_to_refined": lambda q: assay.crude_to_refined_mass(q),
        "refined_to_crude": lambda q: assay.refined_to_crude_mass(q),
    }
    if m not in dispatch:
        raise ValueError(f"Unsupported olive oil mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    return dispatch[m](quantity).to(tgt)


def create_cotton_seed_lint_context(
    assay: CottonSeedLintAssay,
    *,
    name: str | None = None,
) -> pint.Context:
    """Bridges bulk whole cottonseed volume/mass, not lint or raw seed cotton.

    Ginning turnout and lint-bale conversion require explicit methods/modes.
    """
    ctx = pint.Context(name or assay.name.strip().lower().replace(" ", "_"))
    ctx.add_transformation("[volume]", "[mass]", lambda _u, q: assay.cottonseed_volume_to_mass(q))
    ctx.add_transformation("[mass]", "[volume]", lambda _u, q: assay.cottonseed_mass_to_volume(q))
    return ctx


def convert_cotton_seed_lint(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: CottonSeedLintAssay,
    mode: str = "auto",
) -> pint.Quantity:
    """Converts seed cotton, lint, and whole cottonseed across mass and bulk volume."""
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, CottonSeedLintAssay):
        raise TypeError(f"Expected CottonSeedLintAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        return convert_units(quantity, tgt, assay=assay)

    dispatch = {
        "seed_cotton_to_lint": lambda q: assay.seed_cotton_to_lint_mass(q),
        "lint_to_seed_cotton": lambda q: assay.lint_to_seed_cotton_mass(q),
        "seed_cotton_to_seed": lambda q: assay.seed_cotton_to_cottonseed_mass(q),
        "seed_to_seed_cotton": lambda q: assay.cottonseed_to_seed_cotton_mass(q),
        "seed_vol_to_mass": lambda q: assay.cottonseed_volume_to_mass(q),
        "seed_mass_to_vol": lambda q: assay.cottonseed_mass_to_volume(q),
    }
    if m not in dispatch:
        raise ValueError(
            f"Unsupported cotton seed/lint mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    return dispatch[m](quantity).to(tgt)
