"""Commodex agricultural_specs implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final
import math
import pint



_LEGACY_EXPORTS: Final[tuple[str, ...]] = (
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
# Pulp & Lumber Specification (Chemical/Mechanical Pulp & Sawn Timber)
# (NBSK, BEK, ADMT 90% dry / BDMT 100% dry, Softwood Lumber MBF <-> m^3)
# ---------------------------------------------------------------------------

_ADMT_DRY_FRACTION: Final[float] = 0.90  # Standard ISO/TAPPI air-dry pulp: 90% bone-dry fibre, 10% moisture



_CUBIC_FEET_PER_MBF: Final[float] = 1000.0 / 12.0  # Exact: 83.33333333333333 ft^3 per MBF



_CUBIC_METRES_PER_MBF: Final[float] = 2.359737216  # Exact: 1000 * (0.3048 m)^2 * 0.0254 m





# ---------------------------------------------------------------------------
# Natural Rubber Specification (TSR20, RSS3, Field Latex & Cup Lump)
# ---------------------------------------------------------------------------

_KG_PER_MT: Final[float] = 1000.0



_STANDARD_TSR_BALE_KG: Final[float] = 1000.0 / 30.0  # 33.3333 kg (30 bales per metric tonne)





# ---------------------------------------------------------------------------
# Agricultural Ethanol Feedstock & DDGS Co-Product Specification
# (Corn 56 lb/bu, Wheat 60 lb/bu, Sorghum 56 lb/bu, Fuel Ethanol & DDGS)
# ---------------------------------------------------------------------------

_LITRES_PER_US_GALLON: Final[float] = 3.785411784



_KG_PER_SHORT_TON: Final[float] = 907.18474



_ETHANOL_DENSITY_KG_PER_M3: Final[float] = 789.0  # Denatured fuel ethanol (~20°C)





# ---------------------------------------------------------------------------
# Rice Milling & Export Specification (Paddy, Brown, Milled White & Bags)
# ---------------------------------------------------------------------------

_RICE_BUSHEL_LB: Final[float] = 45.0        # Standard US rough rice bushel weight (45 lb/bu)



_RICE_EXPORT_BAG_KG: Final[float] = 50.0   # Standard international milled rice export bag (50 kg)





# ---------------------------------------------------------------------------
# Wool & Natural Fibre Specification (Greasy, Scoured Clean, Tops & Bales)
# ---------------------------------------------------------------------------

_STANDARD_GREASY_WOOL_BALE_KG: Final[float] = 175.0  # Typical Australian/NZ greasy wool bale (~170-180 kg)





# ---------------------------------------------------------------------------
# Poultry & Egg Specification (Broilers, Turkeys, Shell, Liquid & Dried Eggs)
# ---------------------------------------------------------------------------

_EGGS_PER_DOZEN: Final[float] = 12.0



_DOZENS_PER_CASE: Final[float] = 30.0



_EGGS_PER_CASE: Final[float] = _EGGS_PER_DOZEN * _DOZENS_PER_CASE  # 360 eggs per standard wholesale case



_DEFAULT_LARGE_EGG_GRAMS: Final[float] = 56.7  # USDA Large egg minimum (~24 oz/dozen = 56.7 g/egg)





# ---------------------------------------------------------------------------
# Fishmeal & Marine Oil Reduction Specification
# (Anchoveta, Menhaden, Capelin, Sardine, Meal & Crude Marine Oil)
# ---------------------------------------------------------------------------

_FISHMEAL_BAG_KG: Final[float] = 50.0          # Standard export fishmeal bag (50 kg)



_CRUDE_MARINE_OIL_KG_PER_M3: Final[float] = 925.0  # Typical crude fish/marine oil density (~920-930 kg/m^3)



_BULK_FISHMEAL_KG_PER_M3: Final[float] = 640.0     # Typical stowed bulk/bagged fishmeal density (~600-680 kg/m^3)





# ---------------------------------------------------------------------------
# Honey & Apiculture Specification (Bulk Honey, Drums, IBCs & Beeswax)
# ---------------------------------------------------------------------------

_HONEY_EXPORT_DRUM_KG: Final[float] = 300.0  # Standard 200 L steel export honey drum (~290-300 kg net)



_HONEY_IBC_TOTE_KG: Final[float] = 1400.0    # Standard 1,000 L IBC tote (~1,400-1,420 kg net)



_HONEY_PAIL_60LB_KG: Final[float] = 27.2155422  # 27.2155 kg (US 5-gallon / 60 lb pail)





# ---------------------------------------------------------------------------
# Olive Oil Specification (EVOO, Virgin, Lampante, Refined & Pomace)
# ---------------------------------------------------------------------------

_OLIVE_OIL_DENSITY_KG_PER_M3_20C: Final[float] = 915.0  # Typical olive oil density at 20°C (~913-916 kg/m^3)



_OLIVE_OIL_EXPANSION_PER_C: Final[float] = 0.00072    # Approximate volumetric expansion coefficient (/°C)


from core.units.validation import validated_assay







# ---------------------------------------------------------------------------
# Cotton Ginning & Seed-Lint Specification (Seed Cotton, Lint Bales & Seed)
# ---------------------------------------------------------------------------

_US_COTTON_BALE_KG: Final[float] = 217.7243376  # 217.7243376 kg (standard 480 lb US bale)



_INDIAN_COTTON_BALE_KG: Final[float] = 170.0           # Standard 170 kg Indian press bale



_AU_BR_COTTON_BALE_KG: Final[float] = 227.0            # Standard 227 kg (~500 lb) Australian/Brazilian bale





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


from core.units.provenance import ConversionClassification


def _finalize_agricultural_mode(
    converted: pint.Quantity,
    original: pint.Quantity,
    *,
    assay: Any,
    mode: str,
    classification: ConversionClassification = ConversionClassification.YIELD_DEPENDENT,
    return_metadata: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
    rerun: Any = None,
    extra_base_params: Mapping[str, Decimal] | None = None,
) -> Any:
    from core.units.precision import ensure_quantity_ieee754
    from core.units.provenance import build_conversion_result
    out_q = ensure_quantity_ieee754(converted)
    if not return_metadata:
        return out_q
    base_params: dict[str, Decimal] = {"quantity": Decimal(str(original.magnitude))}
    if extra_base_params:
        base_params.update(extra_base_params)
    return build_conversion_result(
        out_q,
        original,
        assay=assay,
        conversion_path=(str(original.units), mode, str(out_q.units)),
        method=f"{type(assay).__name__}.{mode}",
        classification=classification,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=rerun,
        base_parameters=base_params,
    )





def convert_softs(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: SoftsAssay,
    mode: str = "auto",
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    from core.units.assays.agriculture.softs import SoftsAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, SoftsAssay):
        raise TypeError(f"Expected SoftsAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )
    dispatch = {"raw_to_processed": lambda q, a=assay: a.raw_to_processed_mass(q), "processed_to_raw": lambda q, a=assay: a.processed_to_raw_mass(q)}
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported softs mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](
            ureg.Quantity(p["quantity"], quantity.units),
            replace(assay, processing_yield_pct=float(p["processing_yield_pct"])),
        ).to(tgt),
        extra_base_params={"processing_yield_pct": Decimal(str(assay.processing_yield_pct))},
    )





def convert_oilseed_crush(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: OilseedCrushAssay,
    mode: str = "auto",
    temp_c: float | None = None,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    from core.units.assays.agriculture.feedstocks import OilseedCrushAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, OilseedCrushAssay):
        raise TypeError(f"Expected OilseedCrushAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay, observed_temp=temp_c)
        return convert(
            quantity, tgt, assay=assay, observed_temp=temp_c,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )
    dispatch = {
        "seed_to_oil": lambda q, a=assay: a.seed_mass_to_oil_mass(q),
        "seed_to_meal": lambda q, a=assay: a.seed_mass_to_meal_mass(q),
        "oil_volume_to_mass": lambda q, a=assay: a.oil_volume_to_mass(q, temp_c=temp_c),
        "oil_mass_to_volume": lambda q, a=assay: a.oil_mass_to_volume(q, temp_c=temp_c),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported crush mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    cls_tag = ConversionClassification.YIELD_DEPENDENT if m in {"seed_to_oil", "seed_to_meal"} else ConversionClassification.ASSAY_DEPENDENT
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, classification=cls_tag, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](
            ureg.Quantity(p["quantity"], quantity.units),
            replace(assay, oil_yield_pct=float(p["oil_yield_pct"]), meal_yield_pct=float(p["meal_yield_pct"]), oil_density_kg_per_m3=float(p["oil_density_kg_per_m3"])),
        ).to(tgt),
        extra_base_params={
            "oil_yield_pct": Decimal(str(assay.oil_yield_pct)),
            "meal_yield_pct": Decimal(str(assay.meal_yield_pct)),
            "oil_density_kg_per_m3": Decimal(str(assay.oil_density_kg_per_m3)),
        },
    )





def convert_livestock_dairy(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    assay: LivestockDairyAssay,
    mode: str = "auto",
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    from core.units.assays.agriculture.livestock import LivestockDairyAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, LivestockDairyAssay):
        raise TypeError(f"Expected LivestockDairyAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )
    dispatch = {
        "live_to_carcass": lambda q, a=assay: a.live_to_carcass_mass(q),
        "carcass_to_live": lambda q, a=assay: a.carcass_to_live_mass(q),
        "milk_to_solids": lambda q, a=assay: a.milk_to_solids_kg(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported livestock/dairy mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    extra: dict[str, Decimal] = {}
    if assay.dressing_pct is not None:
        extra["dressing_pct"] = Decimal(str(assay.dressing_pct))
        extra["yield_fraction"] = Decimal(str(assay.dressing_pct)) / Decimal(100)
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m,
        classification=ConversionClassification.YIELD_DEPENDENT,
        return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](
            ureg.Quantity(p["quantity"], quantity.units),
            replace(
                assay,
                dressing_pct=(
                    float(p["yield_fraction"] * Decimal(100))
                    if "yield_fraction" in p and p["yield_fraction"] != extra.get("yield_fraction")
                    else (float(p["dressing_pct"]) if "dressing_pct" in p else assay.dressing_pct)
                ),
            ),
        ).to(tgt),
        extra_base_params=extra,
    )




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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    from core.units.assays.agriculture.grains import GrainAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, GrainAssay):
        raise TypeError(f"Expected GrainAssay, got {type(assay).__name__}")
    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )
    dispatch = {
        "apply_moisture": lambda q, mo=actual_moisture_pct: assay.apply_moisture_deduction(q, actual_moisture_pct=mo),
        "reverse_moisture": lambda q, mo=actual_moisture_pct: assay.reverse_moisture_deduction(q, actual_moisture_pct=mo),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported grain mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    eff_m = Decimal(str(actual_moisture_pct if actual_moisture_pct is not None else assay._effective_actual_moisture()))
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m,
        classification=ConversionClassification.ASSAY_DEPENDENT,
        return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units), p["actual_moisture_pct"]).to(tgt),
        extra_base_params={"actual_moisture_pct": eff_m},
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts physical units or an explicit seafood processing/glaze mode."""
    from core.units.assays.agriculture.seafood import SeafoodAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, SeafoodAssay):
        raise TypeError(f"Expected SeafoodAssay, got {type(assay).__name__}")
    target = parse_unit(target_unit)
    key = mode.strip().lower()
    if key == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, target, assay=assay)
        return convert(
            quantity, target, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )
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
        raise UnsupportedConversionError(f"Unsupported seafood mode {mode!r}. Expected {sorted(dispatch)} or 'auto'.")
    res = dispatch[key](quantity).to(target)
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=key, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[key](ureg.Quantity(p["quantity"], quantity.units)).to(target),
    )






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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts pulp or lumber quantities across volume, mass, BDMT, and ADMT."""
    from core.units.assays.agriculture.forestry import PulpLumberAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, PulpLumberAssay):
        raise TypeError(f"Expected PulpLumberAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "to_bdmt": lambda q: assay.mass_to_bdmt(q),
        "from_bdmt": lambda q: assay.bdmt_to_mass(q),
        "to_admt": lambda q: assay.mass_to_admt(q),
        "from_admt": lambda q: assay.admt_to_mass(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported pulp/lumber mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m,
        classification=ConversionClassification.ASSAY_DEPENDENT,
        return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts rubber physical units or explicit wet/dry processing modes.

    auto changes units only; drc_pct applies to wet_to_dry and dry_to_wet.
    Bale weights are nominal specification inputs, not a moisture deduction.
    """
    from core.units.assays.agriculture.rubber import RubberAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, RubberAssay):
        raise TypeError(f"Expected RubberAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "wet_to_dry": lambda q, d=drc_pct: assay.wet_to_dry_mass(q, drc_pct=d),
        "dry_to_wet": lambda q, d=drc_pct: assay.dry_to_wet_mass(q, drc_pct=d),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported rubber mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    eff_drc = Decimal(str(assay.dry_rubber_content_pct if drc_pct is None else drc_pct))
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m,
        classification=ConversionClassification.ASSAY_DEPENDENT,
        return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units), p["drc_pct"]).to(tgt),
        extra_base_params={"drc_pct": eff_drc},
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts feedstock, ethanol, and DDGS quantities across mass and volume."""
    from core.units.assays.agriculture.feedstocks import EthanolFeedstockAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, EthanolFeedstockAssay):
        raise TypeError(f"Expected EthanolFeedstockAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "feedstock_to_ethanol": lambda q: assay.feedstock_mass_to_ethanol_volume(q),
        "ethanol_to_feedstock": lambda q: assay.ethanol_volume_to_feedstock_mass(q),
        "feedstock_to_ddgs": lambda q: assay.feedstock_mass_to_ddgs_mass(q),
        "ddgs_to_feedstock": lambda q: assay.ddgs_mass_to_feedstock_mass(q),
        "ethanol_vol_to_mass": lambda q: assay.ethanol_volume_to_mass(q),
        "ethanol_mass_to_vol": lambda q: assay.ethanol_mass_to_volume(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(
            f"Unsupported ethanol feedstock mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    res = dispatch[m](quantity).to(tgt)
    cls_tag = ConversionClassification.ASSAY_DEPENDENT if m in {"ethanol_vol_to_mass", "ethanol_mass_to_vol"} else ConversionClassification.YIELD_DEPENDENT
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, classification=cls_tag, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts rice quantities across volume, mass, paddy, brown, and milled basis."""
    from core.units.assays.agriculture.rice import RiceAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, RiceAssay):
        raise TypeError(f"Expected RiceAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "paddy_to_brown": lambda q: assay.paddy_to_brown_mass(q),
        "brown_to_paddy": lambda q: assay.brown_to_paddy_mass(q),
        "brown_to_milled": lambda q: assay.brown_to_milled_mass(q),
        "milled_to_brown": lambda q: assay.milled_to_brown_mass(q),
        "paddy_to_milled": lambda q: assay.paddy_to_milled_mass(q),
        "milled_to_paddy": lambda q: assay.milled_to_paddy_mass(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported rice mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts wool quantities across volume, greasy mass, clean scoured mass, and tops."""
    from core.units.assays.agriculture.fibres import WoolAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, WoolAssay):
        raise TypeError(f"Expected WoolAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "greasy_to_clean": lambda q, y=clean_yield_pct: assay.greasy_to_clean_mass(q, clean_yield_pct=y),
        "clean_to_greasy": lambda q, y=clean_yield_pct: assay.clean_to_greasy_mass(q, clean_yield_pct=y),
        "clean_to_top": lambda q, _y=None: assay.clean_to_top_mass(q),
        "top_to_clean": lambda q, _y=None: assay.top_to_clean_mass(q),
        "top_to_greasy": lambda q, y=clean_yield_pct: assay.top_to_greasy_mass(q, clean_yield_pct=y),
        "greasy_to_top": lambda q, y=clean_yield_pct: assay.greasy_to_top_mass(q, clean_yield_pct=y),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported wool mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    eff_y = Decimal(str(assay.clean_yield_pct if clean_yield_pct is None else clean_yield_pct))
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units), p["clean_yield_pct"]).to(tgt),
        extra_base_params={"clean_yield_pct": eff_y},
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts poultry and egg quantities across volume, live, RTC, boneless, shell, liquid, and dried mass."""
    from core.units.assays.agriculture.poultry import PoultryEggAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, PoultryEggAssay):
        raise TypeError(f"Expected PoultryEggAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

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
        raise UnsupportedConversionError(f"Unsupported poultry/egg mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )





def create_fishmeal_marine_oil_context(
    assay: FishmealMarineOilAssay,
    *,
    name: str | None = None,
    product: str = "meal",
) -> pint.Context:
    """Bridges [volume] <-> [mass] for either fishmeal ('meal') or marine oil ('oil')."""
    from core.units.exceptions import UnsupportedConversionError
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
        raise UnsupportedConversionError(f"Unsupported product {product!r}; expected 'meal' or 'oil'.")

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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts whole catch, fishmeal, and crude marine oil across volume and mass."""
    from core.units.assays.agriculture.seafood import FishmealMarineOilAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, FishmealMarineOilAssay):
        raise TypeError(f"Expected FishmealMarineOilAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay, product=product)
        return convert(
            quantity, tgt, assay=assay, product=product,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

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
        raise UnsupportedConversionError(
            f"Unsupported fishmeal/marine oil mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    res = dispatch[m](quantity).to(tgt)
    cls_tag = ConversionClassification.ASSAY_DEPENDENT if "vol" in m else ConversionClassification.YIELD_DEPENDENT
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, classification=cls_tag, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts bulk honey and beeswax quantities across volume, mass, and moisture basis."""
    from core.units.assays.agriculture.softs import HoneyApicultureAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, HoneyApicultureAssay):
        raise TypeError(f"Expected HoneyApicultureAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "normalise_moisture": lambda q: assay.normalise_to_base_moisture(
            q, actual_moisture_pct=actual_moisture_pct
        ),
        "reverse_moisture": lambda q: assay.reverse_moisture_normalisation(q, actual_moisture_pct=actual_moisture_pct),
        "honey_to_beeswax": lambda q: assay.honey_to_beeswax_mass(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(
            f"Unsupported honey/apiculture mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    res = dispatch[m](quantity).to(tgt)
    cls_tag = ConversionClassification.YIELD_DEPENDENT if m == "honey_to_beeswax" else ConversionClassification.ASSAY_DEPENDENT
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, classification=cls_tag, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts olive oil quantities across volume, mass, and refining yield."""
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, OliveOilAssay):
        raise TypeError(f"Expected OliveOilAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay, observed_temp=temp_c)
        return convert(
            quantity, tgt, assay=assay, observed_temp=temp_c,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "vol_to_mass": lambda q: assay.volume_to_mass(q, temp_c=temp_c),
        "mass_to_vol": lambda q: assay.mass_to_volume(q, temp_c=temp_c),
        "crude_to_refined": lambda q: assay.crude_to_refined_mass(q),
        "refined_to_crude": lambda q: assay.refined_to_crude_mass(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(f"Unsupported olive oil mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'.")
    res = dispatch[m](quantity).to(tgt)
    cls_tag = ConversionClassification.YIELD_DEPENDENT if "refined" in m else ConversionClassification.ASSAY_DEPENDENT
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, classification=cls_tag, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )





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
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Converts seed cotton, lint, and whole cottonseed across mass and bulk volume."""
    from core.units.assays.agriculture.fibres import CottonSeedLintAssay
    from core.units.conversions import convert
    from core.units.exceptions import UnsupportedConversionError
    from core.units.provenance import ConversionClassification
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"Expected pint.Quantity, got {type(quantity).__name__}")
    if not isinstance(assay, CottonSeedLintAssay):
        raise TypeError(f"Expected CottonSeedLintAssay, got {type(assay).__name__}")

    tgt = parse_unit(target_unit)
    m = mode.strip().lower()
    if m == "auto":
        if not return_metadata and not strict and uncertainties is None and covariance is None and correlations is None:
            return convert(quantity, tgt, assay=assay)
        return convert(
            quantity, tgt, assay=assay,
            return_metadata=return_metadata, strict=strict,
            uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        )

    dispatch = {
        "seed_cotton_to_lint": lambda q: assay.seed_cotton_to_lint_mass(q),
        "lint_to_seed_cotton": lambda q: assay.lint_to_seed_cotton_mass(q),
        "seed_cotton_to_seed": lambda q: assay.seed_cotton_to_cottonseed_mass(q),
        "seed_to_seed_cotton": lambda q: assay.cottonseed_to_seed_cotton_mass(q),
        "seed_vol_to_mass": lambda q: assay.cottonseed_volume_to_mass(q),
        "seed_mass_to_vol": lambda q: assay.cottonseed_mass_to_volume(q),
    }
    if m not in dispatch:
        raise UnsupportedConversionError(
            f"Unsupported cotton seed/lint mode {mode!r}. Expected one of {sorted(dispatch)} or 'auto'."
        )
    res = dispatch[m](quantity).to(tgt)
    cls_tag = ConversionClassification.ASSAY_DEPENDENT if "vol" in m else ConversionClassification.YIELD_DEPENDENT
    return _finalize_agricultural_mode(
        res, quantity, assay=assay, mode=m, classification=cls_tag, return_metadata=return_metadata,
        uncertainties=uncertainties, covariance=covariance, correlations=correlations,
        rerun=lambda p: dispatch[m](ureg.Quantity(p["quantity"], quantity.units)).to(tgt),
    )


__all__ = ('AgriculturalResidueAssay', 'AnimalFeedAssay', 'CottonSeedLintAssay', 'DairyProductAssay', 'EthanolFeedstockAssay', 'FishmealMarineOilAssay', 'FruitConcentrateAssay', 'GrainAssay', 'HoneyApicultureAssay', 'LivestockDairyAssay', 'MeatProductAssay', 'NutAssay', 'OilseedAssay', 'OilseedCrushAssay', 'OliveOilAssay', 'PaperAssay', 'PoultryEggAssay', 'PulpLumberAssay', 'PulseAssay', 'RiceAssay', 'RootCropAssay', 'RubberAssay', 'SeafoodAssay', 'SilkAndNaturalFibreAssay', 'SoftsAssay', 'SpecialtyGrainAssay', 'SpiceAssay', 'StarchAssay', 'TeaAssay', 'TimberAssay', 'TobaccoAssay', 'VegetableOilAssay', 'WoolAssay', 'convert_agricultural_residue', 'convert_animal_feed', 'convert_cotton_seed_lint', 'convert_dairy_product', 'convert_ethanol_feedstock', 'convert_fishmeal_marine_oil', 'convert_fruit_concentrate', 'convert_grain', 'convert_honey_apiculture', 'convert_livestock_dairy', 'convert_meat_product', 'convert_nut', 'convert_oilseed', 'convert_oilseed_crush', 'convert_olive_oil', 'convert_paper', 'convert_poultry_egg', 'convert_pulp_lumber', 'convert_pulse', 'convert_rice', 'convert_root_crop', 'convert_rubber', 'convert_seafood', 'convert_silk_and_natural_fibre', 'convert_softs', 'convert_specialty_grain', 'convert_spice', 'convert_starch', 'convert_tea', 'convert_timber', 'convert_tobacco', 'convert_vegetable_oil', 'convert_wool', 'create_agricultural_residue_context', 'create_animal_feed_context', 'create_cotton_seed_lint_context', 'create_dairy_product_context', 'create_ethanol_feedstock_context', 'create_fishmeal_marine_oil_context', 'create_fruit_concentrate_context', 'create_grain_context', 'create_honey_apiculture_context', 'create_livestock_dairy_context', 'create_meat_product_context', 'create_nut_context', 'create_oilseed_context', 'create_oilseed_crush_context', 'create_olive_oil_context', 'create_paper_context', 'create_poultry_egg_context', 'create_pulp_lumber_context', 'create_pulse_context', 'create_rice_context', 'create_root_crop_context', 'create_rubber_context', 'create_seafood_context', 'create_silk_and_natural_fibre_context', 'create_softs_context', 'create_specialty_grain_context', 'create_spice_context', 'create_starch_context', 'create_tea_context', 'create_timber_context', 'create_tobacco_context', 'create_vegetable_oil_context', 'create_wool_context')
_PUBLIC_ALIASES = {'convert_fruit_concentrate': ('assays.agriculture.fruit', 'convert_fruit_concentrate'), 'RootCropAssay': ('assays.agriculture.fruit', 'RootCropAssay'), 'SpiceAssay': ('assays.agriculture.specialty', 'SpiceAssay'), 'create_tobacco_context': ('assays.agriculture.specialty', 'create_tobacco_context'), 'create_animal_feed_context': ('assays.agriculture.feedstocks', 'create_animal_feed_context'), 'convert_tobacco': ('assays.agriculture.specialty', 'convert_tobacco'), '_universal_moisture_level_ratio': ('analytical_basis', 'moisture_level_ratio'), 'create_fruit_concentrate_context': ('assays.agriculture.fruit', 'create_fruit_concentrate_context'), 'StarchAssay': ('assays.agriculture.feedstocks', 'StarchAssay'), 'EthanolFeedstockAssay': ('assays.agriculture.feedstocks', 'EthanolFeedstockAssay'), 'create_paper_context': ('assays.agriculture.forestry', 'create_paper_context'), 'create_timber_context': ('assays.agriculture.forestry', 'create_timber_context'), 'PoultryEggAssay': ('assays.agriculture.poultry', 'PoultryEggAssay'), 'convert_root_crop': ('assays.agriculture.fruit', 'convert_root_crop'), 'PulpLumberAssay': ('assays.agriculture.forestry', 'PulpLumberAssay'), 'MissingAssayError': ('exceptions', 'MissingAssayError'), 'create_dairy_product_context': ('assays.agriculture.dairy', 'create_dairy_product_context'), 'GrainAssay': ('assays.agriculture.grains', 'GrainAssay'), 'convert_specialty_grain': ('assays.agriculture.specialty', 'convert_specialty_grain'), 'convert_agricultural_residue': ('assays.agriculture.feedstocks', 'convert_agricultural_residue'), 'LivestockDairyAssay': ('assays.agriculture.livestock', 'LivestockDairyAssay'), 'ensure_quantity_ieee754': ('precision', 'ensure_quantity_ieee754'), 'TeaAssay': ('assays.agriculture.specialty', 'TeaAssay'), 'AnimalFeedAssay': ('assays.agriculture.feedstocks', 'AnimalFeedAssay'), 'create_starch_context': ('assays.agriculture.feedstocks', 'create_starch_context'), 'PulseAssay': ('assays.agriculture.pulses', 'PulseAssay'), 'SeafoodAssay': ('assays.agriculture.seafood', 'SeafoodAssay'), '_universal_wet_to_dry_mass': ('analytical_basis', 'wet_to_dry_mass'), 'convert_animal_feed': ('assays.agriculture.feedstocks', 'convert_animal_feed'), 'RubberAssay': ('assays.agriculture.rubber', 'RubberAssay'), 'create_tea_context': ('assays.agriculture.specialty', 'create_tea_context'), '_universal_dry_matter_fraction': ('analytical_basis', 'dry_matter_fraction'), 'TobaccoAssay': ('assays.agriculture.specialty', 'TobaccoAssay'), 'WoolAssay': ('assays.agriculture.fibres', 'WoolAssay'), 'to_ieee754_float': ('precision', 'to_ieee754_float'), 'convert_starch': ('assays.agriculture.feedstocks', 'convert_starch'), 'FishmealMarineOilAssay': ('assays.agriculture.seafood', 'FishmealMarineOilAssay'), 'create_meat_product_context': ('assays.agriculture.livestock', 'create_meat_product_context'), 'create_agricultural_residue_context': ('assays.agriculture.feedstocks', 'create_agricultural_residue_context'), 'DairyProductAssay': ('assays.agriculture.dairy', 'DairyProductAssay'), '_universal_dry_to_wet_mass': ('analytical_basis', 'dry_to_wet_mass'), 'SilkAndNaturalFibreAssay': ('assays.agriculture.fibres', 'SilkAndNaturalFibreAssay'), 'FruitConcentrateAssay': ('assays.agriculture.fruit', 'FruitConcentrateAssay'), 'MeatProductAssay': ('assays.agriculture.livestock', 'MeatProductAssay'), 'pad_ieee754_decimal': ('precision', 'pad_ieee754_decimal'), 'IncompatibleCommodityError': ('exceptions', 'IncompatibleCommodityError'), 'OilseedAssay': ('assays.agriculture.oilseeds', 'OilseedAssay'), 'PaperAssay': ('assays.agriculture.forestry', 'PaperAssay'), 'create_spice_context': ('assays.agriculture.specialty', 'create_spice_context'), '_universal_convert_moisture_mass': ('analytical_basis', 'convert_moisture_mass'), 'convert_spice': ('assays.agriculture.specialty', 'convert_spice'), 'validated_assay': ('validation', 'validated_assay'), 'OilseedCrushAssay': ('assays.agriculture.feedstocks', 'OilseedCrushAssay'), 'convert_timber': ('assays.agriculture.forestry', 'convert_timber'), 'CottonSeedLintAssay': ('assays.agriculture.fibres', 'CottonSeedLintAssay'), 'create_root_crop_context': ('assays.agriculture.fruit', 'create_root_crop_context'), 'create_vegetable_oil_context': ('assays.agriculture.vegetable_oils', 'create_vegetable_oil_context'), 'convert_tea': ('assays.agriculture.specialty', 'convert_tea'), 'StandardsApplicabilityError': ('exceptions', 'StandardsApplicabilityError'), 'create_nut_context': ('assays.agriculture.specialty', 'create_nut_context'), 'convert_units': ('conversions', 'convert'), 'TimberAssay': ('assays.agriculture.forestry', 'TimberAssay'), 'create_oilseed_context': ('assays.agriculture.oilseeds', 'create_oilseed_context'), 'NutAssay': ('assays.agriculture.specialty', 'NutAssay'), 'convert_paper': ('assays.agriculture.forestry', 'convert_paper'), '_universal_rebase_moisture_intensive': ('analytical_basis', 'rebase_moisture_intensive'), 'create_specialty_grain_context': ('assays.agriculture.specialty', 'create_specialty_grain_context'), 'RiceAssay': ('assays.agriculture.rice', 'RiceAssay'), 'create_pulse_context': ('assays.agriculture.pulses', 'create_pulse_context'), 'HoneyApicultureAssay': ('assays.agriculture.softs', 'HoneyApicultureAssay'), 'SoftsAssay': ('assays.agriculture.softs', 'SoftsAssay'), '_scale_factor': ('conversions', '_scale_factor'), 'ConversionClassification': ('provenance', 'ConversionClassification'), 'AgriculturalResidueAssay': ('assays.agriculture.feedstocks', 'AgriculturalResidueAssay'), 'build_conversion_result': ('provenance', 'build_conversion_result'), 'ureg': ('registry', 'ureg'), 'convert_oilseed': ('assays.agriculture.oilseeds', 'convert_oilseed'), 'SpecialtyGrainAssay': ('assays.agriculture.specialty', 'SpecialtyGrainAssay'), 'preserve_ieee754_scalar': ('precision', 'preserve_ieee754_scalar'), 'convert_silk_and_natural_fibre': ('assays.agriculture.fibres', 'convert_silk_and_natural_fibre'), 'parse_unit': ('registry', 'parse_unit'), 'convert_meat_product': ('assays.agriculture.livestock', 'convert_meat_product'), 'convert_dairy_product': ('assays.agriculture.dairy', 'convert_dairy_product'), 'UnsupportedConversionError': ('exceptions', 'UnsupportedConversionError'), 'create_silk_and_natural_fibre_context': ('assays.agriculture.fibres', 'create_silk_and_natural_fibre_context'), 'VegetableOilAssay': ('assays.agriculture.vegetable_oils', 'VegetableOilAssay'), 'convert_vegetable_oil': ('assays.agriculture.vegetable_oils', 'convert_vegetable_oil'), 'convert_nut': ('assays.agriculture.specialty', 'convert_nut'), 'convert_pulse': ('assays.agriculture.pulses', 'convert_pulse'), 'OliveOilAssay': ('assays.agriculture.vegetable_oils', 'OliveOilAssay')}

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
    from core.units.assays.agriculture.feedstocks import AgriculturalResidueAssay
    from core.units.assays.agriculture.feedstocks import AnimalFeedAssay
    from core.units.assays.agriculture.fibres import CottonSeedLintAssay
    from core.units.assays.agriculture.dairy import DairyProductAssay
    from core.units.assays.agriculture.feedstocks import EthanolFeedstockAssay
    from core.units.assays.agriculture.seafood import FishmealMarineOilAssay
    from core.units.assays.agriculture.fruit import FruitConcentrateAssay
    from core.units.assays.agriculture.grains import GrainAssay
    from core.units.assays.agriculture.softs import HoneyApicultureAssay
    from core.units.assays.agriculture.livestock import LivestockDairyAssay
    from core.units.assays.agriculture.livestock import MeatProductAssay
    from core.units.assays.agriculture.specialty import NutAssay
    from core.units.assays.agriculture.oilseeds import OilseedAssay
    from core.units.assays.agriculture.feedstocks import OilseedCrushAssay
    from core.units.assays.agriculture.vegetable_oils import OliveOilAssay
    from core.units.assays.agriculture.forestry import PaperAssay
    from core.units.assays.agriculture.poultry import PoultryEggAssay
    from core.units.assays.agriculture.forestry import PulpLumberAssay
    from core.units.assays.agriculture.pulses import PulseAssay
    from core.units.assays.agriculture.rice import RiceAssay
    from core.units.assays.agriculture.fruit import RootCropAssay
    from core.units.assays.agriculture.rubber import RubberAssay
    from core.units.assays.agriculture.seafood import SeafoodAssay
    from core.units.assays.agriculture.fibres import SilkAndNaturalFibreAssay
    from core.units.assays.agriculture.softs import SoftsAssay
    from core.units.assays.agriculture.specialty import SpecialtyGrainAssay
    from core.units.assays.agriculture.specialty import SpiceAssay
    from core.units.assays.agriculture.feedstocks import StarchAssay
    from core.units.assays.agriculture.specialty import TeaAssay
    from core.units.assays.agriculture.forestry import TimberAssay
    from core.units.assays.agriculture.specialty import TobaccoAssay
    from core.units.assays.agriculture.vegetable_oils import VegetableOilAssay
    from core.units.assays.agriculture.fibres import WoolAssay
    from core.units.assays.agriculture.feedstocks import convert_agricultural_residue
    from core.units.assays.agriculture.feedstocks import convert_animal_feed
    from core.units.assays.agriculture.dairy import convert_dairy_product
    from core.units.assays.agriculture.fruit import convert_fruit_concentrate
    from core.units.assays.agriculture.livestock import convert_meat_product
    from core.units.assays.agriculture.specialty import convert_nut
    from core.units.assays.agriculture.oilseeds import convert_oilseed
    from core.units.assays.agriculture.forestry import convert_paper
    from core.units.assays.agriculture.pulses import convert_pulse
    from core.units.assays.agriculture.fruit import convert_root_crop
    from core.units.assays.agriculture.fibres import convert_silk_and_natural_fibre
    from core.units.assays.agriculture.specialty import convert_specialty_grain
    from core.units.assays.agriculture.specialty import convert_spice
    from core.units.assays.agriculture.feedstocks import convert_starch
    from core.units.assays.agriculture.specialty import convert_tea
    from core.units.assays.agriculture.forestry import convert_timber
    from core.units.assays.agriculture.specialty import convert_tobacco
    from core.units.assays.agriculture.vegetable_oils import convert_vegetable_oil
    from core.units.assays.agriculture.feedstocks import create_agricultural_residue_context
    from core.units.assays.agriculture.feedstocks import create_animal_feed_context
    from core.units.assays.agriculture.dairy import create_dairy_product_context
    from core.units.assays.agriculture.fruit import create_fruit_concentrate_context
    from core.units.assays.agriculture.livestock import create_meat_product_context
    from core.units.assays.agriculture.specialty import create_nut_context
    from core.units.assays.agriculture.oilseeds import create_oilseed_context
    from core.units.assays.agriculture.forestry import create_paper_context
    from core.units.assays.agriculture.pulses import create_pulse_context
    from core.units.assays.agriculture.fruit import create_root_crop_context
    from core.units.assays.agriculture.fibres import create_silk_and_natural_fibre_context
    from core.units.assays.agriculture.specialty import create_specialty_grain_context
    from core.units.assays.agriculture.specialty import create_spice_context
    from core.units.assays.agriculture.feedstocks import create_starch_context
    from core.units.assays.agriculture.specialty import create_tea_context
    from core.units.assays.agriculture.forestry import create_timber_context
    from core.units.assays.agriculture.specialty import create_tobacco_context
    from core.units.assays.agriculture.vegetable_oils import create_vegetable_oil_context
