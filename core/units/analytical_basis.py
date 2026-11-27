"""Commodex analytical_basis implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, fields
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final, Literal, Self
from typing import Any, Final, overload
from typing import ClassVar
import math
import pint




def _agricultural_families_yield(value, label):
    from core.units.assays.energy.petroleum import _converters_petroleum_fraction
    return None if value is None else _converters_petroleum_fraction(value, label)





def _agricultural_families_forward(quantity, fraction, label):
    from core.units.assays.base import scale_mass
    if fraction is None:
        raise ValueError(f"Explicit {label} is required.")
    return scale_mass(quantity, fraction)





def _agricultural_families_inverse(quantity, fraction, label):
    from core.units.assays.base import mass
    from core.units.assays.energy.petroleum import _converters_petroleum_output
    if fraction is None or fraction == 0:
        raise ValueError(f"Positive {label} is required for the inverse conversion.")
    return _converters_petroleum_output(mass(quantity) / fraction, "kg", quantity)





# ---------------------------------------------------------------------------
# Physical Constants & Validation Helpers
# ---------------------------------------------------------------------------

# Exact bulk density conversion: 1 lb / US Winchester bushel (2150.42 in^3) in kg/hL
_agricultural_quality_WINCHESTER_BUSHEL_M3: Final[Decimal] = Decimal("2150.42") * (Decimal("0.0254") ** 3)



_agricultural_quality_LB_KG: Final[Decimal] = Decimal("0.45359237")



_agricultural_quality_KG_PER_HL_PER_LB_WINCHESTER_BU: Final[Decimal] = (_agricultural_quality_LB_KG / _agricultural_quality_WINCHESTER_BUSHEL_M3) / Decimal(10)



_agricultural_quality_CUBIC_METRES_PER_MBF: Final[Decimal] = Decimal("2.359737216")  # 1,000 board feet = 1000 * 144 in^3




# Fatty acid molar masses (g/mol) and KOH molar mass (g/mol) for Acid Value <-> FFA %
_agricultural_quality_KOH_MG_PER_MOL: Final[Decimal] = Decimal("56105.6")



_agricultural_quality_FFA_MOLAR_MASS_G_MOL: Final[Mapping[str, Decimal]] = MappingProxyType({
    "oleic": Decimal("282.46"),
    "palmitic": Decimal("256.42"),
    "lauric": Decimal("200.32"),
})




# DOC PS 20 American Softwood Lumber Standard actual surfaced dimensions (inches)
# Maps nominal thickness/width in inches -> (S4S_DRY_in, S4S_GREEN_in)
_agricultural_quality_PS20_DIMENSION_TABLE_IN: Final[Mapping[Decimal, tuple[Decimal, Decimal]]] = MappingProxyType({
    Decimal("1"): (Decimal("0.75"), Decimal("0.78125")),
    Decimal("2"): (Decimal("1.5"), Decimal("1.5625")),
    Decimal("3"): (Decimal("2.5"), Decimal("2.5625")),
    Decimal("4"): (Decimal("3.5"), Decimal("3.5625")),
    Decimal("6"): (Decimal("5.5"), Decimal("5.625")),
    Decimal("8"): (Decimal("7.25"), Decimal("7.5")),
    Decimal("10"): (Decimal("9.25"), Decimal("9.5")),
    Decimal("12"): (Decimal("11.25"), Decimal("11.5")),
})




# ISO 2000 physical impurity ceilings for Technically Specified Rubber (% by mass)
_agricultural_quality_TSR_ISO2000_LIMITS: Final[Mapping[str, Mapping[str, Decimal]]] = MappingProxyType({
    "TSR5": MappingProxyType({
        "dirt_pct": Decimal("0.04"),
        "ash_pct": Decimal("0.60"),
        "volatile_matter_pct": Decimal("0.80"),
        "nitrogen_pct": Decimal("0.60"),
    }),
    "TSR10": MappingProxyType({
        "dirt_pct": Decimal("0.08"),
        "ash_pct": Decimal("0.75"),
        "volatile_matter_pct": Decimal("0.80"),
        "nitrogen_pct": Decimal("0.60"),
    }),
    "TSR20": MappingProxyType({
        "dirt_pct": Decimal("0.16"),
        "ash_pct": Decimal("1.00"),
        "volatile_matter_pct": Decimal("0.80"),
        "nitrogen_pct": Decimal("0.60"),
    }),
})





def _agricultural_quality_dec(value: Any, label: str) -> Decimal:
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise TypeError(f"{label} must be int, float, or Decimal, got {type(value).__name__}.")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{label} must be finite, got {value!r}.")
    d = value if isinstance(value, Decimal) else Decimal(str(value))
    if not d.is_finite():
        raise ValueError(f"{label} must be finite, got {value!r}.")
    return d





def _agricultural_quality_pct(
    value: Any,
    label: str,
    *,
    low: Decimal = Decimal(0),
    high: Decimal = Decimal(100),
    allow_low: bool = True,
    allow_high: bool = True,
) -> Decimal:
    d = _agricultural_quality_dec(value, label)
    bad_low = d < low if allow_low else d <= low
    bad_high = d > high if allow_high else d >= high
    if bad_low or bad_high:
        lb = "[" if allow_low else "("
        rb = "]" if allow_high else ")"
        raise ValueError(f"{label} must lie in {lb}{low}, {high}{rb}, got {value!r}.")
    return d





def _agricultural_quality_nonneg(value: Any, label: str, *, strict_positive: bool = False) -> Decimal:
    d = _agricultural_quality_dec(value, label)
    if (d <= Decimal(0)) if strict_positive else (d < Decimal(0)):
        op = "> 0" if strict_positive else ">= 0"
        raise ValueError(f"{label} must be {op}, got {value!r}.")
    return d





def _agricultural_quality_clean_str(value: str, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string.")
    return value.strip()





def _agricultural_quality_blend_weights(w: Decimal | float) -> tuple[Decimal, Decimal]:
    w1 = _agricultural_quality_pct(w, "self_mass_fraction", low=Decimal(0), high=Decimal(1))
    return w1, Decimal(1) - w1





def _agricultural_quality_preserve_type(dec_val: Decimal, *sources: Any) -> Decimal | float:
    for src in sources:
        mag = src.magnitude if isinstance(src, pint.Quantity) else src
        if isinstance(mag, Decimal):
            return dec_val
    for src in sources:
        mag = src.magnitude if isinstance(src, pint.Quantity) else src
        if isinstance(mag, float):
            return float(dec_val)
    return dec_val


from core.units.validation import authorized_conversion


@authorized_conversion
def _agricultural_quality_mass_kg(qty: pint.Quantity, label: str = "mass") -> Decimal:
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(qty, pint.Quantity):
        raise TypeError(f"{label} must be a pint.Quantity, got {type(qty).__name__}.")
    validate_registry_quantity(qty)
    kg = qty.to(ureg.kilogram)
    d = _agricultural_quality_nonneg(kg.magnitude, label)
    return d


from core.units.validation import authorized_conversion


@authorized_conversion
def _agricultural_quality_kg_out(kg_val: Decimal, target_unit: str | Unit, *sources: Any) -> pint.Quantity:
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    dst = parse_unit(target_unit)
    converted = ureg.Quantity(kg_val, ureg.kilogram).to(dst)
    return ureg.Quantity(_agricultural_quality_preserve_type(_agricultural_quality_dec(converted.magnitude, "mass"), *sources), dst)





def rebase_constituent_pct(
    constituent_pct: Decimal | float,
    *,
    source_moisture_pct: Decimal | float,
    target_moisture_pct: Decimal | float,
) -> Decimal:
    """Re-base a constituent percentage (protein, oil, gluten, fat) across moisture bases.

    Conserves constituent-to-dry-matter ratio:
    ``pct_target = pct_source * (100 - target_moisture) / (100 - source_moisture)``.
    Pass ``target_moisture_pct=0`` for dry-matter basis, ``12.0`` for US FGIS wheat
    basis, ``13.5`` for Canadian CGC wheat basis, or ``13.0`` for standard soybean basis.
    """
    val = _agricultural_quality_pct(constituent_pct, "constituent_pct")
    m_src = _agricultural_quality_pct(source_moisture_pct, "source_moisture_pct", allow_high=False)
    m_dst = _agricultural_quality_pct(target_moisture_pct, "target_moisture_pct", allow_high=False)
    res = rebase_moisture_intensive(
        val,
        source_moisture_pct=m_src,
        target_moisture_pct=m_dst,
        allow_negative=False,
    )
    return Decimal(str(res)) if not isinstance(res, Decimal) else res





def normalize_mass_to_moisture(
    mass: pint.Quantity,
    *,
    actual_moisture_pct: Decimal | float,
    base_moisture_pct: Decimal | float,
    allow_moisture_gain: bool = False,
    target_unit: str | Unit | None = None,
) -> pint.Quantity:
    """Normalise delivered agricultural mass to a contractual moisture basis.

    Uses exact dry-matter conservation ``(100 - M_actual) / (100 - M_base)``.
    When ``allow_moisture_gain=False`` (standard commercial elevator/export rule),
    lots drier than ``base_moisture_pct`` are not credited with phantom water weight.
    """
    return convert_moisture_mass(
        mass,
        source_moisture_pct=actual_moisture_pct,
        target_moisture_pct=base_moisture_pct,
        target_unit=target_unit if target_unit is not None else mass.units,
        allow_moisture_gain=allow_moisture_gain,
    )





# ---------------------------------------------------------------------------
# Physical Specification Screening & Benchmark Delta Containers
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class QualityThreshold:
    """Physical acceptance bounds for a single quality attribute (no price schedule)."""

    min_value: Decimal | None = None
    max_value: Decimal | None = None

    def __post_init__(self) -> None:
        lo = _agricultural_quality_dec(self.min_value, "min_value") if self.min_value is not None else None
        hi = _agricultural_quality_dec(self.max_value, "max_value") if self.max_value is not None else None
        if lo is not None and hi is not None and lo > hi:
            raise ValueError(f"min_value ({lo}) cannot exceed max_value ({hi}).")
        object.__setattr__(self, "min_value", lo)
        object.__setattr__(self, "max_value", hi)

    def deviation(self, actual: Decimal | float) -> Decimal:
        """Signed physical excursion outside ``[min_value, max_value]`` (0 inside band)."""
        val = _agricultural_quality_dec(actual, "actual")
        if self.min_value is not None and val < self.min_value:
            return val - self.min_value
        if self.max_value is not None and val > self.max_value:
            return val - self.max_value
        return Decimal(0)

    def passes(self, actual: Decimal | float) -> bool:
        return self.deviation(actual) == Decimal(0)





@dataclass(frozen=True, slots=True)
class QualityAuditResult:
    """Physical quality compliance and signed excursion report for downstream pricing layers."""

    passed: bool
    deviations: Mapping[str, Decimal]
    failed_attributes: tuple[str, ...]

    def __post_init__(self) -> None:
        clean_dev = {str(k): _agricultural_quality_dec(v, k) for k, v in self.deviations.items()}
        object.__setattr__(self, "deviations", MappingProxyType(clean_dev))
        object.__setattr__(self, "failed_attributes", tuple(self.failed_attributes))





class AgriculturalQualityBase:
    """Shared physical audit and benchmark comparison interface for agricultural specs."""

    def evaluate_thresholds(self, thresholds: Mapping[str, QualityThreshold]) -> QualityAuditResult:
        """Screen this lot's physical attributes against ``QualityThreshold`` bounds."""
        deviations: dict[str, Decimal] = {}
        failed: list[str] = []
        for attr, rule in thresholds.items():
            if not isinstance(rule, QualityThreshold):
                raise TypeError(f"Threshold for {attr!r} must be a QualityThreshold.")
            if not hasattr(self, attr):
                raise KeyError(f"{type(self).__name__} has no quality attribute {attr!r}.")
            raw = getattr(self, attr)
            if raw is None:
                raise ValueError(f"Attribute {attr!r} is None on {type(self).__name__}; cannot screen threshold.")
            val = raw.magnitude if isinstance(raw, pint.Quantity) else raw
            dev = rule.deviation(val)
            deviations[attr] = dev
            if dev != Decimal(0):
                failed.append(attr)
        return QualityAuditResult(
            passed=not failed,
            deviations=deviations,
            failed_attributes=tuple(failed),
        )

    def compare_to_benchmark(self, benchmark: Self) -> dict[str, Decimal]:
        """Compute signed physical attribute deltas (``actual - benchmark``) for numeric fields.

        Subclasses override when analytical moisture-basis alignment is required first.
        """
        if type(self) is not type(benchmark):
            raise TypeError(f"Expected {type(self).__name__}, got {type(benchmark).__name__}.")
        deltas: dict[str, Decimal] = {}
        for f in fields(self):  # type: ignore[arg-type]
            v_self = getattr(self, f.name)
            v_bench = getattr(benchmark, f.name)
            if isinstance(v_self, Decimal) and isinstance(v_bench, Decimal):
                deltas[f.name] = v_self - v_bench
        return deltas





# ---------------------------------------------------------------------------
# Non-Destructive Adapter from Existing Agricultural Assays
# ---------------------------------------------------------------------------


def quality_spec_from_assay(assay: Any) -> AgriculturalQualityBase:
    """Build a baseline physical ``*QualitySpec`` from an existing agricultural ``*Assay`` where applicable."""
    from core.units.assays.agriculture.dairy import MilkQualitySpec
    from core.units.assays.agriculture.fibres import CottonQualitySpec
    from core.units.assays.agriculture.fibres import WoolQualitySpec
    from core.units.assays.agriculture.forestry import TimberQualitySpec
    from core.units.assays.agriculture.grains import BarleyQualitySpec
    from core.units.assays.agriculture.grains import CornQualitySpec
    from core.units.assays.agriculture.grains import WheatQualitySpec
    from core.units.assays.agriculture.oilseeds import SoybeanQualitySpec
    from core.units.assays.agriculture.rice import RiceQualitySpec
    from core.units.assays.agriculture.rubber import RubberQualitySpec
    from core.units.assays.agriculture.softs import CocoaQualitySpec
    from core.units.assays.agriculture.softs import CoffeeQualitySpec
    from core.units.assays.agriculture.softs import SugarQualitySpec
    from core.units.assays.agriculture.vegetable_oils import VegetableOilQualitySpec
    from core.units import agricultural_specs as ag

    if isinstance(assay, AgriculturalQualityBase):
        return assay

    if isinstance(assay, ag.GrainAssay):
        moist = Decimal(str(assay.actual_moisture_pct if assay.actual_moisture_pct is not None else assay.base_moisture_pct))
        if assay.grain == "wheat":
            return WheatQualitySpec.from_winchester_test_weight(
                moisture_pct=moist,
                protein_pct=Decimal("11.5"),
                test_weight_lb_per_bu=Decimal(str(assay.bushel_weight_lb)),
            )
        if assay.grain in {"corn", "maize"}:
            return CornQualitySpec(moisture_pct=moist, test_weight_lb_per_bu=Decimal(str(assay.bushel_weight_lb)))
        if assay.grain in {"soybeans", "soybean"}:
            return SoybeanQualitySpec(moisture_pct=moist, protein_pct=Decimal("35.0"), oil_content_pct=Decimal("19.0"))
        if assay.grain == "barley":
            return BarleyQualitySpec(moisture_pct=moist, protein_pct=Decimal("11.0"), germination_energy_pct=Decimal("96.0"))

    if isinstance(assay, ag.RiceAssay):
        return RiceQualitySpec(
            milling_recovery_pct=Decimal(str(assay.milling_yield_pct)),
            head_rice_yield_pct=Decimal(str(assay.head_rice_pct)),
            head_rice_basis="of_milled",
        )

    if isinstance(assay, ag.SoftsAssay):
        if assay.commodity == "sugar":
            pol = Decimal(str(assay.actual_polarisation_deg or assay.reference_polarisation_deg or 96.0))
            return SugarQualitySpec(
                polarisation_deg_z=pol,
                icumsa_colour_iu=Decimal("800") if pol < Decimal(99) else Decimal("45"),
                refining_yield_pct=Decimal(str(assay.processing_yield_pct)),
            )
        if assay.commodity == "coffee":
            return CoffeeQualitySpec(
                screen_size=Decimal("16"),
                moisture_pct=Decimal("11.5"),
                processing_loss_pct=Decimal(100) - Decimal(str(assay.processing_yield_pct)),
            )
        if assay.commodity == "cocoa":
            return CocoaQualitySpec(
                bean_count_per_100g=Decimal("95"),
                fermentation_pct=Decimal("80"),
                fat_content_pct=Decimal("52"),
            )
        if assay.commodity == "cotton":
            return CottonQualitySpec.from_staple_32nds(
                staple_32nds_inch=Decimal("36"),
                micronaire=Decimal("4.2"),
                strength_gpt=Decimal("29.0"),
            )

    if isinstance(assay, ag.LivestockDairyAssay) and assay.category == "milk":
        fat = Decimal(str(assay.fat_pct if assay.fat_pct is not None else 4.0))
        prot = Decimal(str(assay.protein_pct if assay.protein_pct is not None else 3.3))
        return MilkQualitySpec(
            butterfat_pct=fat,
            protein_pct=prot,
            lactose_pct=Decimal("4.75"),
            density_kg_per_m3=Decimal(str(assay.milk_density_kg_per_m3)),
        )

    if isinstance(assay, ag.OliveOilAssay):
        return VegetableOilQualitySpec(
            free_fatty_acid_pct=Decimal(str(assay.free_fatty_acid_pct)),
            peroxide_value_meq_per_kg=Decimal("8.0"),
            iodine_value_g_per_100g=Decimal("82.0"),
            ffa_basis="oleic",
        )

    if isinstance(assay, ag.RubberAssay):
        return RubberQualitySpec(
            dirt_pct=Decimal("0.12"),
            ash_pct=Decimal("0.70"),
            volatile_matter_pct=Decimal("0.60"),
            mooney_viscosity_ml1_4_100c=Decimal("75.0"),
            dry_rubber_content_pct=Decimal(str(assay.dry_rubber_content_pct)),
        )

    if isinstance(assay, ag.WoolAssay):
        return WoolQualitySpec(
            micron_diameter_um=Decimal(str(assay.micron if assay.micron is not None else 21.0)),
            vegetable_matter_pct=Decimal("1.5"),
            clean_yield_pct=Decimal(str(assay.clean_yield_pct)),
        )

    if isinstance(assay, ag.PulpLumberAssay) and assay.product_type == "lumber":
        return TimberQualitySpec.from_north_american_nominal(
            2,
            4,
            species="SPF",
            moisture_content_pct=Decimal(str(assay.moisture_pct)),
        )

    raise TypeError(f"No agricultural quality adapter defined for {type(assay).__name__}.")


from core.units.registry import ureg


_moisture_MASS_FAMILY_DIMS: Final = frozenset({
    ureg.kilogram.dimensionality,
    ureg.wmt.dimensionality,
    ureg.dmt.dimensionality,
})





def _moisture_to_decimal(value: Any, label: str) -> Decimal:
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise TypeError(f"{label} must be int, float, or Decimal, got {type(value).__name__}.")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{label} must be finite, got {value!r}.")
    dec = value if isinstance(value, Decimal) else Decimal(str(value))
    if not dec.is_finite():
        raise ValueError(f"{label} must be finite, got {value!r}.")
    return dec





def _moisture_preserve_scalar_type(dec_val: Decimal, *sources: Any) -> Decimal | float:
    from core.units.precision import pad_ieee754_decimal
    for src in sources:
        mag = src.magnitude if isinstance(src, pint.Quantity) else src
        if isinstance(mag, Decimal):
            return pad_ieee754_decimal(dec_val)
    return float(dec_val)





def validate_moisture_fraction(
    moisture: float | Decimal,
    *,
    label: str = "moisture_fraction",
    allow_zero: bool = True,
) -> Decimal:
    """Validate and convert a wet-basis moisture fraction ``w`` in ``[0, 1)`` to ``Decimal``."""
    from core.units.exceptions import InvalidPhysicalStateError
    w = _moisture_to_decimal(moisture, label)
    too_low = w < Decimal(0) if allow_zero else w <= Decimal(0)
    if too_low or w >= Decimal(1):
        left = "[" if allow_zero else "("
        raise InvalidPhysicalStateError(f"{label} must lie in {left}0, 1), got {moisture!r}.")
    return w





def moisture_pct_to_fraction(
    moisture_pct: float | Decimal,
    *,
    label: str = "moisture_pct",
    allow_zero: bool = True,
) -> Decimal:
    """Validate a moisture percentage in ``[0, 100)`` and return fraction ``w = pct / 100``."""
    from core.units.exceptions import InvalidPhysicalStateError
    pct = _moisture_to_decimal(moisture_pct, label)
    too_low = pct < Decimal(0) if allow_zero else pct <= Decimal(0)
    if too_low or pct >= Decimal(100):
        left = "[" if allow_zero else "("
        raise InvalidPhysicalStateError(f"{label} must lie in {left}0, 100), got {moisture_pct!r}.")
    return pct / Decimal(100)





def _moisture_resolve_moisture(
    moisture: float | Decimal | None,
    moisture_pct: float | Decimal | None,
    *,
    frac_label: str,
    pct_label: str,
) -> Decimal:
    if (moisture is None) == (moisture_pct is None):
        raise ValueError(f"Supply exactly one of {frac_label!r} or {pct_label!r}.")
    if moisture is not None:
        return validate_moisture_fraction(moisture, label=frac_label)
    assert moisture_pct is not None
    return moisture_pct_to_fraction(moisture_pct, label=pct_label)





def dry_matter_fraction(
    moisture: float | Decimal | None = None,
    *,
    moisture_pct: float | Decimal | None = None,
) -> Decimal:
    """Return exact dry-matter fraction ``1 - w`` in ``(0, 1]``."""
    w = _moisture_resolve_moisture(
        moisture,
        moisture_pct,
        frac_label="moisture",
        pct_label="moisture_pct",
    )
    return Decimal(1) - w





def moisture_level_ratio(
    source_moisture: float | Decimal | None = None,
    target_moisture: float | Decimal | None = None,
    *,
    source_moisture_pct: float | Decimal | None = None,
    target_moisture_pct: float | Decimal | None = None,
    allow_moisture_gain: bool = True,
    handling_loss_fraction: float | Decimal = Decimal(0),
) -> Decimal:
    """Return universal mass multiplier ``(1 - w_1) / (1 - w_2)`` between two moisture levels."""
    w1 = _moisture_resolve_moisture(
        source_moisture,
        source_moisture_pct,
        frac_label="source_moisture",
        pct_label="source_moisture_pct",
    )
    w2 = _moisture_resolve_moisture(
        target_moisture,
        target_moisture_pct,
        frac_label="target_moisture",
        pct_label="target_moisture_pct",
    )
    loss = _moisture_to_decimal(handling_loss_fraction, "handling_loss_fraction")
    if loss < Decimal(0) or loss >= Decimal(1):
        raise ValueError(f"handling_loss_fraction must lie in [0, 1), got {handling_loss_fraction!r}.")

    if w1 == w2 or (w1 < w2 and not allow_moisture_gain):
        return Decimal(1)

    ratio = (Decimal(1) - w1) / (Decimal(1) - w2)
    if loss > Decimal(0):
        return max(Decimal(0), ratio - loss)
    return ratio





@overload
def convert_moisture_mass(
    mass: pint.Quantity,
    *,
    source_moisture: float | Decimal | None = ...,
    target_moisture: float | Decimal | None = ...,
    source_moisture_pct: float | Decimal | None = ...,
    target_moisture_pct: float | Decimal | None = ...,
    target_unit: str | Unit | None = ...,
    allow_moisture_gain: bool = ...,
    handling_loss_fraction: float | Decimal = ...,
    allow_negative: bool = ...,
) -> pint.Quantity: ...





@overload
def convert_moisture_mass(
    mass: float | Decimal,
    *,
    source_moisture: float | Decimal | None = ...,
    target_moisture: float | Decimal | None = ...,
    source_moisture_pct: float | Decimal | None = ...,
    target_moisture_pct: float | Decimal | None = ...,
    target_unit: None = ...,
    allow_moisture_gain: bool = ...,
    handling_loss_fraction: float | Decimal = ...,
    allow_negative: bool = ...,
) -> Decimal | float: ...


from core.units.validation import authorized_conversion


@authorized_conversion
def convert_moisture_mass(
    mass: pint.Quantity | float | Decimal,
    *,
    source_moisture: float | Decimal | None = None,
    target_moisture: float | Decimal | None = None,
    source_moisture_pct: float | Decimal | None = None,
    target_moisture_pct: float | Decimal | None = None,
    target_unit: str | Unit | None = None,
    allow_moisture_gain: bool = True,
    handling_loss_fraction: float | Decimal = Decimal(0),
    allow_negative: bool = False,
) -> pint.Quantity | Decimal | float:
    """Universal moisture mass conversion ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``.

    Supports both ``pint.Quantity`` (preserving or converting mass units, including
    semantic ``wmt``, ``dmt``, ``bdmt``, and ``admt`` units) and raw numeric scalars,
    preserving ``Decimal`` precision whenever ``mass`` is ``Decimal``.
    """
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    from core.units.exceptions import InvalidPhysicalStateError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    w1 = _moisture_resolve_moisture(
        source_moisture,
        source_moisture_pct,
        frac_label="source_moisture",
        pct_label="source_moisture_pct",
    )
    w2 = _moisture_resolve_moisture(
        target_moisture,
        target_moisture_pct,
        frac_label="target_moisture",
        pct_label="target_moisture_pct",
    )
    loss = _moisture_to_decimal(handling_loss_fraction, "handling_loss_fraction")
    if loss < Decimal(0) or loss >= Decimal(1):
        raise ValueError(f"handling_loss_fraction must lie in [0, 1), got {handling_loss_fraction!r}.")

    def _apply_exact(m1: Decimal) -> Decimal:
        if w1 == w2 or (w1 < w2 and not allow_moisture_gain):
            return m1
        if loss > Decimal(0):
            factor = max(Decimal(0), ((Decimal(1) - w1) / (Decimal(1) - w2)) - loss)
            return m1 * factor
        if w2 == Decimal(0):
            # Exact wet -> dry: m_dry = m_wet * (1 - w)
            return m1 * (Decimal(1) - w1)
        if w1 == Decimal(0):
            # Exact dry -> wet: m_wet = m_dry / (1 - w)
            return m1 / (Decimal(1) - w2)
        # Exact two-moisture conversion: m_2 = m_1 * (1 - w_1) / (1 - w_2)
        return (m1 * (Decimal(1) - w1)) / (Decimal(1) - w2)

    if isinstance(mass, pint.Quantity):
        validate_registry_quantity(mass)
        if mass.dimensionality not in _moisture_MASS_FAMILY_DIMS:
            raise pint.DimensionalityError(
                mass.units,
                ureg.kilogram,
                mass.dimensionality,
                ureg.kilogram.dimensionality,
            )
        src_tags = set(ureg._semantic_tags(mass.units._units))
        if src_tags & {"contained", "payable"}:
            raise ValueError(f"Moisture conversion requires bulk mass, not {src_tags}-basis unit {mass.units}.")
        if w1 > Decimal(0) and w2 == Decimal(0) and (src_tags & {"dry", "bone_dry"}):
            raise pint.DimensionalityError(
                mass.units,
                ureg.wmt,
                mass.dimensionality,
                ureg.wmt.dimensionality,
            )
        if w1 == Decimal(0) and w2 > Decimal(0) and (src_tags & {"wet"}):
            raise pint.DimensionalityError(
                mass.units,
                ureg.dmt,
                mass.dimensionality,
                ureg.dmt.dimensionality,
            )

        dst = parse_unit(target_unit) if target_unit is not None else mass.units
        if dst.dimensionality not in _moisture_MASS_FAMILY_DIMS:
            raise pint.DimensionalityError(
                ureg.kilogram,
                dst,
                ureg.kilogram.dimensionality,
                dst.dimensionality,
            )
        dst_tags = set(ureg._semantic_tags(dst._units))
        if dst_tags & {"contained", "payable"}:
            raise IncompatibleAnalyticalBasisError(f"Moisture conversion target_unit cannot be {dst_tags}-basis unit {dst}.")
        if w1 > Decimal(0) and w2 == Decimal(0) and ("wet" in dst_tags):
            raise IncompatibleAnalyticalBasisError(f"Wet-to-dry conversion target_unit cannot be wet-basis unit {dst}.")
        if w1 == Decimal(0) and w2 > Decimal(0) and (dst_tags & {"dry", "bone_dry"}):
            raise IncompatibleAnalyticalBasisError(f"Dry-to-wet conversion target_unit cannot be dry-basis unit {dst}.")

        mag_dec = _moisture_to_decimal(mass.magnitude, "mass")
        if not allow_negative and mag_dec < Decimal(0):
            raise InvalidPhysicalStateError(f"Mass must be non-negative, got {mass}.")

        scaled_in_dst = _moisture_to_decimal(
            ureg._convert_mass_scale(mag_dec, mass.units._units, dst._units),
            "scaled_mass",
        )
        converted_dec = _apply_exact(scaled_in_dst)
        return ureg.Quantity(_moisture_preserve_scalar_type(converted_dec, mass), dst)

    if target_unit is not None:
        raise TypeError("target_unit requires mass to be a pint.Quantity.")
    m_dec = _moisture_to_decimal(mass, "mass")
    if not allow_negative and m_dec < Decimal(0):
        raise InvalidPhysicalStateError(f"Mass must be non-negative, got {mass!r}.")
    return _moisture_preserve_scalar_type(_apply_exact(m_dec), mass)





@overload
def wet_to_dry_mass(
    wet_mass: pint.Quantity,
    moisture: float | Decimal | None = ...,
    *,
    moisture_pct: float | Decimal | None = ...,
    target_unit: str | Unit | None = ...,
    allow_negative: bool = ...,
) -> pint.Quantity: ...





@overload
def wet_to_dry_mass(
    wet_mass: float | Decimal,
    moisture: float | Decimal | None = ...,
    *,
    moisture_pct: float | Decimal | None = ...,
    target_unit: None = ...,
    allow_negative: bool = ...,
) -> Decimal | float: ...





def wet_to_dry_mass(
    wet_mass: pint.Quantity | float | Decimal,
    moisture: float | Decimal | None = None,
    *,
    moisture_pct: float | Decimal | None = None,
    target_unit: str | Unit | None = None,
    allow_negative: bool = False,
) -> pint.Quantity | Decimal | float:
    """Universal wet-to-dry mass conversion ``m_dry = m_wet * (1 - w)``."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    resolved_dst = target_unit
    if isinstance(wet_mass, pint.Quantity):
        validate_registry_quantity(wet_mass)
        if resolved_dst is None and wet_mass.units == ureg.wmt:
            resolved_dst = ureg.dmt
    return convert_moisture_mass(
        wet_mass,
        source_moisture=moisture,
        source_moisture_pct=moisture_pct,
        target_moisture=Decimal(0),
        target_unit=resolved_dst,  # type: ignore[arg-type]
        allow_negative=allow_negative,
    )





@overload
def dry_to_wet_mass(
    dry_mass: pint.Quantity,
    moisture: float | Decimal | None = ...,
    *,
    moisture_pct: float | Decimal | None = ...,
    target_unit: str | Unit | None = ...,
    allow_negative: bool = ...,
) -> pint.Quantity: ...





@overload
def dry_to_wet_mass(
    dry_mass: float | Decimal,
    moisture: float | Decimal | None = ...,
    *,
    moisture_pct: float | Decimal | None = ...,
    target_unit: None = ...,
    allow_negative: bool = ...,
) -> Decimal | float: ...





def dry_to_wet_mass(
    dry_mass: pint.Quantity | float | Decimal,
    moisture: float | Decimal | None = None,
    *,
    moisture_pct: float | Decimal | None = None,
    target_unit: str | Unit | None = None,
    allow_negative: bool = False,
) -> pint.Quantity | Decimal | float:
    """Universal dry-to-wet mass conversion ``m_wet = m_dry / (1 - w)``."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    resolved_dst = target_unit
    if isinstance(dry_mass, pint.Quantity):
        validate_registry_quantity(dry_mass)
        if resolved_dst is None and dry_mass.units == ureg.dmt:
            resolved_dst = ureg.wmt
    return convert_moisture_mass(
        dry_mass,
        source_moisture=Decimal(0),
        target_moisture=moisture,
        target_moisture_pct=moisture_pct,
        target_unit=resolved_dst,  # type: ignore[arg-type]
        allow_negative=allow_negative,
    )





def rebase_moisture_intensive(
    value: float | Decimal,
    *,
    source_moisture: float | Decimal | None = None,
    target_moisture: float | Decimal | None = None,
    source_moisture_pct: float | Decimal | None = None,
    target_moisture_pct: float | Decimal | None = None,
    allow_negative: bool = True,
) -> Decimal | float:
    """Re-base an intensive per-unit-mass property (price, GCV, constituent %) across moisture levels.

    Conserves ``value_1 * m_1 = value_2 * m_2``, yielding ``value_2 = value_1 * (1 - w_2) / (1 - w_1)``.
    """
    val_dec = _moisture_to_decimal(value, "value")
    if not allow_negative and val_dec < Decimal(0):
        raise ValueError(f"value must be non-negative, got {value!r}.")
    w1 = _moisture_resolve_moisture(
        source_moisture,
        source_moisture_pct,
        frac_label="source_moisture",
        pct_label="source_moisture_pct",
    )
    w2 = _moisture_resolve_moisture(
        target_moisture,
        target_moisture_pct,
        frac_label="target_moisture",
        pct_label="target_moisture_pct",
    )
    if w1 == w2:
        out = val_dec
    elif w1 == Decimal(0):
        out = val_dec * (Decimal(1) - w2)
    elif w2 == Decimal(0):
        out = val_dec / (Decimal(1) - w1)
    else:
        out = (val_dec * (Decimal(1) - w2)) / (Decimal(1) - w1)
    return _moisture_preserve_scalar_type(out, value)


__all__ = ('AgriculturalQualityBase', 'ProximateAnalysis', 'QualityAuditResult', 'QualityThreshold', 'ThermodynamicBasis', 'UltimateAnalysis', 'convert_moisture_mass', 'dry_matter_fraction', 'dry_to_wet_mass', 'moisture_level_ratio', 'moisture_pct_to_fraction', 'normalize_mass_to_moisture', 'quality_spec_from_assay', 'rebase_constituent_pct', 'rebase_moisture_intensive', 'validate_moisture_fraction', 'wet_to_dry_mass')
_PUBLIC_ALIASES = {'ProximateAnalysis': ('heat_content', 'ProximateAnalysis'), 'ThermodynamicBasis': ('heat_content', 'ThermodynamicBasis'), 'UltimateAnalysis': ('heat_content', 'UltimateAnalysis')}

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
    from core.units.heat_content import ProximateAnalysis
    from core.units.heat_content import ThermodynamicBasis
    from core.units.heat_content import UltimateAnalysis
    from core.units.registry import Unit
