"""Commodex assays energy petroleum implementation and public interfaces."""
from __future__ import annotations

from bisect import bisect_right
from collections.abc import Callable, Iterable
from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, field
from datetime import date
from decimal import Decimal
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from enum import Enum
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final
import inspect
import math
import pint
import warnings

from core.units.validation import validated_assay

# ---------------------------------------------------------------------------
# 2. Hydrocarbon Physical Assay Specification (Crude Oil & Refined Products)
# ---------------------------------------------------------------------------

@validated_assay
@dataclass(frozen=True, slots=True)
class CrudeAssay:
    """Immutable physical specification for a crude oil or refined liquid hydrocarbon.

    Solves API Gravity, Specific Gravity (60°F/60°F), volumetric densities,
    barrel-to-tonne factors and a legacy temperature-only VCF estimate.
    Complete standard algorithms require a supplied validated correction method.
    """

    name: str
    api_gravity: float
    specific_gravity: float
    density_kg_per_m3: float
    density_lb_per_gal: float
    bbl_per_mt: float
    mt_per_bbl: float
    product_group: str = "crude"
    sulfur_percentage: float | None = None
    tan: float | None = None

    def __init__(
        self,
        name: str = "custom_assay",
        *,
        api_gravity: float | None = None,
        specific_gravity: float | None = None,
        bbl_per_mt: float | None = None,
        density_kg_per_m3: float | None = None,
        product_group: str = "crude",
        sulfur_percentage: float | None = None,
        tan: float | None = None,
    ) -> None:
        """Constructs an immutable assay from one density/gravity parameter."""
        from core.units.conversions import MT_PER_BBL_WATER
        from core.units.conversions import WATER_DENSITY_60F_KG_M3
        from core.units.conversions import WATER_DENSITY_60F_LB_GAL
        group = product_group.strip().lower()
        if group not in {"crude", "gasoline", "distillate"}:
            raise ValueError(f"product_group must be 'crude', 'gasoline', or 'distillate', got {product_group!r}")

        provided = sum(p is not None for p in (api_gravity, specific_gravity, density_kg_per_m3, bbl_per_mt))
        if provided != 1:
            raise ValueError(
                "Provide exactly one density parameter: 'api_gravity', 'specific_gravity', 'density_kg_per_m3', or 'bbl_per_mt'."
            )

        resolved_sg: float
        if specific_gravity is not None:
            resolved_sg = float(specific_gravity)
        elif api_gravity is not None:
            api_value = float(api_gravity)
            if not math.isfinite(api_value) or api_value <= -131.5:
                raise ValueError(f"API gravity must be finite and greater than -131.5, got {api_value}")
            resolved_sg = 141.5 / (api_value + 131.5)
        elif density_kg_per_m3 is not None:
            resolved_sg = float(density_kg_per_m3) / WATER_DENSITY_60F_KG_M3
        elif bbl_per_mt is not None:
            bbl_mt = float(bbl_per_mt)
            if not math.isfinite(bbl_mt) or bbl_mt <= 0:
                raise ValueError(f"bbl_per_mt must be finite and strictly positive, got {bbl_mt}")
            resolved_sg = (1.0 / bbl_mt) / MT_PER_BBL_WATER
        else:
            raise ValueError(
                "Must provide at least one density parameter: "
                "'api_gravity', 'specific_gravity', 'density_kg_per_m3', or 'bbl_per_mt'."
            )

        if not math.isfinite(resolved_sg) or resolved_sg <= 0:
            raise ValueError(f"Calculated Specific Gravity must be finite and strictly positive, got {resolved_sg}")

        calc_api = (141.5 / resolved_sg) - 131.5
        calc_rho_kg_m3 = resolved_sg * WATER_DENSITY_60F_KG_M3
        calc_rho_lb_gal = resolved_sg * WATER_DENSITY_60F_LB_GAL
        calc_mt_per_bbl = resolved_sg * MT_PER_BBL_WATER
        calc_bbl_per_mt = 1.0 / calc_mt_per_bbl

        if sulfur_percentage is not None and (sulfur_percentage < 0 or sulfur_percentage > 100):
            raise ValueError(f"Sulfur percentage must be between 0 and 100, got {sulfur_percentage}")
        if tan is not None and tan < 0:
            raise ValueError(f"Total Acid Number (TAN) cannot be negative, got {tan}")

        object.__setattr__(self, "name", str(name))
        object.__setattr__(self, "specific_gravity", (resolved_sg))
        object.__setattr__(self, "api_gravity", (calc_api))
        object.__setattr__(self, "density_kg_per_m3", (calc_rho_kg_m3))
        object.__setattr__(self, "density_lb_per_gal", (calc_rho_lb_gal))
        object.__setattr__(self, "bbl_per_mt", (calc_bbl_per_mt))
        object.__setattr__(self, "mt_per_bbl", (calc_mt_per_bbl))
        object.__setattr__(self, "product_group", group)
        object.__setattr__(self, "sulfur_percentage", sulfur_percentage)
        object.__setattr__(self, "tan", tan)

    @property
    def density_in_air_kg_per_m3(self) -> float:
        """Conventional apparent density using air 1.2 and weights 8000 kg/m³.

        These are declared weighing assumptions, not measured ambient conditions.
        Use physical_to_apparent_mass for explicit air/weight densities.
        """
        return (self.density_kg_per_m3 - 1.2) / (1 - 1.2 / 8000)

    @property
    def gravity_class(self) -> str:
        """Market API gravity classification ('light', 'medium', 'heavy', or 'extra_heavy')."""
        if self.api_gravity > 31.1:
            return "light"
        if self.api_gravity >= 22.3:
            return "medium"
        if self.api_gravity >= 10.0:
            return "heavy"
        return "extra_heavy"

    @property
    def is_sweet(self) -> bool | None:
        """Returns True if sulfur_percentage <= 0.5%, False if > 0.5%, or None if unspecified."""
        if self.sulfur_percentage is None:
            return None
        return self.sulfur_percentage <= 0.5

    def as_context(self, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume] and [mass] for this assay."""
        from core.units.conversions import create_hydrocarbon_context
        return create_hydrocarbon_context(self, name=context_name)

    def as_liquid_spec(self, *, product: LiquidProduct | str | None = None) -> LiquidProductSpec:
        """Expose the assay's 60°F, one-atmosphere density basis explicitly."""
        from core.units.conversions import _as_decimal
        from core.units.registry import ureg
        from core.units.assays.energy.petroleum import LiquidProduct
        from core.units.assays.energy.petroleum import LiquidProductSpec
        from core.units.assays.energy.petroleum import LiquidState
        if product is None:
            if self.product_group == "distillate":
                raise ValueError("Distillate is not an exact correction class; select jet or fuel_oil explicitly.")
            product = self.product_group
        product = LiquidProduct(product)
        allowed = {"crude": {LiquidProduct.CRUDE}, "gasoline": {LiquidProduct.GASOLINE},
                   "distillate": {LiquidProduct.JET, LiquidProduct.FUEL_OIL}}[self.product_group]
        if product not in allowed:
            raise ValueError("Correction product classification conflicts with the assay group.")
        sulfur = None if self.sulfur_percentage is None else _as_decimal(self.sulfur_percentage) / 100
        return LiquidProductSpec(self.name, product, LiquidState(ureg.Quantity(60, "degF"), ureg.Quantity(101325, "Pa")),
                                 ureg.Quantity(_as_decimal(self.density_kg_per_m3), "kg/meter**3"), sulfur)

    def correction_factors(self, observed: LiquidState, *, method: LiquidCorrectionMethod,
                           product: LiquidProduct | str | None = None) -> LiquidCorrectionFactors:
        """Product-specific CTL, CPL and CTPL through a bounded supplied method."""
        from core.units.assays.energy.petroleum import LiquidCorrectionMethod
        if not isinstance(method, LiquidCorrectionMethod):
            raise TypeError("method must be a LiquidCorrectionMethod.")
        if product is None and self.product_group == "distillate" and len(method.products) == 1:
            product = next(iter(method.products))
        return method.factors(self.as_liquid_spec(product=product), observed)

    def bbl_to_mt(self, barrels: float | Decimal, *, in_air: bool = False) -> float:
        """Converts volumetric barrels directly to metric tonnes (in vacuum or in air)."""
        from core.units.conversions import BARREL_TO_M3
        if in_air:
            mt_per_bbl_air = (self.density_in_air_kg_per_m3 * BARREL_TO_M3) / 1000.0
            return float(barrels) * mt_per_bbl_air
        return float(barrels) * self.mt_per_bbl

    def mt_to_bbl(self, metric_tonnes: float | Decimal, *, in_air: bool = False) -> float:
        """Converts metric tonnes (in vacuum or in air) directly to volumetric barrels."""
        from core.units.conversions import BARREL_TO_M3
        if in_air:
            mt_per_bbl_air = (self.density_in_air_kg_per_m3 * BARREL_TO_M3) / 1000.0
            return float(metric_tonnes) / mt_per_bbl_air
        return float(metric_tonnes) * self.bbl_per_mt

    def volume_correction_factor(self, observed_temp_f: float | pint.Quantity) -> float:
        """Legacy temperature-only engineering estimate of VCF to 60°F.

        VCF = exp(-alpha_60 * delta_T * (1 + 0.8 * alpha_60 * delta_T))
        Multiply observed gross volume at `observed_temp_f` (float °F or temperature Quantity)
        by VCF to estimate gross standard volume before water/sediment deductions.
        This omits pressure, full product algorithms, standard rounding and domain
        validation. Use correction_factors with a bounded LiquidCorrectionMethod
        for a declared product method; this function claims no API/ASTM conformity.
        """
        from core.units.conversions import API_VCF_CRUDE_K0
        from core.units.conversions import API_VCF_CRUDE_K1
        from core.units.conversions import API_VCF_CRUDE_K2
        from core.units.conversions import API_VCF_DISTILLATE_K0
        from core.units.conversions import API_VCF_DISTILLATE_K1
        from core.units.conversions import API_VCF_GASOLINE_K0
        from core.units.conversions import API_VCF_GASOLINE_K1
        from core.units.conversions import _absolute_reference_temperature
        if isinstance(observed_temp_f, pint.Quantity):
            _absolute_reference_temperature(observed_temp_f, "Observed liquid temperature")
            temp_f = float(observed_temp_f.to("degF").magnitude)
        else:
            if isinstance(observed_temp_f, bool):
                raise TypeError("Observed temperature must be numeric, not boolean.")
            temp_f = float(observed_temp_f)

        if not math.isfinite(temp_f) or temp_f <= -459.67:
            raise ValueError(f"Observed temperature must be finite and above absolute zero, got {observed_temp_f!r}")

        rho = self.density_kg_per_m3
        if self.product_group == "gasoline":
            k0, k1, k2 = API_VCF_GASOLINE_K0, API_VCF_GASOLINE_K1, 0.0
        elif self.product_group == "distillate":
            k0, k1, k2 = API_VCF_DISTILLATE_K0, API_VCF_DISTILLATE_K1, 0.0
        else:
            k0, k1, k2 = API_VCF_CRUDE_K0, API_VCF_CRUDE_K1, API_VCF_CRUDE_K2

        alpha_60 = (k0 / (rho ** 2)) + (k1 / rho) + k2
        delta_t = temp_f - 60.0
        factor = math.exp(-alpha_60 * delta_t * (1.0 + 0.8 * alpha_60 * delta_t))
        if not math.isfinite(factor) or factor <= 0:
            raise ValueError("Temperature-only VCF estimate is not finite and positive.")
        return factor

    def correct_volume_to_60f(self, observed_volume: pint.Quantity, observed_temp_f: float | pint.Quantity, *, bsw_percentage: float = 0.0,
                              method: LiquidCorrectionMethod | None = None, observed_pressure: pint.Quantity | None = None) -> pint.Quantity:
        """Legacy gross-volume helper; BS&W is a standard-volume percentage.

        Input is GOV, not TOV: free water must already be removed. A supplied
        method requires observed_pressure and uses CTPL; otherwise this remains
        a temperature-only estimate. Use calculate_petroleum_inventory for
        explicit TOV/GOV/GSV/NSV accounting and required water measurements.
        """
        from core.units.conversions import _as_decimal
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        from core.units.registry import validate_dimensionality
        if not isinstance(observed_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(observed_volume).__name__}")
        if not validate_dimensionality(observed_volume, "barrel"):
            raise pint.DimensionalityError(observed_volume.units, ureg.oil_barrel, observed_volume.dimensionality, ureg.oil_barrel.dimensionality,)
        if isinstance(bsw_percentage, bool):
            raise TypeError("BS&W percentage must be numeric, not boolean.")
        bsw = _as_decimal(bsw_percentage)
        if not bsw.is_finite() or not (0 <= bsw < 100):
            raise ValueError(f"bsw_percentage must be in [0, 100), got {bsw_percentage}")
        if (method is None) != (observed_pressure is None):
            raise ValueError("Supply both method and observed_pressure for CTPL correction.")
        if method is not None:
            from core.units.assays.energy.petroleum import LiquidState
            temperature = observed_temp_f if isinstance(observed_temp_f, pint.Quantity) else ureg.Quantity(_as_decimal(observed_temp_f), "degF")
            vcf = self.correction_factors(LiquidState(temperature, observed_pressure), method=method).ctpl
        else:
            vcf = _as_decimal(self.volume_correction_factor(observed_temp_f))
        if not _as_decimal(observed_volume.magnitude).is_finite() or observed_volume.magnitude < 0:
            raise ValueError("Observed gross volume must be finite and non-negative.")
        net_factor = vcf * (1 - bsw / 100)
        scale = _scale_factor(observed_volume.magnitude, net_factor)
        return ureg.Quantity(observed_volume.magnitude * scale, observed_volume.units)

    def net_standard_volume(self, observed_volume: pint.Quantity, observed_temp_f: float | pint.Quantity, bsw_pct: float = 0.0,) -> pint.Quantity:
        """Computes Net Standard Volume (NSV) at 60°F after VCF and BS&W deduction."""
        return self.correct_volume_to_60f(observed_volume, observed_temp_f, bsw_percentage=bsw_pct,)

    def blend_with(self, other: CrudeAssay, *, self_volume_fraction: float, blend_name: str | None = None,) -> CrudeAssay:
        """Computes ideal volumetric blend assay (SG, API, mass-weighted sulfur & TAN)."""
        if not isinstance(other, CrudeAssay):
            raise TypeError(f"Expected CrudeAssay, got {type(other).__name__}")
        if self.product_group != other.product_group:
            raise ValueError("Cross-product blending requires an explicit blend volume model and resulting product classification.")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")
        w_other = 1.0 - self_volume_fraction
        blended_sg = (self.specific_gravity * self_volume_fraction) + (other.specific_gravity * w_other)

        # Sulfur and TAN are mass-fraction properties
        mass_self = self.specific_gravity * self_volume_fraction
        mass_other = other.specific_gravity * w_other
        total_mass = mass_self + mass_other

        blended_sulfur: float | None = None
        if self.sulfur_percentage is not None and other.sulfur_percentage is not None:
            blended_sulfur = (((self.sulfur_percentage * mass_self) + (other.sulfur_percentage * mass_other)) / total_mass)

        blended_tan: float | None = None
        if self.tan is not None and other.tan is not None:
            blended_tan = (((self.tan * mass_self) + (other.tan * mass_other)) / total_mass)

        return CrudeAssay(name=blend_name or f"{self.name}/{other.name} Blend", specific_gravity=blended_sg, product_group=self.product_group, sulfur_percentage=blended_sulfur, tan=blended_tan,)

    # --- Crude Oil Benchmark Presets ---
    @classmethod
    def brent(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.brent")

    @classmethod
    def wti(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.wti")

    @classmethod
    def dubai(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.dubai")

    @classmethod
    def arab_light(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.arab_light")

    @classmethod
    def urals(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.urals")

    @classmethod
    def maya(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.maya")

    @classmethod
    def wcs(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.wcs")

    @classmethod
    def bonny_light(cls) -> CrudeAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "crude_oil", "CrudeAssay.bonny_light")

    # --- Refined Products & Marine Bunker Presets ---
    @classmethod
    def rbob_gasoline(cls) -> CrudeAssay:
        """NYMEX RBOB Reformulated Blendstock (~60.0 API, ~8.45 bbl/mt)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "refined_products", "CrudeAssay.rbob_gasoline")

    @classmethod
    def ulsd_diesel(cls) -> CrudeAssay:
        """Ultra-Low Sulfur Diesel / ICE Gasoil (~36.5 API, ~7.45 bbl/mt)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "refined_products", "CrudeAssay.ulsd_diesel")

    @classmethod
    def jet_a1(cls) -> CrudeAssay:
        """Aviation Turbine Fuel Jet A-1 (~43.5 API, ~7.88 bbl/mt)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "refined_products", "CrudeAssay.jet_a1")

    @classmethod
    def vlsfo_bunker(cls) -> CrudeAssay:
        """IMO 2020 Very Low Sulfur Fuel Oil 0.5% S (~18.5 API, ~6.35 bbl/mt)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "refined_products", "CrudeAssay.vlsfo_bunker")





def _converters_petroleum_number(value, label: str, *, minimum: Decimal | int | None = None) -> Decimal:
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    result = _converters__interpolation_validate_real_scalar(value, label)
    if minimum is not None and result < minimum:
        raise ValueError(f"{label} must be at least {minimum}.")
    return result





def _converters_petroleum_quantity_number(value: pint.Quantity, unit: str, label: str, *, minimum=None) -> Decimal:
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(value, pint.Quantity):
        raise TypeError(f"{label} requires a quantity with explicit units.")
    validate_registry_quantity(value)
    number = _converters_petroleum_number(value.magnitude, label)
    converted = ureg.Quantity(number, value.units).to(unit)
    return _converters_petroleum_number(converted.magnitude, label, minimum=minimum)





def _converters_petroleum_positive_quantity(value: pint.Quantity, unit: str, label: str) -> Decimal:
    result = _converters_petroleum_quantity_number(value, unit, label)
    if result <= 0:
        raise ValueError(f"{label} must be strictly positive.")
    return result





def _converters_petroleum_temperature(value: pint.Quantity) -> pint.Quantity:
    from core.units.conversions import _absolute_reference_temperature
    if not isinstance(value, pint.Quantity):
        raise TypeError("Temperature requires explicit units.")
    _converters_petroleum_number(value.magnitude, "temperature")
    return _absolute_reference_temperature(value, "Liquid temperature")





def _converters_petroleum_pressure(value: pint.Quantity) -> pint.Quantity:
    from core.units.exceptions import InvalidPhysicalStateError
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.conversions import _classify_pressure_datum
    if not isinstance(value, pint.Quantity):
        raise TypeError("Pressure requires explicit units.")
    validate_registry_quantity(value)
    if _classify_pressure_datum(value.units) == "gauge":
        raise InvalidPhysicalStateError("Liquid states require absolute pressure; convert gauge datum explicitly first.")
    return ureg.Quantity(_converters_petroleum_positive_quantity(value, "Pa", "absolute pressure"), "Pa")





def _converters_petroleum_output(number: Decimal, units, source: pint.Quantity) -> pint.Quantity:
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    return ureg.Quantity(_converters__interpolation_match_scalar_type(number, source), units)





def _converters_petroleum_fraction(value, label: str, *, allow_one: bool = True) -> Decimal:
    result = _converters_petroleum_number(value, label, minimum=0)
    if result > 1 or not allow_one and result == 1:
        interval = "[0, 1]" if allow_one else "[0, 1)"
        raise ValueError(f"{label} must be in {interval}.")
    return result





def _converters_petroleum_text(value: str, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} is required.")
    return value.strip()





class LiquidProduct(str, Enum):
    CRUDE = "crude"
    GASOLINE = "gasoline"
    JET = "jet"
    FUEL_OIL = "fuel_oil"
    LUBRICATING_OIL = "lubricating_oil"
    LPG = "lpg"
    NGL = "ngl"
    DENATURED_ETHANOL = "denatured_ethanol"
    ETHANOL_GASOLINE = "ethanol_gasoline"


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LiquidState(_safety_QuantitySnapshot):
    """Single-phase liquid state; neutral pressure units explicitly mean absolute."""
    temperature: pint.Quantity
    absolute_pressure: pint.Quantity
    saturation_pressure: pint.Quantity | None = None

    def __post_init__(self) -> None:
        from core.units.exceptions import InvalidPhysicalStateError
        object.__setattr__(self, "temperature", _converters_petroleum_temperature(self.temperature))
        object.__setattr__(self, "absolute_pressure", _converters_petroleum_pressure(self.absolute_pressure))
        if self.saturation_pressure is not None:
            saturation = _converters_petroleum_pressure(self.saturation_pressure)
            if self.absolute_pressure < saturation:
                raise InvalidPhysicalStateError("Liquid pressure is below the supplied saturation pressure.")
            object.__setattr__(self, "saturation_pressure", saturation)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LiquidProductSpec(_safety_QuantitySnapshot):
    """Density and composition are attached to an explicit reference state."""
    name: str
    product: LiquidProduct
    reference: LiquidState
    reference_density: pint.Quantity
    sulfur_mass_fraction: Decimal | None = None
    ethanol_volume_fraction: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.exceptions import MissingReferenceConditionsError
        from core.units.registry import ureg
        object.__setattr__(self, "name", _converters_petroleum_text(self.name, "product name"))
        object.__setattr__(self, "product", LiquidProduct(self.product))
        if not isinstance(self.reference, LiquidState):
            raise TypeError("reference must be a LiquidState.")
        density = _converters_petroleum_positive_quantity(self.reference_density, "kg/meter**3", "reference density")
        object.__setattr__(self, "reference_density", ureg.Quantity(density, "kg/meter**3"))
        if self.product in {LiquidProduct.LPG, LiquidProduct.NGL} and self.reference.saturation_pressure is None:
            raise MissingReferenceConditionsError("LPG/NGL reference states require explicit saturation pressure.")
        for name in ("sulfur_mass_fraction", "ethanol_volume_fraction"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, _converters_petroleum_fraction(getattr(self, name), name))
        if self.product in {LiquidProduct.DENATURED_ETHANOL, LiquidProduct.ETHANOL_GASOLINE}:
            if self.ethanol_volume_fraction is None:
                raise ValueError("Alcohol corrections require explicit ethanol volume fraction at the reference state.")
        elif self.ethanol_volume_fraction not in {None, Decimal(0)}:
            raise ValueError("Alcohol-containing products require a dedicated ethanol product classification.")

    def as_context(self) -> pint.Context:
        """Mass/volume bridge at this specification's reference state only."""
        context = pint.Context("liquid_reference_density")
        density = self.reference_density.magnitude

        def to_mass(registry, quantity):
            volume = _converters_petroleum_quantity_number(quantity, "meter**3", "reference volume", minimum=0)
            return _converters_petroleum_output(volume * density, "kg", quantity)

        def to_volume(registry, quantity):
            mass = _converters_petroleum_quantity_number(quantity, "kg", "physical mass", minimum=0)
            return _converters_petroleum_output(mass / density, "meter**3", quantity)

        context.add_transformation("[length] ** 3", "[mass]", to_mass)
        context.add_transformation("[mass]", "[length] ** 3", to_volume)
        return context


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LiquidApplicability(_safety_QuantitySnapshot):
    """Inclusive method limits in temperature, absolute pressure and base density."""
    temperature_range: tuple[pint.Quantity, pint.Quantity]
    pressure_range: tuple[pint.Quantity, pint.Quantity]
    density_range: tuple[pint.Quantity, pint.Quantity]
    ethanol_fraction_range: tuple[Decimal, Decimal] | None = None

    def __post_init__(self) -> None:
        from core.units.registry import ureg
        for name in ("temperature_range", "pressure_range", "density_range"):
            if len(getattr(self, name)) != 2:
                raise ValueError(f"{name} requires lower and upper bounds.")
        temperatures = tuple(_converters_petroleum_temperature(value) for value in self.temperature_range)
        pressures = tuple(_converters_petroleum_pressure(value) for value in self.pressure_range)
        densities = tuple(ureg.Quantity(_converters_petroleum_positive_quantity(value, "kg/meter**3", "density bound"), "kg/meter**3") for value in self.density_range)
        for name, values in (("temperature_range", temperatures), ("pressure_range", pressures), ("density_range", densities)):
            if values[0] >= values[1]:
                raise ValueError(f"{name} must be strictly increasing.")
            object.__setattr__(self, name, values)
        if self.ethanol_fraction_range is not None:
            if len(self.ethanol_fraction_range) != 2:
                raise ValueError("ethanol_fraction_range requires two bounds.")
            fractions = tuple(_converters_petroleum_fraction(value, "ethanol fraction bound") for value in self.ethanol_fraction_range)
            if fractions[0] > fractions[1]:
                raise ValueError("Ethanol fraction bounds must be increasing.")
            object.__setattr__(self, "ethanol_fraction_range", fractions)

    def validate(self, spec: LiquidProductSpec, observed: LiquidState) -> None:
        from core.units.exceptions import StandardsApplicabilityError
        for state in (spec.reference, observed):
            if not self.temperature_range[0] <= state.temperature <= self.temperature_range[1]:
                raise StandardsApplicabilityError("Temperature is outside correction-method applicability.")
            if not self.pressure_range[0] <= state.absolute_pressure <= self.pressure_range[1]:
                raise StandardsApplicabilityError("Pressure is outside correction-method applicability.")
        if not self.density_range[0] <= spec.reference_density <= self.density_range[1]:
            raise StandardsApplicabilityError("Reference density is outside correction-method applicability.")
        if spec.product in {LiquidProduct.DENATURED_ETHANOL, LiquidProduct.ETHANOL_GASOLINE}:
            if self.ethanol_fraction_range is None:
                raise StandardsApplicabilityError("Alcohol methods must declare ethanol composition applicability.")
            if not self.ethanol_fraction_range[0] <= spec.ethanol_volume_fraction <= self.ethanol_fraction_range[1]:
                raise StandardsApplicabilityError("Ethanol fraction is outside correction-method applicability.")





@dataclass(frozen=True, slots=True)
class LiquidCorrectionFactors:
    """V_reference = V_observed * CTL * CPL; factors do not remove water/sediment."""
    ctl: Decimal
    cpl: Decimal

    def __post_init__(self) -> None:
        for name in ("ctl", "cpl"):
            number = _converters_petroleum_number(getattr(self, name), name)
            if number <= 0:
                raise ValueError(f"{name} must be strictly positive.")
            object.__setattr__(self, name, number)

    @property
    def ctpl(self) -> Decimal:
        return self.ctl * self.cpl





@dataclass(frozen=True, slots=True)
class LiquidCorrectionMethod:
    """Product-scoped backend, with limits, edition and validation provenance.

    A standard identifier requires validation_reference; this records supplied
    evidence and does not independently certify a third-party algorithm.
    """
    method_id: str
    version: str
    source: str
    products: frozenset[LiquidProduct]
    applicability: LiquidApplicability
    calculator: Callable[[LiquidProductSpec, LiquidState], LiquidCorrectionFactors]
    standard: str | None = None
    validation_reference: str | None = None

    def __post_init__(self) -> None:
        for name in ("method_id", "version", "source"):
            object.__setattr__(self, name, _converters_petroleum_text(getattr(self, name), name))
        products = frozenset(LiquidProduct(value) for value in self.products)
        if not products:
            raise ValueError("Correction method must declare supported products.")
        object.__setattr__(self, "products", products)
        if not isinstance(self.applicability, LiquidApplicability) or not callable(self.calculator):
            raise TypeError("A correction method requires applicability and a callable calculator.")
        if self.standard is not None:
            object.__setattr__(self, "standard", _converters_petroleum_text(self.standard, "standard and edition"))
            object.__setattr__(self, "validation_reference", _converters_petroleum_text(self.validation_reference, "standard validation reference"))

    def factors(self, spec: LiquidProductSpec, observed: LiquidState) -> LiquidCorrectionFactors:
        from core.units.exceptions import IncompatibleCommodityError
        from core.units.exceptions import MissingReferenceConditionsError
        if not isinstance(spec, LiquidProductSpec) or not isinstance(observed, LiquidState):
            raise TypeError("Correction requires a LiquidProductSpec and LiquidState.")
        if spec.product not in self.products:
            raise IncompatibleCommodityError(f"Method {self.method_id!r} does not cover product {spec.product.value!r}.")
        self.applicability.validate(spec, observed)
        if spec.product in {LiquidProduct.LPG, LiquidProduct.NGL}:
            if observed.saturation_pressure is None or spec.reference.saturation_pressure is None:
                raise MissingReferenceConditionsError("LPG/NGL corrections require saturation pressure at both states.")
        factors = self.calculator(spec, observed)
        if not isinstance(factors, LiquidCorrectionFactors):
            raise TypeError("Correction backend must return LiquidCorrectionFactors.")
        same_temperature = observed.temperature == spec.reference.temperature
        same_pressure = observed.absolute_pressure == spec.reference.absolute_pressure
        if same_temperature and same_pressure and abs(factors.cpl - 1) > Decimal("1e-12"):
            raise ValueError("Backend returned non-unity CPL at the reference pressure.")
        if same_temperature and abs(factors.ctl - 1) > Decimal("1e-12"):
            raise ValueError("Backend returned non-unity CTL at the reference temperature.")
        return factors





def constant_coefficient_liquid_method(*, product: LiquidProduct, thermal_expansion: pint.Quantity,
                                      compressibility: pint.Quantity, applicability: LiquidApplicability,
                                      source: str, method_id: str = "constant-coefficient-engineering-v1") -> LiquidCorrectionMethod:
    """Measured-coefficient approximation, not an API/ASTM algorithm or LPG/alcohol fallback."""
    from core.units.exceptions import IncompatibleCommodityError
    product = LiquidProduct(product)
    if product in {LiquidProduct.LPG, LiquidProduct.NGL, LiquidProduct.DENATURED_ETHANOL, LiquidProduct.ETHANOL_GASOLINE}:
        raise IncompatibleCommodityError("Use a dedicated product-specific LPG/NGL or alcohol correction backend.")
    alpha = _converters_petroleum_quantity_number(thermal_expansion, "1/kelvin", "thermal expansion", minimum=0)
    beta = _converters_petroleum_quantity_number(compressibility, "1/Pa", "compressibility", minimum=0)

    def calculate(spec, observed):
        delta_t = observed.temperature.magnitude - spec.reference.temperature.magnitude
        delta_p = observed.absolute_pressure.magnitude - spec.reference.absolute_pressure.magnitude
        return LiquidCorrectionFactors((-alpha * delta_t).exp(), (beta * delta_p).exp())

    return LiquidCorrectionMethod(method_id, "1", source, frozenset({product}), applicability, calculate)





def correct_liquid_volume(volume: pint.Quantity, *, spec: LiquidProductSpec, observed: LiquidState,
                          method: LiquidCorrectionMethod, to_reference: bool = True) -> pint.Quantity:
    """Correct gross liquid volume between the observed and declared reference states."""
    from core.units.registry import ureg
    if not isinstance(to_reference, bool):
        raise TypeError("to_reference must be a boolean.")
    if not isinstance(method, LiquidCorrectionMethod):
        raise TypeError("method must be a LiquidCorrectionMethod.")
    number = _converters_petroleum_quantity_number(volume, "meter**3", "liquid volume", minimum=0)
    factors = method.factors(spec, observed)
    number = number * factors.ctpl if to_reference else number / factors.ctpl
    result = ureg.Quantity(number, "meter**3").to(volume.units)
    return _converters_petroleum_output(result.magnitude, volume.units, volume)





def liquid_density_at_state(spec: LiquidProductSpec, observed: LiquidState, *, method: LiquidCorrectionMethod) -> pint.Quantity:
    """Mass conservation: rho_observed = rho_reference * CTPL."""
    if not isinstance(method, LiquidCorrectionMethod):
        raise TypeError("method must be a LiquidCorrectionMethod.")
    return spec.reference_density * method.factors(spec, observed).ctpl





def infer_liquid_reference_density(measured_density: pint.Quantity, *, name: str, product: LiquidProduct,
                                   reference: LiquidState, observed: LiquidState, method: LiquidCorrectionMethod,
                                   ethanol_volume_fraction: Decimal | None = None,
                                   relative_tolerance: Decimal = Decimal("1e-12"), max_iterations: int = 120) -> pint.Quantity:
    """Invert the backend's density relationship by bounded bisection; never extrapolate."""
    from core.units.exceptions import ConversionConvergenceError
    from core.units.registry import ureg
    if not isinstance(method, LiquidCorrectionMethod):
        raise TypeError("method must be a LiquidCorrectionMethod.")
    measured = _converters_petroleum_positive_quantity(measured_density, "kg/meter**3", "measured density")
    tolerance = _converters_petroleum_number(relative_tolerance, "relative_tolerance")
    if not 0 < tolerance < 1 or isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("Require 0 < relative_tolerance < 1 and positive max_iterations.")
    low, high = (value.magnitude for value in method.applicability.density_range)

    def residual(density):
        from core.units.registry import ureg
        spec = LiquidProductSpec(name, product, reference, ureg.Quantity(density, "kg/meter**3"), ethanol_volume_fraction=ethanol_volume_fraction)
        return liquid_density_at_state(spec, observed, method=method).magnitude - measured

    left, right = residual(low), residual(high)
    for candidate, error in ((low, left), (high, right)):
        if abs(error) <= measured * tolerance:
            result = ureg.Quantity(candidate, "kg/meter**3").to(measured_density.units)
            return _converters_petroleum_output(result.magnitude, result.units, measured_density)
    if left * right >= 0:
        raise ConversionConvergenceError("Measured density cannot be bracketed inside method applicability.")
    for _ in range(max_iterations):
        midpoint = (low + high) / 2
        error = residual(midpoint)
        if abs(error) <= measured * tolerance:
            result = ureg.Quantity(midpoint, "kg/meter**3").to(measured_density.units)
            return _converters_petroleum_output(result.magnitude, result.units, measured_density)
        if left * error > 0:
            low, left = midpoint, error
        else:
            high = midpoint
    raise ConversionConvergenceError("Reference-density solver did not establish convergence.")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class PetroleumInventoryResult(_safety_QuantitySnapshot):
    tov: pint.Quantity
    free_water: pint.Quantity
    gov: pint.Quantity
    gsv: pint.Quantity
    nsv: pint.Quantity
    standard_sediment_and_water: pint.Quantity
    net_physical_mass: pint.Quantity
    factors: LiquidCorrectionFactors
    reference: LiquidState
    observed: LiquidState
    method_id: str
    method_version: str
    method_source: str
    standard: str | None
    validation_reference: str | None
    bsw_fraction: Decimal
    bsw_basis: str





def calculate_petroleum_inventory(total_observed_volume: pint.Quantity, *, free_water_volume: pint.Quantity,
                                  bsw_fraction: Decimal, bsw_basis: str,
                                  spec: LiquidProductSpec, observed: LiquidState,
                                  method: LiquidCorrectionMethod) -> PetroleumInventoryResult:
    """TOV - free water = GOV; GOV*CTPL = GSV; GSV*(1-BS&W) = NSV.

    BS&W is explicitly a fraction of gross standard volume. A laboratory
    observed-volume or mass fraction must be converted separately first.
    """
    from core.units.registry import ureg
    if not isinstance(method, LiquidCorrectionMethod):
        raise TypeError("method must be a LiquidCorrectionMethod.")
    if bsw_basis != "standard_volume":
        raise ValueError("BS&W requires explicit standard_volume basis; other analytical bases need a dedicated conversion.")
    tov = _converters_petroleum_quantity_number(total_observed_volume, "meter**3", "TOV", minimum=0)
    water = _converters_petroleum_quantity_number(free_water_volume, "meter**3", "free water", minimum=0)
    fraction = _converters_petroleum_fraction(bsw_fraction, "BS&W", allow_one=False)
    if water > tov:
        raise ValueError("Free water cannot exceed TOV.")
    factors = method.factors(spec, observed)
    gov = tov - water
    gsv = gov * factors.ctpl
    nsv = gsv * (1 - fraction)

    def volume(value):
        from core.units.registry import ureg
        converted = ureg.Quantity(value, "meter**3").to(total_observed_volume.units)
        return _converters_petroleum_output(converted.magnitude, converted.units, total_observed_volume)

    mass = _converters_petroleum_output(nsv * spec.reference_density.magnitude, "kg", total_observed_volume)
    return PetroleumInventoryResult(volume(tov), volume(water), volume(gov), volume(gsv), volume(nsv),
                                    volume(gsv - nsv), mass, factors, spec.reference, observed,
                                    method.method_id, method.version, method.source, method.standard,
                                    method.validation_reference, fraction, bsw_basis)





def _converters_petroleum_buoyancy_factor(*, liquid_density: pint.Quantity, air_density: pint.Quantity, weight_density: pint.Quantity) -> Decimal:
    liquid = _converters_petroleum_positive_quantity(liquid_density, "kg/meter**3", "liquid density")
    air = _converters_petroleum_quantity_number(air_density, "kg/meter**3", "air density", minimum=0)
    weights = _converters_petroleum_positive_quantity(weight_density, "kg/meter**3", "reference-weight density")
    if air >= liquid or air >= weights:
        raise ValueError("Air density must be below both liquid and reference-weight densities.")
    return (1 - air / liquid) / (1 - air / weights)





def physical_to_apparent_mass(mass: pint.Quantity, *, liquid_density: pint.Quantity,
                              air_density: pint.Quantity, weight_density: pint.Quantity) -> pint.Quantity:
    """Balance-equivalent mass in air, relative to stated calibration weights; not a force."""
    _converters_petroleum_quantity_number(mass, "kg", "physical mass", minimum=0)
    factor = _converters_petroleum_buoyancy_factor(liquid_density=liquid_density, air_density=air_density, weight_density=weight_density)
    return _converters_petroleum_output(_converters_petroleum_number(mass.magnitude, "mass") * factor, mass.units, mass)





def apparent_to_physical_mass(mass: pint.Quantity, *, liquid_density: pint.Quantity,
                              air_density: pint.Quantity, weight_density: pint.Quantity) -> pint.Quantity:
    """Inverse of physical_to_apparent_mass under the same weighing conditions."""
    _converters_petroleum_quantity_number(mass, "kg", "apparent mass", minimum=0)
    factor = _converters_petroleum_buoyancy_factor(liquid_density=liquid_density, air_density=air_density, weight_density=weight_density)
    return _converters_petroleum_output(_converters_petroleum_number(mass.magnitude, "mass") / factor, mass.units, mass)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class TankCalibrationPoint(_safety_QuantitySnapshot):
    level: pint.Quantity
    volume: pint.Quantity
    standard_uncertainty: pint.Quantity

    def __post_init__(self) -> None:
        from core.units.registry import ureg
        for name, unit in (("level", "meter"), ("volume", "meter**3"), ("standard_uncertainty", "meter**3")):
            value = _converters_petroleum_quantity_number(getattr(self, name), unit, name, minimum=0)
            object.__setattr__(self, name, ureg.Quantity(value, unit))


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class TankVolumeReading(_safety_QuantitySnapshot):
    volume: pint.Quantity
    standard_uncertainty: pint.Quantity
    calibration_temperature: pint.Quantity
    tank_id: str
    datum_id: str
    certificate_id: str


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class TankCalibrationTable(_safety_QuantitySnapshot):
    """Certified-datum, piecewise-linear capacity table with explicit covariance.

    Uncertainty is a standard uncertainty. This table does not apply tank-shell
    thermal expansion, roof displacement or an observed-to-standard liquid factor.
    """
    tank_id: str
    datum_id: str
    certificate_id: str
    calibration_temperature: pint.Quantity
    points: tuple[TankCalibrationPoint, ...]
    endpoint_correlation: Decimal
    valid_from: date | None = None
    valid_to: date | None = None

    def __post_init__(self) -> None:
        from core.units.currencies import _date
        for name in ("tank_id", "datum_id", "certificate_id"):
            object.__setattr__(self, name, _converters_petroleum_text(getattr(self, name), name))
        object.__setattr__(self, "calibration_temperature", _converters_petroleum_temperature(self.calibration_temperature))
        points = tuple(self.points)
        if len(points) < 2 or not all(isinstance(point, TankCalibrationPoint) for point in points):
            raise ValueError("Tank calibration requires at least two explicit calibration points.")
        for left, right in zip(points, points[1:]):
            if right.level <= left.level or right.volume <= left.volume:
                raise ValueError("Calibration levels and volumes must be strictly increasing.")
        object.__setattr__(self, "points", points)
        correlation = _converters_petroleum_number(self.endpoint_correlation, "endpoint correlation")
        if not -1 <= correlation <= 1:
            raise ValueError("Endpoint correlation must be in [-1, 1].")
        object.__setattr__(self, "endpoint_correlation", correlation)
        for name in ("valid_from", "valid_to"):
            if getattr(self, name) is not None:
                object.__setattr__(self, name, _date(getattr(self, name)))
        if self.valid_from and self.valid_to and self.valid_from >= self.valid_to:
            raise ValueError("Calibration validity must be a non-empty half-open interval.")

    def volume_at(self, level: pint.Quantity, *, datum_id: str, level_standard_uncertainty: pint.Quantity,
                  measurement_date: date | str | None = None,
                  calibration_level_correlation: Decimal = Decimal(0)) -> TankVolumeReading:
        from core.units.currencies import _date
        from core.units.registry import ureg
        if datum_id != self.datum_id:
            raise ValueError("Measured level datum does not match the calibration certificate.")
        if measurement_date is not None:
            measurement_date = _date(measurement_date)
        if self.valid_from or self.valid_to:
            if measurement_date is None:
                raise ValueError("A dated calibration table requires measurement_date.")
            when = _date(measurement_date)
            if self.valid_from and when < self.valid_from or self.valid_to and when >= self.valid_to:
                raise ValueError("Calibration table is not effective on the measurement date.")
        measured = _converters_petroleum_quantity_number(level, "meter", "measured level", minimum=0)
        level_u = _converters_petroleum_quantity_number(level_standard_uncertainty, "meter", "level standard uncertainty", minimum=0)
        correlation = _converters_petroleum_number(calibration_level_correlation, "calibration/level correlation")
        if not -1 <= correlation <= 1:
            raise ValueError("Calibration/level correlation must be in [-1, 1].")
        levels = tuple(point.level.magnitude for point in self.points)
        if measured < levels[0] or measured > levels[-1]:
            raise ValueError("Measured level is outside the calibration table; extrapolation is prohibited.")
        index = min(max(bisect_right(levels, measured) - 1, 0), len(levels) - 2)
        left, right = self.points[index:index + 2]
        weight = (measured - left.level.magnitude) / (right.level.magnitude - left.level.magnitude)
        slope = (right.volume.magnitude - left.volume.magnitude) / (right.level.magnitude - left.level.magnitude)
        if measured == levels[index] and index > 0:
            previous = self.points[index - 1]
            slope = max(slope, (left.volume.magnitude - previous.volume.magnitude) / (left.level.magnitude - previous.level.magnitude))
        volume = left.volume.magnitude + weight * (right.volume.magnitude - left.volume.magnitude)
        a = (1 - weight) * left.standard_uncertainty.magnitude
        b = weight * right.standard_uncertainty.magnitude
        calibration_u = ((a - b) ** 2 + 2 * a * b * (1 + self.endpoint_correlation)).sqrt()
        level_volume_u = slope * level_u
        combined_u = ((calibration_u - level_volume_u) ** 2 + 2 * calibration_u * level_volume_u * (1 + correlation)).sqrt()
        return TankVolumeReading(ureg.Quantity(volume, "meter**3"), ureg.Quantity(combined_u, "meter**3"),
                                 self.calibration_temperature, self.tank_id, self.datum_id, self.certificate_id)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class LiquidBlendComponent(_safety_QuantitySnapshot):
    spec: LiquidProductSpec
    amount: pint.Quantity
    kinematic_viscosity: pint.Quantity | None = None
    viscosity_temperature: pint.Quantity | None = None

    def __post_init__(self) -> None:
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        from core.units.validation import _safety_copy_quantities
        if not isinstance(self.spec, LiquidProductSpec):
            raise TypeError("Blend components require a LiquidProductSpec.")
        if not isinstance(self.amount, pint.Quantity):
            raise TypeError("Blend amount requires mass or reference-volume units.")
        validate_registry_quantity(self.amount)
        if self.amount.dimensionality not in (ureg.kilogram.dimensionality, (ureg.meter**3).dimensionality):
            raise ValueError("Blend amounts must be physical mass or volume at the product reference state.")
        if _converters_petroleum_number(self.amount.magnitude, "blend amount") <= 0:
            raise ValueError("Blend component amounts must be strictly positive.")
        object.__setattr__(self, "amount", _safety_copy_quantities(self.amount))
        if (self.kinematic_viscosity is None) != (self.viscosity_temperature is None):
            raise ValueError("Kinematic viscosity and its measurement temperature must be supplied together.")
        if self.kinematic_viscosity is not None:
            viscosity = _converters_petroleum_positive_quantity(self.kinematic_viscosity, "centistokes", "kinematic viscosity")
            object.__setattr__(self, "kinematic_viscosity", ureg.Quantity(viscosity, "centistokes"))
            object.__setattr__(self, "viscosity_temperature", _converters_petroleum_temperature(self.viscosity_temperature))


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class ResolvedLiquidBlendComponent(_safety_QuantitySnapshot):
    spec: LiquidProductSpec
    mass: pint.Quantity
    reference_volume: pint.Quantity





@dataclass(frozen=True, slots=True)
class BlendVolumeModel:
    method_id: str
    source: str
    calculator: Callable[[tuple[ResolvedLiquidBlendComponent, ...]], pint.Quantity]

    def __post_init__(self) -> None:
        object.__setattr__(self, "method_id", _converters_petroleum_text(self.method_id, "blend volume method"))
        object.__setattr__(self, "source", _converters_petroleum_text(self.source, "blend model source"))
        if not callable(self.calculator):
            raise TypeError("Blend volume calculator must be callable.")





def ideal_blend_volume(components: tuple[ResolvedLiquidBlendComponent, ...]) -> pint.Quantity:
    """Explicit engineering assumption of additive reference volumes, without contraction."""
    from core.units.registry import ureg
    return ureg.Quantity(sum((item.reference_volume.magnitude for item in components), Decimal(0)), "meter**3")


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class PetroleumBlendResult(_safety_QuantitySnapshot):
    reference: LiquidState
    total_mass: pint.Quantity
    additive_volume: pint.Quantity
    blend_volume: pint.Quantity
    reference_density: pint.Quantity
    sulfur_mass_fraction: Decimal | None
    kinematic_viscosity: pint.Quantity | None
    viscosity_temperature: pint.Quantity | None
    volume_method_id: str
    volume_method_source: str
    viscosity_method: str | None
    components: tuple[ResolvedLiquidBlendComponent, ...]

    @property
    def api_gravity(self) -> Decimal:
        """API is derived nonlinearly from density, only on its 60°F basis."""
        from core.units.exceptions import StandardsApplicabilityError
        from core.units.registry import ureg
        reference_k = ureg.Quantity(Decimal(60), "degF").to("kelvin").magnitude
        if abs(self.reference.temperature.magnitude - reference_k) > Decimal("1e-9") or self.reference.absolute_pressure.magnitude != 101325:
            raise StandardsApplicabilityError("API gravity requires density at 60°F and one standard atmosphere.")
        return Decimal("141.5") / (self.reference_density.magnitude / Decimal("999.016")) - Decimal("131.5")





def blend_petroleum(components: Iterable[LiquidBlendComponent], *, volume_model: BlendVolumeModel,
                     viscosity_method: str | None = None) -> PetroleumBlendResult:
    """Conserve mass and sulfur; require explicit volume and nonlinear-viscosity models."""
    from core.units.registry import ureg
    components = tuple(components)
    if not components or not all(isinstance(item, LiquidBlendComponent) for item in components):
        raise ValueError("Supply at least one explicit liquid blend component.")
    if not isinstance(volume_model, BlendVolumeModel):
        raise TypeError("Supply an explicit BlendVolumeModel; no additive-volume assumption is inferred.")
    initial_reference = components[0].spec.reference
    if any(item.spec.reference.temperature != initial_reference.temperature or item.spec.reference.absolute_pressure != initial_reference.absolute_pressure for item in components[1:]):
        raise ValueError("Normalize blend components to the same temperature and pressure reference first.")
    # Saturation pressure belongs to each composition. A blend does not inherit
    # the first component's bubble point; a separate mixture model must supply it.
    reference = LiquidState(initial_reference.temperature, initial_reference.absolute_pressure)
    if any(item.spec.product in {LiquidProduct.DENATURED_ETHANOL, LiquidProduct.ETHANOL_GASOLINE} for item in components) and volume_model.calculator is ideal_blend_volume:
        raise ValueError("Alcohol blends require a dedicated non-ideal volume model.")
    resolved = []
    for item in components:
        density = item.spec.reference_density.magnitude
        if item.amount.dimensionality == ureg.kilogram.dimensionality:
            mass = _converters_petroleum_quantity_number(item.amount, "kg", "blend mass")
            volume = mass / density
        else:
            volume = _converters_petroleum_quantity_number(item.amount, "meter**3", "reference blend volume")
            mass = volume * density
        resolved.append(ResolvedLiquidBlendComponent(item.spec, ureg.Quantity(mass, "kg"), ureg.Quantity(volume, "meter**3")))
    resolved = tuple(resolved)
    total_mass = sum((item.mass.magnitude for item in resolved), Decimal(0))
    additive = ideal_blend_volume(resolved)
    volume = _converters_petroleum_positive_quantity(volume_model.calculator(resolved), "meter**3", "blended reference volume")
    sulfur = None
    if all(item.spec.sulfur_mass_fraction is not None for item in resolved):
        sulfur = sum((item.mass.magnitude * item.spec.sulfur_mass_fraction for item in resolved), Decimal(0)) / total_mass
    viscosity, viscosity_temperature = None, None
    viscosities = [item.kinematic_viscosity for item in components]
    if any(value is not None for value in viscosities):
        if not all(value is not None for value in viscosities) or viscosity_method != "refutas_mass":
            raise ValueError("Complete viscosity data and explicit refutas_mass selection are required.")
        viscosity_temperature = components[0].viscosity_temperature
        if any(item.viscosity_temperature != viscosity_temperature for item in components):
            raise ValueError("Viscosity blending requires a common measurement temperature.")
        vbn = Decimal(0)
        for item, liquid in zip(components, resolved):
            nu = item.kinematic_viscosity.magnitude
            if nu <= Decimal(".2"):
                raise ValueError("Refutas viscosity requires more than 0.2 cSt.")
            index = Decimal("14.534") * (nu + Decimal(".8")).ln().ln() + Decimal("10.975")
            vbn += index * liquid.mass.magnitude / total_mass
        nu = ((vbn - Decimal("10.975")) / Decimal("14.534")).exp().exp() - Decimal(".8")
        viscosity = ureg.Quantity(nu, "centistokes")
    elif viscosity_method is not None:
        raise ValueError("A viscosity method requires measured component viscosities.")
    return PetroleumBlendResult(reference, ureg.Quantity(total_mass, "kg"), additive, ureg.Quantity(volume, "meter**3"),
                                ureg.Quantity(total_mass / volume, "kg/meter**3"), sulfur, viscosity,
                                viscosity_temperature, volume_model.method_id, volume_model.source, viscosity_method, resolved)


__all__ = ('BlendVolumeModel', 'CrudeAssay', 'DistillationCut', 'LiquidApplicability', 'LiquidBlendComponent', 'LiquidCorrectionFactors', 'LiquidCorrectionMethod', 'LiquidProduct', 'LiquidProductSpec', 'LiquidState', 'PetroleumBlendResult', 'PetroleumInventoryResult', 'RefineryCrudeAssay', 'ResolvedLiquidBlendComponent', 'TankCalibrationPoint', 'TankCalibrationTable', 'TankVolumeReading', 'apparent_to_physical_mass', 'blend_petroleum', 'calculate_petroleum_inventory', 'constant_coefficient_liquid_method', 'correct_liquid_volume', 'ideal_blend_volume', 'infer_liquid_reference_density', 'liquid_density_at_state', 'physical_to_apparent_mass')
_PUBLIC_ALIASES = {'RefineryCrudeAssay': ('assays.energy.refining', 'RefineryCrudeAssay'), 'DistillationCut': ('assays.energy.refining', 'DistillationCut')}

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
    from core.units.assays.energy.refining import DistillationCut
    from core.units.assays.energy.refining import RefineryCrudeAssay
