"""Commodex assays agriculture vegetable_oils implementation and public interfaces."""
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
from core.units.validation import validated_assay

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class VegetableOilAssay(_agricultural_families_AgriculturalFamilyAssay):
    refining_recovery_fraction: Decimal | None = None
    thermal_expansion_per_kelvin: Decimal | None = None
    temperature_range: tuple[pint.Quantity, pint.Quantity] | None = None
    _commodities = frozenset({"palm", "soybean", "sunflower", "rapeseed", "coconut", "groundnut"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "refining_recovery_fraction", _agricultural_families_yield(self.refining_recovery_fraction, "refining recovery"))
        if self.thermal_expansion_per_kelvin is not None:
            alpha = _converters_petroleum_number(self.thermal_expansion_per_kelvin, "thermal expansion", minimum=0)
            object.__setattr__(self, "thermal_expansion_per_kelvin", alpha)
            if self.reference_temperature is None or self.temperature_range is None:
                raise ValueError("Expansion correction requires reference temperature and explicit applicability bounds.")
        if self.temperature_range is not None:
            bounds = tuple(_converters_petroleum_temperature(q) for q in self.temperature_range)
            if len(bounds) != 2 or bounds[0] >= bounds[1]:
                raise ValueError("Temperature bounds must be an increasing pair.")
            if self.reference_temperature is not None and not bounds[0] <= self.reference_temperature <= bounds[1]:
                raise ValueError("Reference temperature is outside applicability bounds.")
            object.__setattr__(self, "temperature_range", bounds)

    def refined_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.refining_recovery_fraction, "refining recovery")

    def crude_oil_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.refining_recovery_fraction, "refining recovery")

    def volume_correction_factor(self, observed_temp):
        """Declared constant-expansion approximation; numeric temperatures are C."""
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        from core.units.registry import ureg
        if self.thermal_expansion_per_kelvin is None:
            raise ValueError("Declare thermal expansion and its temperature applicability.")
        observed = _converters_petroleum_temperature(observed_temp if isinstance(observed_temp, pint.Quantity) else ureg.Quantity(_converters_petroleum_number(observed_temp, "temperature"), "degC"))
        if not self.temperature_range[0] <= observed <= self.temperature_range[1]:
            raise ValueError("Observed temperature is outside declared applicability.")
        denominator = 1 + self.thermal_expansion_per_kelvin * (observed.magnitude - self.reference_temperature.magnitude)
        if denominator <= 0:
            raise ValueError("Expansion model gives nonpositive volume.")
        return 1 / denominator


from core.units.analytical_basis import AgriculturalQualityBase


# ---------------------------------------------------------------------------
# 12. Vegetable Oil Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class VegetableOilQualitySpec(AgriculturalQualityBase):
    """Physical quality specification for crude and refined edible vegetable oils."""

    free_fatty_acid_pct: Decimal
    peroxide_value_meq_per_kg: Decimal
    iodine_value_g_per_100g: Decimal
    moisture_and_impurities_pct: Decimal = Decimal("0.10")  # MIU %
    ffa_basis: Literal["oleic", "palmitic", "lauric"] = "oleic"

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_FFA_MOLAR_MASS_G_MOL
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        if self.ffa_basis not in _agricultural_quality_FFA_MOLAR_MASS_G_MOL:
            raise ValueError(f"ffa_basis must be 'oleic', 'palmitic', or 'lauric', got {self.ffa_basis!r}.")
        ffa = _agricultural_quality_pct(self.free_fatty_acid_pct, "free_fatty_acid_pct", allow_high=False)
        pv = _agricultural_quality_nonneg(self.peroxide_value_meq_per_kg, "peroxide_value_meq_per_kg")
        iv = _agricultural_quality_pct(self.iodine_value_g_per_100g, "iodine_value_g_per_100g", high=Decimal(260))
        miu = _agricultural_quality_pct(self.moisture_and_impurities_pct, "moisture_and_impurities_pct", high=Decimal(25))
        if ffa + miu >= Decimal(100):
            raise ValueError("free_fatty_acid_pct + moisture_and_impurities_pct must be < 100%.")

        object.__setattr__(self, "free_fatty_acid_pct", ffa)
        object.__setattr__(self, "peroxide_value_meq_per_kg", pv)
        object.__setattr__(self, "iodine_value_g_per_100g", iv)
        object.__setattr__(self, "moisture_and_impurities_pct", miu)

    @property
    def acid_value_mg_koh_per_g(self) -> Decimal:
        """Exact stoichiometric Acid Value (mg KOH / g oil) from ``free_fatty_acid_pct`` and ``ffa_basis``."""
        from core.units.analytical_basis import _agricultural_quality_FFA_MOLAR_MASS_G_MOL
        from core.units.analytical_basis import _agricultural_quality_KOH_MG_PER_MOL
        mw = _agricultural_quality_FFA_MOLAR_MASS_G_MOL[self.ffa_basis]
        return (self.free_fatty_acid_pct / Decimal(100)) * (_agricultural_quality_KOH_MG_PER_MOL / mw)

    def ffa_on_basis(self, target_basis: Literal["oleic", "palmitic", "lauric"]) -> Decimal:
        from core.units.analytical_basis import _agricultural_quality_FFA_MOLAR_MASS_G_MOL
        from core.units.analytical_basis import _agricultural_quality_KOH_MG_PER_MOL
        if target_basis not in _agricultural_quality_FFA_MOLAR_MASS_G_MOL:
            raise ValueError(f"Unsupported target_basis: {target_basis!r}.")
        return self.acid_value_mg_koh_per_g * _agricultural_quality_FFA_MOLAR_MASS_G_MOL[target_basis] / _agricultural_quality_KOH_MG_PER_MOL * Decimal(100)

    def neutral_oil_yield_pct(self, *, refining_loss_factor: Decimal | float = Decimal("1.4")) -> Decimal:
        """Estimated refined neutral oil recovery ``100 - MIU - k_ref * FFA``."""
        from core.units.analytical_basis import _agricultural_quality_nonneg
        k = _agricultural_quality_nonneg(refining_loss_factor, "refining_loss_factor")
        return max(Decimal(0), Decimal(100) - self.moisture_and_impurities_pct - k * self.free_fatty_acid_pct)

    def net_pure_oil_mass(
        self,
        gross_mass: pint.Quantity,
        *,
        deduct_refining_ffa_loss: bool = False,
        refining_loss_factor: Decimal | float = Decimal("1.4"),
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_kg_out
        from core.units.analytical_basis import _agricultural_quality_mass_kg
        kg = _agricultural_quality_mass_kg(gross_mass, "gross_mass")
        if deduct_refining_ffa_loss:
            frac = self.neutral_oil_yield_pct(refining_loss_factor=refining_loss_factor) / Decimal(100)
        else:
            frac = (Decimal(100) - self.moisture_and_impurities_pct) / Decimal(100)
        return _agricultural_quality_kg_out(kg * frac, target_unit or gross_mass.units, gross_mass)

    def blend_with(self, other: VegetableOilQualitySpec, *, self_mass_fraction: Decimal | float) -> VegetableOilQualitySpec:
        from core.units.analytical_basis import _agricultural_quality_blend_weights
        if not isinstance(other, VegetableOilQualitySpec):
            raise TypeError(f"Expected VegetableOilQualitySpec, got {type(other).__name__}.")
        w1, w2 = _agricultural_quality_blend_weights(self_mass_fraction)
        other_ffa_aligned = other.ffa_on_basis(self.ffa_basis)
        return VegetableOilQualitySpec(
            free_fatty_acid_pct=w1 * self.free_fatty_acid_pct + w2 * other_ffa_aligned,
            peroxide_value_meq_per_kg=w1 * self.peroxide_value_meq_per_kg + w2 * other.peroxide_value_meq_per_kg,
            iodine_value_g_per_100g=w1 * self.iodine_value_g_per_100g + w2 * other.iodine_value_g_per_100g,
            moisture_and_impurities_pct=w1 * self.moisture_and_impurities_pct + w2 * other.moisture_and_impurities_pct,
            ffa_basis=self.ffa_basis,
        )

    def to_material_spec(self, name: str = "Vegetable Oil") -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import CompositionBasis
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        miu_f = self.moisture_and_impurities_pct / Decimal(100)
        ffa_f = self.free_fatty_acid_pct / Decimal(100)
        pur_f = Decimal(1) - miu_f
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(purity=pur_f, moisture=miu_f),
            composition=MaterialComposition(
                components={"neutral_triglycerides": max(Decimal(0), pur_f - ffa_f)},
                basis=CompositionBasis.AS_RECEIVED,
                moisture_fraction=miu_f,
                impurities={"FFA": ffa_f},
            ),
        )


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (VegetableOilAssay,))

from core.units.agricultural_specs import Decimal, _LITRES_PER_US_GALLON, _OLIVE_OIL_DENSITY_KG_PER_M3_20C, _OLIVE_OIL_EXPANSION_PER_C, _validate_finite_float, create_olive_oil_context, pint
@validated_assay
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
        from core.units.exceptions import StandardsApplicabilityError
        from core.units.registry import ureg
        value = observed_temp.to(ureg.degC).magnitude if isinstance(observed_temp, pint.Quantity) else observed_temp
        temperature = Decimal(str(value))
        if not temperature.is_finite() or not (-20 <= temperature <= 100):
            raise StandardsApplicabilityError("Observed temperature must be finite and in [-20, 100] degrees C.")
        delta = temperature - Decimal(str(self.reference_temp_c))
        return 1 - Decimal(str(_OLIVE_OIL_EXPANSION_PER_C)) * delta

    def _effective_density_kg_per_m3(
        self, temp_c: float | Decimal | pint.Quantity | None = None,
    ) -> Decimal:
        factor = Decimal(1) if temp_c is None else self.volume_correction_factor(temp_c)
        return Decimal(str(self.density_kg_per_m3)) * factor

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
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
        from core.units.registry import ureg
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
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self._effective_density_kg_per_m3(temp_c) / 1000
        scale = mt_per_m3 if isinstance(mt_q.magnitude, Decimal) else float(mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    # --- Refining Yield Conversions ---
    def crude_to_refined_mass(self, crude_mass: pint.Quantity) -> pint.Quantity:
        """Converts crude/lampante olive oil mass to refined oil mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(crude_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(crude_mass).__name__}")
        mt_q = self._validate_mass(crude_mass)
        fraction = Decimal(str(self.refining_yield_pct)) / 100
        scale = fraction if isinstance(mt_q.magnitude, Decimal) else float(fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def refined_to_crude_mass(self, refined_mass: pint.Quantity) -> pint.Quantity:
        """Converts refined olive oil mass back to required crude/lampante mass (metric_ton)."""
        from core.units.registry import ureg
        if not isinstance(refined_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(refined_mass).__name__}")
        mt_q = self._validate_mass(refined_mass)
        fraction = Decimal(str(self.refining_yield_pct)) / 100
        scale = fraction if isinstance(mt_q.magnitude, Decimal) else float(fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    # --- Pricing Parity Helpers (same currency; FX is handled by the routers) ---
    @staticmethod
    def _price_result(price: float | Decimal, factor: Decimal, *, divide: bool = False) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        value = Decimal(str(price))
        if not value.is_finite():
            raise ValueError("Olive oil price must be finite.")
        value = value / factor if divide else value * factor
        return preserve_ieee754_scalar(value, price)

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
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "OliveOilAssay.extra_virgin_evoo")

    @classmethod
    def lampante_crude(cls, *, free_fatty_acid_pct: float = 3.0, refining_yield_pct: float = 95.0) -> OliveOilAssay:
        """Lampante crude olive oil requiring physical/chemical refining."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "OliveOilAssay.lampante_crude", overrides={'name': f'Lampante Olive Oil ({free_fatty_acid_pct:g}% FFA)', 'free_fatty_acid_pct': free_fatty_acid_pct, 'refining_yield_pct': refining_yield_pct})

    @classmethod
    def refined_olive_oil(cls) -> OliveOilAssay:
        """Refined olive oil benchmark (low FFA, ~914 kg/m^3 at 20°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "OliveOilAssay.refined_olive_oil")

__all__ = ('VegetableOilAssay', 'VegetableOilQualitySpec', 'convert_vegetable_oil', 'create_vegetable_oil_context', 'OliveOilAssay')
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
    from core.units.assays.agriculture.vegetable_oils import convert_vegetable_oil
    from core.units.assays.agriculture.vegetable_oils import create_vegetable_oil_context
