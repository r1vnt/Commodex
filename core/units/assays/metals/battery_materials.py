"""Commodex assays metals battery_materials implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final
import inspect
import math
import pint
import warnings

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class GraphiteAssay(MaterialSpecification):
    fixed_carbon_fraction: Decimal
    product_form: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.composition import ElementComposition
        fraction = _converters_petroleum_fraction(self.fixed_carbon_fraction, "fixed carbon")
        if self.product_form not in {"flake", "amorphous", "spherical", "purified", "synthetic"}:
            raise ValueError("Unknown graphite product form.")
        if self.composition is None:
            object.__setattr__(self, "composition", ElementComposition({"C": fraction}))
        elif self.composition.fraction("C") != fraction:
            raise ValueError("Fixed carbon conflicts with assayed composition.")
        object.__setattr__(self, "fixed_carbon_fraction", fraction)
        MaterialSpecification.__post_init__(self)

    def purified_mass(self, quantity, *, target_fraction, recovery_fraction, mass_basis="dry"):
        from core.units.assays.base import beneficiated_mass
        return beneficiated_mass(self, quantity, component="C", target_fraction=target_fraction,
                                  recovery_fraction=recovery_fraction, mass_basis=mass_basis)


from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class BatteryMaterialAssay:
    """Immutable physical and stoichiometric specification for battery raw materials.

    Bridges Wet Metric Tonnes (WMT), Dry Metric Tonnes (DMT), contained battery
    metal (Li, Ni, Co, Mn), and Lithium Carbonate Equivalent (LCE) including
    refinery metallurgical recovery rates.
    """

    name: str
    active_metal: str
    active_metal_pct: float
    moisture_pct: float = 0.0
    purity_pct: float | None = None
    lce_factor_per_dmt: float | None = None
    default_recovery_rate: float = 1.0

    def __post_init__(self) -> None:
        if not self.active_metal or not self.active_metal.strip():
            raise ValueError("active_metal must be a non-empty string (e.g. 'Li', 'Ni', 'Co', 'Mn').")
        if not math.isfinite(self.active_metal_pct) or not (0.0 < self.active_metal_pct <= 100.0):
            raise ValueError(f"active_metal_pct must be finite and in (0, 100], got {self.active_metal_pct}")
        if not math.isfinite(self.moisture_pct) or not (0.0 <= self.moisture_pct < 60.0):
            raise ValueError(f"moisture_pct must be finite and in [0, 60), got {self.moisture_pct}")
        if self.purity_pct is not None:
            if not math.isfinite(self.purity_pct) or not (0.0 < self.purity_pct <= 100.0):
                raise ValueError(f"purity_pct must be finite and in (0, 100], got {self.purity_pct}")
        if self.lce_factor_per_dmt is not None:
            if not math.isfinite(self.lce_factor_per_dmt) or self.lce_factor_per_dmt <= 0.0:
                raise ValueError(f"lce_factor_per_dmt must be strictly positive, got {self.lce_factor_per_dmt}")
        if not math.isfinite(self.default_recovery_rate) or not (0.0 < self.default_recovery_rate <= 1.0):
            raise ValueError(f"default_recovery_rate must be in (0, 1], got {self.default_recovery_rate}")

    # --- Core Stoichiometric & Dry-Basis Properties ---
    @property
    def dry_mass_fraction(self) -> float:
        """Ratio of Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT)."""
        from core.units.analytical_basis import dry_matter_fraction
        return float(dry_matter_fraction(moisture_pct=self.moisture_pct))

    @property
    def dmt_per_mt_contained_metal(self) -> float:
        """Dry Metric Tonnes of material required per 1 mt of contained active metal."""
        return (100.0 / self.active_metal_pct)

    @property
    def dmt_per_mt_lce_stoichiometric(self) -> float | None:
        """Stoichiometric DMT required per 1 mt LCE (100% recovery), or None if non-lithium."""
        if self.lce_factor_per_dmt is None:
            return None
        return (1.0 / self.lce_factor_per_dmt)

    @property
    def dmt_per_mt_lce_recoverable(self) -> float | None:
        """Effective DMT required per 1 mt LCE at default_recovery_rate, or None if non-lithium."""
        if self.lce_factor_per_dmt is None:
            return None
        return (1.0 / (self.lce_factor_per_dmt * self.default_recovery_rate))

    @property
    def li2o_pct(self) -> float | None:
        """Equivalent Li2O grade (%) on a dry basis for lithium-bearing materials."""
        from core.units.conversions import _LI2O_TO_LCE_FACTOR
        if self.lce_factor_per_dmt is None:
            return None
        return ((self.lce_factor_per_dmt / _LI2O_TO_LCE_FACTOR) * 100.0)

    def as_context(self, *, basis: str = "contained", context_name: str | None = None) -> pint.Context:
        """Returns a named context; use explicit methods for metal and LCE masses."""
        from core.units.conversions import create_battery_material_context
        return create_battery_material_context(self, basis=basis, name=context_name)

    # --- Physical Mass Conversions (WMT <-> DMT <-> Contained Metal <-> LCE) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def wmt_to_dmt(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) to Dry Metric Tonnes (DMT)."""
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(wet_mass)
        return wet_to_dry_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.dmt, allow_negative=True)

    def dmt_to_wmt(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to Wet Metric Tonnes (WMT)."""
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(dry_mass)
        return dry_to_wet_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.wmt, allow_negative=True)

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture_pct: float | Decimal,
        *,
        source_moisture_pct: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity:
        """Converts ore mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=self.moisture_pct if source_moisture_pct is None else source_moisture_pct,
            target_moisture_pct=target_moisture_pct,
            target_unit=target_unit or ureg.metric_ton,
        )

    def dmt_to_contained_metal(self, dry_mass: pint.Quantity) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) to contained active metal mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(dry_mass)
        frac = self.active_metal_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.contained_metal_tonne)

    def contained_metal_to_dmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained active metal mass to required Dry Metric Tonnes (DMT)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(metal_mass)
        frac = self.active_metal_pct / 100.0
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.dmt)

    def wmt_to_contained_metal(self, wet_mass: pint.Quantity) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to contained active metal mass (metric_ton)."""
        return self.dmt_to_contained_metal(self.wmt_to_dmt(wet_mass))

    def contained_metal_to_wmt(self, metal_mass: pint.Quantity) -> pint.Quantity:
        """Converts contained active metal mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.contained_metal_to_dmt(metal_mass))

    def dmt_to_lce(
        self,
        dry_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts Dry Metric Tonnes (DMT) of lithium feedstock/chemical to LCE (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.exceptions import IncompatibleCommodityError
        from core.units.registry import ureg
        if self.lce_factor_per_dmt is None:
            raise IncompatibleCommodityError(f"{self.name!r} is not a lithium-bearing material (lce_factor_per_dmt is None).")
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError(f"recovery_rate must be in (0, 1], got {rec}")
        mt_q = self._validate_mass(dry_mass)
        factor = self.lce_factor_per_dmt * rec
        scale = _scale_factor(mt_q.magnitude, factor)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def lce_to_dmt(
        self,
        lce_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts target LCE mass (metric_ton) to required Dry Metric Tonnes (DMT) of feedstock/chemical."""
        from core.units.conversions import _scale_factor
        from core.units.exceptions import IncompatibleCommodityError
        from core.units.registry import ureg
        if self.lce_factor_per_dmt is None:
            raise IncompatibleCommodityError(f"{self.name!r} is not a lithium-bearing material (lce_factor_per_dmt is None).")
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError(f"recovery_rate must be in (0, 1], got {rec}")
        mt_q = self._validate_mass(lce_mass)
        factor = self.lce_factor_per_dmt * rec
        scale = _scale_factor(mt_q.magnitude, factor)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def wmt_to_lce(
        self,
        wet_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts Wet Metric Tonnes (WMT) directly to LCE (metric_ton)."""
        return self.dmt_to_lce(self.wmt_to_dmt(wet_mass), recovery_rate=recovery_rate)

    def lce_to_wmt(
        self,
        lce_mass: pint.Quantity,
        *,
        recovery_rate: float | None = None,
    ) -> pint.Quantity:
        """Converts target LCE mass directly to required Wet Metric Tonnes (WMT)."""
        return self.dmt_to_wmt(self.lce_to_dmt(lce_mass, recovery_rate=recovery_rate))

    # --- Commercial Pricing Parity ($/dmt Product <-> $/mt Contained Metal or LCE) ---
    def dmt_price_to_metal_price(self, price_per_dmt: float) -> float:
        """Converts a $/dmt product quote (e.g. NiSO4, CoSO4) to $/mt contained active metal."""
        p = float(price_per_dmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt must be finite, got {price_per_dmt}")
        return (p / (self.active_metal_pct / 100.0))

    def metal_price_to_dmt_price(self, price_per_mt_metal: float) -> float:
        """Converts a $/mt contained metal quote to $/dmt product equivalent."""
        p = float(price_per_mt_metal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_metal must be finite, got {price_per_mt_metal}")
        return (p * (self.active_metal_pct / 100.0))

    def dmt_price_to_lce_price(
        self,
        price_per_dmt: float,
        *,
        recovery_rate: float | None = None,
        conversion_cost_per_mt_lce: float = 0.0,
    ) -> float:
        """Converts feedstock/chemical $/dmt (e.g. SC6) to implied $/mt LCE parity."""
        from core.units.exceptions import IncompatibleCommodityError
        if self.lce_factor_per_dmt is None:
            raise IncompatibleCommodityError(f"{self.name!r} is not a lithium-bearing material.")
        p = float(price_per_dmt)
        cost = float(conversion_cost_per_mt_lce)
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(p) or not math.isfinite(cost) or not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError("Inputs must be finite and recovery_rate must be in (0, 1].")
        return ((p / (self.lce_factor_per_dmt * rec)) + cost)

    def lce_price_to_dmt_price(
        self,
        price_per_mt_lce: float,
        *,
        recovery_rate: float | None = None,
        conversion_cost_per_mt_lce: float = 0.0,
    ) -> float:
        """Computes implied feedstock/chemical $/dmt netback from spot $/mt LCE."""
        from core.units.exceptions import IncompatibleCommodityError
        if self.lce_factor_per_dmt is None:
            raise IncompatibleCommodityError(f"{self.name!r} is not a lithium-bearing material.")
        p_lce = float(price_per_mt_lce)
        cost = float(conversion_cost_per_mt_lce)
        rec = self.default_recovery_rate if recovery_rate is None else float(recovery_rate)
        if not math.isfinite(p_lce) or not math.isfinite(cost) or not math.isfinite(rec) or not (0.0 < rec <= 1.0):
            raise ValueError("Inputs must be finite and recovery_rate must be in (0, 1].")
        return ((p_lce - cost) * (self.lce_factor_per_dmt * rec))

    def blend_with(
        self,
        other: BatteryMaterialAssay,
        *,
        self_wet_mass_fraction: float,
        blend_name: str | None = None,
    ) -> BatteryMaterialAssay:
        """Blends two battery material lots of the same active metal.

        - Moisture blends linearly by Wet Metric Tonnes (WMT).
        - Active metal %, LCE factor, purity, and recovery blend by Dry Metric Tonnes (DMT).
        """
        from core.units.exceptions import IncompatibleCommodityError
        if not isinstance(other, BatteryMaterialAssay):
            raise TypeError(f"Expected BatteryMaterialAssay, got {type(other).__name__}")
        if self.active_metal != other.active_metal:
            raise IncompatibleCommodityError(
                f"Cannot blend materials with different active metals: {self.active_metal!r} vs {other.active_metal!r}"
            )
        if not math.isfinite(self_wet_mass_fraction) or not (0.0 <= self_wet_mass_fraction <= 1.0):
            raise ValueError("self_wet_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_wet_mass_fraction)
        w2 = 1.0 - w1
        b_moisture = (self.moisture_pct * w1) + (other.moisture_pct * w2)

        d1 = w1 * self.dry_mass_fraction
        d2 = w2 * other.dry_mass_fraction
        total_d = d1 + d2
        dw1, dw2 = d1 / total_d, d2 / total_d

        b_metal = (self.active_metal_pct * dw1) + (other.active_metal_pct * dw2)
        b_rec = (self.default_recovery_rate * dw1) + (other.default_recovery_rate * dw2)

        b_purity: float | None = None
        if self.purity_pct is not None and other.purity_pct is not None:
            b_purity = ((self.purity_pct * dw1) + (other.purity_pct * dw2))

        b_lce: float | None = None
        if self.lce_factor_per_dmt is not None and other.lce_factor_per_dmt is not None:
            b_lce = ((self.lce_factor_per_dmt * dw1) + (other.lce_factor_per_dmt * dw2))

        return BatteryMaterialAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            active_metal=self.active_metal,
            active_metal_pct=(b_metal),
            moisture_pct=(b_moisture),
            purity_pct=b_purity,
            lce_factor_per_dmt=b_lce,
            default_recovery_rate=(b_rec),
        )

    # --- Global Battery Material Benchmark Presets & Factories ---
    @classmethod
    def from_spodumene_grade(
        cls,
        li2o_pct: float,
        *,
        moisture_pct: float = 8.0,
        recovery_rate: float = 0.85,
        name: str | None = None,
    ) -> BatteryMaterialAssay:
        """Constructs a custom Spodumene concentrate assay directly from % Li2O."""
        from core.units.conversions import _LI2O_TO_LCE_FACTOR
        from core.units.conversions import _LI_TO_LCE_FACTOR
        grade = float(li2o_pct)
        if not math.isfinite(grade) or not (0.5 <= grade <= 8.03):
            raise ValueError(f"Spodumene li2o_pct must be in [0.5, 8.03] (pure spodumene max), got {li2o_pct}")
        li_pct = grade * (_LI2O_TO_LCE_FACTOR / _LI_TO_LCE_FACTOR)
        return cls(
            name=name or f"Spodumene Concentrate ({grade:.2f}% Li2O)",
            active_metal="Li",
            active_metal_pct=(li_pct),
            moisture_pct=moisture_pct,
            lce_factor_per_dmt=((grade / 100.0) * _LI2O_TO_LCE_FACTOR),
            default_recovery_rate=recovery_rate,
        )
    @classmethod
    def spodumene_sc6(cls) -> BatteryMaterialAssay:
        """Spodumene Concentrate 6.0% Li2O (SC6 CIF China, ~85% typical refinery recovery)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.spodumene_sc6")

    @classmethod
    def lithium_carbonate_bg(cls) -> BatteryMaterialAssay:
        """Battery-Grade Lithium Carbonate (Li2CO3, 99.5% min purity)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.lithium_carbonate_bg")

    @classmethod
    def lithium_hydroxide_bg(cls) -> BatteryMaterialAssay:
        """Battery-Grade Lithium Hydroxide Monohydrate (LiOH·H2O, 56.5% min LiOH)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.lithium_hydroxide_bg")

    @classmethod
    def nickel_sulphate_hexahydrate(cls) -> BatteryMaterialAssay:
        """Battery-Grade Nickel Sulphate Hexahydrate (NiSO4·6H2O, 22.2% Ni)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.nickel_sulphate_hexahydrate")

    @classmethod
    def cobalt_sulphate_heptahydrate(cls) -> BatteryMaterialAssay:
        """Battery-Grade Cobalt Sulphate Heptahydrate (CoSO4·7H2O, 20.5% Co)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.cobalt_sulphate_heptahydrate")

    @classmethod
    def manganese_sulphate_monohydrate(cls) -> BatteryMaterialAssay:
        """Battery-Grade Manganese Sulphate Monohydrate (MnSO4·H2O, 32.0% Mn)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.manganese_sulphate_monohydrate")

    @classmethod
    def mhp_nickel(cls) -> BatteryMaterialAssay:
        """Mixed Hydroxide Precipitate (MHP — 38.0% Ni dry basis, ~30% moisture)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.mhp_nickel")

    @classmethod
    def cobalt_hydroxide(cls) -> BatteryMaterialAssay:
        """DRC Cobalt Hydroxide intermediate (30.0% Co dry basis, ~35% moisture)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "metals", "BatteryMaterialAssay.cobalt_hydroxide")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (GraphiteAssay,))
__all__ = ('BatteryMaterialAssay', 'GraphiteAssay', 'convert_graphite', 'create_graphite_context')
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
    from core.units.assays.metals.battery_materials import convert_graphite
    from core.units.assays.metals.battery_materials import create_graphite_context
