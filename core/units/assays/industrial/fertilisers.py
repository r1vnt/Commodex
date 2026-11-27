"""Commodex assays industrial fertilisers implementation and public interfaces."""
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

from core.units.assays.industrial.chemicals import ChemicalSolutionAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class AmmoniaAssay(ChemicalSolutionAssay):
    solute_formula: str = "NH3"
    state: str

    def __post_init__(self):
        from core.units.assays.industrial.chemicals import ChemicalSolutionAssay
        ChemicalSolutionAssay.__post_init__(self)
        if self.solute_formula != "NH3" or self.state not in {"refrigerated", "pressurised", "aqueous"}:
            raise ValueError("Declare NH3 and refrigerated, pressurised or aqueous ammonia.")

    def ammonia_mass(self, quantity):
        return self.active_mass(quantity)


from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class FertilizerAssay:
    """Immutable nutrient-grade and physical specification for fertilisers.

    Bridges physical product mass (mt, short_ton, lb), nutrient-equivalent mass
    (N, P2O5, P, K2O, K, S), liquid/bulk volume (gal, m^3), NOLA $/st <-> $/mt
    pricing parity, and NPK blending.
    """

    name: str
    n_pct: float = 0.0
    p2o5_pct: float = 0.0
    k2o_pct: float = 0.0
    s_pct: float = 0.0
    density_kg_per_m3: float | None = None
    is_liquid: bool = False

    def __post_init__(self) -> None:
        for field_name, val in (
            ("n_pct", self.n_pct),
            ("p2o5_pct", self.p2o5_pct),
            ("k2o_pct", self.k2o_pct),
            ("s_pct", self.s_pct),
        ):
            if not math.isfinite(val) or not (0.0 <= val <= 100.0):
                raise ValueError(f"{field_name} must be finite and in [0, 100], got {val}")
        if (self.n_pct + self.p2o5_pct + self.k2o_pct + self.s_pct) <= 0.0:
            raise ValueError("At least one nutrient percentage (N, P2O5, K2O, S) must be strictly positive.")
        if (self.n_pct + self.p2o5_pct + self.k2o_pct + self.s_pct) > 100.0 + 1e-9:
            raise ValueError("Sum of nutrient percentages cannot exceed 100%.")
        if self.density_kg_per_m3 is not None:
            if not math.isfinite(self.density_kg_per_m3) or not (300.0 <= self.density_kg_per_m3 <= 2500.0):
                raise ValueError(
                    f"density_kg_per_m3 must be finite and in [300, 2500] kg/m^3, got {self.density_kg_per_m3}"
                )

    # --- Nutrient Grade & Volumetric Parity Properties ---
    @property
    def npk_grade(self) -> str:
        """Standard N-P2O5-K2O (and optional -S) grade string (e.g. '46-0-0', '18-46-0')."""
        base = f"{self.n_pct:g}-{self.p2o5_pct:g}-{self.k2o_pct:g}"
        return f"{base}-{self.s_pct:g}S" if self.s_pct > 0.0 else base

    @property
    def elemental_p_pct(self) -> float:
        """Elemental phosphorus content (wt% P) converted from P2O5."""
        from core.units.conversions import _P_IN_P2O5
        return (self.p2o5_pct * _P_IN_P2O5)

    @property
    def elemental_k_pct(self) -> float:
        """Elemental potassium content (wt% K) converted from K2O."""
        from core.units.conversions import _K_IN_K2O
        return (self.k2o_pct * _K_IN_K2O)

    @property
    def gal_per_mt(self) -> float | None:
        """US gallons per metric tonne (primarily for liquid UAN), or None if density unset."""
        from core.units.conversions import _LITRES_PER_US_GAL
        if self.density_kg_per_m3 is None:
            return None
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def lb_per_gal(self) -> float | None:
        """US pounds per US gallon (US liquid fertiliser standard), or None if density unset."""
        from core.units.conversions import _LB_PER_MT
        if self.gal_per_mt is None:
            return None
        return (_LB_PER_MT / self.gal_per_mt)

    @property
    def lb_n_per_gal(self) -> float | None:
        """Pounds of elemental Nitrogen per US gallon (UAN application benchmark)."""
        if self.lb_per_gal is None:
            return None
        return (self.lb_per_gal * (self.n_pct / 100.0))

    @property
    def stowage_factor_m3_per_mt(self) -> float | None:
        """Bulk stowage factor (m^3 per mt), or None if density unset."""
        if self.density_kg_per_m3 is None:
            return None
        return (1000.0 / self.density_kg_per_m3)

    def _nutrient_fraction(self, nutrient: str) -> float:
        from core.units.conversions import _K_IN_K2O
        from core.units.conversions import _P_IN_P2O5
        key = nutrient.strip().upper()
        mapping = {
            "N": self.n_pct / 100.0,
            "P2O5": self.p2o5_pct / 100.0,
            "P": (self.p2o5_pct * _P_IN_P2O5) / 100.0,
            "K2O": self.k2o_pct / 100.0,
            "K": (self.k2o_pct * _K_IN_K2O) / 100.0,
            "S": self.s_pct / 100.0,
        }
        if key not in mapping:
            raise ValueError(f"nutrient must be one of {sorted(mapping)}, got {nutrient!r}")
        frac = mapping[key]
        if frac <= 0.0:
            raise ValueError(f"{self.name!r} contains 0% {key}.")
        return frac

    def as_context(self, *, nutrient: str = "N", context_name: str | None = None) -> pint.Context:
        """Bridges volume and product mass; use explicit methods for nutrient masses."""
        from core.units.conversions import create_fertilizer_context
        return create_fertilizer_context(self, nutrient=nutrient, name=context_name)

    # --- Physical Mass & Volume Conversions ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def product_to_nutrient(self, product_mass: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts physical fertiliser mass (e.g. mt, short_ton, lb) to nutrient mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(product_mass)
        frac = self._nutrient_fraction(nutrient)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def nutrient_to_product(self, nutrient_mass: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts target nutrient mass (e.g. mt N, mt P2O5, mt K2O) to required product mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(nutrient_mass)
        frac = self._nutrient_fraction(nutrient)
        scale = _scale_factor(mt_q.magnitude, frac)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid or bulk stowage volume (e.g. gal, m^3) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts mass (e.g. mt, short_ton) to liquid or bulk stowage volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def volume_to_nutrient(self, volume: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts liquid/bulk volume (e.g. gal, m^3) directly to nutrient mass (metric_ton)."""
        return self.product_to_nutrient(self.volume_to_mass(volume), nutrient)

    def nutrient_to_volume(self, nutrient_mass: pint.Quantity, nutrient: str) -> pint.Quantity:
        """Converts target nutrient mass (e.g. mt N, lb N) directly to liquid/bulk volume (m^3)."""
        return self.mass_to_volume(self.nutrient_to_product(nutrient_mass, nutrient))

    def equivalent_product_mass(
        self,
        product_mass: pint.Quantity,
        other: FertilizerAssay,
        nutrient: str,
    ) -> pint.Quantity:
        """Computes the mass of `other` required to deliver the same `nutrient` tonnage as `product_mass`."""
        if not isinstance(other, FertilizerAssay):
            raise TypeError(f"Expected FertilizerAssay, got {type(other).__name__}")
        target_nutrient = self.product_to_nutrient(product_mass, nutrient)
        return other.nutrient_to_product(target_nutrient, nutrient)

    # --- Commercial Pricing Parity ($/st <-> $/mt, $/product <-> $/nutrient, $/gal <-> $/mt) ---
    @staticmethod
    def st_price_to_mt_price(price_per_st: float) -> float:
        """Converts US NOLA barge $/short_ton quote to international $/metric_ton."""
        from core.units.conversions import _MT_PER_SHORT_TON
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / _MT_PER_SHORT_TON)

    @staticmethod
    def mt_price_to_st_price(price_per_mt: float) -> float:
        """Converts international $/metric_ton quote to US NOLA barge $/short_ton."""
        from core.units.conversions import _MT_PER_SHORT_TON
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p * _MT_PER_SHORT_TON)

    def product_price_to_nutrient_price(self, price_per_mt_product: float, nutrient: str) -> float:
        """Converts a $/mt product quote to implied $/mt contained nutrient (e.g. $/mt N)."""
        p = float(price_per_mt_product)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_product must be finite, got {price_per_mt_product}")
        return (p / self._nutrient_fraction(nutrient))

    def nutrient_price_to_product_price(self, price_per_mt_nutrient: float, nutrient: str) -> float:
        """Converts a $/mt nutrient quote to implied $/mt product equivalent."""
        p = float(price_per_mt_nutrient)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt_nutrient must be finite, got {price_per_mt_nutrient}")
        return (p * self._nutrient_fraction(nutrient))

    def gal_price_to_mt_price(self, price_per_gal: float) -> float:
        """Converts liquid fertiliser $/gal quote to $/mt equivalent."""
        if self.gal_per_mt is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        p = float(price_per_gal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gal must be finite, got {price_per_gal}")
        return (p * self.gal_per_mt)

    def mt_price_to_gal_price(self, price_per_mt: float) -> float:
        """Converts liquid fertiliser $/mt quote to $/gal equivalent."""
        if self.gal_per_mt is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_m3.")
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.gal_per_mt)

    def st_price_to_lb_nutrient_price(self, price_per_st: float, nutrient: str) -> float:
        """Converts a US $/short_ton product quote to $/lb contained nutrient (e.g. $/lb N)."""
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / (2000.0 * self._nutrient_fraction(nutrient)))

    def st_price_to_unit_price(self, price_per_st: float, nutrient: str) -> float:
        """Converts a US $/short_ton product quote to $/short-ton unit (1 unit = 20 lb nutrient)."""
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / (self._nutrient_fraction(nutrient) * 100.0))

    def blend_with(
        self,
        other: FertilizerAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> FertilizerAssay:
        """Computes a mass-weighted NPK/S fertiliser blend (e.g. bulk-blended NPK)."""
        if not isinstance(other, FertilizerAssay):
            raise TypeError(f"Expected FertilizerAssay, got {type(other).__name__}")
        if self.is_liquid != other.is_liquid:
            raise ValueError("Cannot directly blend liquid and solid fertilisers.")
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1

        b_rho: float | None = None
        if self.density_kg_per_m3 is not None and other.density_kg_per_m3 is not None:
            vol_per_kg = (w1 / self.density_kg_per_m3) + (w2 / other.density_kg_per_m3)
            b_rho = (1.0 / vol_per_kg)

        return FertilizerAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            n_pct=((self.n_pct * w1) + (other.n_pct * w2)),
            p2o5_pct=((self.p2o5_pct * w1) + (other.p2o5_pct * w2)),
            k2o_pct=((self.k2o_pct * w1) + (other.k2o_pct * w2)),
            s_pct=((self.s_pct * w1) + (other.s_pct * w2)),
            density_kg_per_m3=b_rho,
            is_liquid=self.is_liquid,
        )

    # --- Global Fertiliser Benchmark Presets ---
    @classmethod
    def urea_granular(cls) -> FertilizerAssay:
        """Granular Urea (46-0-0 — Middle East FOB / NOLA Barge / Brazil CFR)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.urea_granular")

    @classmethod
    def uan_32(cls) -> FertilizerAssay:
        """Urea Ammonium Nitrate 32% N liquid solution (32-0-0, 11.06 lb/gal)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.uan_32")

    @classmethod
    def uan_28(cls) -> FertilizerAssay:
        """Urea Ammonium Nitrate 28% N liquid solution (28-0-0, 10.67 lb/gal)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.uan_28")

    @classmethod
    def dap(cls) -> FertilizerAssay:
        """Diammonium Phosphate (DAP, 18-46-0 — Tampa / Morocco / India CFR)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.dap")

    @classmethod
    def map_fertilizer(cls) -> FertilizerAssay:
        """Monoammonium Phosphate (MAP, 11-52-0 — Brazil CFR / NOLA)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.map_fertilizer")

    @classmethod
    def potash_mop(cls) -> FertilizerAssay:
        """Muriate of Potash (MOP / KCl, 0-0-60 — Vancouver / Brazil CFR / SE Asia)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.potash_mop")

    @classmethod
    def ammonium_sulphate(cls) -> FertilizerAssay:
        """Granular Ammonium Sulphate (21-0-0-24S)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.ammonium_sulphate")

    @classmethod
    def anhydrous_ammonia(cls) -> FertilizerAssay:
        """Anhydrous Ammonia agricultural grade (82-0-0, refrigerated liquid at -33°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.anhydrous_ammonia")

    @classmethod
    def tsp(cls) -> FertilizerAssay:
        """Triple Superphosphate (TSP, 0-46-0 — Morocco / Brazil CFR)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.tsp")

    @classmethod
    def can_27(cls) -> FertilizerAssay:
        """Calcium Ammonium Nitrate (CAN 27% N — European benchmark)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.can_27")

    @classmethod
    def potash_sop(cls) -> FertilizerAssay:
        """Sulphate of Potash (SOP / K2SO4, 0-0-50-18S — specialty chloride-free potash)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "FertilizerAssay.potash_sop")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (AmmoniaAssay,))
__all__ = ('AmmoniaAssay', 'FertilizerAssay', 'convert_ammonia', 'create_ammonia_context')
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
    from core.units.assays.industrial.fertilisers import convert_ammonia
    from core.units.assays.industrial.fertilisers import create_ammonia_context
