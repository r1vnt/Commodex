"""Commodex assays industrial petrochemicals implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final
import inspect
import math
import pint
import warnings

from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class PetrochemAssay:
    """Immutable physical and stoichiometric specification for petrochemicals.

    Bridges liquid/cryogenic volume (gal, bbl, m^3), mass (mt, lb, kg), molar
    quantities (kmol), ASTM D1555 thermal volume correction, and stoichiometric
    feedstock-to-product conversion yields.
    """

    name: str
    formula: str
    molecular_weight_g_per_mol: float
    density_kg_per_m3: float
    reference_temp_c: float = 15.0
    purity_pct: float = 99.5
    thermal_expansion_coeff_per_c: float | None = None
    freezing_point_c: float | None = None

    def __post_init__(self) -> None:
        if not self.formula or not self.formula.strip():
            raise ValueError("formula must be a non-empty string (e.g. 'C2H4', 'C6H6').")
        if not math.isfinite(self.molecular_weight_g_per_mol) or self.molecular_weight_g_per_mol <= 0.0:
            raise ValueError(
                f"molecular_weight_g_per_mol must be strictly positive, got {self.molecular_weight_g_per_mol}"
            )
        if not math.isfinite(self.density_kg_per_m3) or not (300.0 <= self.density_kg_per_m3 <= 2000.0):
            raise ValueError(
                f"density_kg_per_m3 must be finite and in [300, 2000] kg/m^3, got {self.density_kg_per_m3}"
            )
        if not math.isfinite(self.reference_temp_c) or not (-200.0 <= self.reference_temp_c <= 100.0):
            raise ValueError(f"reference_temp_c must be in [-200, 100] °C, got {self.reference_temp_c}")
        if not math.isfinite(self.purity_pct) or not (0.0 < self.purity_pct <= 100.0):
            raise ValueError(f"purity_pct must be in (0, 100], got {self.purity_pct}")
        if self.thermal_expansion_coeff_per_c is not None:
            if not math.isfinite(self.thermal_expansion_coeff_per_c) or self.thermal_expansion_coeff_per_c <= 0.0:
                raise ValueError(
                    f"thermal_expansion_coeff_per_c must be strictly positive, got {self.thermal_expansion_coeff_per_c}"
                )
        if self.freezing_point_c is not None and not math.isfinite(self.freezing_point_c):
            raise ValueError(f"freezing_point_c must be finite, got {self.freezing_point_c}")

    # --- Volumetric & Gravimetric Parity Properties ---
    @property
    def specific_gravity(self) -> float:
        """Specific gravity relative to water at 15°C (999.102 kg/m^3)."""
        return (self.density_kg_per_m3 / 999.102)

    @property
    def gal_per_mt(self) -> float:
        """US gallons per metric tonne at reference_temp_c."""
        from core.units.conversions import _LITRES_PER_US_GAL
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def bbl_per_mt(self) -> float:
        """Barrels per metric tonne at reference_temp_c."""
        return (self.gal_per_mt / 42.0)

    @property
    def lb_per_gal(self) -> float:
        """US pounds per US gallon at reference_temp_c (USGC aromatics standard)."""
        from core.units.conversions import _LB_PER_MT
        return (_LB_PER_MT / self.gal_per_mt)

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [substance]."""
        from core.units.conversions import create_petrochem_context
        return create_petrochem_context(self, name=context_name)

    # --- Physical Volume <-> Mass <-> Moles Conversions ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid/cryogenic volume (e.g. gal, bbl, m^3) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts mass (e.g. mt, lb, kg) to liquid/cryogenic volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def mass_to_moles(self, mass: pint.Quantity, *, active_only: bool = False) -> pint.Quantity:
        """Converts mass (e.g. mt, kg) to molar amount (kilomole).

        If active_only=True, scales mass by purity_pct before dividing by molecular weight.
        """
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        kg_q = mass.to(ureg.kilogram)
        purity_mult = (self.purity_pct / 100.0) if active_only else 1.0
        # 1 kg / (MW g/mol) = 1 / MW kmol
        kmol_per_kg = purity_mult / self.molecular_weight_g_per_mol
        scale = _scale_factor(kg_q.magnitude, kmol_per_kg)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilomole)

    def moles_to_mass(self, moles: pint.Quantity, *, active_only: bool = False) -> pint.Quantity:
        """Converts molar amount (e.g. kmol, mol) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(moles, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(moles).__name__}")
        kmol_q = moles.to(ureg.kilomole)
        purity_mult = (self.purity_pct / 100.0) if active_only else 1.0
        # 1 kmol * MW (kg/kmol) / 1000 = MW / 1000 mt
        mt_per_kmol = (self.molecular_weight_g_per_mol / 1000.0) / purity_mult
        scale = _scale_factor(kmol_q.magnitude, mt_per_kmol)
        return ureg.Quantity(kmol_q.magnitude * scale, ureg.metric_ton)

    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float | None = None,
    ) -> float:
        """Computes ASTM D1555 thermal Volume Correction Factor (VCF) to reference temperature."""
        from core.units.registry import ureg
        if self.thermal_expansion_coeff_per_c is None:
            raise ValueError(f"{self.name!r} does not define thermal_expansion_coeff_per_c.")
        t_ref = self.reference_temp_c if ref_temp_c is None else float(ref_temp_c)
        if isinstance(observed_temp, pint.Quantity):
            t_obs_c = float(observed_temp.to(ureg.degC).magnitude)
        else:
            t_obs_c = float(observed_temp)
        if not math.isfinite(t_obs_c) or not math.isfinite(t_ref):
            raise ValueError("Temperatures must be finite.")
        alpha = self.thermal_expansion_coeff_per_c
        dt = t_obs_c - t_ref
        return (math.exp(-alpha * dt * (1.0 + 0.8 * alpha * dt)))

    def density_at_temp(self, observed_temp: float | pint.Quantity) -> float:
        """Computes hot-loaded liquid density (kg/m^3) at observed_temp via ASTM D4311."""
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        return (self.density_kg_per_m3 * vcf)

    def observed_volume_to_mass(
        self,
        observed_volume: pint.Quantity,
        observed_temp: float | pint.Quantity,
    ) -> pint.Quantity:
        """Converts observed liquid volume at observed_temp directly to mass (metric_ton) via ASTM D1555 VCF."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(observed_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(observed_volume).__name__}")
        vcf = self.volume_correction_factor(observed_temp, ref_temp_c=self.reference_temp_c)
        scale = _scale_factor(observed_volume.magnitude, vcf)
        std_vol = ureg.Quantity(observed_volume.magnitude * scale, observed_volume.units)
        return self.volume_to_mass(std_vol)

    # --- Stoichiometric Feedstock-to-Product Yield Parity ---
    def stoichiometric_product_mass(
        self,
        feedstock_mass: pint.Quantity,
        product_assay: PetrochemAssay,
        *,
        molar_ratio: float = 1.0,
        yield_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes downstream product mass (metric_ton) from this feedstock mass.

        Args:
            feedstock_mass: Input mass of this feedstock (e.g. Q_(1000, 'metric_ton')).
            product_assay: Target downstream PetrochemAssay.
            molar_ratio: Moles of product formed per 1 mole of this feedstock.
            yield_efficiency: Plant stoichiometric yield efficiency in (0, 1].
        """
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(product_assay, PetrochemAssay):
            raise TypeError(f"Expected PetrochemAssay, got {type(product_assay).__name__}")
        if not math.isfinite(molar_ratio) or molar_ratio <= 0.0:
            raise ValueError(f"molar_ratio must be strictly positive, got {molar_ratio}")
        if not math.isfinite(yield_efficiency) or not (0.0 < yield_efficiency <= 1.0):
            raise ValueError(f"yield_efficiency must be in (0, 1], got {yield_efficiency}")

        feed_kmol = self.mass_to_moles(feedstock_mass, active_only=True)
        scale = _scale_factor(feed_kmol.magnitude, molar_ratio * yield_efficiency)
        prod_kmol = ureg.Quantity(feed_kmol.magnitude * scale, ureg.kilomole)
        return product_assay.moles_to_mass(prod_kmol, active_only=True)

    def stoichiometric_feedstock_mass(
        self,
        product_mass: pint.Quantity,
        feedstock_assay: PetrochemAssay,
        *,
        molar_ratio: float = 1.0,
        yield_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes required feedstock mass (metric_ton) to produce a target mass of this product.

        Args:
            product_mass: Target mass of this product (e.g. Q_(1000, 'metric_ton')).
            feedstock_assay: Upstream feedstock PetrochemAssay.
            molar_ratio: Moles of feedstock consumed per 1 mole of this product.
            yield_efficiency: Plant stoichiometric yield efficiency in (0, 1].
        """
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(feedstock_assay, PetrochemAssay):
            raise TypeError(f"Expected PetrochemAssay, got {type(feedstock_assay).__name__}")
        if not math.isfinite(molar_ratio) or molar_ratio <= 0.0:
            raise ValueError(f"molar_ratio must be strictly positive, got {molar_ratio}")
        if not math.isfinite(yield_efficiency) or not (0.0 < yield_efficiency <= 1.0):
            raise ValueError(f"yield_efficiency must be in (0, 1], got {yield_efficiency}")

        prod_kmol = self.mass_to_moles(product_mass, active_only=True)
        scale = _scale_factor(prod_kmol.magnitude, molar_ratio / yield_efficiency)
        feed_kmol = ureg.Quantity(prod_kmol.magnitude * scale, ureg.kilomole)
        return feedstock_assay.moles_to_mass(feed_kmol, active_only=True)

    # --- Commercial Regional Arbitrage Pricing ($/gal & $/lb <-> $/mt) ---
    def gal_price_to_mt_price(self, price_per_gal: float) -> float:
        """Converts USGC $/gal quote to Europe/Asia $/mt equivalent."""
        p = float(price_per_gal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gal must be finite, got {price_per_gal}")
        return (p * self.gal_per_mt)

    def mt_price_to_gal_price(self, price_per_mt: float) -> float:
        """Converts Europe/Asia $/mt quote to USGC $/gal equivalent."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.gal_per_mt)

    @staticmethod
    def lb_price_to_mt_price(price_per_lb: float) -> float:
        """Converts USGC $/lb quote (e.g. Ethylene/Propylene) to $/mt."""
        from core.units.conversions import _LB_PER_MT
        p = float(price_per_lb)
        if not math.isfinite(p):
            raise ValueError(f"price_per_lb must be finite, got {price_per_lb}")
        return (p * _LB_PER_MT)

    @staticmethod
    def mt_price_to_lb_price(price_per_mt: float) -> float:
        """Converts $/mt quote to USGC $/lb equivalent."""
        from core.units.conversions import _LB_PER_MT
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / _LB_PER_MT)

    def blend_with(
        self,
        other: PetrochemAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> PetrochemAssay:
        """Blends two miscible liquid petrochemical streams (e.g. BTX / Mixed Xylenes).

        - Density & thermal expansion: volume-weighted
        - Molecular weight: mole-weighted
        - Purity: mass-weighted
        """
        if not isinstance(other, PetrochemAssay):
            raise TypeError(f"Expected PetrochemAssay, got {type(other).__name__}")
        if abs(self.reference_temp_c - other.reference_temp_c) > 1e-6:
            raise ValueError(
                f"Cannot blend streams with different reference temperatures: "
                f"{self.reference_temp_c}°C vs {other.reference_temp_c}°C."
            )
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        v1 = float(self_volume_fraction)
        v2 = 1.0 - v1

        b_rho = (self.density_kg_per_m3 * v1) + (other.density_kg_per_m3 * v2)
        m1 = self.density_kg_per_m3 * v1
        m2 = other.density_kg_per_m3 * v2
        total_m = m1 + m2
        w1, w2 = m1 / total_m, m2 / total_m

        n1 = m1 / self.molecular_weight_g_per_mol
        n2 = m2 / other.molecular_weight_g_per_mol
        x1, x2 = n1 / (n1 + n2), n2 / (n1 + n2)
        b_mw = (self.molecular_weight_g_per_mol * x1) + (other.molecular_weight_g_per_mol * x2)
        b_purity = (self.purity_pct * w1) + (other.purity_pct * w2)

        b_alpha: float | None = None
        if self.thermal_expansion_coeff_per_c is not None and other.thermal_expansion_coeff_per_c is not None:
            b_alpha = ((self.thermal_expansion_coeff_per_c * v1) + (other.thermal_expansion_coeff_per_c * v2))

        return PetrochemAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            formula=self.formula if self.formula == other.formula else "Blend",
            molecular_weight_g_per_mol=(b_mw),
            density_kg_per_m3=(b_rho),
            reference_temp_c=self.reference_temp_c,
            purity_pct=(b_purity),
            thermal_expansion_coeff_per_c=b_alpha,
            freezing_point_c=None,
        )

    # --- Global Petrochemical Benchmark Presets ---
    @classmethod
    def ethylene(cls) -> PetrochemAssay:
        """Polymer-Grade Ethylene (C2H4 — cryogenic liquid at -103.7°C, 1 atm)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.ethylene")

    @classmethod
    def propylene_pgp(cls) -> PetrochemAssay:
        """Polymer-Grade Propylene (C3H6 — pressurised liquid at 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.propylene_pgp")

    @classmethod
    def benzene(cls) -> PetrochemAssay:
        """Nitration-Grade Benzene (C6H6 — ASTM D2359, 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.benzene")

    @classmethod
    def toluene(cls) -> PetrochemAssay:
        """Nitration-Grade Toluene (C7H8 — ASTM D841, 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.toluene")

    @classmethod
    def paraxylene(cls) -> PetrochemAssay:
        """Paraxylene (PX, C8H10 — polyester chain benchmark, 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.paraxylene")

    @classmethod
    def styrene_monomer(cls) -> PetrochemAssay:
        """Styrene Monomer (SM, C8H8 — 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.styrene_monomer")

    @classmethod
    def meg(cls) -> PetrochemAssay:
        """Fibre-Grade Monoethylene Glycol (MEG, C2H6O2 — 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.meg")

    @classmethod
    def mixed_xylenes(cls) -> PetrochemAssay:
        """Isomer-Grade Mixed Xylenes (MX, C8H10 — USGC / Asian aromatics benchmark, 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.mixed_xylenes")

    @classmethod
    def butadiene(cls) -> PetrochemAssay:
        """1,3-Butadiene (BD, C4H6 — pressurised liquid at 15°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "PetrochemAssay.butadiene")


__all__ = ('PetrochemAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
