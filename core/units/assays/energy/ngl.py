"""Commodex assays energy ngl implementation and public interfaces."""
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
class NGLAssay:
    """Immutable physical and thermodynamic specification for NGLs and LPG (GPA 2145).

    Bridges Mont Belvieu liquid volume (gal, bbl), international cargo mass (mt),
    energy (MMBtu, GJ), and gas-plant shrinkage volume (scf) for frac-spread and
    LPG arbitrage workflows.
    """

    name: str
    specific_gravity_60f: float
    hhv_btu_per_gal: float
    lhv_btu_per_gal: float
    gas_equivalent_scf_per_gal: float
    molecular_weight_g_per_mol: float | None = None
    vapor_pressure_psia_100f: float | None = None
    thermal_expansion_coeff_per_c: float | None = None
    refrigerated_density_kg_per_m3: float | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.specific_gravity_60f) or not (0.30 <= self.specific_gravity_60f <= 0.78):
            raise ValueError(
                f"specific_gravity_60f for NGL/LPG must be in [0.30, 0.78], got {self.specific_gravity_60f}"
            )
        if not math.isfinite(self.hhv_btu_per_gal) or self.hhv_btu_per_gal <= 0.0:
            raise ValueError(f"hhv_btu_per_gal must be finite and strictly positive, got {self.hhv_btu_per_gal}")
        if not math.isfinite(self.lhv_btu_per_gal) or not (0.0 < self.lhv_btu_per_gal <= self.hhv_btu_per_gal):
            raise ValueError(
                f"lhv_btu_per_gal ({self.lhv_btu_per_gal}) must be in (0, hhv_btu_per_gal ({self.hhv_btu_per_gal})]."
            )
        if not math.isfinite(self.gas_equivalent_scf_per_gal) or self.gas_equivalent_scf_per_gal <= 0.0:
            raise ValueError(
                f"gas_equivalent_scf_per_gal must be strictly positive, got {self.gas_equivalent_scf_per_gal}"
            )
        if self.molecular_weight_g_per_mol is not None:
            if not math.isfinite(self.molecular_weight_g_per_mol) or self.molecular_weight_g_per_mol <= 0.0:
                raise ValueError(f"molecular_weight_g_per_mol must be strictly positive, got {self.molecular_weight_g_per_mol}")
        if self.vapor_pressure_psia_100f is not None:
            if not math.isfinite(self.vapor_pressure_psia_100f) or self.vapor_pressure_psia_100f < 0.0:
                raise ValueError(f"vapor_pressure_psia_100f cannot be negative, got {self.vapor_pressure_psia_100f}")
        if self.thermal_expansion_coeff_per_c is not None:
            if not math.isfinite(self.thermal_expansion_coeff_per_c) or self.thermal_expansion_coeff_per_c <= 0.0:
                raise ValueError(
                    f"thermal_expansion_coeff_per_c must be strictly positive, got {self.thermal_expansion_coeff_per_c}"
                )
        if self.refrigerated_density_kg_per_m3 is not None:
            if not math.isfinite(self.refrigerated_density_kg_per_m3) or self.refrigerated_density_kg_per_m3 <= 0.0:
                raise ValueError(
                    f"refrigerated_density_kg_per_m3 must be strictly positive, got {self.refrigerated_density_kg_per_m3}"
                )

    # --- Mass, Volume & Energy Parity Properties ---
    @property
    def density_kg_per_m3(self) -> float:
        """Liquid density at 60°F (kg/m^3)."""
        from core.units.conversions import _WATER_DENSITY_60F_KG_M3
        return (self.specific_gravity_60f * _WATER_DENSITY_60F_KG_M3)

    @property
    def gal_per_mt(self) -> float:
        """US gallons per metric tonne at 60°F (Mont Belvieu <-> FEI/ARA factor)."""
        from core.units.conversions import _LITRES_PER_US_GAL
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def bbl_per_mt(self) -> float:
        """Barrels per metric tonne at 60°F."""
        return (self.gal_per_mt / 42.0)

    @property
    def hhv_mmbtu_per_bbl(self) -> float:
        """Gross heating value per barrel (MMBtu_IT/bbl)."""
        return ((self.hhv_btu_per_gal * 42.0) / 1e6)

    @property
    def lhv_mmbtu_per_bbl(self) -> float:
        """Net heating value per barrel (MMBtu_IT/bbl)."""
        return ((self.lhv_btu_per_gal * 42.0) / 1e6)

    @property
    def hhv_mmbtu_per_mt(self) -> float:
        """Gross heating value per metric tonne (MMBtu_IT/mt)."""
        return ((self.hhv_btu_per_gal * self.gal_per_mt) / 1e6)

    @property
    def lhv_mmbtu_per_mt(self) -> float:
        """Net heating value per metric tonne (MMBtu_IT/mt)."""
        return ((self.lhv_btu_per_gal * self.gal_per_mt) / 1e6)

    @property
    def hhv_mj_per_kg(self) -> float:
        """Gross heating value (MJ/kg = GJ/mt)."""
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.hhv_mmbtu_per_mt * MJ_PER_BTU_IT * 1000.0)

    @property
    def lhv_mj_per_kg(self) -> float:
        """Net heating value (MJ/kg = GJ/mt)."""
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.lhv_mmbtu_per_mt * MJ_PER_BTU_IT * 1000.0)

    def _btu_per_gal(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"hhv", "gcv"}:
            return self.hhv_btu_per_gal
        if key in {"lhv", "ncv"}:
            return self.lhv_btu_per_gal
        raise ValueError(f"basis must be 'hhv'/'gcv' or 'lhv'/'ncv', got {basis!r}")

    def as_context(self, *, basis: str = "hhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [energy]."""
        from core.units.conversions import create_ngl_context
        return create_ngl_context(self, basis=basis, name=context_name)

    # --- Physical Volume <-> Mass <-> Energy <-> Shrinkage Conversions ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid NGL/LPG volume (e.g. gal, bbl, m^3) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        gal_q = volume.to(ureg.gallon)
        scale = _scale_factor(gal_q.magnitude, self.gal_per_mt)
        return ureg.Quantity(gal_q.magnitude / scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts NGL/LPG mass (e.g. mt, kg) to liquid volume (gallon at 60°F)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, self.gal_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.gallon)

    def volume_to_energy(self, volume: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts liquid NGL/LPG volume (e.g. gal, bbl) to energy (MMBtu_it)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        mmbtu_per_gal = self._btu_per_gal(basis) / 1e6
        gal_q = volume.to(ureg.gallon)
        scale = _scale_factor(gal_q.magnitude, mmbtu_per_gal)
        return ureg.Quantity(gal_q.magnitude * scale, ureg.MMBtu_it)

    def energy_to_volume(self, energy: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts energy (e.g. MMBtu_it, GJ) to liquid NGL/LPG volume (gallon)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        mmbtu_per_gal = self._btu_per_gal(basis) / 1e6
        mmbtu_q = energy.to(ureg.MMBtu_it)
        scale = _scale_factor(mmbtu_q.magnitude, mmbtu_per_gal)
        return ureg.Quantity(mmbtu_q.magnitude / scale, ureg.gallon)

    def liquid_to_gas_shrinkage(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts extracted liquid NGL volume (gal/bbl) to gas-plant shrinkage volume (scf)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        gal_q = volume.to(ureg.gallon)
        scale = _scale_factor(gal_q.magnitude, self.gas_equivalent_scf_per_gal)
        return ureg.Quantity(gal_q.magnitude * scale, ureg.scf_60F_14696psi_dry)

    def mass_to_energy(self, mass: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts NGL/LPG mass (e.g. mt, kg) directly to energy (MMBtu_it)."""
        return self.volume_to_energy(self.mass_to_volume(mass), basis=basis)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "hhv") -> pint.Quantity:
        """Converts energy (e.g. MMBtu_it, GJ) directly to NGL/LPG mass (metric_ton)."""
        return self.volume_to_mass(self.energy_to_volume(energy, basis=basis))

    def gas_shrinkage_to_liquid(self, gas_volume: pint.Quantity) -> pint.Quantity:
        """Converts gas-plant shrinkage volume (e.g. scf, Mcf) to liquid NGL volume (gallon)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(gas_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(gas_volume).__name__}")
        scf_q = gas_volume.to(ureg.scf_60F_14696psi_dry)
        scale = _scale_factor(scf_q.magnitude, self.gas_equivalent_scf_per_gal)
        return ureg.Quantity(scf_q.magnitude / scale, ureg.gallon)

    def refrigerated_volume_to_mass(self, liquid_volume: pint.Quantity) -> pint.Quantity:
        """Converts refrigerated VLGC cargo volume (m^3 at boiling point) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(liquid_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(liquid_volume).__name__}")
        if self.refrigerated_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define refrigerated_density_kg_per_m3.")
        m3_q = liquid_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.refrigerated_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_refrigerated_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts LPG mass (e.g. mt) to refrigerated VLGC cargo volume (m^3 at boiling point)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.refrigerated_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define refrigerated_density_kg_per_m3.")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.refrigerated_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float = 15.5556,
    ) -> float:
        """Computes light-end thermal Volume Correction Factor (VCF) to 60°F (or 15°C)."""
        from core.units.registry import ureg
        if self.thermal_expansion_coeff_per_c is None:
            raise ValueError(f"{self.name!r} does not define thermal_expansion_coeff_per_c.")
        if isinstance(observed_temp, pint.Quantity):
            t_obs_c = float(observed_temp.to(ureg.degC).magnitude)
        else:
            t_obs_c = float(observed_temp)
        if not math.isfinite(t_obs_c) or not math.isfinite(ref_temp_c):
            raise ValueError("Temperatures must be finite.")
        alpha = self.thermal_expansion_coeff_per_c
        dt = t_obs_c - ref_temp_c
        return (math.exp(-alpha * dt * (1.0 + 0.8 * alpha * dt)))

    def blend_with(
        self,
        other: NGLAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> NGLAssay:
        """Blends two NGL/LPG streams volumetrically (e.g. Propane/Butane autogas or Y-grade).

        - SG, HHV/LHV (Btu/gal), gas equivalent, thermal expansion: volume-weighted
        - Vapour pressure: mole-weighted (Raoult's Law) when molecular weights are known
        """
        if not isinstance(other, NGLAssay):
            raise TypeError(f"Expected NGLAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        v1 = float(self_volume_fraction)
        v2 = 1.0 - v1

        b_sg = (self.specific_gravity_60f * v1) + (other.specific_gravity_60f * v2)
        b_hhv = (self.hhv_btu_per_gal * v1) + (other.hhv_btu_per_gal * v2)
        b_lhv = (self.lhv_btu_per_gal * v1) + (other.lhv_btu_per_gal * v2)
        b_scf = (self.gas_equivalent_scf_per_gal * v1) + (other.gas_equivalent_scf_per_gal * v2)

        b_mw: float | None = None
        b_vp: float | None = None
        if self.molecular_weight_g_per_mol is not None and other.molecular_weight_g_per_mol is not None:
            n1 = (v1 * self.specific_gravity_60f) / self.molecular_weight_g_per_mol
            n2 = (v2 * other.specific_gravity_60f) / other.molecular_weight_g_per_mol
            x1, x2 = n1 / (n1 + n2), n2 / (n1 + n2)
            b_mw = ((x1 * self.molecular_weight_g_per_mol) + (x2 * other.molecular_weight_g_per_mol))
            if self.vapor_pressure_psia_100f is not None and other.vapor_pressure_psia_100f is not None:
                b_vp = ((x1 * self.vapor_pressure_psia_100f) + (x2 * other.vapor_pressure_psia_100f))
        elif self.vapor_pressure_psia_100f is not None and other.vapor_pressure_psia_100f is not None:
            b_vp = ((v1 * self.vapor_pressure_psia_100f) + (v2 * other.vapor_pressure_psia_100f))

        b_alpha: float | None = None
        if self.thermal_expansion_coeff_per_c is not None and other.thermal_expansion_coeff_per_c is not None:
            b_alpha = ((v1 * self.thermal_expansion_coeff_per_c) + (v2 * other.thermal_expansion_coeff_per_c))

        b_refrig_rho: float | None = None
        if self.refrigerated_density_kg_per_m3 is not None and other.refrigerated_density_kg_per_m3 is not None:
            b_refrig_rho = ((v1 * self.refrigerated_density_kg_per_m3) + (v2 * other.refrigerated_density_kg_per_m3))

        return NGLAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            specific_gravity_60f=(b_sg),
            hhv_btu_per_gal=(b_hhv),
            lhv_btu_per_gal=(b_lhv),
            gas_equivalent_scf_per_gal=(b_scf),
            molecular_weight_g_per_mol=b_mw,
            vapor_pressure_psia_100f=b_vp,
            thermal_expansion_coeff_per_c=b_alpha,
            refrigerated_density_kg_per_m3=b_refrig_rho,
        )

    # --- GPA 2145-16 Pure Component & Commercial LPG Presets ---
    @classmethod
    def ethane(cls) -> NGLAssay:
        """Pure Ethane (C2H6) — Mont Belvieu purity ethane."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "NGLAssay.ethane")

    @classmethod
    def propane(cls) -> NGLAssay:
        """Commercial Propane (C3H8) — Mont Belvieu / Argus FEI / Saudi CP."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "NGLAssay.propane")

    @classmethod
    def isobutane(cls) -> NGLAssay:
        """Isobutane (i-C4H10) — Alkylation feedstock."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "NGLAssay.isobutane")

    @classmethod
    def n_butane(cls) -> NGLAssay:
        """Normal Butane (n-C4H10) — Gasoline blending & LPG benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "NGLAssay.n_butane")

    @classmethod
    def natural_gasoline(cls) -> NGLAssay:
        """Natural Gasoline (C5+ / Pentanes Plus) — Diluent & naphtha blendstock."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "NGLAssay.natural_gasoline")

    @classmethod
    def y_grade(cls) -> NGLAssay:
        """US Y-Grade Raw NGL Mix (~42% C2, 28% C3, 7% iC4, 13% nC4, 10% C5+ by volume)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "NGLAssay.y_grade")


__all__ = ('NGLAssay',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
