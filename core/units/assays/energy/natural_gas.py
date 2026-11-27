"""Commodex assays energy natural_gas implementation and public interfaces."""
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

#-------------------------------------------------------------------------------
# Gas Assay
#--------------------------------------------------------------------------------

@validated_assay
@dataclass(frozen=True, slots=True)
class GasAssay:
    """Immutable thermodynamic specification for natural gas and LNG streams.

    Resolves calorific, Wobbe Index, and cryogenic LNG cargo conversion factors
    across international reference standards (FERC, Texas RRC, ISO 13443, DIN 1343, Groningen).
    """

    name: str
    gcv_btu_per_scf: float
    ncv_btu_per_scf: float
    gcv_mj_per_sm3: float
    gcv_mj_per_nm3: float
    mmbtu_per_mt_lng: float | None = None
    relative_density: float | None = None
    wobbe_index_mj_per_sm3: float | None = None
    molar_composition: Mapping[str, float] | None = None

    def __init__(
        self,
        name: str = "custom_gas",
        *,
        gcv_btu_per_scf: float | None = None,
        ncv_btu_per_scf: float | None = None,
        gcv_mj_per_sm3: float | None = None,
        gcv_mj_per_nm3: float | None = None,
        ncv_to_gcv_ratio: float = 0.9009,
        relative_density: float | None = None,
        wobbe_index: float | None = None,
        molar_composition: Mapping[str, float] | None = None,
        normalize_composition: bool = False,
    ) -> None:
        """Constructs an immutable gas assay from any single standard calorific benchmark."""
        from core.units.composition import GAS_MOLAR_MASSES
        from core.units.conversions import DRY_AIR_DENSITY_15C_KG_M3
        from core.units.conversions import DRY_AIR_MOLAR_MASS_KG_KMOL
        from core.units.conversions import FT3_PER_M3
        from core.units.conversions import MJ_PER_BTU_IT
        from core.units.conversions import _DENSITY_RATIO_DIN0
        from core.units.conversions import _DENSITY_RATIO_ISO15
        from core.units.conversions import _ideal_gas_density_per_scf
        from core.units.conversions import _validated_composition
        provided = sum(p is not None for p in (gcv_btu_per_scf, gcv_mj_per_sm3, gcv_mj_per_nm3))
        if provided != 1:
            raise ValueError(
                "Provide exactly one calorific parameter: "
                "'gcv_btu_per_scf', 'gcv_mj_per_sm3', or 'gcv_mj_per_nm3'."
            )

        resolved_gcv_scf: float
        if gcv_btu_per_scf is not None:
            resolved_gcv_scf = float(gcv_btu_per_scf)
        elif gcv_mj_per_sm3 is not None:
            val_sm3 = float(gcv_mj_per_sm3)
            if not math.isfinite(val_sm3):
                raise ValueError(f"GCV must be finite, got {gcv_mj_per_sm3}")
            resolved_gcv_scf = val_sm3 / (MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_ISO15)
        else:
            val_nm3 = float(gcv_mj_per_nm3)  # type: ignore[arg-type]
            if not math.isfinite(val_nm3):
                raise ValueError(f"GCV must be finite, got {gcv_mj_per_nm3}")
            resolved_gcv_scf = val_nm3 / (MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_DIN0)

        if not math.isfinite(resolved_gcv_scf) or resolved_gcv_scf <= 0:
            raise ValueError(f"GCV must be finite and strictly positive, got {resolved_gcv_scf}")

        resolved_ncv_scf: float
        if ncv_btu_per_scf is not None:
            resolved_ncv_scf = float(ncv_btu_per_scf)
            if not math.isfinite(resolved_ncv_scf) or resolved_ncv_scf <= 0 or resolved_ncv_scf > resolved_gcv_scf:
                raise ValueError("NCV must be finite, strictly positive, and cannot exceed GCV.")
        else:
            if not math.isfinite(ncv_to_gcv_ratio) or not (0.5 < ncv_to_gcv_ratio <= 1.0):
                raise ValueError(f"ncv_to_gcv_ratio must be finite and in (0.5, 1.0], got {ncv_to_gcv_ratio}")
            resolved_ncv_scf = resolved_gcv_scf * ncv_to_gcv_ratio

        calc_mj_sm3 = resolved_gcv_scf * MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_ISO15
        calc_mj_nm3 = resolved_gcv_scf * MJ_PER_BTU_IT * FT3_PER_M3 * _DENSITY_RATIO_DIN0

        calc_wobbe: float | None = float(wobbe_index) if wobbe_index is not None else None
        calc_rel_dens: float | None = float(relative_density) if relative_density is not None else None
        if calc_rel_dens is not None:
            if not math.isfinite(calc_rel_dens) or calc_rel_dens <= 0:
                raise ValueError("Relative density must be finite and strictly positive.")
            if calc_wobbe is not None and (not math.isfinite(calc_wobbe) or calc_wobbe <= 0):
                raise ValueError("Wobbe index must be finite and strictly positive.")
            if calc_wobbe is None:
                calc_wobbe = calc_mj_sm3 / math.sqrt(calc_rel_dens)
        elif calc_wobbe is not None:
            if not math.isfinite(calc_wobbe) or calc_wobbe <= 0:
                raise ValueError("Wobbe index must be finite and strictly positive.")
            calc_rel_dens = (calc_mj_sm3 / calc_wobbe) ** 2

        clean_comp: Mapping[str, float] | None = None
        if molar_composition is not None:
            clean_comp = MappingProxyType(_validated_composition(molar_composition, normalize=normalize_composition))
            if calc_rel_dens is None:
                molar_mass = sum(float(GAS_MOLAR_MASSES[name]) * fraction for name, fraction in clean_comp.items())
                calc_rel_dens = molar_mass / DRY_AIR_MOLAR_MASS_KG_KMOL
                calc_wobbe = calc_mj_sm3 / math.sqrt(calc_rel_dens)

        # Derive LNG gravimetric heating value (MMBtu/mt) if relative density is known
        calc_mmbtu_mt: float | None = None
        if clean_comp is not None:
            molar_mass = sum(float(GAS_MOLAR_MASSES[name]) * fraction for name, fraction in clean_comp.items())
            density_per_scf = _ideal_gas_density_per_scf(molar_mass)
            calc_mmbtu_mt = resolved_gcv_scf / (density_per_scf * 1000)
        elif calc_rel_dens is not None:
            gas_density_15c = calc_rel_dens * DRY_AIR_DENSITY_15C_KG_M3
            mj_per_kg = calc_mj_sm3 / gas_density_15c
            calc_mmbtu_mt = (mj_per_kg * 1000.0) / (MJ_PER_BTU_IT * 1e6)

        object.__setattr__(self, "name", str(name))
        object.__setattr__(self, "gcv_btu_per_scf", (resolved_gcv_scf))
        object.__setattr__(self, "ncv_btu_per_scf", (resolved_ncv_scf))
        object.__setattr__(self, "gcv_mj_per_sm3", (calc_mj_sm3))
        object.__setattr__(self, "gcv_mj_per_nm3", (calc_mj_nm3))
        object.__setattr__(self, "mmbtu_per_mt_lng", (calc_mmbtu_mt) if calc_mmbtu_mt is not None else None)
        object.__setattr__(self, "relative_density", (calc_rel_dens) if calc_rel_dens is not None else None)
        object.__setattr__(
            self,
            "wobbe_index_mj_per_sm3",
            (calc_wobbe) if calc_wobbe is not None else None,
        )
        object.__setattr__(self, "molar_composition", clean_comp)

    @classmethod
    def from_composition(cls, molar_composition: Mapping[str, float | Decimal], *, name: str = "composition gas",
                         normalize_composition: bool = False) -> GasAssay:
        """Derive both calorific values and molar mass on the compatibility ideal-reference basis."""
        from core.units.composition import GasComposition
        result = GasComposition(molar_composition, normalize=normalize_composition).to_gas_assay(name=name)
        if cls is GasAssay:
            return result
        return cls(name=name, gcv_btu_per_scf=result.gcv_btu_per_scf, ncv_btu_per_scf=result.ncv_btu_per_scf,
                   relative_density=result.relative_density, molar_composition=result.molar_composition)

    def __hash__(self) -> int:
        comp_tuple = (
            tuple(sorted(self.molar_composition.items()))
            if self.molar_composition is not None
            else None
        )
        return hash(
            (
                self.name,
                self.gcv_btu_per_scf,
                self.ncv_btu_per_scf,
                self.gcv_mj_per_sm3,
                self.gcv_mj_per_nm3,
                self.mmbtu_per_mt_lng,
                self.relative_density,
                self.wobbe_index_mj_per_sm3,
                comp_tuple,
            )
        )

    @property
    def wobbe_index_btu_per_scf(self) -> float | None:
        """US FERC Wobbe Index in Btu/scf (GCV_Btu_per_scf / sqrt(SG))."""
        if self.relative_density is None:
            return None
        return (self.gcv_btu_per_scf / math.sqrt(self.relative_density))

    @property
    def ncv_mj_per_sm3(self) -> float:
        """Net Calorific Value (LHV) in MJ/Sm³ at 15°C, 101.325 kPa."""
        return (self.gcv_mj_per_sm3 * (self.ncv_btu_per_scf / self.gcv_btu_per_scf))

    @property
    def ncv_mj_per_nm3(self) -> float:
        """Net Calorific Value (LHV) in MJ/Nm³ at 0°C, 101.325 kPa."""
        return (self.gcv_mj_per_nm3 * (self.ncv_btu_per_scf / self.gcv_btu_per_scf))

    @property
    def molar_mass_g_per_mol(self) -> float | None:
        """Composition-derived conventional mass, or unknown when no composition exists."""
        from core.units.composition import GAS_MOLAR_MASSES
        if self.molar_composition is None:
            return None
        return sum(float(GAS_MOLAR_MASSES[name]) * fraction for name, fraction in self.molar_composition.items())

    @property
    def density_kg_per_scf(self) -> float | None:
        """Density matching this compatibility assay's ideal reference-volume basis."""
        from core.units.conversions import DRY_AIR_DENSITY_15C_KG_M3
        from core.units.conversions import FT3_PER_M3
        from core.units.conversions import _DENSITY_RATIO_ISO15
        from core.units.conversions import _ideal_gas_density_per_scf
        molar_mass = self.molar_mass_g_per_mol
        if molar_mass is not None:
            return _ideal_gas_density_per_scf(molar_mass)
        if self.relative_density is None:
            return None
        return self.relative_density * DRY_AIR_DENSITY_15C_KG_M3 / (FT3_PER_M3 * _DENSITY_RATIO_ISO15)

    @property
    def net_wobbe_index_mj_per_sm3(self) -> float | None:
        if self.relative_density is None:
            return None
        return self.ncv_mj_per_sm3 / math.sqrt(self.relative_density)

    @property
    def is_h_gas(self) -> bool | None:
        """Returns True for High-Calorific H-Gas (Wobbe >= 45.0 MJ/Sm³), False for L-Gas, or None."""
        if self.wobbe_index_mj_per_sm3 is None:
            return None
        return self.wobbe_index_mj_per_sm3 >= 45.0

    @property
    def gas_group(self) -> str | None:
        """European EN 16726 pipeline interchangeability group ('H-Gas', 'L-Gas', or None)."""
        if self.is_h_gas is None:
            return None
        return "H-Gas" if self.is_h_gas else "L-Gas"

    def _resolve_basis_ratio(self, basis: str) -> float:
        b = basis.strip().lower()
        if b not in {"gcv", "hhv", "ncv", "lhv"}:
            raise ValueError(f"Calorific basis must be 'gcv'/'hhv' or 'ncv'/'lhv', got {basis!r}")
        return (self.ncv_btu_per_scf / self.gcv_btu_per_scf) if b in {"ncv", "lhv"} else 1.0

    def as_context(self, *, basis: str = "gcv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging gas reference volumes and energy."""
        from core.units.conversions import create_gas_heating_context
        return create_gas_heating_context(self, basis=basis, name=context_name)

    def solve_z_factor(
        self,
        pressure: pint.Quantity,
        temperature: pint.Quantity,
        *,
        method: str = "auto",
        co2_mol_frac: float | None = None,
        h2s_mol_frac: float | None = None,
        tol: float = 1e-10,
        max_iter: int = 50,
        raise_on_nonconvergence: bool = True,
    ) -> SolverResult:
        """Solves real-gas compressibility factor (Z) for this assay and returns SolverResult."""
        from core.units.thermodynamics import SolverDomainError
        from core.units.thermodynamics import solve_z_factor
        if self.relative_density is None and method.strip().lower() in {"dak", "auto"}:
            raise SolverDomainError(
                "DAK requires measured or composition-derived relative density; no default gas is assumed."
            )
        rel_dens = self.relative_density if self.relative_density is not None else 1.0
        return solve_z_factor(
            pressure,
            temperature,
            relative_density=rel_dens,
            method=method,
            molar_composition=self.molar_composition,
            co2_mol_frac=co2_mol_frac,
            h2s_mol_frac=h2s_mol_frac,
            tol=tol,
            max_iter=max_iter,
            raise_on_nonconvergence=raise_on_nonconvergence,
        )

    def z_factor(
        self,
        pressure: pint.Quantity,
        temperature: pint.Quantity,
        *,
        method: str = "auto",
        co2_mol_frac: float | None = None,
        h2s_mol_frac: float | None = None,
        tol: float = 1e-10,
        max_iter: int = 50,
        return_result: bool = False,
    ) -> float | SolverResult:
        """Computes real-gas compressibility factor (Z) at (pressure, temperature) for this assay."""
        if not isinstance(return_result, bool):
            raise TypeError("return_result must be a bool.")
        solved = self.solve_z_factor(
            pressure,
            temperature,
            method=method,
            co2_mol_frac=co2_mol_frac,
            h2s_mol_frac=h2s_mol_frac,
            tol=tol,
            max_iter=max_iter,
            raise_on_nonconvergence=True,
        )
        return solved if return_result else solved.value

    def lng_liquid_to_energy(
        self,
        liquid_volume: pint.Quantity,
        *,
        liquid_density_kg_m3: float | Decimal | None = None,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Legacy ideal-reference parity using explicit cargo density; prefer LNGCargo."""
        from core.units.conversions import _as_decimal
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        if liquid_density_kg_m3 is None:
            raise ValueError("Supply cargo-specific liquid_density_kg_m3 or use LNGCargo; no universal LNG density is assumed.")
        liquid_density_kg_m3 = _as_decimal(liquid_density_kg_m3)
        if not liquid_density_kg_m3.is_finite() or liquid_density_kg_m3 <= 0:
            raise ValueError(f"liquid_density_kg_m3 must be finite and strictly positive, got {liquid_density_kg_m3}")
        ratio = self._resolve_basis_ratio(basis)
        cbm_q = liquid_volume.to(ureg.cbm_lng)
        if not _as_decimal(cbm_q.magnitude).is_finite() or cbm_q.magnitude < 0:
            raise ValueError("LNG liquid volume must be finite and non-negative.")
        mt_per_cbm = liquid_density_kg_m3 / 1000
        mmbtu_per_cbm = mt_per_cbm * _as_decimal(self.mmbtu_per_mt_lng) * _as_decimal(ratio)
        scale = _scale_factor(cbm_q.magnitude, mmbtu_per_cbm)
        return ureg.Quantity(cbm_q.magnitude * scale, ureg.MMBtu_it)

    def energy_to_lng_liquid(
        self,
        energy: pint.Quantity,
        *,
        liquid_density_kg_m3: float | Decimal | None = None,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Legacy inverse parity requiring explicit cargo density; prefer LNGCargo."""
        from core.units.conversions import _as_decimal
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        if liquid_density_kg_m3 is None:
            raise ValueError("Supply cargo-specific liquid_density_kg_m3 or use LNGCargo; no universal LNG density is assumed.")
        liquid_density_kg_m3 = _as_decimal(liquid_density_kg_m3)
        if not liquid_density_kg_m3.is_finite() or liquid_density_kg_m3 <= 0:
            raise ValueError(f"liquid_density_kg_m3 must be finite and strictly positive, got {liquid_density_kg_m3}")
        ratio = self._resolve_basis_ratio(basis)
        mmbtu_q = energy.to(ureg.MMBtu_it)
        if not _as_decimal(mmbtu_q.magnitude).is_finite() or mmbtu_q.magnitude < 0:
            raise ValueError("LNG cargo energy must be finite and non-negative.")
        mt_per_cbm = liquid_density_kg_m3 / 1000
        mmbtu_per_cbm = mt_per_cbm * _as_decimal(self.mmbtu_per_mt_lng) * _as_decimal(ratio)
        scale = _scale_factor(mmbtu_q.magnitude, mmbtu_per_cbm)
        return ureg.Quantity(mmbtu_q.magnitude / scale, ureg.cbm_lng)

    lng_volume_to_energy = lng_liquid_to_energy
    energy_to_lng_volume = energy_to_lng_liquid

    def lng_mass_to_energy(
        self,
        mass: pint.Quantity,
        *,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Converts LNG cargo mass (e.g. mt, tonne) directly to cargo energy (MMBtu_it)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        ratio = self._resolve_basis_ratio(basis)
        mt_q = mass.to(ureg.metric_tonne)
        mmbtu_per_mt = self.mmbtu_per_mt_lng * ratio
        scale = _scale_factor(mt_q.magnitude, mmbtu_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.MMBtu_it)

    def energy_to_lng_mass(
        self,
        energy: pint.Quantity,
        *,
        basis: str = "gcv",
    ) -> pint.Quantity:
        """Converts cargo energy (e.g. MMBtu_it, MWh, GJ) directly to LNG cargo mass (metric_tonne)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.mmbtu_per_mt_lng is None:
            raise ValueError("GasAssay requires relative_density to compute LNG mass-energy parity.")
        ratio = self._resolve_basis_ratio(basis)
        mmbtu_q = energy.to(ureg.MMBtu_it)
        mmbtu_per_mt = self.mmbtu_per_mt_lng * ratio
        scale = _scale_factor(mmbtu_q.magnitude, mmbtu_per_mt)
        return ureg.Quantity(mmbtu_q.magnitude / scale, ureg.metric_tonne)

    def lng_liquid_to_gas_volume(
        self,
        liquid_volume: pint.Quantity,
        *,
        target_standard: str = "iso15",
        liquid_density_kg_m3: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Legacy ideal-reference expansion with explicit cargo density; prefer LNGCargo."""
        from core.units.conversions import FT3_PER_M3
        from core.units.conversions import _DENSITY_RATIO_DIN0
        from core.units.conversions import _DENSITY_RATIO_FERC
        from core.units.conversions import _DENSITY_RATIO_GOST20
        from core.units.conversions import _DENSITY_RATIO_ISO15
        from core.units.conversions import _DENSITY_RATIO_TX
        from core.units.conversions import _as_decimal
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.relative_density is None:
            raise ValueError("GasAssay requires relative_density to compute regasified gas expansion.")
        if liquid_density_kg_m3 is None:
            raise ValueError("Supply cargo-specific liquid_density_kg_m3 or use LNGCargo; no universal LNG density is assumed.")
        liquid_density_kg_m3 = _as_decimal(liquid_density_kg_m3)
        if not liquid_density_kg_m3.is_finite() or liquid_density_kg_m3 <= 0:
            raise ValueError(f"liquid_density_kg_m3 must be finite and strictly positive, got {liquid_density_kg_m3}")

        std_key = target_standard.strip().lower()
        gost_unit = ureg.m3_20C_101325kPa_dry
        gas_density_15c = _as_decimal(self.density_kg_per_scf) * _as_decimal(FT3_PER_M3) * _as_decimal(_DENSITY_RATIO_ISO15)
        sm3_per_cbm = liquid_density_kg_m3 / gas_density_15c

        std_targets: dict[str, tuple[Unit, float]] = {
            "iso15": (ureg.Sm3_15C_101325kPa_dry, sm3_per_cbm),
            "din0": (ureg.Nm3_0C_101325kPa_dry, sm3_per_cbm * _as_decimal(_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_DIN0)),
            "gost20": (gost_unit, sm3_per_cbm * _as_decimal(_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_GOST20)),
            "ferc": (ureg.scf_60F_14_73psia, sm3_per_cbm * _as_decimal(FT3_PER_M3) * _as_decimal(_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_FERC)),
            "tx": (ureg.scf_60F_14_65psia, sm3_per_cbm * _as_decimal(FT3_PER_M3) * _as_decimal(_DENSITY_RATIO_ISO15 / _DENSITY_RATIO_TX)),
        }
        if std_key not in std_targets:
            raise ValueError(f"Unknown target_standard {target_standard!r}. Choose from {tuple(std_targets)}.")

        dst_unit, factor = std_targets[std_key]
        cbm_q = liquid_volume.to(ureg.cbm_lng)
        if not _as_decimal(cbm_q.magnitude).is_finite() or cbm_q.magnitude < 0:
            raise ValueError("LNG liquid volume must be finite and non-negative.")
        scale = _scale_factor(cbm_q.magnitude, factor)
        return ureg.Quantity(cbm_q.magnitude * scale, dst_unit)

    def blend_with(
        self,
        other: GasAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> GasAssay:
        """Computes ideal volumetric blend of two gas streams (GCV, NCV, SG, Wobbe, and molar composition)."""
        if not isinstance(other, GasAssay):
            raise TypeError(f"Expected GasAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        w_other = 1.0 - self_volume_fraction
        blended_gcv = (self.gcv_btu_per_scf * self_volume_fraction) + (other.gcv_btu_per_scf * w_other)
        blended_ncv = (self.ncv_btu_per_scf * self_volume_fraction) + (other.ncv_btu_per_scf * w_other)

        blended_sg: float | None = None
        if self.relative_density is not None and other.relative_density is not None:
            blended_sg = (self.relative_density * self_volume_fraction) + (other.relative_density * w_other)

        blended_comp: dict[str, float] | None = None
        if self.molar_composition is not None and other.molar_composition is not None:
            all_keys = set(self.molar_composition) | set(other.molar_composition)
            blended_comp = {k: (self.molar_composition.get(k, 0.0) * self_volume_fraction) + (other.molar_composition.get(k, 0.0) * w_other) for k in all_keys}

        return GasAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            gcv_btu_per_scf=blended_gcv,
            ncv_btu_per_scf=blended_ncv,
            relative_density=blended_sg,
            molar_composition=blended_comp,
        )

    @classmethod
    def henry_hub(cls) -> GasAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.henry_hub")

    @classmethod
    def uk_nbp(cls) -> GasAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.uk_nbp")

    @classmethod
    def dutch_ttf(cls) -> GasAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.dutch_ttf")

    @classmethod
    def groningen_lgas(cls) -> GasAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.groningen_lgas")

    @classmethod
    def lng_lean(cls) -> GasAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.lng_lean")

    @classmethod
    def lng_rich(cls) -> GasAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.lng_rich")

    @classmethod
    def jkm_lng(cls) -> GasAssay:
        """Platts JKM (Japan/Korea Marker) rich seaborne LNG benchmark (~1,100 Btu/scf)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.jkm_lng")

    @classmethod
    def gost_russian(cls) -> GasAssay:
        """Russian/CIS export pipeline gas benchmark under GOST 2939-63 (20°C, 101.325 kPa)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.gost_russian")

    @classmethod
    def pure_methane(cls) -> GasAssay:
        from core.units.standards import construct_reference
        return construct_reference(cls, "natural_gas", "GasAssay.pure_methane")


__all__ = ('GasAssay', 'GasReferenceConditions')
_PUBLIC_ALIASES = {'GasReferenceConditions': ('reference_conditions', 'GasReferenceConditions')}

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
    from core.units.reference_conditions import GasReferenceConditions
    from core.units.thermodynamics import SolverResult
    from core.units.registry import Unit
