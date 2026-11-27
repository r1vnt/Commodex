"""Commodex assays energy biofuels implementation and public interfaces."""
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
class CleanFuelAssay:
    """Immutable thermodynamic specification for zero/low-carbon energy carriers.

    Supports gravimetric (LHV/HHV), gaseous volumetric (Nm3), cryogenic/liquid
    shipping parity (m3), and hydrogen-carrier cracking equivalency for H2, NH3,
    and Methanol.
    """

    name: str
    lhv_mj_per_kg: float
    hhv_mj_per_kg: float
    density_kg_per_nm3: float | None = None
    liquid_density_kg_per_m3: float | None = None
    hydrogen_mass_fraction: float | None = None
    carbon_intensity_kg_co2e_per_kg: float | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.lhv_mj_per_kg) or self.lhv_mj_per_kg <= 0.0:
            raise ValueError(f"lhv_mj_per_kg must be finite and strictly positive, got {self.lhv_mj_per_kg}")
        if not math.isfinite(self.hhv_mj_per_kg) or self.hhv_mj_per_kg < self.lhv_mj_per_kg:
            raise ValueError(
                f"hhv_mj_per_kg ({self.hhv_mj_per_kg}) must be finite and >= lhv_mj_per_kg ({self.lhv_mj_per_kg})."
            )
        if self.density_kg_per_nm3 is not None:
            if not math.isfinite(self.density_kg_per_nm3) or self.density_kg_per_nm3 <= 0.0:
                raise ValueError(f"density_kg_per_nm3 must be strictly positive, got {self.density_kg_per_nm3}")
        if self.liquid_density_kg_per_m3 is not None:
            if not math.isfinite(self.liquid_density_kg_per_m3) or self.liquid_density_kg_per_m3 <= 0.0:
                raise ValueError(
                    f"liquid_density_kg_per_m3 must be strictly positive, got {self.liquid_density_kg_per_m3}"
                )
        if self.hydrogen_mass_fraction is not None:
            if not math.isfinite(self.hydrogen_mass_fraction) or not (0.0 < self.hydrogen_mass_fraction <= 1.0):
                raise ValueError(
                    f"hydrogen_mass_fraction must be in (0.0, 1.0], got {self.hydrogen_mass_fraction}"
                )
        if self.carbon_intensity_kg_co2e_per_kg is not None:
            if not math.isfinite(self.carbon_intensity_kg_co2e_per_kg) or self.carbon_intensity_kg_co2e_per_kg < 0.0:
                raise ValueError(
                    f"carbon_intensity_kg_co2e_per_kg cannot be negative, got {self.carbon_intensity_kg_co2e_per_kg}"
                )

    # --- Gravimetric & Volumetric Energy Density Properties ---
    @property
    def lhv_kwh_per_kg(self) -> float:
        return (self.lhv_mj_per_kg / 3.6)

    @property
    def hhv_kwh_per_kg(self) -> float:
        return (self.hhv_mj_per_kg / 3.6)

    @property
    def lhv_mwh_per_mt(self) -> float:
        return (self.lhv_mj_per_kg / 3.6)

    @property
    def hhv_mwh_per_mt(self) -> float:
        return (self.hhv_mj_per_kg / 3.6)

    @property
    def lhv_mmbtu_per_mt(self) -> float:
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.lhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mmbtu_per_mt(self) -> float:
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.hhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    @property
    def lhv_mj_per_nm3(self) -> float | None:
        if self.density_kg_per_nm3 is None:
            return None
        return (self.lhv_mj_per_kg * self.density_kg_per_nm3)

    @property
    def hhv_mj_per_nm3(self) -> float | None:
        if self.density_kg_per_nm3 is None:
            return None
        return (self.hhv_mj_per_kg * self.density_kg_per_nm3)

    @property
    def lhv_gj_per_m3_liquid(self) -> float | None:
        if self.liquid_density_kg_per_m3 is None:
            return None
        return ((self.lhv_mj_per_kg * self.liquid_density_kg_per_m3) / 1000.0)

    @property
    def carbon_intensity_g_co2e_per_mj_lhv(self) -> float | None:
        """Well-to-gate carbon intensity in gCO2e/MJ (LHV) for EU RFNBO / FuelEU Maritime."""
        if self.carbon_intensity_kg_co2e_per_kg is None:
            return None
        # (kg CO2e / kg fuel) * 1000 g/kg / (MJ / kg fuel) = gCO2e / MJ
        return ((self.carbon_intensity_kg_co2e_per_kg * 1000.0) / self.lhv_mj_per_kg)

    def _heating_value_mj_per_kg(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"lhv", "ncv"}:
            return self.lhv_mj_per_kg
        if key in {"hhv", "gcv"}:
            return self.hhv_mj_per_kg
        raise ValueError(f"basis must be 'lhv'/'ncv' or 'hhv'/'gcv', got {basis!r}")

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [mass] and [energy] for this carrier."""
        from core.units.conversions import create_clean_fuel_context
        return create_clean_fuel_context(self, basis=basis, name=context_name)

    # --- Mass <-> Energy & Liquid Shipping Conversions ---
    def mass_to_energy(self, mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts carrier mass (e.g. kg, mt) to total energy (MWh)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        hv_mwh_per_mt = self._heating_value_mj_per_kg(basis) / 3.6
        mt_q = mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, hv_mwh_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.MWh)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. MWh, GJ, MMBtu) to required carrier mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        hv_mwh_per_mt = self._heating_value_mj_per_kg(basis) / 3.6
        mwh_q = energy.to(ureg.MWh)
        scale = _scale_factor(mwh_q.magnitude, hv_mwh_per_mt)
        return ureg.Quantity(mwh_q.magnitude / scale, ureg.metric_ton)

    def gas_volume_to_mass(self, gas_volume: pint.Quantity) -> pint.Quantity:
        """Converts normal gaseous volume (DIN 1343 Nm3) to mass (kg)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(gas_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(gas_volume).__name__}")
        if self.density_kg_per_nm3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_nm3.")
        nm3_q = gas_volume.to(ureg.Nm3_0C_101325kPa_dry)
        scale = _scale_factor(nm3_q.magnitude, self.density_kg_per_nm3)
        return ureg.Quantity(nm3_q.magnitude * scale, ureg.kilogram)

    def mass_to_gas_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts carrier mass (e.g. kg, mt) to normal gaseous volume (DIN 1343 Nm3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.density_kg_per_nm3 is None:
            raise ValueError(f"{self.name!r} does not define density_kg_per_nm3.")
        kg_q = mass.to(ureg.kilogram)
        scale = _scale_factor(kg_q.magnitude, self.density_kg_per_nm3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.Nm3_0C_101325kPa_dry)

    def liquid_volume_to_mass(self, liquid_volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid/cryogenic cargo volume (e.g. m^3, L) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(liquid_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(liquid_volume).__name__}")
        if self.liquid_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define liquid_density_kg_per_m3.")
        m3_q = liquid_volume.to(ureg.meter ** 3)
        mt_per_m3 = self.liquid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_liquid_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts carrier mass (e.g. mt, kg) to liquid/cryogenic cargo volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.liquid_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define liquid_density_kg_per_m3.")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.liquid_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def liquid_volume_to_energy(self, liquid_volume: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts liquid/cryogenic cargo volume (e.g. m^3) directly to energy (MWh)."""
        return self.mass_to_energy(self.liquid_volume_to_mass(liquid_volume), basis=basis)

    def energy_to_liquid_volume(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. MWh, GJ, MMBtu) directly to liquid/cryogenic cargo volume (m^3)."""
        return self.mass_to_liquid_volume(self.energy_to_mass(energy, basis=basis))

    # --- Hydrogen-Carrier Equivalency (Ammonia / Methanol Cracking) ---
    def contained_hydrogen_mass(
        self,
        carrier_mass: pint.Quantity,
        *,
        cracking_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes recoverable H2 mass from a carrier cargo given cracking efficiency."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        from core.units.registry import validate_dimensionality
        if not isinstance(carrier_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(carrier_mass).__name__}")
        if not validate_dimensionality(carrier_mass, "kilogram"):
            raise pint.DimensionalityError(
                carrier_mass.units,
                ureg.kilogram,
                carrier_mass.dimensionality,
                ureg.kilogram.dimensionality,
            )
        if self.hydrogen_mass_fraction is None:
            raise ValueError(f"{self.name!r} does not define hydrogen_mass_fraction.")
        if not math.isfinite(cracking_efficiency) or not (0.0 < cracking_efficiency <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0.0, 1.0], got {cracking_efficiency}")
        factor = self.hydrogen_mass_fraction * cracking_efficiency
        scale = _scale_factor(carrier_mass.magnitude, factor)
        return ureg.Quantity(carrier_mass.magnitude * scale, carrier_mass.units)

    def carrier_mass_for_hydrogen(
        self,
        h2_mass: pint.Quantity,
        *,
        cracking_efficiency: float = 1.0,
    ) -> pint.Quantity:
        """Computes required carrier mass to deliver a target H2 mass after cracking."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        from core.units.registry import validate_dimensionality
        if not isinstance(h2_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(h2_mass).__name__}")
        if not validate_dimensionality(h2_mass, "kilogram"):
            raise pint.DimensionalityError(
                h2_mass.units,
                ureg.kilogram,
                h2_mass.dimensionality,
                ureg.kilogram.dimensionality,
            )
        if self.hydrogen_mass_fraction is None:
            raise ValueError(f"{self.name!r} does not define hydrogen_mass_fraction.")
        if not math.isfinite(cracking_efficiency) or not (0.0 < cracking_efficiency <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0.0, 1.0], got {cracking_efficiency}")
        factor = self.hydrogen_mass_fraction * cracking_efficiency
        scale = _scale_factor(h2_mass.magnitude, factor)
        return ureg.Quantity(h2_mass.magnitude / scale, h2_mass.units)

    # --- Clean Fuel Carrier Presets ---
    @classmethod
    def hydrogen(cls) -> CleanFuelAssay:
        """Pure Hydrogen (H2) ISO 6976 / NIST specification."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "CleanFuelAssay.hydrogen")

    @classmethod
    def ammonia(cls) -> CleanFuelAssay:
        """Anhydrous Ammonia (NH3) maritime & power carrier specification."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "CleanFuelAssay.ammonia")

    @classmethod
    def methanol(cls) -> CleanFuelAssay:
        """Pure Methanol (CH3OH) maritime & chemical carrier specification."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "CleanFuelAssay.methanol")


from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class BiofuelAssay:
    """Immutable physical and regulatory specification for liquid biofuels and SAF.

    Bridges volumetric (gal, m^3, bbl), gravimetric (mt, kg), and energy (MJ, GJ,
    MWh, MMBtu) conversions alongside US EPA RIN equivalence values and LCFS/RED II
    carbon-intensity (gCO2e/MJ) credit calculations.
    """

    name: str
    density_kg_per_m3: float
    lhv_mj_per_kg: float
    hhv_mj_per_kg: float
    epa_rin_ev: float | None = None
    rin_d_code: str | None = None
    carbon_intensity_g_co2e_per_mj: float | None = None
    thermal_expansion_coeff_per_c: float | None = None
    red_energy_multiplier: float = 1.0

    def __post_init__(self) -> None:
        if not math.isfinite(self.density_kg_per_m3) or not (500.0 <= self.density_kg_per_m3 <= 1200.0):
            raise ValueError(f"density_kg_per_m3 must be finite and in [500, 1200] kg/m^3, got {self.density_kg_per_m3}")
        if not math.isfinite(self.lhv_mj_per_kg) or self.lhv_mj_per_kg <= 0.0:
            raise ValueError(f"lhv_mj_per_kg must be finite and strictly positive, got {self.lhv_mj_per_kg}")
        if not math.isfinite(self.hhv_mj_per_kg) or self.hhv_mj_per_kg < self.lhv_mj_per_kg:
            raise ValueError(f"hhv_mj_per_kg ({self.hhv_mj_per_kg}) must be finite and >= lhv_mj_per_kg ({self.lhv_mj_per_kg}).")
        if self.epa_rin_ev is not None:
            if not math.isfinite(self.epa_rin_ev) or not (0.0 < self.epa_rin_ev <= 3.0):
                raise ValueError(f"epa_rin_ev must be in (0.0, 3.0], got {self.epa_rin_ev}")
        if self.carbon_intensity_g_co2e_per_mj is not None:
            if not math.isfinite(self.carbon_intensity_g_co2e_per_mj):
                raise ValueError(
                    f"carbon_intensity_g_co2e_per_mj must be finite, got {self.carbon_intensity_g_co2e_per_mj}"
                )
        if self.thermal_expansion_coeff_per_c is not None:
            if not math.isfinite(self.thermal_expansion_coeff_per_c) or self.thermal_expansion_coeff_per_c <= 0.0:
                raise ValueError(
                    f"thermal_expansion_coeff_per_c must be finite and strictly positive, got {self.thermal_expansion_coeff_per_c}"
                )
        if not math.isfinite(self.red_energy_multiplier) or self.red_energy_multiplier < 1.0:
            raise ValueError(f"red_energy_multiplier must be finite and >= 1.0, got {self.red_energy_multiplier}")

    # --- Volumetric & Gravimetric Parity Properties ---
    @property
    def specific_gravity(self) -> float:
        """Specific gravity relative to water at 15°C (999.102 kg/m^3)."""
        return (self.density_kg_per_m3 / 999.102)
    @property
    def bbl_per_mt(self) -> float:
        """Barrels per metric tonne at 15°C / 60°F."""
        from core.units.conversions import _LITRES_PER_BBL
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_BBL / 1000.0)))

    @property
    def gal_per_mt(self) -> float:
        """US gallons per metric tonne at 15°C / 60°F."""
        from core.units.conversions import _LITRES_PER_US_GAL
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def lhv_mj_per_l(self) -> float:
        """EU RED II volumetric Lower Heating Value (MJ/L)."""
        return ((self.lhv_mj_per_kg * self.density_kg_per_m3) / 1000.0)

    @property
    def hhv_mj_per_l(self) -> float:
        """Volumetric Higher Heating Value (MJ/L)."""
        return ((self.hhv_mj_per_kg * self.density_kg_per_m3) / 1000.0)

    @property
    def lhv_btu_per_gal(self) -> float:
        """US volumetric Lower Heating Value (Btu_IT/gal)."""
        from core.units.conversions import MJ_PER_BTU_IT
        from core.units.conversions import _LITRES_PER_US_GAL
        return ((self.lhv_mj_per_l * _LITRES_PER_US_GAL) / MJ_PER_BTU_IT)

    @property
    def lhv_mwh_per_mt(self) -> float:
        return (self.lhv_mj_per_kg / 3.6)

    @property
    def lhv_mmbtu_per_mt(self) -> float:
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.lhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    @property
    def hhv_mwh_per_mt(self) -> float:
        return (self.hhv_mj_per_kg / 3.6)

    @property
    def hhv_mmbtu_per_mt(self) -> float:
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.hhv_mj_per_kg / (MJ_PER_BTU_IT * 1000.0))

    def volume_correction_factor(
        self,
        observed_temp: float | pint.Quantity,
        *,
        ref_temp_c: float = 15.0,
    ) -> float:
        """Computes thermal Volume Correction Factor (VCF) from observed temperature to reference.

        Accepts observed_temp as float (°C) or pint.Quantity (e.g. Q_(68, 'degF'), Q_(20, 'degC')).
        Use ref_temp_c=15.0 for ISO/ARA, 15.5556 (60°F) for US, or 20.0 for Brazilian ANP ethanol.
        """
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
        vcf = math.exp(-alpha * dt * (1.0 + 0.8 * alpha * dt))
        return (vcf)

    def ghg_savings_pct(self, fossil_comparator_g_co2e_per_mj: float = 94.0) -> float:
        """Computes EU RED II/III (94.0) or ICAO CORSIA (89.0) GHG emissions savings (%)."""
        if self.carbon_intensity_g_co2e_per_mj is None:
            raise ValueError(f"{self.name!r} does not define carbon_intensity_g_co2e_per_mj.")
        if not math.isfinite(fossil_comparator_g_co2e_per_mj) or fossil_comparator_g_co2e_per_mj <= 0.0:
            raise ValueError("fossil_comparator_g_co2e_per_mj must be strictly positive.")
        return ((1.0 - (self.carbon_intensity_g_co2e_per_mj / fossil_comparator_g_co2e_per_mj)) * 100.0)

    def _heating_value_mj_per_kg(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"lhv", "ncv"}:
            return self.lhv_mj_per_kg
        if key in {"hhv", "gcv"}:
            return self.hhv_mj_per_kg
        raise ValueError(f"basis must be 'lhv'/'ncv' or 'hhv'/'gcv', got {basis!r}")

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [energy]."""
        from core.units.conversions import create_biofuel_context
        return create_biofuel_context(self, basis=basis, name=context_name)

    # --- Physical Volume <-> Mass <-> Energy Conversions ---
    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts liquid volume (e.g. gal, bbl, m^3) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts mass (e.g. mt, kg) to liquid volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        mt_q = mass.to(ureg.metric_ton)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def mass_to_energy(self, mass: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts mass (e.g. mt, kg) to total energy (GJ)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        # 1 MJ/kg = 1 GJ/mt
        hv_gj_per_mt = self._heating_value_mj_per_kg(basis)
        mt_q = mass.to(ureg.metric_ton)
        scale = _scale_factor(mt_q.magnitude, hv_gj_per_mt)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.GJ)

    def energy_to_mass(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. GJ, MWh, MMBtu) to mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        hv_gj_per_mt = self._heating_value_mj_per_kg(basis)
        gj_q = energy.to(ureg.GJ)
        scale = _scale_factor(gj_q.magnitude, hv_gj_per_mt)
        return ureg.Quantity(gj_q.magnitude / scale, ureg.metric_ton)

    def volume_to_energy(self, volume: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts liquid volume (e.g. gal, m^3, bbl) directly to energy (GJ)."""
        return self.mass_to_energy(self.volume_to_mass(volume), basis=basis)

    def energy_to_volume(self, energy: pint.Quantity, *, basis: str = "lhv") -> pint.Quantity:
        """Converts energy (e.g. GJ, MWh, MMBtu) directly to liquid volume (m^3)."""
        return self.mass_to_volume(self.energy_to_mass(energy, basis=basis))

    # --- Regulatory Credit Calculations (EPA RFS RINs & LCFS / CORSIA) ---
    def rins_generated(self, volume: pint.Quantity) -> float:
        """Computes US EPA RFS Renewable Identification Numbers (RINs) generated by a fuel volume."""
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.epa_rin_ev is None:
            raise ValueError(f"{self.name!r} does not define an epa_rin_ev equivalence value.")
        gal = float(volume.to(ureg.gallon).magnitude)
        return (gal * self.epa_rin_ev)

    def lcfs_credits_mt(
        self,
        volume: pint.Quantity,
        benchmark_ci_g_co2e_per_mj: float,
        *,
        eer: float = 1.0,
    ) -> float:
        """Computes LCFS / Clean Fuel Standard credits (or deficits if negative) in mt CO2e.

        Formula:
            Credits (mt CO2e) = (CI_benchmark * EER - CI_fuel) * Energy_LHV_MJ / 1e6
        """
        from core.units.registry import ureg
        if self.carbon_intensity_g_co2e_per_mj is None:
            raise ValueError(f"{self.name!r} does not define carbon_intensity_g_co2e_per_mj.")
        if not math.isfinite(benchmark_ci_g_co2e_per_mj) or benchmark_ci_g_co2e_per_mj < 0.0:
            raise ValueError(f"benchmark_ci_g_co2e_per_mj must be non-negative, got {benchmark_ci_g_co2e_per_mj}")
        if not math.isfinite(eer) or eer <= 0.0:
            raise ValueError(f"eer (Energy Economy Ratio) must be strictly positive, got {eer}")

        energy_mj = float(self.volume_to_energy(volume, basis="lhv").to(ureg.MJ).magnitude)
        ci_delta = (benchmark_ci_g_co2e_per_mj * eer) - self.carbon_intensity_g_co2e_per_mj
        return ((ci_delta * energy_mj) / 1e6)

    def red_compliance_energy(self, volume: pint.Quantity) -> pint.Quantity:
        """Computes EU RED II/III / UK RTFO compliance energy (GJ LHV) including double-counting."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        phys_gj = self.volume_to_energy(volume, basis="lhv")
        scale = _scale_factor(phys_gj.magnitude, self.red_energy_multiplier)
        return ureg.Quantity(phys_gj.magnitude * scale, phys_gj.units)

    def blend_with(
        self,
        other: BiofuelAssay,
        *,
        self_volume_fraction: float,
        blend_name: str | None = None,
    ) -> BiofuelAssay:
        """Computes a physically consistent liquid biofuel blend.

        - Density & RIN EV: volume-weighted
        - LHV & HHV (MJ/kg): mass-weighted
        - Carbon Intensity (gCO2e/MJ): LHV-energy-weighted
        """
        if not isinstance(other, BiofuelAssay):
            raise TypeError(f"Expected BiofuelAssay, got {type(other).__name__}")
        if not math.isfinite(self_volume_fraction) or not (0.0 <= self_volume_fraction <= 1.0):
            raise ValueError("self_volume_fraction must be between 0.0 and 1.0.")

        v1 = float(self_volume_fraction)
        v2 = 1.0 - v1

        blended_density = (self.density_kg_per_m3 * v1) + (other.density_kg_per_m3 * v2)
        m1 = self.density_kg_per_m3 * v1
        m2 = other.density_kg_per_m3 * v2
        total_m = m1 + m2
        w1, w2 = m1 / total_m, m2 / total_m

        blended_lhv = (self.lhv_mj_per_kg * w1) + (other.lhv_mj_per_kg * w2)
        blended_hhv = (self.hhv_mj_per_kg * w1) + (other.hhv_mj_per_kg * w2)

        blended_rin: float | None = None
        if self.epa_rin_ev is not None and other.epa_rin_ev is not None:
            blended_rin = ((self.epa_rin_ev * v1) + (other.epa_rin_ev * v2))

        e1 = m1 * self.lhv_mj_per_kg
        e2 = m2 * other.lhv_mj_per_kg

        blended_ci: float | None = None
        if self.carbon_intensity_g_co2e_per_mj is not None and other.carbon_intensity_g_co2e_per_mj is not None:
            blended_ci = ((self.carbon_intensity_g_co2e_per_mj * e1 + other.carbon_intensity_g_co2e_per_mj * e2) / (e1 + e2))

        blended_alpha: float | None = None
        if self.thermal_expansion_coeff_per_c is not None and other.thermal_expansion_coeff_per_c is not None:
            blended_alpha = ((self.thermal_expansion_coeff_per_c * v1) + (other.thermal_expansion_coeff_per_c * v2))

        blended_red_mult = ((self.red_energy_multiplier * e1 + other.red_energy_multiplier * e2) / (e1 + e2))

        return BiofuelAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            density_kg_per_m3=(blended_density),
            lhv_mj_per_kg=(blended_lhv),
            hhv_mj_per_kg=(blended_hhv),
            epa_rin_ev=blended_rin,
            rin_d_code=self.rin_d_code if self.rin_d_code == other.rin_d_code else None,
            carbon_intensity_g_co2e_per_mj=blended_ci,
            thermal_expansion_coeff_per_c=blended_alpha,
            red_energy_multiplier=blended_red_mult,
        )

    # --- Global Biofuel & SAF Benchmark Presets ---
    @classmethod
    def ethanol(cls) -> BiofuelAssay:
        """Anhydrous Fuel Ethanol (CBOT / Argus T2 ARA / US RFS D6)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiofuelAssay.ethanol")

    @classmethod
    def hydrous_ethanol(cls) -> BiofuelAssay:
        """Brazilian B3 / ANP Hydrous Ethanol (E100, ~95.5% v/v ethanol)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiofuelAssay.hydrous_ethanol")

    @classmethod
    def fame_biodiesel(cls) -> BiofuelAssay:
        """FAME Biodiesel / B100 (ARA FAME 0 / US RFS D4)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiofuelAssay.fame_biodiesel")

    @classmethod
    def ucome(cls) -> BiofuelAssay:
        """Used Cooking Oil Methyl Ester (ARA UCOME, EU RED II Annex IX-B double-counting)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiofuelAssay.ucome")

    @classmethod
    def hvo_renewable_diesel(cls) -> BiofuelAssay:
        """Hydrotreated Vegetable Oil / Renewable Diesel (HVO / RD100, US RFS D4)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiofuelAssay.hvo_renewable_diesel")

    @classmethod
    def saf_hefa(cls) -> BiofuelAssay:
        """HEFA-SPK Sustainable Aviation Fuel (ASTM D7566 / CORSIA / US RFS D4)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "BiofuelAssay.saf_hefa")


__all__ = ('BiofuelAssay', 'CleanFuelAssay')
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
