"""Commodex assays energy hydrogen implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, field
from decimal import Decimal
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from types import MappingProxyType
from typing import TYPE_CHECKING, Any, Final
import inspect
import math
import pint
import warnings

from core.units.conversions import H2_HHV_MJ_PER_KG
from core.units.conversions import H2_LHV_MJ_PER_KG
from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class HydrogenAssay:
    """Immutable physical, energetic, and regulatory specification for hydrogen and carriers.

    Bridges storage/transport volume (m^3, Nm^3, L), carrier and contained H2 mass
    (kg, metric_ton), energy (MWh, GJ, MMBtu), electrolyser power demand, ammonia/LOHC
    cracking parity, and US 45V / EU RFNBO carbon-intensity screening.

    Density describes storage volume at the stated temperature and pressure,
    not an equation of state. Purity is metadata; h2_mass_fraction controls yield.
    Context masses always represent total product/carrier mass.
    """

    name: str
    state: str
    density_kg_per_m3: float
    h2_mass_fraction: float = 1.0
    reference_temp_c: float = 0.0
    reference_pressure_bar: float = 1.01325
    lhv_mj_per_kg_h2: float = H2_LHV_MJ_PER_KG
    hhv_mj_per_kg_h2: float = H2_HHV_MJ_PER_KG
    default_cracking_efficiency: float = 1.0
    purity_pct: float = 99.97
    carbon_intensity_kg_co2e_per_kg_h2: float | None = None

    def __post_init__(self) -> None:
        if self.state not in {"gas", "liquid", "carrier"}:
            raise ValueError(f"state must be 'gas', 'liquid', or 'carrier', got {self.state!r}")
        if not math.isfinite(self.density_kg_per_m3) or not (0.01 <= self.density_kg_per_m3 <= 2000.0):
            raise ValueError(
                f"density_kg_per_m3 must be finite and in [0.01, 2000] kg/m^3, got {self.density_kg_per_m3}"
            )
        if not math.isfinite(self.h2_mass_fraction) or not (0.0 < self.h2_mass_fraction <= 1.0):
            raise ValueError(f"h2_mass_fraction must be in (0, 1], got {self.h2_mass_fraction}")
        if not math.isfinite(self.reference_temp_c) or not (-273.15 < self.reference_temp_c <= 200.0):
            raise ValueError(f"reference_temp_c must be in (-273.15, 200] °C, got {self.reference_temp_c}")
        if not math.isfinite(self.reference_pressure_bar) or self.reference_pressure_bar <= 0.0:
            raise ValueError(f"reference_pressure_bar must be strictly positive, got {self.reference_pressure_bar}")
        if not math.isfinite(self.lhv_mj_per_kg_h2) or self.lhv_mj_per_kg_h2 <= 0.0:
            raise ValueError(f"lhv_mj_per_kg_h2 must be strictly positive, got {self.lhv_mj_per_kg_h2}")
        if not math.isfinite(self.hhv_mj_per_kg_h2) or self.hhv_mj_per_kg_h2 < self.lhv_mj_per_kg_h2:
            raise ValueError("hhv_mj_per_kg_h2 must be finite and >= lhv_mj_per_kg_h2.")
        if not math.isfinite(self.default_cracking_efficiency) or not (0.0 < self.default_cracking_efficiency <= 1.0):
            raise ValueError(
                f"default_cracking_efficiency must be in (0, 1], got {self.default_cracking_efficiency}"
            )
        if not math.isfinite(self.purity_pct) or not (0.0 < self.purity_pct <= 100.0):
            raise ValueError(f"purity_pct must be in (0, 100], got {self.purity_pct}")
        if self.carbon_intensity_kg_co2e_per_kg_h2 is not None:
            if not math.isfinite(self.carbon_intensity_kg_co2e_per_kg_h2) or self.carbon_intensity_kg_co2e_per_kg_h2 < 0.0:
                raise ValueError(
                    f"carbon_intensity_kg_co2e_per_kg_h2 must be non-negative, got {self.carbon_intensity_kg_co2e_per_kg_h2}"
                )

    # --- Energy Density, Volumetric Storage & Regulatory Properties ---
    @property
    def lhv_kwh_per_kg_h2(self) -> float:
        """Lower Heating Value in kWh per kg of contained H2 (~33.32 kWh/kg)."""
        return (self.lhv_mj_per_kg_h2 / 3.6)

    @property
    def hhv_kwh_per_kg_h2(self) -> float:
        """Higher Heating Value in kWh per kg of contained H2 (~39.41 kWh/kg)."""
        return (self.hhv_mj_per_kg_h2 / 3.6)

    @property
    def lhv_mmbtu_per_kg_h2(self) -> float:
        """Lower Heating Value in MMBtu per kg of contained H2."""
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.lhv_mj_per_kg_h2 / (MJ_PER_BTU_IT * 1_000_000.0))

    @property
    def hhv_mmbtu_per_kg_h2(self) -> float:
        """Higher Heating Value in MMBtu per kg of contained H2."""
        from core.units.conversions import MJ_PER_BTU_IT
        return (self.hhv_mj_per_kg_h2 / (MJ_PER_BTU_IT * 1_000_000.0))

    @property
    def h2_kg_per_m3(self) -> float:
        """Volumetric H2 storage density (kg contained H2 per m^3 of gas/liquid/carrier)."""
        return (self.density_kg_per_m3 * self.h2_mass_fraction)

    @property
    def is_eu_rfnbo_compliant(self) -> bool | None:
        """CI-only RFNBO screen (28.2 g CO2e/MJ LHV); None if CI is unset.

        Assumes supplied CI includes the required lifecycle/delivery boundary.
        This does not verify renewable origin, additionality, or electricity
        temporal/geographical correlation and cannot certify full compliance.
        """
        from core.units.conversions import _EU_RFNBO_MAX_CI_G_PER_MJ
        if self.carbon_intensity_kg_co2e_per_kg_h2 is None:
            return None
        return self.carbon_intensity_kg_co2e_per_kg_h2 <= (_EU_RFNBO_MAX_CI_G_PER_MJ * self.lhv_mj_per_kg_h2 / 1000.0)

    @property
    def us_45v_credit_usd_per_kg_h2(self) -> float | None:
        """Unadjusted statutory 45V tier reference, assuming the 5x PWA rate.

        Uses a $0.60 base and 33.4% for the 0.45-to-1.5 tier ($1.002/kg).
        Not a current-year credit or eligibility determination. Use
        us_45v_credit_for_year() for supported inflation-adjusted amounts.
        https://uscode.house.gov/us/usc/t26/s45V/
        """
        ci = self.carbon_intensity_kg_co2e_per_kg_h2
        if ci is None:
            return None
        if ci < 0.45:
            return 3.00
        if ci < 1.50:
            return 1.002
        if ci < 2.50:
            return 0.75
        if ci <= 4.00:
            return 0.60
        return 0.0

    def us_45v_credit_for_year(
        self,
        production_year: int,
        *,
        prevailing_wage_and_apprenticeship: bool = False,
    ) -> float | None:
        """Returns the 2026 CI-tier rate, conditional on all other eligibility rules.

        PWA must be explicitly selected for the 5x increased credit. Other
        years require an updated rate table rather than silently reusing 2026.
        Source: IRS Notice 2026-41, section 3.02.
        https://www.irs.gov/irb/2026-29_irb
        """
        if production_year != 2026:
            raise ValueError("Only production_year=2026 has a verified rate table.")
        ci = self.carbon_intensity_kg_co2e_per_kg_h2
        if ci is None:
            return None
        if ci < 0.45:
            rate = 0.656
        elif ci < 1.5:
            rate = 0.219
        elif ci < 2.5:
            rate = 0.164
        elif ci <= 4.0:
            rate = 0.131
        else:
            rate = 0.0
        return float(Decimal(str(rate)) * (5 if prevailing_wage_and_apprenticeship else 1))

    def _mj_per_kg_h2(self, basis: str) -> float:
        key = basis.strip().lower()
        if key in {"lhv", "ncv"}:
            return self.lhv_mj_per_kg_h2
        if key in {"hhv", "gcv"}:
            return self.hhv_mj_per_kg_h2
        raise ValueError(f"basis must be 'lhv'/'ncv' or 'hhv'/'gcv', got {basis!r}")

    def as_context(self, *, basis: str = "lhv", context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [energy]."""
        from core.units.conversions import create_hydrogen_context
        return create_hydrogen_context(self, basis=basis, name=context_name)

    # --- Physical Conversions (Volume <-> Carrier Mass <-> H2 Mass <-> Nm^3 <-> Energy) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.kilogram)

    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts storage/transport volume (m^3, L) to total product/carrier mass (kg)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        scale = _scale_factor(m3_q.magnitude, self.density_kg_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.kilogram)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts total product/carrier mass (kg, metric_ton) to storage/transport volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        kg_q = self._validate_mass(mass)
        scale = _scale_factor(kg_q.magnitude, self.density_kg_per_m3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.meter ** 3)

    def carrier_to_h2_mass(
        self,
        carrier_mass: pint.Quantity,
        *,
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts carrier mass (e.g. mt NH3, mt LOHC) to delivered pure H2 mass (kg)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0, 1], got {eff}")
        kg_q = self._validate_mass(carrier_mass)
        factor = self.h2_mass_fraction * eff
        scale = _scale_factor(kg_q.magnitude, factor)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilogram)

    def h2_to_carrier_mass(
        self,
        h2_mass: pint.Quantity,
        *,
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts target delivered H2 mass (kg, mt) to required carrier mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError(f"cracking_efficiency must be in (0, 1], got {eff}")
        kg_q = self._validate_mass(h2_mass)
        factor = self.h2_mass_fraction * eff
        scale = _scale_factor(kg_q.magnitude, factor)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.kilogram).to(ureg.metric_ton)

    @staticmethod
    def h2_mass_to_nm3(h2_mass: pint.Quantity) -> pint.Quantity:
        """Converts pure H2 mass (kg, mt) to Normal cubic metres (Nm^3 at 0°C, 1 atm)."""
        from core.units.conversions import H2_DENSITY_0C_KG_NM3
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(h2_mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(h2_mass).__name__}")
        kg_q = h2_mass.to(ureg.kilogram)
        scale = _scale_factor(kg_q.magnitude, H2_DENSITY_0C_KG_NM3)
        return ureg.Quantity(kg_q.magnitude / scale, ureg.Nm3_H2)

    @staticmethod
    def nm3_to_h2_mass(nm3_volume: pint.Quantity) -> pint.Quantity:
        """Converts Normal cubic metres (Nm^3 at 0°C, 1 atm) to pure H2 mass (kg)."""
        from core.units.conversions import H2_DENSITY_0C_KG_NM3
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(nm3_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(nm3_volume).__name__}")
        if nm3_volume.dimensionality == ureg.Nm3_H2.dimensionality:
            m3_q = nm3_volume.to(ureg.Nm3_H2)
        else:
            m3_q = nm3_volume.to(ureg.meter ** 3)
        scale = _scale_factor(m3_q.magnitude, H2_DENSITY_0C_KG_NM3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.kilogram)

    def mass_to_energy(
        self,
        mass: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts product/carrier mass to delivered H2 chemical energy (MWh)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        h2_kg = self.carrier_to_h2_mass(mass, cracking_efficiency=cracking_efficiency)
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        scale = _scale_factor(h2_kg.magnitude, mwh_per_kg)
        return ureg.Quantity(h2_kg.magnitude * scale, ureg.MWh)

    def energy_to_mass(
        self,
        energy: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts target delivered H2 energy (MWh, GJ, MMBtu) to required product/carrier mass (kg)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(energy, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(energy).__name__}")
        mwh_q = energy.to(ureg.MWh)
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        scale = _scale_factor(mwh_q.magnitude, mwh_per_kg)
        req_h2_kg = ureg.Quantity(mwh_q.magnitude / scale, ureg.kilogram)
        return self.h2_to_carrier_mass(req_h2_kg, cracking_efficiency=cracking_efficiency).to(ureg.kilogram)

    def volume_to_energy(
        self,
        volume: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts storage/transport volume (m^3) directly to delivered H2 energy (MWh)."""
        return self.mass_to_energy(
            self.volume_to_mass(volume),
            basis=basis,
            cracking_efficiency=cracking_efficiency,
        )

    def energy_to_volume(
        self,
        energy: pint.Quantity,
        *,
        basis: str = "lhv",
        cracking_efficiency: float | None = None,
    ) -> pint.Quantity:
        """Converts target delivered H2 energy (MWh, GJ) directly to storage/transport volume (m^3)."""
        return self.mass_to_volume(
            self.energy_to_mass(energy, basis=basis, cracking_efficiency=cracking_efficiency)
        )

    def electrolysis_power_required(
        self,
        h2_mass: pint.Quantity,
        *,
        electrolyser_efficiency_lhv: float = 0.65,
    ) -> pint.Quantity:
        """Computes electrical energy input (MWh_e) required to produce h2_mass via electrolysis.

        At 65% LHV efficiency, 1 kg H2 requires ~51.26 kWh_e (0.05126 MWh_e).
        """
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not math.isfinite(electrolyser_efficiency_lhv) or not (0.0 < electrolyser_efficiency_lhv <= 1.0):
            raise ValueError(
                f"electrolyser_efficiency_lhv must be in (0, 1], got {electrolyser_efficiency_lhv}"
            )
        kg_q = self._validate_mass(h2_mass)
        mwh_e_per_kg = (self.lhv_mj_per_kg_h2 / 3600.0) / electrolyser_efficiency_lhv
        scale = _scale_factor(kg_q.magnitude, mwh_e_per_kg)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.MWh)

    # --- Commercial Pricing Parity ($/kg H2 <-> $/MWh, $/MMBtu, $/mt Carrier) ---
    def kg_price_to_mwh_price(self, price_per_kg_h2: float, *, basis: str = "lhv") -> float:
        """Converts $/kg H2 quote to $/MWh energy parity."""
        p = float(price_per_kg_h2)
        if not math.isfinite(p):
            raise ValueError(f"price_per_kg_h2 must be finite, got {price_per_kg_h2}")
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        return (p / mwh_per_kg)

    def mwh_price_to_kg_price(self, price_per_mwh: float, *, basis: str = "lhv") -> float:
        """Converts $/MWh energy quote to $/kg H2 parity."""
        p = float(price_per_mwh)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mwh must be finite, got {price_per_mwh}")
        mwh_per_kg = self._mj_per_kg_h2(basis) / 3600.0
        return (p * mwh_per_kg)

    def kg_price_to_mmbtu_price(self, price_per_kg_h2: float, *, basis: str = "hhv") -> float:
        """Converts $/kg H2 quote to $/MMBtu natural gas switching parity."""
        from core.units.conversions import MJ_PER_BTU_IT
        p = float(price_per_kg_h2)
        if not math.isfinite(p):
            raise ValueError(f"price_per_kg_h2 must be finite, got {price_per_kg_h2}")
        mmbtu_per_kg = self._mj_per_kg_h2(basis) / (MJ_PER_BTU_IT * 1_000_000.0)
        return (p / mmbtu_per_kg)

    def mmbtu_price_to_kg_price(self, price_per_mmbtu: float, *, basis: str = "hhv") -> float:
        """Converts $/MMBtu quote to $/kg H2 parity."""
        from core.units.conversions import MJ_PER_BTU_IT
        p = float(price_per_mmbtu)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mmbtu must be finite, got {price_per_mmbtu}")
        mmbtu_per_kg = self._mj_per_kg_h2(basis) / (MJ_PER_BTU_IT * 1_000_000.0)
        return (p * mmbtu_per_kg)

    def carrier_mt_price_to_h2_kg_price(
        self,
        price_per_mt_carrier: float,
        *,
        cracking_efficiency: float | None = None,
        cracking_cost_per_kg_h2: float = 0.0,
    ) -> float:
        """Converts carrier $/mt (e.g. $/mt NH3) to implied delivered $/kg H2."""
        p = float(price_per_mt_carrier)
        cost = float(cracking_cost_per_kg_h2)
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(p) or not math.isfinite(cost) or not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError("Inputs must be finite and cracking_efficiency must be in (0, 1].")
        delivered_kg_h2_per_mt = 1000.0 * self.h2_mass_fraction * eff
        return ((p / delivered_kg_h2_per_mt) + cost)

    def h2_kg_price_to_carrier_mt_price(
        self,
        price_per_kg_h2: float,
        *,
        cracking_efficiency: float | None = None,
        cracking_cost_per_kg_h2: float = 0.0,
    ) -> float:
        """Computes implied carrier $/mt netback from delivered $/kg H2."""
        p = float(price_per_kg_h2)
        cost = float(cracking_cost_per_kg_h2)
        eff = self.default_cracking_efficiency if cracking_efficiency is None else float(cracking_efficiency)
        if not math.isfinite(p) or not math.isfinite(cost) or not math.isfinite(eff) or not (0.0 < eff <= 1.0):
            raise ValueError("Inputs must be finite and cracking_efficiency must be in (0, 1].")
        delivered_kg_h2_per_mt = 1000.0 * self.h2_mass_fraction * eff
        return ((p - cost) * delivered_kg_h2_per_mt)

    def blend_with(
        self,
        other: HydrogenAssay,
        *,
        self_h2_mass_fraction: float,
        blend_name: str | None = None,
    ) -> HydrogenAssay:
        """Blends two compatible H2 streams (e.g. Green + Blue H2 mass-balancing).

        Requires identical physical state and carrier fraction; weights carbon
        intensity and purity linearly by H2 mass.
        """
        if not isinstance(other, HydrogenAssay):
            raise TypeError(f"Expected HydrogenAssay, got {type(other).__name__}")
        if self.state != other.state or abs(self.h2_mass_fraction - other.h2_mass_fraction) > 1e-6:
            raise ValueError("Can only blend hydrogen streams with matching state and h2_mass_fraction.")
        for attribute in (
            "density_kg_per_m3", "reference_temp_c", "reference_pressure_bar",
            "lhv_mj_per_kg_h2", "hhv_mj_per_kg_h2",
        ):
            if not math.isclose(getattr(self, attribute), getattr(other, attribute), rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"Cannot blend streams with different {attribute}; normalise their physical basis first.")
        if not math.isfinite(self_h2_mass_fraction) or not (0.0 <= self_h2_mass_fraction <= 1.0):
            raise ValueError("self_h2_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_h2_mass_fraction)
        w2 = 1.0 - w1
        if w1 == 1.0:
            b_ci = self.carbon_intensity_kg_co2e_per_kg_h2
        elif w1 == 0.0:
            b_ci = other.carbon_intensity_kg_co2e_per_kg_h2
        elif self.carbon_intensity_kg_co2e_per_kg_h2 is None or other.carbon_intensity_kg_co2e_per_kg_h2 is None:
            b_ci = None
        else:
            b_ci = self.carbon_intensity_kg_co2e_per_kg_h2 * w1 + other.carbon_intensity_kg_co2e_per_kg_h2 * w2
        return HydrogenAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            state=self.state,
            density_kg_per_m3=self.density_kg_per_m3,
            h2_mass_fraction=self.h2_mass_fraction,
            reference_temp_c=self.reference_temp_c,
            reference_pressure_bar=self.reference_pressure_bar,
            lhv_mj_per_kg_h2=self.lhv_mj_per_kg_h2,
            hhv_mj_per_kg_h2=self.hhv_mj_per_kg_h2,
            default_cracking_efficiency=(self.default_cracking_efficiency * w1) + (other.default_cracking_efficiency * w2),
            purity_pct=((self.purity_pct * w1) + (other.purity_pct * w2)),
            carbon_intensity_kg_co2e_per_kg_h2=b_ci,
        )

    # --- Global Hydrogen & Carrier Benchmark Presets ---
    @classmethod
    def gaseous_h2_ntp(cls) -> HydrogenAssay:
        """Pure Gaseous Hydrogen at NTP (0°C, 1 atm — Nm^3 metering basis)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.gaseous_h2_ntp")

    @classmethod
    def compressed_h2_350bar(cls) -> HydrogenAssay:
        """Compressed Hydrogen at 350 bar, 15°C (heavy-duty mobility & tube trailers)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.compressed_h2_350bar")

    @classmethod
    def compressed_h2_700bar(cls) -> HydrogenAssay:
        """Compressed Hydrogen at 700 bar, 15°C (light-duty FCEV standard)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.compressed_h2_700bar")

    @classmethod
    def liquid_h2(cls) -> HydrogenAssay:
        """Cryogenic Liquid Hydrogen (LH2 at -252.87°C, 1 atm)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.liquid_h2")

    @classmethod
    def ammonia_carrier(cls) -> HydrogenAssay:
        """Liquid Ammonia as H2 Carrier (NH3 at -33°C, 17.755 wt% H2, ~85% cracking efficiency)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.ammonia_carrier")

    @classmethod
    def lohc_dbt(cls) -> HydrogenAssay:
        """Liquid Organic Hydrogen Carrier (H18-Dibenzyltoluene, 6.2 wt% H2, ~90% release)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.lohc_dbt")

    @classmethod
    def green_h2_rfnbo(cls) -> HydrogenAssay:
        """Electrolytic Green Hydrogen (illustrative low-CI stream, 0.4 kg CO2e/kg H2)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.green_h2_rfnbo")

    @classmethod
    def blue_h2_atr(cls) -> HydrogenAssay:
        """Autothermal Reforming Blue Hydrogen with CCS (~2.0 kg CO2e/kg H2)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "alternative_fuels", "HydrogenAssay.blue_h2_atr")


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class HydrogenCarrierAssay(MaterialSpecification):
    carrier: str
    h2_mass_fraction: Decimal
    cracking_recovery_fraction: Decimal
    hydrogen_heat_content: object

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        if self.carrier not in {"ammonia", "LOHC", "methanol", "other"}:
            raise ValueError("Declare ammonia, LOHC, methanol or other carrier.")
        for name in ("h2_mass_fraction", "cracking_recovery_fraction"):
            object.__setattr__(self, name, _converters_petroleum_fraction(getattr(self, name), name))
        from core.units.heat_content import HeatContentAssay
        if not isinstance(self.hydrogen_heat_content, HeatContentAssay):
            raise TypeError("Supply an explicit pure-H2 HeatContentAssay.")
        self.hydrogen_heat_content.energy_density("MJ/kg", kind="gcv")
        self.hydrogen_heat_content.energy_density("MJ/kg", kind="ncv")

    def delivered_hydrogen(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, self.h2_mass_fraction * self.cracking_recovery_fraction)

    def carrier_required(self, quantity):
        from core.units.assays.base import scale_mass
        factor = self.h2_mass_fraction * self.cracking_recovery_fraction
        if factor == 0:
            raise ValueError("Cannot invert zero delivered-hydrogen recovery.")
        return scale_mass(quantity, 1 / factor)

    def as_heat_content(self, **overrides):
        from core.units.heat_content import HeatContentAssay
        factor = self.h2_mass_fraction * self.cracking_recovery_fraction
        specification = HeatContentAssay(name=self.name, gcv=self.hydrogen_heat_content.energy_density("MJ/kg", kind="gcv") * factor,
                                ncv=self.hydrogen_heat_content.energy_density("MJ/kg", kind="ncv") * factor,
                                density=self.density, energy_role="recovered_hydrogen")
        from dataclasses import replace
        return replace(specification, **overrides)

    def delivered_energy(self, quantity, *, basis="lhv", target_unit="MWh"):
        from core.units.heat_content import convert_heat_content
        return convert_heat_content(quantity, target_unit, assay=self.as_heat_content(), basis=basis)

    def combustion_energy(self, quantity, *, basis="lhv", target_unit="MWh"):
        """Carrier combustion uses its separately declared heat_content field."""
        from core.units.heat_content import convert_heat_content
        if self.heat_content is None or self.heat_content.energy_role != "combustion":
            raise ValueError("Declare a separate carrier combustion heat_content specification.")
        return convert_heat_content(quantity, target_unit, assay=self.heat_content, basis=basis)

    def cracking_balance(self, quantity, *, cracking_energy_per_kg_carrier, basis="lhv"):
        from core.units.processing import HydrogenCrackingSpecification
        return HydrogenCrackingSpecification(self.h2_mass_fraction, self.cracking_recovery_fraction,
            self.hydrogen_heat_content, cracking_energy_per_kg_carrier, self.heat_content).balance(quantity, basis=basis)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (HydrogenCarrierAssay,))
__all__ = ('HydrogenAssay', 'HydrogenCarrierAssay', 'convert_hydrogen_carrier', 'create_hydrogen_carrier_context')
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
    from core.units.assays.energy.hydrogen import convert_hydrogen_carrier
    from core.units.assays.energy.hydrogen import create_hydrogen_carrier_context
