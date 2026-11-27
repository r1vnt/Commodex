"""Commodex assays industrial chemicals implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from pathlib import Path
from types import MappingProxyType
from typing import Final, Literal
from typing import TYPE_CHECKING, Any, Final
import inspect
import json
import math
import pint
import warnings




def dilution_addition(quantity, source_fraction, target_fraction, diluent_fraction):
    from core.units.assays.base import scale_mass
    from core.units.assays.energy.petroleum import _converters_petroleum_fraction
    source = _converters_petroleum_fraction(source_fraction, "source concentration")
    target = _converters_petroleum_fraction(target_fraction, "target concentration")
    diluent = _converters_petroleum_fraction(diluent_fraction, "diluent concentration")
    if target == source:
        return scale_mass(quantity, 0)
    if target == diluent or not min(source, diluent) < target < max(source, diluent):
        raise ValueError("Target concentration must lie between feed and diluent concentrations.")
    return scale_mass(quantity, (source - target) / (target - diluent))


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class ChemicalSolutionAssay(MaterialSpecification):
    solute_formula: str
    active_mass_fraction: Decimal
    solute_molar_mass: pint.Quantity | None = None

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.composition import StoichiometricCompound
        from core.units.registry import ureg
        MaterialSpecification.__post_init__(self)
        StoichiometricCompound(self.solute_formula)
        fraction = _converters_petroleum_fraction(self.active_mass_fraction, "active mass fraction")
        if fraction == 0:
            raise ValueError("Declare positive active ingredient concentration.")
        object.__setattr__(self, "active_mass_fraction", fraction)
        if self.solute_molar_mass is not None:
            object.__setattr__(self, "solute_molar_mass", ureg.Quantity(_converters_petroleum_positive_quantity(self.solute_molar_mass, "g/mol", "declared solute molar mass"), "g/mol"))

    @property
    def molar_mass(self):
        from core.units.composition import StoichiometricCompound
        return self.solute_molar_mass if self.solute_molar_mass is not None else StoichiometricCompound(self.solute_formula).molar_mass

    def as_context(self, *, name=None):
        # Moles denote active solute; the mass endpoint denotes total solution.
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.base import add_molar_context
        return add_molar_context(MaterialSpecification.as_context(self, name=name), self.molar_mass / self.active_mass_fraction)

    def as_heat_content(self, **overrides):
        from core.units.assays.base import MaterialSpecification
        from dataclasses import replace
        heat = MaterialSpecification.as_heat_content(self, **overrides)
        if heat.molar_mass is None:
            heat = replace(heat, molar_mass=self.molar_mass / self.active_mass_fraction)
        return heat

    def active_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, self.active_mass_fraction)

    def solution_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, 1 / self.active_mass_fraction)

    def active_to_moles(self, quantity):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        return _converters_petroleum_output(mass(quantity) / self.molar_mass.to("kg/mol").magnitude, "mol", quantity)

    def moles_to_active(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        amount = _converters_petroleum_quantity_number(quantity, "mol", "active moles", minimum=0)
        return _converters_petroleum_output(amount * self.molar_mass.to("kg/mol").magnitude, "kg", quantity)

    def solution_to_moles(self, quantity):
        return self.active_to_moles(self.active_mass(quantity))

    @property
    def molarity(self):
        if self.density is None:
            raise ValueError("Molarity requires measured solution density at the declared temperature.")
        return (self.density * self.active_mass_fraction / self.molar_mass).to("mol/liter")

    def diluent_required(self, quantity, *, target_fraction, diluent_fraction=Decimal(0)):
        return dilution_addition(quantity, self.active_mass_fraction, target_fraction, diluent_fraction)

    def diluted_mass(self, quantity, **conditions):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        addition = self.diluent_required(quantity, **conditions)
        return _converters_petroleum_output(mass(quantity) + mass(addition), "kg", quantity).to(quantity.units)

    def to_chemical_assay(self):
        from core.units.assays.industrial.chemicals import ChemicalAssay
        if self.density is None or self.reference_temperature is None:
            raise ValueError("Legacy chemical adaptation requires density and reference temperature.")
        return ChemicalAssay(name=self.name, formula=self.solute_formula,
                             molecular_weight_g_per_mol=float(self.molar_mass.magnitude),
                             concentration_pct=float(self.active_mass_fraction * 100),
                             density_kg_per_m3=float(self.density.magnitude),
                             reference_temp_c=float(self.reference_temperature.to("degC").magnitude))

    @classmethod
    def from_chemical_assay(cls, assay):
        from core.units.registry import ureg
        from core.units.assays.industrial.chemicals import ChemicalAssay
        if not isinstance(assay, ChemicalAssay):
            raise TypeError("Expected ChemicalAssay.")
        result = cls(name=assay.name, solute_formula=assay.formula,
                     active_mass_fraction=Decimal(str(assay.concentration_pct)) / 100,
                     solute_molar_mass=ureg.Quantity(Decimal(str(assay.molecular_weight_g_per_mol)), "g/mol"),
                     density=ureg.Quantity(assay.density_kg_per_m3, "kg/m^3"),
                     reference_temperature=ureg.Quantity(assay.reference_temp_c, "degC"))
        return result

    def to_petrochem_assay(self):
        from core.units.assays.industrial.petrochemicals import PetrochemAssay
        if self.density is None or self.reference_temperature is None:
            raise ValueError("Legacy petrochem adaptation requires density and reference temperature.")
        return PetrochemAssay(name=self.name, formula=self.solute_formula,
                             molecular_weight_g_per_mol=float(self.molar_mass.magnitude),
                             purity_pct=float(self.active_mass_fraction * 100),
                             density_kg_per_m3=float(self.density.magnitude),
                             reference_temp_c=float(self.reference_temperature.to("degC").magnitude))





@dataclass(frozen=True, slots=True, kw_only=True)
class SolventAssay(ChemicalSolutionAssay):
    def pure_solvent_mass(self, quantity):
        return self.active_mass(quantity)

    @classmethod
    def from_petrochem_assay(cls, assay):
        from core.units.registry import ureg
        from core.units.assays.industrial.petrochemicals import PetrochemAssay
        if not isinstance(assay, PetrochemAssay):
            raise TypeError("Expected PetrochemAssay.")
        return cls(name=assay.name, solute_formula=assay.formula,
                   active_mass_fraction=Decimal(str(assay.purity_pct)) / 100,
                   solute_molar_mass=ureg.Quantity(Decimal(str(assay.molecular_weight_g_per_mol)), "g/mol"),
                   density=ureg.Quantity(assay.density_kg_per_m3, "kg/m^3"),
                   reference_temperature=ureg.Quantity(assay.reference_temp_c, "degC"))


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class SulphurProductAssay(MaterialSpecification):
    state: str
    sulphur_mass_fraction: Decimal

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        if self.state not in {"solid", "molten", "acid"}:
            raise ValueError("Sulphur state must be solid, molten or acid.")
        object.__setattr__(self, "sulphur_mass_fraction", _converters_petroleum_fraction(self.sulphur_mass_fraction, "contained sulphur"))

    def sulphur_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, self.sulphur_mass_fraction)

    def acid_equivalent(self, quantity, *, acid_mass_fraction):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.composition import StoichiometricCompound
        strength = _converters_petroleum_fraction(acid_mass_fraction, "acid strength")
        if strength == 0:
            raise ValueError("Acid strength must be positive.")
        pure_acid = StoichiometricCompound("S").equivalent_mass(self.sulphur_mass(quantity), "H2SO4", element="S")
        return scale_mass(pure_acid, 1 / strength)


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class ChlorAlkaliAssay(MaterialSpecification):
    salt_purity_fraction: Decimal
    salt_conversion_fraction: Decimal

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        for name in ("salt_purity_fraction", "salt_conversion_fraction"):
            object.__setattr__(self, name, _converters_petroleum_fraction(getattr(self, name), name))

    def _reacted_salt(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, self.salt_purity_fraction * self.salt_conversion_fraction)

    def co_product_mass(self, quantity, *, product, recovery_fraction):
        # 2 NaCl + 2 H2O -> Cl2 + 2 NaOH + H2 (Euro Chlor BREF).
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.composition import StoichiometricCompound
        coefficients = {"Cl2": Decimal(1), "NaOH": Decimal(2), "H2": Decimal(1)}
        if product not in coefficients:
            raise ValueError("Chlor-alkali product must be Cl2, NaOH or H2.")
        ratio = coefficients[product] * StoichiometricCompound(product).molar_mass.magnitude / (2 * StoichiometricCompound("NaCl").molar_mass.magnitude)
        return scale_mass(self._reacted_salt(quantity), ratio * _converters_petroleum_fraction(recovery_fraction, "co-product recovery"))

    def water_consumed(self, quantity):
        from core.units.assays.base import scale_mass
        from core.units.composition import StoichiometricCompound
        ratio = StoichiometricCompound("H2O").molar_mass.magnitude / StoichiometricCompound("NaCl").molar_mass.magnitude
        return scale_mass(self._reacted_salt(quantity), ratio)


from core.units.validation import validated_assay



# ---------------------------------------------------------------------------
# Industrial & Inorganic Chemicals Specification
# (Caustic Soda, Sulphuric/HCl Acids, Methanol, Soda Ash, Sulphur, Chlorine)
# ---------------------------------------------------------------------------

@validated_assay
@dataclass(frozen=True, slots=True)
class ChemicalAssay:
    """Immutable physical, concentration, and molar specification for industrial chemicals.

    Bridges commercial/solution mass (WMT), 100% active solute mass (DMT),
    liquid/bulk volume (m^3, gal), molar quantities (kmol), Baumé gravity,
    dilution water requirements, and regional pricing parity.
    """

    name: str
    formula: str
    molecular_weight_g_per_mol: float
    concentration_pct: float
    density_kg_per_m3: float
    reference_temp_c: float = 15.0
    is_liquid: bool = True

    def __post_init__(self) -> None:
        if not self.formula or not self.formula.strip():
            raise ValueError("formula must be a non-empty string (e.g. 'NaOH', 'H2SO4', 'CH3OH').")
        if not math.isfinite(self.molecular_weight_g_per_mol) or self.molecular_weight_g_per_mol <= 0.0:
            raise ValueError(
                f"molecular_weight_g_per_mol must be strictly positive, got {self.molecular_weight_g_per_mol}"
            )
        if not math.isfinite(self.concentration_pct) or not (0.0 < self.concentration_pct <= 100.0):
            raise ValueError(f"concentration_pct must be finite and in (0, 100], got {self.concentration_pct}")
        if not math.isfinite(self.density_kg_per_m3) or not (300.0 <= self.density_kg_per_m3 <= 3000.0):
            raise ValueError(
                f"density_kg_per_m3 must be finite and in [300, 3000] kg/m^3, got {self.density_kg_per_m3}"
            )
        if not math.isfinite(self.reference_temp_c) or not (-100.0 <= self.reference_temp_c <= 200.0):
            raise ValueError(f"reference_temp_c must be in [-100, 200] °C, got {self.reference_temp_c}")

    # --- Concentration, Gravity & Volumetric Properties ---
    @property
    def active_mass_fraction(self) -> float:
        """Mass fraction of active chemical solute (e.g. 0.50 for 50% Caustic Soda)."""
        return self.concentration_pct / 100.0

    @property
    def specific_gravity(self) -> float:
        """Specific gravity relative to water at 15°C (999.102 kg/m^3)."""
        from core.units.conversions import WATER_DENSITY_15C_KG_M3
        return (self.density_kg_per_m3 / WATER_DENSITY_15C_KG_M3)

    @property
    def baume_degrees(self) -> float | None:
        """Baumé gravity (°Bé) for liquid chemicals; None for bulk solids.

        Uses the standard US/industrial scales:
        - Heavy liquids (SG >= 1.0): 145 - (145 / SG)
        - Light liquids (SG < 1.0):  (140 / SG) - 130
        """
        from core.units.conversions import WATER_DENSITY_15C_KG_M3
        if not self.is_liquid:
            return None
        sg = self.density_kg_per_m3 / WATER_DENSITY_15C_KG_M3
        if sg >= 1.0:
            return (145.0 - (145.0 / sg))
        return ((140.0 / sg) - 130.0)

    @property
    def gal_per_mt(self) -> float:
        """US gallons per commercial metric tonne at reference_temp_c."""
        from core.units.conversions import _LITRES_PER_US_GAL
        return (1000.0 / (self.density_kg_per_m3 * (_LITRES_PER_US_GAL / 1000.0)))

    @property
    def lb_per_gal(self) -> float:
        """US pounds per US gallon at reference_temp_c."""
        from core.units.conversions import _LB_PER_MT
        return (_LB_PER_MT / self.gal_per_mt)

    @property
    def active_kg_per_m3(self) -> float:
        """Active solute mass (kg) contained in 1 m^3 of commercial solution/bulk."""
        return (self.density_kg_per_m3 * self.active_mass_fraction)

    def as_context(self, *, active_only: bool = True, context_name: str | None = None) -> pint.Context:
        """Builds a Pint transformation context bridging [volume], [mass], and [substance]."""
        from core.units.conversions import create_chemical_context
        return create_chemical_context(self, active_only=active_only, name=context_name)

    # --- Physical Conversions (Volume <-> Solution Mass <-> Active Mass <-> Moles) ---
    def _validate_mass(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        return mass.to(ureg.metric_ton)

    def volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts commercial liquid/bulk volume (m^3, gal, L) to total mass (metric_ton)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        m3_q = volume.to(ureg.meter ** 3)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_volume(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts commercial total mass (mt, short_ton, lb) to volume (m^3)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def solution_to_active_mass(self, solution_mass: pint.Quantity) -> pint.Quantity:
        """Converts commercial/wet solution mass (WMT) to 100% active solute mass (DMT)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(solution_mass)
        scale = _scale_factor(mt_q.magnitude, self.active_mass_fraction)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def active_to_solution_mass(self, active_mass: pint.Quantity) -> pint.Quantity:
        """Converts 100% active solute mass (DMT) to required commercial solution mass (WMT)."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        mt_q = self._validate_mass(active_mass)
        scale = _scale_factor(mt_q.magnitude, self.active_mass_fraction)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.metric_ton)

    def volume_to_active_mass(self, volume: pint.Quantity) -> pint.Quantity:
        """Converts commercial volume (e.g. m^3, gal) directly to 100% active solute mass (metric_ton)."""
        return self.solution_to_active_mass(self.volume_to_mass(volume))

    def active_mass_to_volume(self, active_mass: pint.Quantity) -> pint.Quantity:
        """Converts target 100% active solute mass directly to required commercial volume (m^3)."""
        return self.mass_to_volume(self.active_to_solution_mass(active_mass))

    def mass_to_moles(self, mass: pint.Quantity, *, active_only: bool = True) -> pint.Quantity:
        """Converts mass to kilomoles; active_only applies the solution concentration."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        kg_q = mass.to(ureg.kilogram)
        conc_mult = self.active_mass_fraction if active_only else 1.0
        kmol_per_kg = conc_mult / self.molecular_weight_g_per_mol
        scale = _scale_factor(kg_q.magnitude, kmol_per_kg)
        return ureg.Quantity(kg_q.magnitude * scale, ureg.kilomole)

    def moles_to_mass(self, moles: pint.Quantity, *, active_only: bool = True) -> pint.Quantity:
        """Converts moles to tonnes; active_only returns the required solution mass."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(moles, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(moles).__name__}")
        kmol_q = moles.to(ureg.kilomole)
        conc_mult = self.active_mass_fraction if active_only else 1.0
        mt_per_kmol = (self.molecular_weight_g_per_mol / 1000.0) / conc_mult
        scale = _scale_factor(kmol_q.magnitude, mt_per_kmol)
        return ureg.Quantity(kmol_q.magnitude * scale, ureg.metric_ton)

    # --- Dilution & Strength Blending ---
    def water_required_to_dilute(
        self,
        solution_mass: pint.Quantity,
        target_concentration_pct: float,
    ) -> pint.Quantity:
        """Computes mass of water (metric_ton) required to dilute solution_mass to target_concentration_pct."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        c_target = float(target_concentration_pct)
        if not math.isfinite(c_target) or not (0.0 < c_target < self.concentration_pct):
            raise ValueError(
                f"target_concentration_pct must be strictly in (0, {self.concentration_pct}), got {target_concentration_pct}"
            )
        mt_q = self._validate_mass(solution_mass)
        water_ratio = (self.concentration_pct / c_target) - 1.0
        scale = _scale_factor(mt_q.magnitude, water_ratio)
        return ureg.Quantity(mt_q.magnitude * scale, ureg.metric_ton)

    def dilute_to(
        self,
        target_concentration_pct: float,
        *,
        new_density_kg_per_m3: float | None = None,
        name: str | None = None,
    ) -> ChemicalAssay:
        """Returns a diluted ChemicalAssay at target_concentration_pct.

        If new_density_kg_per_m3 is omitted, estimates density via ideal volume
        additivity with water at 15 degrees C. This approximation does not model
        thermal effects or volume contraction; supply a measured density when needed.
        """
        from core.units.conversions import WATER_DENSITY_15C_KG_M3
        c_target = float(target_concentration_pct)
        if not math.isfinite(c_target) or not (0.0 < c_target < self.concentration_pct):
            raise ValueError(
                f"target_concentration_pct must be strictly in (0, {self.concentration_pct}), got {target_concentration_pct}"
            )
        if new_density_kg_per_m3 is None:
            sol_frac = c_target / self.concentration_pct
            water_frac = 1.0 - sol_frac
            vol_per_kg = (sol_frac / self.density_kg_per_m3) + (water_frac / WATER_DENSITY_15C_KG_M3)
            est_rho = (1.0 / vol_per_kg)
        else:
            est_rho = float(new_density_kg_per_m3)

        return ChemicalAssay(
            name=name or f"{self.formula} ({c_target:g}% Solution)",
            formula=self.formula,
            molecular_weight_g_per_mol=self.molecular_weight_g_per_mol,
            concentration_pct=c_target,
            density_kg_per_m3=est_rho,
            reference_temp_c=self.reference_temp_c,
            is_liquid=True,
        )

    def blend_with(
        self,
        other: ChemicalAssay,
        *,
        self_mass_fraction: float,
        blend_name: str | None = None,
    ) -> ChemicalAssay:
        """Blends two lots/strengths of the same chemical by total mass fraction."""
        if not isinstance(other, ChemicalAssay):
            raise TypeError(f"Expected ChemicalAssay, got {type(other).__name__}")
        if self.formula != other.formula:
            raise ValueError(f"Cannot blend different chemicals: {self.formula!r} vs {other.formula!r}")
        if self.is_liquid != other.is_liquid:
            raise ValueError("Cannot directly blend a liquid solution with a solid grade.")
        if not math.isfinite(self_mass_fraction) or not (0.0 <= self_mass_fraction <= 1.0):
            raise ValueError("self_mass_fraction must be between 0.0 and 1.0.")

        w1 = float(self_mass_fraction)
        w2 = 1.0 - w1
        b_conc = (self.concentration_pct * w1) + (other.concentration_pct * w2)
        vol_per_kg = (w1 / self.density_kg_per_m3) + (w2 / other.density_kg_per_m3)

        return ChemicalAssay(
            name=blend_name or f"{self.name}/{other.name} Blend",
            formula=self.formula,
            molecular_weight_g_per_mol=self.molecular_weight_g_per_mol,
            concentration_pct=b_conc,
            density_kg_per_m3=(1.0 / vol_per_kg),
            reference_temp_c=self.reference_temp_c,
            is_liquid=self.is_liquid,
        )

    # --- Commercial Pricing Parity ($/WMT <-> $/DMT 100%, $/gal <-> $/mt, $/st <-> $/mt) ---
    def solution_price_to_active_price(self, price_per_wmt: float) -> float:
        """Converts a $/WMT commercial solution quote to $/DMT (100% active basis)."""
        p = float(price_per_wmt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_wmt must be finite, got {price_per_wmt}")
        return (p / self.active_mass_fraction)

    def active_price_to_solution_price(self, price_per_dmt_active: float) -> float:
        """Converts a $/DMT (100% active basis) quote to $/WMT commercial solution."""
        p = float(price_per_dmt_active)
        if not math.isfinite(p):
            raise ValueError(f"price_per_dmt_active must be finite, got {price_per_dmt_active}")
        return (p * self.active_mass_fraction)

    def gal_price_to_mt_price(self, price_per_gal: float) -> float:
        """Converts USGC $/gal quote (e.g. Methanol) to international $/mt equivalent."""
        p = float(price_per_gal)
        if not math.isfinite(p):
            raise ValueError(f"price_per_gal must be finite, got {price_per_gal}")
        return (p * self.gal_per_mt)

    def mt_price_to_gal_price(self, price_per_mt: float) -> float:
        """Converts international $/mt quote to USGC $/gal equivalent."""
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p / self.gal_per_mt)

    @staticmethod
    def st_price_to_mt_price(price_per_st: float) -> float:
        """Converts US $/short_ton quote (e.g. US Caustic/Soda Ash) to $/metric_ton."""
        from core.units.conversions import _MT_PER_SHORT_TON
        p = float(price_per_st)
        if not math.isfinite(p):
            raise ValueError(f"price_per_st must be finite, got {price_per_st}")
        return (p / _MT_PER_SHORT_TON)

    @staticmethod
    def mt_price_to_st_price(price_per_mt: float) -> float:
        """Converts $/metric_ton quote to US $/short_ton."""
        from core.units.conversions import _MT_PER_SHORT_TON
        p = float(price_per_mt)
        if not math.isfinite(p):
            raise ValueError(f"price_per_mt must be finite, got {price_per_mt}")
        return (p * _MT_PER_SHORT_TON)

    # --- Global Industrial Chemical Benchmark Presets ---
    @classmethod
    def caustic_soda_50pct(cls) -> ChemicalAssay:
        """Liquid Caustic Soda (50% NaOH aqueous — global chlor-alkali benchmark)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.caustic_soda_50pct")

    @classmethod
    def caustic_soda_flakes(cls) -> ChemicalAssay:
        """Solid Caustic Soda Flakes (99% NaOH)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.caustic_soda_flakes")

    @classmethod
    def sulphuric_acid_98pct(cls) -> ChemicalAssay:
        """Commercial Sulphuric Acid (98% H2SO4 — smelter/fertiliser benchmark)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.sulphuric_acid_98pct")

    @classmethod
    def sulphuric_acid_66be(cls) -> ChemicalAssay:
        """66°Bé Sulphuric Acid (93.19% H2SO4 — North American commercial grade)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.sulphuric_acid_66be")

    @classmethod
    def hydrochloric_acid_32pct(cls) -> ChemicalAssay:
        """Commercial Hydrochloric Acid (32% HCl, ~20°Bé)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.hydrochloric_acid_32pct")

    @classmethod
    def methanol(cls) -> ChemicalAssay:
        """IMPCA-Grade Methanol (99.85% CH3OH — USGC $/gal vs. Rotterdam/China $/mt)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.methanol")

    @classmethod
    def soda_ash_dense(cls) -> ChemicalAssay:
        """Dense Soda Ash (99.2% Na2CO3 — glass & lithium refining benchmark)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.soda_ash_dense")

    @classmethod
    def molten_sulphur(cls) -> ChemicalAssay:
        """Molten Elemental Sulphur (99.5% S — pumped at ~135°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.molten_sulphur")

    @classmethod
    def hydrogen_peroxide_50pct(cls) -> ChemicalAssay:
        """Industrial Hydrogen Peroxide (50% H2O2 aqueous)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.hydrogen_peroxide_50pct")

    @classmethod
    def liquid_chlorine(cls) -> ChemicalAssay:
        """Liquid Chlorine (99.5% Cl2 — pressurised at 0°C)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "chemicals", "ChemicalAssay.liquid_chlorine")




_converters_chemical_BRIX_DATA_FILE: Final[Path] = (Path(__file__).resolve().parents[2] / "reference_data") / "units" / "icumsa_brix.json"


from core.units.registry import ureg

_converters_chemical_MOLARITY_DIM: Final = ureg.molarity.dimensionality


from core.units.registry import ureg
_converters_chemical_NORMALITY_DIM: Final = ureg.normality.dimensionality


from core.units.registry import ureg
_converters_chemical_BAUME_HEAVY_DIM: Final = ureg.degree_Baume_heavy.dimensionality


from core.units.registry import ureg
_converters_chemical_BAUME_LIGHT_DIM: Final = ureg.degree_Baume_light.dimensionality


from core.units.registry import ureg
_converters_chemical_PLATO_DIM: Final = ureg.degree_Plato.dimensionality


from core.units.registry import ureg
_converters_chemical_BRIX_DIM: Final = ureg.degree_Brix.dimensionality





def _converters_chemical_load_icumsa_brix(path: Path = _converters_chemical_BRIX_DATA_FILE) -> tuple[MonotonicCubicTable, float]:
    from core.units.thermodynamics import MonotonicCubicTable
    payload = json.loads(path.read_text(encoding="utf-8"))
    if payload.get("schema_version") != 1:
        raise ValueError(f"Unsupported icumsa_brix.json schema_version: {payload.get('schema_version')!r}")
    return MonotonicCubicTable(payload["table_brix_sg"]), float(payload["reference_temp_c"])





_converters_chemical_ICUMSA_BRIX_TABLE, _converters_chemical_ICUMSA_REF_TEMP_C = _converters_chemical_load_icumsa_brix()




# ASBC 20 °C / 20 °C Plato cubic coefficients:
# °P = -616.868 + 1111.14 * SG - 630.272 * SG^2 + 135.997 * SG^3
_converters_chemical_ASBC_C0: Final[float] = -616.868



_converters_chemical_ASBC_C1: Final[float] = 1111.14



_converters_chemical_ASBC_C2: Final[float] = -630.272



_converters_chemical_ASBC_C3: Final[float] = 135.997





def _converters_chemical_plato_poly(sg: float) -> float:
    return _converters_chemical_ASBC_C0 + sg * (_converters_chemical_ASBC_C1 + sg * (_converters_chemical_ASBC_C2 + sg * _converters_chemical_ASBC_C3))





def _converters_chemical_plato_poly_deriv(sg: float) -> float:
    return _converters_chemical_ASBC_C1 + sg * (2.0 * _converters_chemical_ASBC_C2 + 3.0 * _converters_chemical_ASBC_C3 * sg)





def _converters_chemical_solve_plato_sg(plato: float) -> float:
    sg = 1.0 + 0.004 * plato
    for _ in range(12):
        f_val = _converters_chemical_plato_poly(sg) - plato
        if abs(f_val) <= 1e-15:
            break
        sg -= f_val / _converters_chemical_plato_poly_deriv(sg)
    return sg





_converters_chemical_PLATO_MAX_DEG: Final[float] = 30.0



_converters_chemical_PLATO_MAX_SG: Final[float] = _converters_chemical_solve_plato_sg(_converters_chemical_PLATO_MAX_DEG)





def _converters_chemical_validate_20c_reference(ref_temp_c: float | int | Decimal, scale_name: str) -> None:
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    t_dec = _converters__interpolation_validate_real_scalar(ref_temp_c, "ref_temp_c")
    if abs(t_dec - Decimal(20)) > Decimal("1e-9"):
        raise StandardsApplicabilityError(
            f"{scale_name} calibration requires hydrometer specific gravity at 20 °C / 20 °C "
            f"(ref_temp_c=20.0), got {ref_temp_c}°C. Correct readings to 20 °C first."
        )





def _converters_chemical_extract_dimensionless_sg(sg: pint.Quantity | float | int | Decimal, label: str = "sg") -> tuple[Decimal, bool]:
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if isinstance(sg, pint.Quantity):
        validate_registry_quantity(sg)
        if not sg.dimensionless:
            raise pint.DimensionalityError(
                sg.units, ureg.dimensionless, sg.dimensionality, ureg.dimensionless.dimensionality
            )
        _converters__interpolation_validate_real_scalar(sg.magnitude, label)
        reduced = sg.to(ureg.dimensionless)
        return _converters__interpolation_validate_real_scalar(reduced.magnitude, label), True
    return _converters__interpolation_validate_real_scalar(sg, label), False


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# 1. Normality <-> Molarity (N = M * z)
# ---------------------------------------------------------------------------


@authorized_conversion
def molarity_to_normality(
    molarity: pint.Quantity | float | int | Decimal,
    equivalence_factor: float | int | Decimal,
) -> pint.Quantity | float | Decimal:
    """Convert molarity (M, mol/L) to normality (N = M * z) using explicit ``equivalence_factor > 0``."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    z_dec = _converters__interpolation_validate_real_scalar(equivalence_factor, "equivalence_factor")
    if z_dec <= 0:
        raise ValueError(f"equivalence_factor must be strictly positive, got {equivalence_factor}.")

    if isinstance(molarity, pint.Quantity):
        validate_registry_quantity(molarity)
        if molarity.dimensionality != _converters_chemical_MOLARITY_DIM:
            raise pint.DimensionalityError(
                molarity.units, ureg.molarity, molarity.dimensionality, _converters_chemical_MOLARITY_DIM
            )
        m_raw = _converters__interpolation_validate_real_scalar(molarity.magnitude, "molarity")
        if m_raw < 0:
            raise ValueError("Molarity must be non-negative.")
        m_mol_l = _converters__interpolation_validate_real_scalar(molarity.to(ureg.molarity).magnitude, "molarity")
        n_dec = m_mol_l * z_dec
        return ureg.Quantity(_converters__interpolation_match_scalar_type(n_dec, molarity, equivalence_factor), ureg.normality)

    m_dec = _converters__interpolation_validate_real_scalar(molarity, "molarity")
    if m_dec < 0:
        raise ValueError("Molarity must be non-negative.")
    return _converters__interpolation_match_scalar_type(m_dec * z_dec, molarity, equivalence_factor)


from core.units.validation import authorized_conversion


@authorized_conversion
def normality_to_molarity(
    normality: pint.Quantity | float | int | Decimal,
    equivalence_factor: float | int | Decimal,
    *,
    target_unit: str | pint.Unit = "molarity",
) -> pint.Quantity | float | Decimal:
    """Convert normality (N, eq/L) to molarity (M = N / z) using explicit ``equivalence_factor > 0``."""
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    z_dec = _converters__interpolation_validate_real_scalar(equivalence_factor, "equivalence_factor")
    if z_dec <= 0:
        raise ValueError(f"equivalence_factor must be strictly positive, got {equivalence_factor}.")

    if isinstance(normality, pint.Quantity):
        validate_registry_quantity(normality)
        if normality.dimensionality != _converters_chemical_NORMALITY_DIM:
            raise pint.DimensionalityError(
                normality.units, ureg.normality, normality.dimensionality, _converters_chemical_NORMALITY_DIM
            )
        n_raw = _converters__interpolation_validate_real_scalar(normality.magnitude, "normality")
        if n_raw < 0:
            raise ValueError("Normality must be non-negative.")
        n_eq_l = _converters__interpolation_validate_real_scalar(normality.to(ureg.normality).magnitude, "normality")
        dst = parse_unit(target_unit)
        if dst.dimensionality != _converters_chemical_MOLARITY_DIM:
            raise pint.DimensionalityError(ureg.molarity, dst, _converters_chemical_MOLARITY_DIM, dst.dimensionality)
        m_q = ureg.Quantity(n_eq_l / z_dec, ureg.molarity).to(dst)
        return ureg.Quantity(_converters__interpolation_match_scalar_type(m_q.magnitude, normality, equivalence_factor), dst)

    n_dec = _converters__interpolation_validate_real_scalar(normality, "normality")
    if n_dec < 0:
        raise ValueError("Normality must be non-negative.")
    return _converters__interpolation_match_scalar_type(n_dec / z_dec, normality, equivalence_factor)





# ---------------------------------------------------------------------------
# 2. Heavy and Light Baumé <-> Specific Gravity (SG)
# ---------------------------------------------------------------------------


def _converters_chemical_resolve_baume_moduli(
    modulus: Literal["us", "rational"] | float | int | Decimal,
) -> tuple[Decimal, Decimal, Decimal]:
    """Return (heavy_modulus, light_numerator, light_offset) as exact Decimals."""
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if isinstance(modulus, str):
        key = modulus.strip().lower()
        if key == "us":
            # US 60 °F / 60 °F standard: heavy = 145 - 145/SG; light = 140/SG - 130
            return Decimal(145), Decimal(140), Decimal(130)
        if key in {"rational", "european"}:
            # Rational / European 144.3 scale: heavy = 144.3 - 144.3/SG; light = 144.3/SG - 134.3
            return Decimal("144.3"), Decimal("144.3"), Decimal("134.3")
        raise ValueError(f"Unknown Baumé modulus {modulus!r}; use 'us', 'rational', or a positive numeric modulus.")

    mod_dec = _converters__interpolation_validate_real_scalar(modulus, "modulus")
    if mod_dec <= Decimal(10):
        raise ValueError(f"Custom Baumé modulus must exceed 10, got {modulus}.")
    return mod_dec, mod_dec, mod_dec - Decimal(10)


from core.units.validation import authorized_conversion


@authorized_conversion
def sg_to_baume(
    sg: pint.Quantity | float | int | Decimal,
    *,
    scale: Literal["heavy", "light"],
    modulus: Literal["us", "rational"] | float | int | Decimal = "us",
) -> pint.Quantity | float | Decimal:
    """Convert Specific Gravity (SG) to heavy or light Degrees Baumé (°Bé)."""
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    if not isinstance(scale, str) or scale.strip().lower() not in {"heavy", "light"}:
        raise ValueError(f"scale must be 'heavy' or 'light', got {scale!r}.")
    branch = scale.strip().lower()
    heavy_mod, light_num, light_off = _converters_chemical_resolve_baume_moduli(modulus)
    sg_dec, is_qty = _converters_chemical_extract_dimensionless_sg(sg, "sg")
    if sg_dec <= 0:
        raise ValueError(f"Specific gravity must be strictly positive, got {sg}.")

    if branch == "heavy":
        if sg_dec < Decimal(1):
            raise ValueError(f"Heavy Baumé requires SG >= 1.0, got {sg_dec}.")
        be_dec = heavy_mod - (heavy_mod / sg_dec)
        out = _converters__interpolation_match_scalar_type(be_dec, sg)
        return ureg.Quantity(out, ureg.degree_Baume_heavy) if is_qty else out

    # branch == "light"
    if sg_dec > Decimal(1):
        raise ValueError(f"Light Baumé requires SG <= 1.0, got {sg_dec}.")
    be_dec = (light_num / sg_dec) - light_off
    out = _converters__interpolation_match_scalar_type(be_dec, sg)
    return ureg.Quantity(out, ureg.degree_Baume_light) if is_qty else out


from core.units.validation import authorized_conversion


@authorized_conversion
def baume_to_sg(
    baume: pint.Quantity | float | int | Decimal,
    *,
    scale: Literal["heavy", "light"],
    modulus: Literal["us", "rational"] | float | int | Decimal = "us",
) -> pint.Quantity | float | Decimal:
    """Convert heavy or light Degrees Baumé (°Bé) to Specific Gravity (SG)."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    if not isinstance(scale, str) or scale.strip().lower() not in {"heavy", "light"}:
        raise ValueError(f"scale must be 'heavy' or 'light', got {scale!r}.")
    branch = scale.strip().lower()
    heavy_mod, light_num, light_off = _converters_chemical_resolve_baume_moduli(modulus)

    is_qty = isinstance(baume, pint.Quantity)
    if is_qty:
        validate_registry_quantity(baume)
        expected_unit = ureg.degree_Baume_heavy if branch == "heavy" else ureg.degree_Baume_light
        expected_dim = _converters_chemical_BAUME_HEAVY_DIM if branch == "heavy" else _converters_chemical_BAUME_LIGHT_DIM
        if baume.dimensionality != expected_dim:
            raise pint.DimensionalityError(baume.units, expected_unit, baume.dimensionality, expected_dim)
        _converters__interpolation_validate_real_scalar(baume.magnitude, "baume")
        be_dec = _converters__interpolation_validate_real_scalar(baume.to(expected_unit).magnitude, "baume")
    else:
        be_dec = _converters__interpolation_validate_real_scalar(baume, "baume")

    if branch == "heavy":
        if be_dec < Decimal(0) or be_dec >= heavy_mod:
            raise ValueError(
                f"Heavy Baumé must satisfy 0 <= °Bé < {heavy_mod}, got {be_dec}."
            )
        sg_dec = heavy_mod / (heavy_mod - be_dec)
    else:
        if be_dec < Decimal(10):
            raise ValueError(
                f"Light Baumé must satisfy °Bé >= 10 (SG <= 1.0), got {be_dec}."
            )
        sg_dec = light_num / (light_off + be_dec)

    out = _converters__interpolation_match_scalar_type(sg_dec, baume)
    return ureg.Quantity(out, ureg.dimensionless) if is_qty else out





def sg_to_baume_heavy(
    sg: pint.Quantity | float | int | Decimal,
    *,
    modulus: Literal["us", "rational"] | float | int | Decimal = "us",
) -> pint.Quantity | float | Decimal:
    """Convert Specific Gravity (SG >= 1.0) to heavy Baumé (°Bé)."""
    return sg_to_baume(sg, scale="heavy", modulus=modulus)





def baume_heavy_to_sg(
    baume: pint.Quantity | float | int | Decimal,
    *,
    modulus: Literal["us", "rational"] | float | int | Decimal = "us",
) -> pint.Quantity | float | Decimal:
    """Convert heavy Baumé (0 <= °Bé < modulus) to Specific Gravity (SG)."""
    return baume_to_sg(baume, scale="heavy", modulus=modulus)





def sg_to_baume_light(
    sg: pint.Quantity | float | int | Decimal,
    *,
    modulus: Literal["us", "rational"] | float | int | Decimal = "us",
) -> pint.Quantity | float | Decimal:
    """Convert Specific Gravity (0 < SG <= 1.0) to light Baumé (°Bé)."""
    return sg_to_baume(sg, scale="light", modulus=modulus)





def baume_light_to_sg(
    baume: pint.Quantity | float | int | Decimal,
    *,
    modulus: Literal["us", "rational"] | float | int | Decimal = "us",
) -> pint.Quantity | float | Decimal:
    """Convert light Baumé (°Bé >= 10) to Specific Gravity (SG)."""
    return baume_to_sg(baume, scale="light", modulus=modulus)


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# 3. Plato <-> Specific Gravity (ASBC 20 °C / 20 °C Cubic, 0–30 °P)
# ---------------------------------------------------------------------------


@authorized_conversion
def sg_to_plato(
    sg: pint.Quantity | float | int | Decimal,
    *,
    ref_temp_c: float = 20.0,
) -> pint.Quantity | float | Decimal:
    """Convert Specific Gravity (20 °C / 20 °C) to Degrees Plato (°P) via ASBC cubic (valid 0–30 °P)."""
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    _converters_chemical_validate_20c_reference(ref_temp_c, "Plato (ASBC)")
    sg_dec, is_qty = _converters_chemical_extract_dimensionless_sg(sg, "sg")
    sg_flt = float(sg_dec)
    if sg_flt < 1.0 or sg_flt > _converters_chemical_PLATO_MAX_SG + 1e-6:
        raise StandardsApplicabilityError(
            f"ASBC Plato polynomial is valid for SG (20/20 °C) in [1.0, {_converters_chemical_PLATO_MAX_SG:.5f}] (0–30 °P), got {sg_dec}."
        )
    plato_val = _converters_chemical_plato_poly(sg_flt)
    out = _converters__interpolation_match_scalar_type(plato_val, sg)
    return ureg.Quantity(out, ureg.degree_Plato) if is_qty else out


from core.units.validation import authorized_conversion


@authorized_conversion
def plato_to_sg(
    plato: pint.Quantity | float | int | Decimal,
    *,
    ref_temp_c: float = 20.0,
) -> pint.Quantity | float | Decimal:
    """Invert Degrees Plato (0–30 °P) to Specific Gravity (20 °C / 20 °C) via Newton-Raphson."""
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    _converters_chemical_validate_20c_reference(ref_temp_c, "Plato (ASBC)")
    is_qty = isinstance(plato, pint.Quantity)
    if is_qty:
        validate_registry_quantity(plato)
        if plato.dimensionality != _converters_chemical_PLATO_DIM:
            raise pint.DimensionalityError(plato.units, ureg.degree_Plato, plato.dimensionality, _converters_chemical_PLATO_DIM)
        _converters__interpolation_validate_real_scalar(plato.magnitude, "plato")
        p_flt = float(plato.to(ureg.degree_Plato).magnitude)
    else:
        p_flt = float(_converters__interpolation_validate_real_scalar(plato, "plato"))

    if p_flt < 0.0 or p_flt > _converters_chemical_PLATO_MAX_DEG:
        raise StandardsApplicabilityError(f"ASBC Plato conversion is valid for 0 <= °P <= 30, got {p_flt}.")

    sg_val = _converters_chemical_solve_plato_sg(p_flt)
    out = _converters__interpolation_match_scalar_type(sg_val, plato)
    return ureg.Quantity(out, ureg.dimensionless) if is_qty else out


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# 4. Brix <-> Specific Gravity (ICUMSA 20 °C / 20 °C Sucrose Table, 0–85 °Bx)
# ---------------------------------------------------------------------------


@authorized_conversion
def brix_to_sg(
    brix: pint.Quantity | float | int | Decimal,
    *,
    ref_temp_c: float = 20.0,
) -> pint.Quantity | float | Decimal:
    """Convert Degrees Brix (0–85 °Bx) to Specific Gravity (20 °C / 20 °C) via ICUMSA SPS-4 table."""
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    _converters_chemical_validate_20c_reference(ref_temp_c, "Brix (ICUMSA)")
    is_qty = isinstance(brix, pint.Quantity)
    if is_qty:
        validate_registry_quantity(brix)
        if brix.dimensionality != _converters_chemical_BRIX_DIM:
            raise pint.DimensionalityError(brix.units, ureg.degree_Brix, brix.dimensionality, _converters_chemical_BRIX_DIM)
        _converters__interpolation_validate_real_scalar(brix.magnitude, "brix")
        bx_flt = float(brix.to(ureg.degree_Brix).magnitude)
    else:
        bx_flt = float(_converters__interpolation_validate_real_scalar(brix, "brix"))

    if bx_flt < _converters_chemical_ICUMSA_BRIX_TABLE.min_x or bx_flt > _converters_chemical_ICUMSA_BRIX_TABLE.max_x:
        raise StandardsApplicabilityError(
            f"ICUMSA Brix table is valid for {_converters_chemical_ICUMSA_BRIX_TABLE.min_x} <= °Bx <= {_converters_chemical_ICUMSA_BRIX_TABLE.max_x}, got {bx_flt}."
        )

    sg_val = _converters_chemical_ICUMSA_BRIX_TABLE.evaluate(bx_flt)
    out = _converters__interpolation_match_scalar_type(sg_val, brix)
    return ureg.Quantity(out, ureg.dimensionless) if is_qty else out


from core.units.validation import authorized_conversion


@authorized_conversion
def sg_to_brix(
    sg: pint.Quantity | float | int | Decimal,
    *,
    ref_temp_c: float = 20.0,
) -> pint.Quantity | float | Decimal:
    """Invert Specific Gravity (20 °C / 20 °C) to Degrees Brix (0–85 °Bx) via ICUMSA SPS-4 table."""
    from core.units.exceptions import StandardsApplicabilityError
    from core.units.registry import ureg
    from core.units.thermodynamics import _converters__interpolation_match_scalar_type
    _converters_chemical_validate_20c_reference(ref_temp_c, "Brix (ICUMSA)")
    sg_dec, is_qty = _converters_chemical_extract_dimensionless_sg(sg, "sg")
    sg_flt = float(sg_dec)
    if sg_flt < _converters_chemical_ICUMSA_BRIX_TABLE.min_y or sg_flt > _converters_chemical_ICUMSA_BRIX_TABLE.max_y:
        raise StandardsApplicabilityError(
            f"ICUMSA Brix SG (20/20 °C) must be in "
            f"[{_converters_chemical_ICUMSA_BRIX_TABLE.min_y}, {_converters_chemical_ICUMSA_BRIX_TABLE.max_y}] (0–85 °Bx), got {sg_dec}."
        )

    bx_val = _converters_chemical_ICUMSA_BRIX_TABLE.invert(sg_flt)
    out = _converters__interpolation_match_scalar_type(bx_val, sg)
    return ureg.Quantity(out, ureg.degree_Brix) if is_qty else out


from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# High-Level Quantity Router for Chemical Strength & Hydrometer Readings
# ---------------------------------------------------------------------------


@authorized_conversion
def convert_chemical_strength(
    quantity: pint.Quantity,
    target_unit: str | pint.Unit,
    *,
    equivalence_factor: float | int | Decimal | None = None,
    modulus: Literal["us", "rational"] | float | int | Decimal = "us",
    ref_temp_c: float = 20.0,
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties: Any = None,
    covariance: Any = None,
    correlations: Any = None,
) -> Any:
    """Explicit router for Normality <-> Molarity and SG <-> Baumé / Brix / Plato Quantities."""
    from core.units.exceptions import MissingAssayError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.thermodynamics import _converters__interpolation_validate_real_scalar
    from core.units.precision import ensure_quantity_ieee754

    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"quantity must be a pint.Quantity, got {type(quantity).__name__}.")
    validate_registry_quantity(quantity)
    _converters__interpolation_validate_real_scalar(quantity.magnitude, "quantity")

    dst = parse_unit(target_unit)
    src_dim = quantity.dimensionality
    dst_dim = dst.dimensionality

    if src_dim == dst_dim:
        out_q = ensure_quantity_ieee754(quantity.to(dst))
        if not return_metadata:
            return out_q
        from core.units.provenance import ConversionClassification
        from core.units.provenance import build_conversion_result
        return build_conversion_result(
            out_q,
            quantity,
            method="convert_chemical_strength",
            classification=ConversionClassification.EXACT,
            uncertainties=uncertainties,
            covariance=covariance,
            correlations=correlations,
        )

    # Molarity <-> Normality
    if {src_dim, dst_dim} == {_converters_chemical_MOLARITY_DIM, _converters_chemical_NORMALITY_DIM}:
        if equivalence_factor is None:
            raise MissingAssayError("Converting between molarity and normality requires an explicit equivalence_factor > 0.")
        if src_dim == _converters_chemical_MOLARITY_DIM:
            res = molarity_to_normality(quantity, equivalence_factor)
            assert isinstance(res, pint.Quantity)
            out_q = ensure_quantity_ieee754(res.to(dst))
        else:
            res = normality_to_molarity(quantity, equivalence_factor, target_unit=dst)
            assert isinstance(res, pint.Quantity)
            out_q = ensure_quantity_ieee754(res)
        if not return_metadata:
            return out_q
        from core.units.provenance import ConversionClassification
        from core.units.reference_conditions import ReferenceConditions
        from core.units.provenance import build_conversion_result
        ref = ReferenceConditions(
            classification=ConversionClassification.COMPOSITION_DEPENDENT.value,
            parameters={"equivalence_factor": equivalence_factor},
        )
        return build_conversion_result(
            out_q,
            quantity,
            reference_conditions=ref,
            method="convert_chemical_strength:equivalence",
            estimated=False,
            uncertainties=uncertainties,
            covariance=covariance,
            correlations=correlations,
        )

    # Hydrometer scales <-> dimensionless SG
    hydrometer_dims = {_converters_chemical_BAUME_HEAVY_DIM, _converters_chemical_BAUME_LIGHT_DIM, _converters_chemical_PLATO_DIM, _converters_chemical_BRIX_DIM}
    std_name: str | None = None
    estimated = False
    if src_dim in hydrometer_dims and dst.dimensionless:
        if src_dim == _converters_chemical_BAUME_HEAVY_DIM:
            sg_q = baume_to_sg(quantity, scale="heavy", modulus=modulus)
            std_name = "NIST Baumé Heavy"
        elif src_dim == _converters_chemical_BAUME_LIGHT_DIM:
            sg_q = baume_to_sg(quantity, scale="light", modulus=modulus)
            std_name = "NIST Baumé Light"
        elif src_dim == _converters_chemical_PLATO_DIM:
            sg_q = plato_to_sg(quantity, ref_temp_c=ref_temp_c)
            std_name = "ASBC Plato"
            estimated = True
        else:
            sg_q = brix_to_sg(quantity, ref_temp_c=ref_temp_c)
            std_name = "ICUMSA SPS-4 Brix"
            estimated = True
        assert isinstance(sg_q, pint.Quantity)
        out_q = ensure_quantity_ieee754(sg_q.to(dst))
    elif quantity.dimensionless and dst_dim in hydrometer_dims:
        if dst_dim == _converters_chemical_BAUME_HEAVY_DIM:
            res_q = sg_to_baume(quantity, scale="heavy", modulus=modulus)
            std_name = "NIST Baumé Heavy"
        elif dst_dim == _converters_chemical_BAUME_LIGHT_DIM:
            res_q = sg_to_baume(quantity, scale="light", modulus=modulus)
            std_name = "NIST Baumé Light"
        elif dst_dim == _converters_chemical_PLATO_DIM:
            res_q = sg_to_plato(quantity, ref_temp_c=ref_temp_c)
            std_name = "ASBC Plato"
            estimated = True
        else:
            res_q = sg_to_brix(quantity, ref_temp_c=ref_temp_c)
            std_name = "ICUMSA SPS-4 Brix"
            estimated = True
        assert isinstance(res_q, pint.Quantity)
        out_q = ensure_quantity_ieee754(res_q.to(dst))
    else:
        raise pint.DimensionalityError(quantity.units, dst, src_dim, dst_dim)

    if not return_metadata:
        return out_q
    from core.units.provenance import ConversionClassification
    from core.units.reference_conditions import ReferenceConditions
    from core.units.provenance import build_conversion_result
    ref = ReferenceConditions(
        temperature=ureg.Quantity(ref_temp_c, "degC"),
        classification=ConversionClassification.STATE_AND_MODEL_DEPENDENT.value if estimated else ConversionClassification.ASSAY_DEPENDENT.value,
        parameters={"modulus": str(modulus), "ref_temp_c": ref_temp_c},
    )
    return build_conversion_result(
        out_q,
        quantity,
        reference_conditions=ref,
        method="convert_chemical_strength:hydrometer",
        standard=std_name,
        estimated=estimated,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
    )


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (ChemicalSolutionAssay, SolventAssay, SulphurProductAssay, ChlorAlkaliAssay,))
__all__ = ('ChemicalAssay', 'ChemicalSolutionAssay', 'ChlorAlkaliAssay', 'SolventAssay', 'SulphurProductAssay', 'baume_heavy_to_sg', 'baume_light_to_sg', 'baume_to_sg', 'brix_to_sg', 'convert_chemical_solution', 'convert_chemical_strength', 'convert_chlor_alkali', 'convert_solvent', 'convert_sulphur_product', 'create_chemical_solution_context', 'create_chlor_alkali_context', 'create_solvent_context', 'create_sulphur_product_context', 'molarity_to_normality', 'normality_to_molarity', 'plato_to_sg', 'sg_to_baume', 'sg_to_baume_heavy', 'sg_to_baume_light', 'sg_to_brix', 'sg_to_plato')
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
    from core.units.thermodynamics import MonotonicCubicTable
    from core.units.assays.industrial.chemicals import convert_chemical_solution
    from core.units.assays.industrial.chemicals import convert_chlor_alkali
    from core.units.assays.industrial.chemicals import convert_solvent
    from core.units.assays.industrial.chemicals import convert_sulphur_product
    from core.units.assays.industrial.chemicals import create_chemical_solution_context
    from core.units.assays.industrial.chemicals import create_chlor_alkali_context
    from core.units.assays.industrial.chemicals import create_solvent_context
    from core.units.assays.industrial.chemicals import create_sulphur_product_context
