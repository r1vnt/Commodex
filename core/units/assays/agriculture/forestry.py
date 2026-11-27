"""Commodex assays agriculture forestry implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, fields
from dataclasses import dataclass, replace
from decimal import Decimal
from types import MappingProxyType
from typing import Any, Final
from typing import Any, Final, Literal, Self
from typing import ClassVar
import math
import pint

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class TimberAssay(_agricultural_families_AgriculturalFamilyAssay):
    """Dry-basis wood MC is water/dry wood; e.g. 1.25 means 125% MC.

    Convention: https://www.fpl.fs.usda.gov/documnts/usda/ah188/chapter01.pdf
    Timber recovery is a volume yield on the declared input/output volume bases.
    """
    moisture_content_dry_basis: Decimal | None = None
    lumber_recovery_fraction: Decimal | None = None
    _commodities = frozenset({"roundwood", "sawlogs", "hardwood", "softwood"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        if self.moisture_content_dry_basis is not None:
            ratio = _converters_petroleum_number(self.moisture_content_dry_basis, "dry-basis wood moisture", minimum=0)
            wet_fraction = ratio / (1 + ratio)
            if self.moisture_fraction is not None and abs(_converters_petroleum_fraction(self.moisture_fraction, "wet-basis moisture") - wet_fraction) > Decimal("1e-20"):
                raise ValueError("Dry- and wet-basis timber moisture declarations conflict.")
            object.__setattr__(self, "moisture_content_dry_basis", ratio)
            object.__setattr__(self, "moisture_fraction", wet_fraction)
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "lumber_recovery_fraction", _agricultural_families_yield(self.lumber_recovery_fraction, "lumber volume recovery"))

    def wet_to_dry_mass(self, quantity):
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        from core.units.assays.base import mass
        from core.units.assays.base import require_mass_basis
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        if self.moisture_content_dry_basis is None:
            return _agricultural_families_AgriculturalFamilyAssay.wet_to_dry_mass(self, quantity)
        require_mass_basis(quantity, "wet")
        return _converters_petroleum_output(mass(quantity) / (1 + self.moisture_content_dry_basis), "kg", quantity)

    def dry_to_wet_mass(self, quantity):
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        from core.units.assays.base import require_mass_basis
        from core.units.assays.base import scale_mass
        if self.moisture_content_dry_basis is None:
            return _agricultural_families_AgriculturalFamilyAssay.dry_to_wet_mass(self, quantity)
        require_mass_basis(quantity, "dry")
        return scale_mass(quantity, 1 + self.moisture_content_dry_basis)

    def sawn_volume(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if self.lumber_recovery_fraction is None:
            raise ValueError("Declare lumber volume recovery.")
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "roundwood volume", minimum=0)
        return _converters_petroleum_output(volume * self.lumber_recovery_fraction, "m^3", quantity)

    def roundwood_required(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if self.lumber_recovery_fraction is None or self.lumber_recovery_fraction == 0:
            raise ValueError("Positive lumber volume recovery is required.")
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "sawn volume", minimum=0)
        return _converters_petroleum_output(volume / self.lumber_recovery_fraction, "m^3", quantity)

    def solid_volume(self, quantity, *, solid_volume_fraction):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        fraction = _converters_petroleum_fraction(solid_volume_fraction, "solid volume fraction")
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "stacked volume", minimum=0)
        return _converters_petroleum_output(volume * fraction, "m^3", quantity)

    def stacked_volume(self, quantity, *, solid_volume_fraction):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        fraction = _converters_petroleum_fraction(solid_volume_fraction, "solid volume fraction")
        if fraction == 0:
            raise ValueError("Positive solid volume fraction is required.")
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "solid volume", minimum=0)
        return _converters_petroleum_output(volume / fraction, "m^3", quantity)


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class PaperAssay(_agricultural_families_AgriculturalFamilyAssay):
    basis_weight: pint.Quantity | None = None
    pulping_recovery_fraction: Decimal | None = None
    _commodities = frozenset({"recovered-paper", "occ", "newsprint", "containerboard"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.registry import ureg
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        if self.basis_weight is not None:
            object.__setattr__(self, "basis_weight", ureg.Quantity(_converters_petroleum_positive_quantity(self.basis_weight, "kg/m^2", "paper basis weight"), "kg/m^2"))
        object.__setattr__(self, "pulping_recovery_fraction", _agricultural_families_yield(self.pulping_recovery_fraction, "pulping recovery"))

    def area_to_mass(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if self.basis_weight is None:
            raise ValueError("Declare paper basis weight.")
        area = _converters_petroleum_quantity_number(quantity, "m^2", "paper area", minimum=0)
        return _converters_petroleum_output(area * self.basis_weight.magnitude, "kg", quantity)

    def mass_to_area(self, quantity):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        if self.basis_weight is None:
            raise ValueError("Declare paper basis weight.")
        return _converters_petroleum_output(mass(quantity) / self.basis_weight.magnitude, "m^2", quantity)

    def as_context(self, *, name=None):
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        context = _agricultural_families_AgriculturalFamilyAssay.as_context(self, name=name)
        if self.basis_weight is not None:
            context.add_transformation("[area]", "[mass]", lambda _r, q: self.area_to_mass(q))
            context.add_transformation("[mass]", "[area]", lambda _r, q: self.mass_to_area(q))
        return context

    def roll_mass(self, quantity, *, width):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        length = _converters_petroleum_quantity_number(quantity, "m", "paper length", minimum=0)
        width_m = _converters_petroleum_positive_quantity(width, "m", "paper width")
        return self.area_to_mass(_converters_petroleum_output(length * width_m, "m^2", quantity))

    def recovered_fibre_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.pulping_recovery_fraction, "pulping recovery")

    def recovered_paper_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.pulping_recovery_fraction, "pulping recovery")


from core.units.analytical_basis import AgriculturalQualityBase
from core.units.validation import authorized_conversion


# ---------------------------------------------------------------------------
# 14. Sawn Timber / Dimensional Lumber Quality Specification
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class TimberQualitySpec(AgriculturalQualityBase):
    """Physical specification for sawn timber species, moisture basis, and nominal vs actual geometry.

    Forestry and ASTM D4442 express timber moisture on an **oven-dry basis**
    (``MC_od = 100 * (m_wet - m_od) / m_od``), which can exceed 100% in green sapwood.
    ``moisture_basis='wet'`` is also accepted and converted seamlessly.
    """

    species: str
    moisture_content_pct: Decimal
    nominal_thickness_mm: Decimal
    nominal_width_mm: Decimal
    actual_thickness_mm: Decimal
    actual_width_mm: Decimal
    moisture_basis: Literal["oven_dry", "wet"] = "oven_dry"
    oven_dry_density_kg_per_m3: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.analytical_basis import _agricultural_quality_clean_str
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_pct
        sp = _agricultural_quality_clean_str(self.species, "species").upper()
        if self.moisture_basis == "oven_dry":
            mc = _agricultural_quality_pct(self.moisture_content_pct, "moisture_content_pct", low=Decimal(0), high=Decimal(300))
        elif self.moisture_basis == "wet":
            mc = _agricultural_quality_pct(self.moisture_content_pct, "moisture_content_pct", low=Decimal(0), high=Decimal(100), allow_high=False)
        else:
            raise ValueError("moisture_basis must be 'oven_dry' or 'wet'.")

        nt = _agricultural_quality_nonneg(self.nominal_thickness_mm, "nominal_thickness_mm", strict_positive=True)
        nw = _agricultural_quality_nonneg(self.nominal_width_mm, "nominal_width_mm", strict_positive=True)
        at = _agricultural_quality_nonneg(self.actual_thickness_mm, "actual_thickness_mm", strict_positive=True)
        aw = _agricultural_quality_nonneg(self.actual_width_mm, "actual_width_mm", strict_positive=True)
        if at > nt * Decimal("1.05") or aw > nw * Decimal("1.05"):
            raise ValueError("Actual surfaced dimensions cannot exceed nominal sawn dimensions (+5% green oversize).")

        rho_od: Decimal | None = None
        if self.oven_dry_density_kg_per_m3 is not None:
            rho_od = _agricultural_quality_pct(self.oven_dry_density_kg_per_m3, "oven_dry_density_kg_per_m3", low=Decimal(100), high=Decimal(1400))

        object.__setattr__(self, "species", sp)
        object.__setattr__(self, "moisture_content_pct", mc)
        object.__setattr__(self, "nominal_thickness_mm", nt)
        object.__setattr__(self, "nominal_width_mm", nw)
        object.__setattr__(self, "actual_thickness_mm", at)
        object.__setattr__(self, "actual_width_mm", aw)
        object.__setattr__(self, "oven_dry_density_kg_per_m3", rho_od)

    @classmethod
    def from_north_american_nominal(
        cls,
        nominal_thickness_in: int | float | Decimal,
        nominal_width_in: int | float | Decimal,
        *,
        species: str = "SPF",
        moisture_content_pct: Decimal | float = Decimal("19.0"),
        surfaced: Literal["S4S_DRY", "S4S_GREEN", "ROUGH"] = "S4S_DRY",
        oven_dry_density_kg_per_m3: Decimal | float | None = Decimal("420.0"),
    ) -> TimberQualitySpec:
        """Construct a DOC PS 20 North American softwood lumber spec (e.g. 2x4, 2x6, 2x10)."""
        from core.units.analytical_basis import _agricultural_quality_PS20_DIMENSION_TABLE_IN
        from core.units.analytical_basis import _agricultural_quality_dec
        from core.units.analytical_basis import _agricultural_quality_nonneg
        t_nom = _agricultural_quality_nonneg(nominal_thickness_in, "nominal_thickness_in", strict_positive=True)
        w_nom = _agricultural_quality_nonneg(nominal_width_in, "nominal_width_in", strict_positive=True)

        def _actual_in(nom_in: Decimal) -> Decimal:
            from core.units.analytical_basis import _agricultural_quality_PS20_DIMENSION_TABLE_IN
            if surfaced == "ROUGH":
                return nom_in
            if nom_in not in _agricultural_quality_PS20_DIMENSION_TABLE_IN:
                raise ValueError(f"Nominal dimension {nom_in}\" is not in standard DOC PS 20 table {sorted(_agricultural_quality_PS20_DIMENSION_TABLE_IN)}.")
            dry_in, green_in = _agricultural_quality_PS20_DIMENSION_TABLE_IN[nom_in]
            return dry_in if surfaced == "S4S_DRY" else green_in

        mm_per_in = Decimal("25.4")
        return cls(
            species=species,
            moisture_content_pct=_agricultural_quality_dec(moisture_content_pct, "moisture_content_pct"),
            nominal_thickness_mm=t_nom * mm_per_in,
            nominal_width_mm=w_nom * mm_per_in,
            actual_thickness_mm=_actual_in(t_nom) * mm_per_in,
            actual_width_mm=_actual_in(w_nom) * mm_per_in,
            moisture_basis="oven_dry",
            oven_dry_density_kg_per_m3=(
                _agricultural_quality_dec(oven_dry_density_kg_per_m3, "oven_dry_density_kg_per_m3")
                if oven_dry_density_kg_per_m3 is not None
                else None
            ),
        )

    @property
    def oven_dry_moisture_pct(self) -> Decimal:
        if self.moisture_basis == "oven_dry":
            return self.moisture_content_pct
        return (self.moisture_content_pct / (Decimal(100) - self.moisture_content_pct)) * Decimal(100)

    @property
    def wet_basis_moisture_pct(self) -> Decimal:
        if self.moisture_basis == "wet":
            return self.moisture_content_pct
        od = self.moisture_content_pct
        return (od / (Decimal(100) + od)) * Decimal(100)

    @property
    def actual_to_nominal_area_ratio(self) -> Decimal:
        """Ratio of actual surfaced cross-section to nominal board-foot cross-section (e.g. 0.65625 for 2x4 S4S)."""
        return (self.actual_thickness_mm * self.actual_width_mm) / (self.nominal_thickness_mm * self.nominal_width_mm)

    @authorized_conversion
    def nominal_mbf_to_actual_volume(
        self,
        mbf: Decimal | float,
        *,
        target_unit: str | Unit = "meter ** 3",
    ) -> pint.Quantity:
        """Convert nominal Thousand Board Feet (MBF) to true physical wood volume at actual dimensions."""
        from core.units.analytical_basis import _agricultural_quality_CUBIC_METRES_PER_MBF
        from core.units.analytical_basis import _agricultural_quality_dec
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_preserve_type
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        count = _agricultural_quality_nonneg(mbf, "mbf")
        actual_m3 = count * _agricultural_quality_CUBIC_METRES_PER_MBF * self.actual_to_nominal_area_ratio
        dst = parse_unit(target_unit)
        q = ureg.Quantity(actual_m3, "meter ** 3").to(dst)
        return ureg.Quantity(_agricultural_quality_preserve_type(_agricultural_quality_dec(q.magnitude, "volume"), mbf), dst)

    @authorized_conversion
    def actual_volume_to_nominal_mbf(self, actual_volume: pint.Quantity) -> Decimal | float:
        """Inverse of ``nominal_mbf_to_actual_volume``: nominal MBF from true physical wood volume."""
        from core.units.analytical_basis import _agricultural_quality_CUBIC_METRES_PER_MBF
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_preserve_type
        from core.units.registry import validate_registry_quantity
        if not isinstance(actual_volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(actual_volume).__name__}.")
        validate_registry_quantity(actual_volume)
        m3 = _agricultural_quality_nonneg(actual_volume.to("meter ** 3").magnitude, "actual_volume")
        mbf = m3 / (_agricultural_quality_CUBIC_METRES_PER_MBF * self.actual_to_nominal_area_ratio)
        return _agricultural_quality_preserve_type(mbf, actual_volume)

    @authorized_conversion
    def piece_nominal_board_feet(self, length: pint.Quantity) -> Decimal | float:
        """Nominal board feet ``(T_nom_in * W_nom_in * L_ft) / 12`` for a single piece."""
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_preserve_type
        from core.units.registry import validate_registry_quantity
        if not isinstance(length, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(length).__name__}.")
        validate_registry_quantity(length)
        len_ft = _agricultural_quality_nonneg(length.to("foot").magnitude, "length")
        t_in = self.nominal_thickness_mm / Decimal("25.4")
        w_in = self.nominal_width_mm / Decimal("25.4")
        bf = (t_in * w_in * len_ft) / Decimal(12)
        return _agricultural_quality_preserve_type(bf, length)

    @authorized_conversion
    def piece_actual_volume(self, length: pint.Quantity, *, target_unit: str | Unit = "meter ** 3") -> pint.Quantity:
        from core.units.analytical_basis import _agricultural_quality_dec
        from core.units.analytical_basis import _agricultural_quality_nonneg
        from core.units.analytical_basis import _agricultural_quality_preserve_type
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        if not isinstance(length, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(length).__name__}.")
        validate_registry_quantity(length)
        len_m = _agricultural_quality_nonneg(length.to("meter").magnitude, "length")
        area_m2 = (self.actual_thickness_mm / Decimal(1000)) * (self.actual_width_mm / Decimal(1000))
        dst = parse_unit(target_unit)
        q = ureg.Quantity(area_m2 * len_m, "meter ** 3").to(dst)
        return ureg.Quantity(_agricultural_quality_preserve_type(_agricultural_quality_dec(q.magnitude, "volume"), length), dst)

    @property
    def density_at_moisture_kg_per_m3(self) -> Decimal | None:
        """Total mass per unit actual volume ``rho_od * (1 + MC_od / 100)``."""
        if self.oven_dry_density_kg_per_m3 is None:
            return None
        return self.oven_dry_density_kg_per_m3 * (Decimal(1) + self.oven_dry_moisture_pct / Decimal(100))

    def to_material_spec(self, name: str | None = None) -> MaterialSpec:
        from core.units.assays.base import MaterialSpec
        from core.units.composition import MaterialProperties
        from core.units.processing import ProcessingSpecification
        from core.units.registry import ureg
        moist_wet_f = self.wet_basis_moisture_pct / Decimal(100)
        rho = self.density_at_moisture_kg_per_m3
        return MaterialSpec(
            name=name or f"{self.species} Timber",
            properties=MaterialProperties(
                bulk_density=ureg.Quantity(rho, "kg / meter ** 3") if rho is not None else None,
                moisture=moist_wet_f,
            ),
            processing=ProcessingSpecification(
                input_material="rough_nominal_timber",
                output_yields={"planed_dressed_timber": self.actual_to_nominal_area_ratio},
                process_losses={"planer_shavings_kerf": max(Decimal(0), Decimal(1) - self.actual_to_nominal_area_ratio)},
            ),
        )


from core.units.validation import validated_assay


@validated_assay
@dataclass(frozen=True, slots=True)
class PulpLumberAssay:
    """Immutable physical, moisture, and volumetric specification for pulp and lumber.

    Bridges pulp moisture bases (ADMT at 90% fibre, BDMT at 100% fibre, or custom
    moisture_pct), sawn timber volume (MBF <-> m^3 with optional dressed recovery),
    bulk stowage density, unit bale/bundle counts, and price parity.
    """

    name: str
    product_type: str  # "pulp" or "lumber"
    moisture_pct: float = 10.0  # 10.0% = standard ADMT (90% bone-dry fibre); 0.0% = BDMT
    dressed_recovery_pct: float = 100.0  # Dressed/finished volume relative to nominal MBF
    bulk_density_kg_per_m3: float | None = None
    unit_package_kg: float | None = None  # e.g. pulp bale (~250 kg) or lumber pack

    def __post_init__(self) -> None:
        from core.units.agricultural_specs import _validate_finite_float
        if not isinstance(self.product_type, str):
            raise ValueError("product_type must be 'pulp' or 'lumber'.")
        prod = self.product_type.strip().lower()
        if prod not in {"pulp", "lumber"}:
            raise ValueError(f"product_type must be 'pulp' or 'lumber', got {self.product_type!r}")

        moist = _validate_finite_float(self.moisture_pct, "moisture_pct")
        if not (0.0 <= moist < 100.0):
            raise ValueError(f"moisture_pct must be in [0, 100), got {self.moisture_pct}")

        dressed_rec = _validate_finite_float(self.dressed_recovery_pct, "dressed_recovery_pct")
        if not (0.0 < dressed_rec <= 100.0):
            raise ValueError(f"dressed_recovery_pct must be in (0, 100], got {self.dressed_recovery_pct}")

        rho: float | None = None
        if self.bulk_density_kg_per_m3 is not None:
            rho = _validate_finite_float(self.bulk_density_kg_per_m3, "bulk_density_kg_per_m3")
            if not (100.0 <= rho <= 2000.0):
                raise ValueError(f"bulk_density_kg_per_m3 must be in [100, 2000], got {self.bulk_density_kg_per_m3}")

        pkg_kg: float | None = None
        if self.unit_package_kg is not None:
            pkg_kg = _validate_finite_float(self.unit_package_kg, "unit_package_kg")
            if pkg_kg <= 0.0:
                raise ValueError(f"unit_package_kg must be strictly positive, got {self.unit_package_kg}")

        object.__setattr__(self, "product_type", prod)
        object.__setattr__(self, "moisture_pct", moist)
        object.__setattr__(self, "dressed_recovery_pct", dressed_rec)
        object.__setattr__(self, "bulk_density_kg_per_m3", rho)
        object.__setattr__(self, "unit_package_kg", pkg_kg)

    @property
    def bone_dry_fraction(self) -> float:
        """Bone-dry fibre fraction in [0, 1] (0.90 for standard 10% moisture ADMT)."""
        from core.units.analytical_basis import dry_matter_fraction
        return float(dry_matter_fraction(moisture_pct=self.moisture_pct))

    @property
    def admt_factor(self) -> float:
        """Multiplier converting mass at this assay's moisture to standard ADMT (90% dry / 10% moisture)."""
        from core.units.analytical_basis import moisture_level_ratio
        return float(moisture_level_ratio(source_moisture_pct=self.moisture_pct, target_moisture_pct=10.0))

    def as_context(self, *, context_name: str | None = None) -> pint.Context:
        from core.units.agricultural_specs import create_pulp_lumber_context
        return create_pulp_lumber_context(self, name=context_name)

    @staticmethod
    def _validate_mass(mass: pint.Quantity) -> pint.Quantity:
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        quantity = mass.to(ureg.metric_ton)
        magnitude = Decimal(str(quantity.magnitude))
        if not magnitude.is_finite() or magnitude < 0:
            raise ValueError("Pulp/lumber mass must be finite and non-negative.")
        return quantity

    @staticmethod
    def _validated_number(value: float | Decimal, label: str, *, nonnegative: bool = True) -> Decimal:
        number = Decimal(str(value))
        if not number.is_finite() or (nonnegative and number < 0):
            raise ValueError(f"{label} must be finite" + (" and non-negative." if nonnegative else "."))
        return number

    def _dry_scale(self, magnitude: Any, *, air_dry: bool = False) -> float | Decimal:
        fraction = 1 - Decimal(str(self.moisture_pct)) / 100
        if air_dry:
            fraction /= Decimal("0.9")
        return fraction if isinstance(magnitude, Decimal) else float(fraction)

    # --- Pulp Moisture Conversions (ADMT <-> BDMT <-> Actual Moisture) ---
    def mass_to_bdmt(self, mass: pint.Quantity) -> pint.Quantity:
        """Converts pulp/fibre mass (mt, kg, short_ton) to Bone-Dry Metric Tonnes (BDMT)."""
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        return wet_to_dry_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.bdmt)

    def bdmt_to_mass(self, bdmt_mass: pint.Quantity) -> pint.Quantity:
        """Converts Bone-Dry Metric Tonnes (BDMT) to mass at this assay's moisture_pct."""
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(bdmt_mass)
        return dry_to_wet_mass(mt_q, moisture_pct=self.moisture_pct, target_unit=ureg.metric_ton)

    def mass_to_admt(
        self,
        mass: pint.Quantity,
        *,
        source_moisture_pct: float | Decimal | None = None,
    ) -> pint.Quantity:
        """Converts mass at this assay's moisture_pct to Air-Dry Metric Tonnes (10% moisture / 90% dry)."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        src_m = self.moisture_pct if source_moisture_pct is None else source_moisture_pct
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=src_m,
            target_moisture_pct=Decimal("10.0"),
            target_unit=ureg.admt,
        )

    def admt_to_mass(self, admt_mass: pint.Quantity) -> pint.Quantity:
        """Converts Air-Dry Metric Tonnes (10% moisture / 90% dry) to mass at this assay's moisture_pct."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(admt_mass)
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=Decimal("10.0"),
            target_moisture_pct=self.moisture_pct,
            target_unit=ureg.metric_ton,
        )

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture_pct: float | Decimal,
        *,
        source_moisture_pct: float | Decimal | None = None,
        target_unit: str | pint.Unit | None = None,
    ) -> pint.Quantity:
        """Converts pulp/lumber mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.registry import ureg
        mt_q = self._validate_mass(mass)
        return convert_moisture_mass(
            mt_q,
            source_moisture_pct=self.moisture_pct if source_moisture_pct is None else source_moisture_pct,
            target_moisture_pct=target_moisture_pct,
            target_unit=target_unit or ureg.metric_ton,
        )

    # --- Lumber Volume Conversions (MBF <-> m^3) ---
    def mbf_to_cubic_metres(self, mbf: float | Decimal, *, apply_dressed_recovery: bool = False) -> pint.Quantity:
        """Converts nominal MBF count to m^3, optionally applying finished recovery."""
        from core.units.agricultural_specs import _CUBIC_METRES_PER_MBF
        from core.units.registry import ureg
        count = self._validated_number(mbf, "mbf")
        recovery = Decimal(str(self.dressed_recovery_pct)) / 100 if apply_dressed_recovery else Decimal(1)
        volume = count * Decimal(str(_CUBIC_METRES_PER_MBF)) * recovery
        magnitude = volume if isinstance(mbf, Decimal) else float(volume)
        return ureg.Quantity(magnitude, ureg.meter ** 3)

    def cubic_metres_to_mbf(self, volume: pint.Quantity, *, apply_dressed_recovery: bool = False) -> float | Decimal:
        """Returns nominal MBF count, reversing any requested dressed recovery."""
        from core.units.agricultural_specs import _CUBIC_METRES_PER_MBF
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        quantity = volume.to(ureg.meter ** 3)
        cubic_metres = self._validated_number(quantity.magnitude, "volume")
        recovery = Decimal(str(self.dressed_recovery_pct)) / 100 if apply_dressed_recovery else Decimal(1)
        count = cubic_metres / (Decimal(str(_CUBIC_METRES_PER_MBF)) * recovery)
        return count if isinstance(quantity.magnitude, Decimal) else float(count)

    # --- Bulk Stowage & Package Conversions ---
    def bulk_volume_to_mass(self, volume: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        m3_q = volume.to(ureg.meter ** 3)
        self._validated_number(m3_q.magnitude, "volume")
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(m3_q.magnitude, mt_per_m3)
        return ureg.Quantity(m3_q.magnitude * scale, ureg.metric_ton)

    def mass_to_bulk_volume(self, mass: pint.Quantity) -> pint.Quantity:
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if not isinstance(mass, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(mass).__name__}")
        if self.bulk_density_kg_per_m3 is None:
            raise ValueError(f"{self.name!r} does not define bulk_density_kg_per_m3.")
        mt_q = self._validate_mass(mass)
        mt_per_m3 = self.bulk_density_kg_per_m3 / 1000.0
        scale = _scale_factor(mt_q.magnitude, mt_per_m3)
        return ureg.Quantity(mt_q.magnitude / scale, ureg.meter ** 3)

    def mass_to_packages(self, mass: pint.Quantity) -> float | Decimal:
        """Returns unrounded bale/bundle count, preserving Decimal magnitudes."""
        from core.units.conversions import _scale_factor
        from core.units.registry import ureg
        if self.unit_package_kg is None:
            raise ValueError(f"{self.name!r} does not define unit_package_kg.")
        kilograms = self._validate_mass(mass).to(ureg.kilogram)
        scale = _scale_factor(kilograms.magnitude, self.unit_package_kg)
        return kilograms.magnitude / scale

    def packages_to_mass(self, count: float | Decimal) -> pint.Quantity:
        from core.units.registry import ureg
        if self.unit_package_kg is None:
            raise ValueError(f"{self.name!r} does not define unit_package_kg.")
        number = self._validated_number(count, "count")
        kilograms = number * Decimal(str(self.unit_package_kg))
        magnitude = kilograms if isinstance(count, Decimal) else float(kilograms)
        return ureg.Quantity(magnitude, ureg.kilogram).to(ureg.metric_ton)

    @staticmethod
    def _price_parity(price: float | Decimal, factor: Decimal, *, divide: bool) -> float | Decimal:
        from core.units.precision import preserve_ieee754_scalar
        number = PulpLumberAssay._validated_number(price, "price", nonnegative=False)
        value = number / factor if divide else number * factor
        return preserve_ieee754_scalar(value, price)

    @staticmethod
    def admt_price_to_bdmt_price(price_per_admt: float | Decimal) -> float | Decimal:
        """Converts nominal $/ADMT to $/BDMT, rounded to four decimals."""
        return PulpLumberAssay._price_parity(price_per_admt, Decimal("0.9"), divide=True)

    @staticmethod
    def bdmt_price_to_admt_price(price_per_bdmt: float | Decimal) -> float | Decimal:
        """Converts nominal $/BDMT to $/ADMT, rounded to four decimals."""
        return PulpLumberAssay._price_parity(price_per_bdmt, Decimal("0.9"), divide=False)

    @staticmethod
    def mbf_price_to_m3_price(price_per_mbf: float | Decimal) -> float | Decimal:
        """Converts nominal $/MBF to $/m^3; no dressed-recovery adjustment."""
        from core.units.agricultural_specs import _CUBIC_METRES_PER_MBF
        return PulpLumberAssay._price_parity(price_per_mbf, Decimal(str(_CUBIC_METRES_PER_MBF)), divide=True)

    @staticmethod
    def m3_price_to_mbf_price(price_per_m3: float | Decimal) -> float | Decimal:
        """Converts $/m^3 to nominal $/MBF; no dressed-recovery adjustment."""
        from core.units.agricultural_specs import _CUBIC_METRES_PER_MBF
        return PulpLumberAssay._price_parity(price_per_m3, Decimal(str(_CUBIC_METRES_PER_MBF)), divide=False)

    # --- Standard Forestry Presets ---
    @classmethod
    def nbsk_softwood_pulp(cls) -> PulpLumberAssay:
        """Northern Bleached Softwood Kraft (NBSK) standard air-dry pulp (90% dry / 10% moisture)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "PulpLumberAssay.nbsk_softwood_pulp")

    @classmethod
    def bek_hardwood_pulp(cls) -> PulpLumberAssay:
        """Bleached Eucalyptus Kraft (BEK) standard air-dry pulp (90% dry / 10% moisture)."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "PulpLumberAssay.bek_hardwood_pulp")

    @classmethod
    def spf_dimensional_lumber(cls, *, dressed_recovery_pct: float = 92.0) -> PulpLumberAssay:
        """Spruce-Pine-Fir (SPF) dimensional lumber benchmark."""
        from core.units.standards import construct_reference
        return construct_reference(cls, "agriculture", "PulpLumberAssay.spf_dimensional_lumber", overrides={'dressed_recovery_pct': dressed_recovery_pct})


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (TimberAssay, PaperAssay,))
__all__ = ('PaperAssay', 'PulpLumberAssay', 'TimberAssay', 'TimberQualitySpec', 'convert_paper', 'convert_timber', 'create_paper_context', 'create_timber_context')
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
    from core.units.assays.agriculture.forestry import convert_paper
    from core.units.assays.agriculture.forestry import convert_timber
    from core.units.assays.agriculture.forestry import create_paper_context
    from core.units.assays.agriculture.forestry import create_timber_context
