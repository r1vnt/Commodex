"""Commodex composition implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Callable, Mapping
from collections.abc import Mapping
from dataclasses import dataclass, field
from dataclasses import dataclass, field, replace
from decimal import Decimal
from enum import Enum
from itertools import permutations
from pathlib import Path
from types import MappingProxyType
from typing import Any
from typing import Any, Final, Literal
import inspect
import json
import math
import pint
import re

from core.units.registry import ureg


# ---------------------------------------------------------------------------
# Strict Scalar & Quantity Validation Helpers
# ---------------------------------------------------------------------------

_materials_MASS_DIM: Final = ureg.kilogram.dimensionality


from core.units.registry import ureg
_materials_VOLUME_DIM: Final = (ureg.meter ** 3).dimensionality


from core.units.registry import ureg
_materials_SUBSTANCE_DIM: Final = ureg.mole.dimensionality


from core.units.registry import ureg
_materials_DENSITY_DIM: Final = (ureg.kilogram / (ureg.meter ** 3)).dimensionality


from core.units.registry import ureg
_materials_MOLAR_MASS_DIM: Final = (ureg.kilogram / ureg.mole).dimensionality


from core.units.registry import ureg
_materials_DYNAMIC_VISC_DIM: Final = ureg.centipoise.dimensionality


from core.units.registry import ureg
_materials_KINEMATIC_VISC_DIM: Final = ureg.centistokes.dimensionality


from core.units.registry import ureg

_materials_MASS_FAMILY_DIMS: Final = frozenset({
    ureg.kilogram.dimensionality,
    ureg.wmt.dimensionality,
    ureg.dmt.dimensionality,
    ureg.contained_metal_tonne.dimensionality,
    ureg.payable_metal_tonne.dimensionality,
})





def _materials_to_decimal(value: Any, label: str) -> Decimal:
    if isinstance(value, bool) or not isinstance(value, (int, float, Decimal)):
        raise TypeError(f"{label} must be a real numeric scalar (int, float, or Decimal), got {type(value).__name__}.")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{label} must be finite, got {value!r}.")
    dec = value if isinstance(value, Decimal) else Decimal(str(value))
    if not dec.is_finite():
        raise ValueError(f"{label} must be finite, got {value!r}.")
    return dec





def _materials_validate_fraction(
    value: Any,
    label: str,
    *,
    allow_zero: bool = True,
    allow_one: bool = True,
) -> Decimal:
    dec = _materials_to_decimal(value, label)
    low_bad = dec < Decimal(0) if allow_zero else dec <= Decimal(0)
    high_bad = dec > Decimal(1) if allow_one else dec >= Decimal(1)
    if low_bad or high_bad:
        left = "[" if allow_zero else "("
        right = "]" if allow_one else ")"
        raise ValueError(f"{label} must lie in {left}0, 1{right}, got {value!r}.")
    return dec





def _materials_clean_name(value: str, label: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{label} must be a non-empty string.")
    return value.strip()





def _materials_validate_positive_quantity(
    qty: pint.Quantity,
    expected_dims: tuple[Any, ...],
    label: str,
) -> pint.Quantity:
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(qty, pint.Quantity):
        raise TypeError(f"{label} must be a pint.Quantity, got {type(qty).__name__}.")
    validate_registry_quantity(qty)
    mag = _materials_to_decimal(qty.magnitude, label)
    if mag <= Decimal(0):
        raise ValueError(f"{label} must be strictly positive, got {qty}.")
    if qty.dimensionality not in expected_dims:
        raise pint.DimensionalityError(
            qty.units,
            ureg.dimensionless,
            qty.dimensionality,
            expected_dims[0],
            extra_msg=f"; invalid dimensionality for {label}",
        )
    return ureg.Quantity(mag, qty.units)





def _materials_preserve_numeric_type(decimal_val: Decimal, *sources: Any) -> Decimal | float:
    from core.units.precision import pad_ieee754_decimal
    for src in sources:
        mag = src.magnitude if isinstance(src, pint.Quantity) else src
        if isinstance(mag, Decimal):
            return pad_ieee754_decimal(decimal_val)
    for src in sources:
        mag = src.magnitude if isinstance(src, pint.Quantity) else src
        if isinstance(mag, float):
            return float(decimal_val)
    return pad_ieee754_decimal(decimal_val)


from core.units.validation import authorized_conversion


@authorized_conversion
def _materials_mass_in_kg(quantity: pint.Quantity, label: str = "mass") -> tuple[Decimal, set[str]]:
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError(f"{label} must be a pint.Quantity, got {type(quantity).__name__}.")
    validate_registry_quantity(quantity)
    if quantity.dimensionality not in _materials_MASS_FAMILY_DIMS:
        raise pint.DimensionalityError(
            quantity.units,
            ureg.kilogram,
            quantity.dimensionality,
            _materials_MASS_DIM,
        )
    mag = _materials_to_decimal(quantity.magnitude, label)
    if mag < Decimal(0):
        raise ValueError(f"{label} must be non-negative, got {quantity}.")
    tags = set(ureg._semantic_tags(quantity.units._units))
    kg_val = ureg._convert_mass_scale(mag, quantity.units._units, ureg.kilogram._units)
    return _materials_to_decimal(kg_val, label), tags


from core.units.validation import authorized_conversion


@authorized_conversion
def _materials_kg_to_target(kg_decimal: Decimal, target_unit: Unit, *sources: Any) -> pint.Quantity:
    from core.units.registry import ureg
    if target_unit.dimensionality not in _materials_MASS_FAMILY_DIMS:
        raise pint.DimensionalityError(
            ureg.kilogram,
            target_unit,
            _materials_MASS_DIM,
            target_unit.dimensionality,
        )
    scaled = ureg._convert_mass_scale(kg_decimal, ureg.kilogram._units, target_unit._units)
    out_mag = _materials_preserve_numeric_type(_materials_to_decimal(scaled, "converted mass"), *sources)
    return ureg.Quantity(out_mag, target_unit)





# ---------------------------------------------------------------------------
# 1. CompositionBasis & MaterialComposition
# ---------------------------------------------------------------------------


class CompositionBasis(str, Enum):
    """Analytical basis for component, ash, and impurity fractions."""

    AS_RECEIVED = "as_received"
    AIR_DRIED = "air_dried"
    DRY = "dry"
    DRY_ASH_FREE = "dry_ash_free"
    MOLAR = "molar"
    VOLUME = "volume"

    @classmethod
    def _missing_(cls, value: object) -> CompositionBasis | None:
        if isinstance(value, str):
            alias_map = {
                "ar": cls.AS_RECEIVED,
                "wet": cls.AS_RECEIVED,
                "wmt": cls.AS_RECEIVED,
                "gross": cls.AS_RECEIVED,
                "as_received": cls.AS_RECEIVED,
                "ad": cls.AIR_DRIED,
                "adb": cls.AIR_DRIED,
                "air_dry": cls.AIR_DRIED,
                "air_dried": cls.AIR_DRIED,
                "db": cls.DRY,
                "dmt": cls.DRY,
                "dry": cls.DRY,
                "bone_dry": cls.DRY,
                "daf": cls.DRY_ASH_FREE,
                "dry_ash_free": cls.DRY_ASH_FREE,
                "mol": cls.MOLAR,
                "mole": cls.MOLAR,
                "molar": cls.MOLAR,
                "vol": cls.VOLUME,
                "volumetric": cls.VOLUME,
                "volume": cls.VOLUME,
            }
            return alias_map.get(value.strip().lower())
        return None





_materials_GRAVIMETRIC_BASES: Final = frozenset({
    CompositionBasis.AS_RECEIVED,
    CompositionBasis.AIR_DRIED,
    CompositionBasis.DRY,
    CompositionBasis.DRY_ASH_FREE,
})





@dataclass(frozen=True, slots=True)
class MaterialComposition:
    """Immutable analytical composition with moisture/ash re-basing, contained mass, and blending.

    Conventions
    -----------
    - ``components`` and ``impurities`` are expressed on ``basis`` with disjoint keys.
    - ``moisture_fraction`` is the wet-sample free/total moisture fraction in ``[0, 1)``.
      When ``basis`` is ``DRY`` or ``DRY_ASH_FREE``, ``moisture_fraction`` may still be
      stored as sample metadata to enable ``DRY <-> AS_RECEIVED`` mass and basis conversions.
    - ``ash_fraction`` is expressed on ``basis`` for ``AS_RECEIVED`` / ``AIR_DRIED`` / ``DRY``,
      and on a ``DRY`` basis when ``basis == DRY_ASH_FREE``.
    """

    components: Mapping[str, Decimal] = field(default_factory=dict)
    basis: CompositionBasis = CompositionBasis.DRY
    moisture_fraction: Decimal | None = None
    ash_fraction: Decimal | None = None
    impurities: Mapping[str, Decimal] = field(default_factory=dict)

    def __post_init__(self) -> None:
        from core.units.exceptions import InvalidPhysicalStateError
        try:
            resolved_basis = CompositionBasis(self.basis)
        except ValueError as error:
            raise ValueError(f"Unsupported composition basis: {self.basis!r}.") from error

        if not isinstance(self.components, Mapping) or not isinstance(self.impurities, Mapping):
            raise TypeError("components and impurities must be mappings.")

        clean_comp: dict[str, Decimal] = {}
        for raw_key, raw_val in self.components.items():
            key = _materials_clean_name(raw_key, "component name")
            if key in clean_comp:
                raise ValueError(f"Duplicate component key after trimming: {key!r}.")
            clean_comp[key] = _materials_validate_fraction(raw_val, f"component {key!r}")

        clean_imp: dict[str, Decimal] = {}
        for raw_key, raw_val in self.impurities.items():
            key = _materials_clean_name(raw_key, "impurity name")
            if key in clean_imp:
                raise ValueError(f"Duplicate impurity key after trimming: {key!r}.")
            if key in clean_comp:
                raise ValueError(f"Identifier {key!r} cannot appear in both components and impurities.")
            clean_imp[key] = _materials_validate_fraction(raw_val, f"impurity {key!r}")

        moist = (
            _materials_validate_fraction(self.moisture_fraction, "moisture_fraction", allow_one=False)
            if self.moisture_fraction is not None
            else None
        )
        ash = (
            _materials_validate_fraction(self.ash_fraction, "ash_fraction", allow_one=False)
            if self.ash_fraction is not None
            else None
        )

        # Closure check: fractions expressed on the same basis cannot exceed unity.
        comp_plus_imp = sum(clean_comp.values(), Decimal(0)) + sum(clean_imp.values(), Decimal(0))
        closure_total = comp_plus_imp
        if resolved_basis in {CompositionBasis.AS_RECEIVED, CompositionBasis.AIR_DRIED}:
            if moist is not None:
                closure_total += moist
            if ash is not None:
                closure_total += ash
        elif resolved_basis == CompositionBasis.DRY:
            if ash is not None:
                closure_total += ash

        # Allow 1e-12 slack only for binary float literals converted via str; reject real > 1 sums.
        if closure_total > Decimal("1.000000000001"):
            raise InvalidPhysicalStateError(
                f"Composition mass/mole/volume closure exceeded 1.0 on {resolved_basis.value} basis (got {closure_total})."
            )

        object.__setattr__(self, "components", MappingProxyType(clean_comp))
        object.__setattr__(self, "basis", resolved_basis)
        object.__setattr__(self, "moisture_fraction", moist)
        object.__setattr__(self, "ash_fraction", ash)
        object.__setattr__(self, "impurities", MappingProxyType(clean_imp))

    def __hash__(self) -> int:
        return hash((
            tuple(sorted(self.components.items())),
            self.basis,
            self.moisture_fraction,
            self.ash_fraction,
            tuple(sorted(self.impurities.items())),
        ))

    @property
    def dry_mass_fraction(self) -> Decimal:
        """Dry-matter mass fraction ``1 - moisture_fraction`` of the wet sample."""
        from core.units.analytical_basis import dry_matter_fraction
        if self.moisture_fraction is None:
            raise ValueError("Explicit moisture_fraction is required to compute wet/dry mass ratio.")
        return dry_matter_fraction(self.moisture_fraction)

    @property
    def total_component_fraction(self) -> Decimal:
        return sum(self.components.values(), Decimal(0))

    @property
    def total_impurity_fraction(self) -> Decimal:
        return sum(self.impurities.values(), Decimal(0))

    def _dry_ash_fraction(self, override_moisture: Decimal | None = None, override_ash: Decimal | None = None) -> Decimal:
        ash = override_ash if override_ash is not None else self.ash_fraction
        if ash is None:
            raise ValueError("ash_fraction is required for DRY_ASH_FREE basis conversions.")
        if self.basis in {CompositionBasis.DRY, CompositionBasis.DRY_ASH_FREE}:
            dry_ash = ash
        elif self.basis in {CompositionBasis.AS_RECEIVED, CompositionBasis.AIR_DRIED}:
            moist = override_moisture if override_moisture is not None else self.moisture_fraction
            if moist is None:
                raise ValueError(f"moisture_fraction is required to convert ash from {self.basis.value} to dry basis.")
            dry_ash = ash / (Decimal(1) - moist)
        else:
            raise ValueError(f"Ash re-basing is undefined for {self.basis.value} compositions.")
        if dry_ash >= Decimal(1):
            raise ValueError(f"Derived dry-basis ash fraction must be < 1.0, got {dry_ash}.")
        return dry_ash

    def to_basis(
        self,
        target_basis: CompositionBasis | str,
        *,
        target_moisture_fraction: Decimal | float | None = None,
        ash_fraction: Decimal | float | None = None,
    ) -> MaterialComposition:
        """Project gravimetric component, impurity, and ash fractions onto ``target_basis``."""
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        dst = CompositionBasis(target_basis)
        if self.basis not in _materials_GRAVIMETRIC_BASES or dst not in _materials_GRAVIMETRIC_BASES:
            if self.basis == dst:
                return self
            raise IncompatibleAnalyticalBasisError(
                f"Analytical basis projection requires gravimetric bases (AR, AD, DRY, DAF), got {self.basis.value} -> {dst.value}."
            )

        tgt_moist = (
            _materials_validate_fraction(target_moisture_fraction, "target_moisture_fraction", allow_one=False)
            if target_moisture_fraction is not None
            else self.moisture_fraction
        )
        ovr_ash = (
            _materials_validate_fraction(ash_fraction, "ash_fraction", allow_one=False)
            if ash_fraction is not None
            else None
        )

        if self.basis == dst and (
            dst in {CompositionBasis.DRY, CompositionBasis.DRY_ASH_FREE}
            or tgt_moist == self.moisture_fraction
        ) and ovr_ash is None:
            return self

        # Step 1: Convert source fractions to DRY basis
        if self.basis == CompositionBasis.DRY:
            to_dry = Decimal(1)
            dry_ash = ovr_ash if ovr_ash is not None else self.ash_fraction
        elif self.basis in {CompositionBasis.AS_RECEIVED, CompositionBasis.AIR_DRIED}:
            if self.moisture_fraction is None:
                raise ValueError(f"source moisture_fraction is required to convert from {self.basis.value} basis.")
            to_dry = Decimal(1) / (Decimal(1) - self.moisture_fraction)
            base_ash = ovr_ash if ovr_ash is not None else self.ash_fraction
            dry_ash = base_ash * to_dry if base_ash is not None else None
        else:  # DRY_ASH_FREE
            dry_ash = self._dry_ash_fraction(override_ash=ovr_ash)
            to_dry = Decimal(1) - dry_ash

        dry_comp = {k: v * to_dry for k, v in self.components.items()}
        dry_imp = {k: v * to_dry for k, v in self.impurities.items()}

        # Step 2: Project DRY fractions onto destination basis
        if dst == CompositionBasis.DRY:
            from_dry = Decimal(1)
            out_ash = dry_ash
            out_moist = tgt_moist
        elif dst in {CompositionBasis.AS_RECEIVED, CompositionBasis.AIR_DRIED}:
            if tgt_moist is None:
                raise ValueError(f"target_moisture_fraction (or stored moisture_fraction) is required for {dst.value} basis.")
            from_dry = Decimal(1) - tgt_moist
            out_ash = dry_ash * from_dry if dry_ash is not None else None
            out_moist = tgt_moist
        else:  # DRY_ASH_FREE
            if dry_ash is None:
                raise ValueError("ash_fraction is required to project onto DRY_ASH_FREE basis.")
            if dry_ash >= Decimal(1):
                raise ValueError("Dry-basis ash fraction must be < 1.0.")
            from_dry = Decimal(1) / (Decimal(1) - dry_ash)
            out_ash = dry_ash
            out_moist = tgt_moist

        return MaterialComposition(
            components={k: v * from_dry for k, v in dry_comp.items()},
            basis=dst,
            moisture_fraction=out_moist,
            ash_fraction=out_ash,
            impurities={k: v * from_dry for k, v in dry_imp.items()},
        )

    def component_fraction(
        self,
        component: str,
        *,
        on_basis: CompositionBasis | str | None = None,
    ) -> Decimal:
        """Return the fraction of ``component`` (from ``components`` or ``impurities``) on ``on_basis``."""
        comp = self if on_basis is None else self.to_basis(on_basis)
        key = _materials_clean_name(component, "component")
        if key in comp.components:
            return comp.components[key]
        if key in comp.impurities:
            return comp.impurities[key]
        raise KeyError(f"Component or impurity {key!r} is not defined in this MaterialComposition.")

    def wet_to_dry_mass(
        self,
        wet_mass: pint.Quantity,
        *,
        moisture_fraction: Decimal | float | None = None,
        target_unit: str | Unit = "dmt",
    ) -> pint.Quantity:
        """Convert wet/as-received mass to bone-dry mass."""
        from core.units.analytical_basis import wet_to_dry_mass
        moist = (
            _materials_validate_fraction(moisture_fraction, "moisture_fraction", allow_one=False)
            if moisture_fraction is not None
            else self.moisture_fraction
        )
        if moist is None:
            raise ValueError("Explicit moisture_fraction is required for wet_to_dry_mass.")
        return wet_to_dry_mass(wet_mass, moist, target_unit=target_unit)

    def dry_to_wet_mass(
        self,
        dry_mass: pint.Quantity,
        *,
        moisture_fraction: Decimal | float | None = None,
        target_unit: str | Unit = "wmt",
    ) -> pint.Quantity:
        """Convert bone-dry mass to wet/as-received mass."""
        from core.units.analytical_basis import dry_to_wet_mass
        moist = (
            _materials_validate_fraction(moisture_fraction, "moisture_fraction", allow_one=False)
            if moisture_fraction is not None
            else self.moisture_fraction
        )
        if moist is None:
            raise ValueError("Explicit moisture_fraction is required for dry_to_wet_mass.")
        return dry_to_wet_mass(dry_mass, moist, target_unit=target_unit)

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture_fraction: Decimal | float,
        *,
        source_moisture_fraction: Decimal | float | None = None,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Convert mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        w1 = (
            _materials_validate_fraction(source_moisture_fraction, "source_moisture_fraction", allow_one=False)
            if source_moisture_fraction is not None
            else self.moisture_fraction
        )
        if w1 is None:
            raise ValueError("Explicit moisture_fraction is required for mass_at_moisture.")
        return convert_moisture_mass(
            mass,
            source_moisture=w1,
            target_moisture=target_moisture_fraction,
            target_unit=target_unit,
        )

    def contained_mass(
        self,
        mass: pint.Quantity,
        component: str,
        *,
        input_basis: CompositionBasis | str | None = None,
        target_unit: str | Unit = "contained_metal_tonne",
    ) -> pint.Quantity:
        """Compute contained mass of ``component`` from a wet, dry, or basis-matched mass quantity."""
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        if self.basis not in _materials_GRAVIMETRIC_BASES:
            raise ValueError(f"contained_mass requires a gravimetric composition basis, not {self.basis.value}.")
        kg_in, tags = _materials_mass_in_kg(mass, "mass")
        if "contained" in tags or "payable" in tags:
            raise ValueError("contained_mass input must be bulk feedstock mass, not already-contained mass.")

        if "wet" in tags:
            inferred = CompositionBasis.AS_RECEIVED
        elif "dry" in tags:
            inferred = CompositionBasis.DRY
        else:
            inferred = None

        if input_basis is not None:
            eff_basis = CompositionBasis(input_basis)
            if inferred is not None and eff_basis != inferred:
                raise IncompatibleAnalyticalBasisError(f"Mass unit {mass.units} ({inferred.value}) contradicts input_basis={eff_basis.value!r}.")
        else:
            eff_basis = inferred if inferred is not None else self.basis

        frac = self.component_fraction(component, on_basis=eff_basis)
        dst = parse_unit(target_unit)
        dst_tags = set(ureg._semantic_tags(dst._units))
        if dst_tags and dst_tags != {"contained"}:
            raise ValueError(f"contained_mass target_unit must be contained-basis or unqualified mass, got {dst}.")
        return _materials_kg_to_target(kg_in * frac, dst, mass)

    def required_mass_for_contained(
        self,
        contained: pint.Quantity,
        component: str,
        *,
        output_basis: CompositionBasis | str | None = None,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Inverse of ``contained_mass``: required feedstock mass to deliver ``contained`` component mass."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        if self.basis not in _materials_GRAVIMETRIC_BASES:
            raise ValueError(f"required_mass_for_contained requires a gravimetric basis, not {self.basis.value}.")
        kg_cont, tags = _materials_mass_in_kg(contained, "contained")
        if tags and tags != {"contained"}:
            raise ValueError(f"required_mass_for_contained expects contained-basis or unqualified mass, got {contained.units}.")

        eff_basis = CompositionBasis(output_basis) if output_basis is not None else self.basis
        frac = self.component_fraction(component, on_basis=eff_basis)
        if frac <= Decimal(0):
            raise ValueError(f"Component {component!r} has zero fraction; required feedstock mass is undefined.")

        if target_unit is None:
            dst = ureg.wmt if eff_basis == CompositionBasis.AS_RECEIVED else ureg.dmt if eff_basis == CompositionBasis.DRY else ureg.metric_tonne
        else:
            dst = parse_unit(target_unit)
        return _materials_kg_to_target(kg_cont / frac, dst, contained)

    def blend_with(
        self,
        other: MaterialComposition,
        *,
        self_weight: Decimal | float,
        other_weight: Decimal | float | None = None,
        weight_basis: CompositionBasis | str | None = None,
    ) -> MaterialComposition:
        """Conserve mass/moles/volume across two MaterialCompositions.

        When ``other_weight`` is omitted, ``self_weight`` must lie in ``[0, 1]`` and
        ``other_weight = 1 - self_weight``. When blending two ``DRY``-basis compositions
        using wet lot weights (``weight_basis='as_received'``), dry-basis grades are
        weighted by each lot's bone-dry mass while ``moisture_fraction`` is weighted by
        wet mass, strictly conserving both water and contained elements.
        """
        if not isinstance(other, MaterialComposition):
            raise TypeError(f"Expected MaterialComposition, got {type(other).__name__}.")

        if other_weight is None:
            w1 = _materials_validate_fraction(self_weight, "self_weight")
            w2 = Decimal(1) - w1
        else:
            w1 = _materials_to_decimal(self_weight, "self_weight")
            w2 = _materials_to_decimal(other_weight, "other_weight")
            if w1 < Decimal(0) or w2 < Decimal(0) or (w1 + w2) <= Decimal(0):
                raise ValueError("Blend weights must be non-negative with positive sum.")

        total_w = w1 + w2
        f1, f2 = w1 / total_w, w2 / total_w

        if self.basis != other.basis:
            if self.basis in _materials_GRAVIMETRIC_BASES and other.basis in _materials_GRAVIMETRIC_BASES:
                other_aligned = other.to_basis(self.basis)
            else:
                raise ValueError(f"Cannot blend incompatible composition bases: {self.basis.value} vs {other.basis.value}.")
        else:
            other_aligned = other

        eff_weight_basis = CompositionBasis(weight_basis) if weight_basis is not None else self.basis
        if eff_weight_basis == CompositionBasis.AS_RECEIVED and self.basis == CompositionBasis.DRY:
            if self.moisture_fraction is None or other.moisture_fraction is None:
                raise ValueError("Wet-weight blending of DRY-basis compositions requires moisture_fraction on both lots.")
            dry1 = f1 * (Decimal(1) - self.moisture_fraction)
            dry2 = f2 * (Decimal(1) - other.moisture_fraction)
            total_dry = dry1 + dry2
            c_w1, c_w2 = dry1 / total_dry, dry2 / total_dry
            blended_moist: Decimal | None = (f1 * self.moisture_fraction) + (f2 * other.moisture_fraction)
        else:
            if eff_weight_basis != self.basis:
                raise ValueError(
                    f"weight_basis={eff_weight_basis.value!r} does not match composition basis {self.basis.value!r}."
                )
            c_w1, c_w2 = f1, f2
            if self.moisture_fraction is not None and other_aligned.moisture_fraction is not None:
                blended_moist = (f1 * self.moisture_fraction) + (f2 * other_aligned.moisture_fraction)
            else:
                blended_moist = None

        all_comp = sorted(set(self.components) | set(other_aligned.components))
        blended_comp = {
            k: (c_w1 * self.components.get(k, Decimal(0))) + (c_w2 * other_aligned.components.get(k, Decimal(0)))
            for k in all_comp
        }
        all_imp = sorted(set(self.impurities) | set(other_aligned.impurities))
        blended_imp = {
            k: (c_w1 * self.impurities.get(k, Decimal(0))) + (c_w2 * other_aligned.impurities.get(k, Decimal(0)))
            for k in all_imp
        }
        blended_ash: Decimal | None = None
        if self.ash_fraction is not None and other_aligned.ash_fraction is not None:
            blended_ash = (c_w1 * self.ash_fraction) + (c_w2 * other_aligned.ash_fraction)

        return MaterialComposition(
            components=blended_comp,
            basis=self.basis,
            moisture_fraction=blended_moist,
            ash_fraction=blended_ash,
            impurities=blended_imp,
        )





# ---------------------------------------------------------------------------
# 2. MaterialProperties
# ---------------------------------------------------------------------------


def _materials_refutas_vbn(cst: float) -> float:
    if not math.isfinite(cst) or cst <= 0.0:
        raise ValueError(f"Kinematic viscosity must be strictly positive for Refutas blending, got {cst}.")
    return 14.534 * math.log(math.log(cst + 0.8)) + 10.975





def _materials_cst_from_refutas_vbn(vbn: float) -> float:
    return math.exp(math.exp((vbn - 10.975) / 14.534)) - 0.8


from core.units.validation import _safety_QuantitySnapshot
from core.units.validation import authorized_conversion


@dataclass(frozen=True, slots=True)
class MaterialProperties(_safety_QuantitySnapshot):
    """Immutable physical property bundle with density, purity, moisture, molar, and price helpers."""

    density: Quantity | None = None
    bulk_density: Quantity | None = None
    molar_mass: Quantity | None = None
    viscosity: Quantity | None = None
    moisture: Decimal | None = None
    purity: Decimal | None = None
    porosity: Decimal | None = None

    def __post_init__(self) -> None:
        from core.units.exceptions import InvalidPhysicalStateError
        dens = (
            _materials_validate_positive_quantity(self.density, (_materials_DENSITY_DIM,), "density")
            if self.density is not None
            else None
        )
        bulk = (
            _materials_validate_positive_quantity(self.bulk_density, (_materials_DENSITY_DIM,), "bulk_density")
            if self.bulk_density is not None
            else None
        )
        mm = (
            _materials_validate_positive_quantity(self.molar_mass, (_materials_MOLAR_MASS_DIM,), "molar_mass")
            if self.molar_mass is not None
            else None
        )
        visc = (
            _materials_validate_positive_quantity(
                self.viscosity, (_materials_DYNAMIC_VISC_DIM, _materials_KINEMATIC_VISC_DIM), "viscosity"
            )
            if self.viscosity is not None
            else None
        )
        moist = (
            _materials_validate_fraction(self.moisture, "moisture", allow_one=False)
            if self.moisture is not None
            else None
        )
        pur = (
            _materials_validate_fraction(self.purity, "purity", allow_zero=False)
            if self.purity is not None
            else None
        )
        por = (
            _materials_validate_fraction(self.porosity, "porosity", allow_one=False)
            if self.porosity is not None
            else None
        )

        if dens is not None and bulk is not None:
            rho_true = dens.to("kg / meter ** 3").magnitude
            rho_bulk = bulk.to("kg / meter ** 3").magnitude
            if rho_bulk > rho_true * Decimal("1.000000000001"):
                raise InvalidPhysicalStateError(
                    f"bulk_density ({rho_bulk} kg/m^3) cannot exceed particle/true density ({rho_true} kg/m^3)."
                )
            if por is not None:
                expected_bulk = rho_true * (Decimal(1) - por)
                if abs(rho_bulk - expected_bulk) > rho_true * Decimal("1e-6"):
                    raise InvalidPhysicalStateError(
                        f"porosity={por} conflicts with bulk_density={rho_bulk} kg/m^3 and density={rho_true} kg/m^3."
                    )

        object.__setattr__(self, "density", dens)
        object.__setattr__(self, "bulk_density", bulk)
        object.__setattr__(self, "molar_mass", mm)
        object.__setattr__(self, "viscosity", visc)
        object.__setattr__(self, "moisture", moist)
        object.__setattr__(self, "purity", pur)
        object.__setattr__(self, "porosity", por)

    @property
    def effective_bulk_density(self) -> Quantity | None:
        """Bulk density if declared, or ``density * (1 - porosity)`` when derivable."""
        from core.units.registry import ureg
        if self.bulk_density is not None:
            return self.bulk_density
        if self.density is not None and self.porosity is not None:
            d = self.density
            return ureg.Quantity(d.magnitude * (Decimal(1) - self.porosity), d.units)
        return None

    @property
    def effective_true_density(self) -> Quantity | None:
        """True/particle density if declared, or ``bulk_density / (1 - porosity)`` when derivable."""
        from core.units.registry import ureg
        if self.density is not None:
            return self.density
        if self.bulk_density is not None and self.porosity is not None:
            b = self.bulk_density
            return ureg.Quantity(b.magnitude / (Decimal(1) - self.porosity), b.units)
        return None

    @property
    def implied_porosity(self) -> Decimal | None:
        """Declared void fraction ``porosity``, or ``1 - bulk_density / density`` when both exist."""
        if self.porosity is not None:
            return self.porosity
        if self.density is not None and self.bulk_density is not None:
            rho_true = self.density.to("kg / meter ** 3").magnitude
            rho_bulk = self.bulk_density.to("kg / meter ** 3").magnitude
            return max(Decimal(0), Decimal(1) - (rho_bulk / rho_true))
        return None

    def _select_density(self, *, use_bulk: bool = False) -> pint.Quantity:
        from core.units.exceptions import MissingAssayError
        if use_bulk:
            chosen = self.effective_bulk_density
            if chosen is None:
                raise MissingAssayError("bulk_density (or density + porosity) is required when use_bulk=True.")
            return chosen
        chosen = self.effective_true_density or self.effective_bulk_density
        if chosen is None:
            raise MissingAssayError("MaterialProperties requires density or bulk_density for volume/mass conversion.")
        return chosen

    def volume_to_mass(
        self,
        volume: pint.Quantity,
        *,
        target_unit: str | Unit = "kilogram",
        use_bulk: bool = False,
    ) -> pint.Quantity:
        """Convert volume to mass using ``density`` (or ``bulk_density`` when ``use_bulk=True``)."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        if not isinstance(volume, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(volume).__name__}.")
        validate_registry_quantity(volume)
        if volume.dimensionality != _materials_VOLUME_DIM:
            raise pint.DimensionalityError(volume.units, ureg.meter ** 3, volume.dimensionality, _materials_VOLUME_DIM)
        vol_dec = _materials_to_decimal(volume.magnitude, "volume")
        if vol_dec < Decimal(0):
            raise ValueError("Volume must be non-negative.")

        rho = self._select_density(use_bulk=use_bulk).to("kg / meter ** 3").magnitude
        vol_m3 = ureg.Quantity(vol_dec, volume.units).to("meter ** 3").magnitude
        dst = parse_unit(target_unit)
        return _materials_kg_to_target(vol_m3 * rho, dst, volume)

    def mass_to_volume(
        self,
        mass: pint.Quantity,
        *,
        target_unit: str | Unit = "meter ** 3",
        use_bulk: bool = False,
    ) -> pint.Quantity:
        """Convert mass to volume using ``density`` (or ``bulk_density`` when ``use_bulk=True``)."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        kg_val, _ = _materials_mass_in_kg(mass, "mass")
        rho = self._select_density(use_bulk=use_bulk).to("kg / meter ** 3").magnitude
        dst = parse_unit(target_unit)
        if dst.dimensionality != _materials_VOLUME_DIM:
            raise pint.DimensionalityError(ureg.meter ** 3, dst, _materials_VOLUME_DIM, dst.dimensionality)
        vol_q = ureg.Quantity(kg_val / rho, "meter ** 3").to(dst)
        return ureg.Quantity(_materials_preserve_numeric_type(vol_q.magnitude, mass), dst)

    def gross_to_pure(
        self,
        gross_quantity: pint.Quantity,
        *,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Scale a gross/commercial mass or volume by ``purity`` to active/fine content."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        if self.purity is None:
            raise ValueError("purity is required for gross_to_pure conversion.")
        if not isinstance(gross_quantity, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(gross_quantity).__name__}.")
        validate_registry_quantity(gross_quantity)
        mag = _materials_to_decimal(gross_quantity.magnitude, "gross_quantity")
        if mag < Decimal(0):
            raise ValueError("Quantity must be non-negative.")
        dst = parse_unit(target_unit) if target_unit is not None else gross_quantity.units
        scaled = ureg.Quantity(mag * self.purity, gross_quantity.units).to(dst)
        return ureg.Quantity(_materials_preserve_numeric_type(scaled.magnitude, gross_quantity), dst)

    def pure_to_gross(
        self,
        pure_quantity: pint.Quantity,
        *,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Inverse of ``gross_to_pure``: required gross/solution quantity to supply ``pure_quantity``."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        if self.purity is None:
            raise ValueError("purity is required for pure_to_gross conversion.")
        if not isinstance(pure_quantity, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(pure_quantity).__name__}.")
        validate_registry_quantity(pure_quantity)
        mag = _materials_to_decimal(pure_quantity.magnitude, "pure_quantity")
        if mag < Decimal(0):
            raise ValueError("Quantity must be non-negative.")
        dst = parse_unit(target_unit) if target_unit is not None else pure_quantity.units
        scaled = ureg.Quantity(mag / self.purity, pure_quantity.units).to(dst)
        return ureg.Quantity(_materials_preserve_numeric_type(scaled.magnitude, pure_quantity), dst)

    def wet_to_dry(
        self,
        wet_mass: pint.Quantity,
        *,
        target_unit: str | Unit = "dmt",
    ) -> pint.Quantity:
        """Convert wet mass to dry mass using ``moisture``."""
        from core.units.analytical_basis import wet_to_dry_mass
        if self.moisture is None:
            raise ValueError("moisture is required for wet_to_dry conversion.")
        return wet_to_dry_mass(wet_mass, self.moisture, target_unit=target_unit)

    def dry_to_wet(
        self,
        dry_mass: pint.Quantity,
        *,
        target_unit: str | Unit = "wmt",
    ) -> pint.Quantity:
        """Convert dry mass to wet mass using ``moisture``."""
        from core.units.analytical_basis import dry_to_wet_mass
        if self.moisture is None:
            raise ValueError("moisture is required for dry_to_wet conversion.")
        return dry_to_wet_mass(dry_mass, self.moisture, target_unit=target_unit)

    def mass_at_moisture(
        self,
        mass: pint.Quantity,
        target_moisture: Decimal | float,
        *,
        source_moisture: Decimal | float | None = None,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        """Convert mass between two moisture levels ``m_2 = m_1 * (1 - w_1) / (1 - w_2)``."""
        from core.units.analytical_basis import convert_moisture_mass
        w1 = self.moisture if source_moisture is None else _materials_validate_fraction(source_moisture, "source_moisture", allow_one=False)
        if w1 is None:
            raise ValueError("moisture is required for mass_at_moisture conversion.")
        return convert_moisture_mass(
            mass,
            source_moisture=w1,
            target_moisture=target_moisture,
            target_unit=target_unit,
        )

    def mass_to_moles(
        self,
        mass: pint.Quantity,
        *,
        target_unit: str | Unit = "mole",
        active_only: bool = False,
    ) -> pint.Quantity:
        """Convert mass to chemical substance amount using ``molar_mass`` (and ``purity`` if ``active_only``)."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        if self.molar_mass is None:
            raise ValueError("molar_mass is required for mass_to_moles conversion.")
        kg_val, _ = _materials_mass_in_kg(mass, "mass")
        if active_only:
            if self.purity is None:
                raise ValueError("purity is required when active_only=True.")
            kg_val *= self.purity
        kg_per_mol = self.molar_mass.to("kg / mole").magnitude
        dst = parse_unit(target_unit)
        if dst.dimensionality != _materials_SUBSTANCE_DIM:
            raise pint.DimensionalityError(ureg.mole, dst, _materials_SUBSTANCE_DIM, dst.dimensionality)
        mol_q = ureg.Quantity(kg_val / kg_per_mol, "mole").to(dst)
        return ureg.Quantity(_materials_preserve_numeric_type(mol_q.magnitude, mass), dst)

    def moles_to_mass(
        self,
        moles: pint.Quantity,
        *,
        target_unit: str | Unit = "kilogram",
        active_only: bool = False,
    ) -> pint.Quantity:
        """Convert chemical substance amount to mass using ``molar_mass`` (dividing by ``purity`` if ``active_only``)."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        if self.molar_mass is None:
            raise ValueError("molar_mass is required for moles_to_mass conversion.")
        if not isinstance(moles, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(moles).__name__}.")
        validate_registry_quantity(moles)
        if moles.dimensionality != _materials_SUBSTANCE_DIM:
            raise pint.DimensionalityError(moles.units, ureg.mole, moles.dimensionality, _materials_SUBSTANCE_DIM)
        mol_dec = _materials_to_decimal(moles.magnitude, "moles")
        if mol_dec < Decimal(0):
            raise ValueError("Moles must be non-negative.")
        kg_per_mol = self.molar_mass.to("kg / mole").magnitude
        kg_val = ureg.Quantity(mol_dec, moles.units).to("mole").magnitude * kg_per_mol
        if active_only:
            if self.purity is None:
                raise ValueError("purity is required when active_only=True.")
            kg_val /= self.purity
        dst = parse_unit(target_unit)
        return _materials_kg_to_target(kg_val, dst, moles)

    def kinematic_viscosity(self, target_unit: str | Unit = "centistokes") -> pint.Quantity:
        """Return kinematic viscosity, converting dynamic viscosity via ``density`` if needed."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        if self.viscosity is None:
            raise ValueError("viscosity is not defined on MaterialProperties.")
        dst = parse_unit(target_unit)
        if dst.dimensionality != _materials_KINEMATIC_VISC_DIM:
            raise pint.DimensionalityError(ureg.centistokes, dst, _materials_KINEMATIC_VISC_DIM, dst.dimensionality)
        v = self.viscosity
        if v.dimensionality == _materials_KINEMATIC_VISC_DIM:
            return v.to(dst)
        rho = self._select_density(use_bulk=False)
        return (v / rho).to(dst)

    def dynamic_viscosity(self, target_unit: str | Unit = "centipoise") -> pint.Quantity:
        """Return dynamic viscosity, converting kinematic viscosity via ``density`` if needed."""
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        if self.viscosity is None:
            raise ValueError("viscosity is not defined on MaterialProperties.")
        dst = parse_unit(target_unit)
        if dst.dimensionality != _materials_DYNAMIC_VISC_DIM:
            raise pint.DimensionalityError(ureg.centipoise, dst, _materials_DYNAMIC_VISC_DIM, dst.dimensionality)
        v = self.viscosity
        if v.dimensionality == _materials_DYNAMIC_VISC_DIM:
            return v.to(dst)
        rho = self._select_density(use_bulk=False)
        return (v * rho).to(dst)

    def as_context(
        self,
        *,
        use_bulk_density: bool = False,
        active_only: bool = False,
        name: str | None = None,
    ) -> pint.Context:
        """Create a Pint Context bridging [volume], [mass], and [substance] for these properties."""
        from core.units.context import create_material_context
        return create_material_context(
            self,
            use_bulk_density=use_bulk_density,
            active_only=active_only,
            name=name,
        )

    @authorized_conversion
    def normalize_price(
        self,
        price: pint.Quantity,
        target_unit: str | Unit | None = None,
        *,
        from_purity_basis: Literal["gross", "pure"] = "gross",
        to_purity_basis: Literal["gross", "pure"] = "gross",
        from_moisture_basis: Literal["wet", "dry"] | None = None,
        to_moisture_basis: Literal["wet", "dry"] | None = None,
        use_bulk_density: bool = False,
        rates: Mapping[str, Decimal | float | int] | None = None,
        base_currency: str = "USD",
    ) -> pint.Quantity:
        """Normalise a commodity price quote across currencies, units, purity bases, and moisture bases.

        Remember inverse price scaling:
        - ``price_per_pure = price_per_gross / purity``
        - ``price_per_dry = price_per_wet / (1 - moisture)``
        """
        from core.units.analytical_basis import rebase_moisture_intensive
        from core.units.exceptions import MissingFXRateError
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        from core.units.conversions import _convert_price_quantity
        from core.units.conversions import _numerator_currency_iso
        from core.units.conversions import _split_ratio_units
        from core.units.conversions import convert_currency

        if not isinstance(price, pint.Quantity):
            raise TypeError(f"Expected pint.Quantity, got {type(price).__name__}.")
        validate_registry_quantity(price)
        if from_purity_basis not in {"gross", "pure"} or to_purity_basis not in {"gross", "pure"}:
            raise ValueError("Purity bases must be 'gross' or 'pure'.")
        if from_moisture_basis not in {None, "wet", "dry"} or to_moisture_basis not in {None, "wet", "dry"}:
            raise ValueError("Moisture bases must be 'wet', 'dry', or None.")

        src_num, src_den = _split_ratio_units(price)
        dst = parse_unit(target_unit) if target_unit is not None else price.units
        dst_num, dst_den = _split_ratio_units(dst)

        src_tags = set(ureg._semantic_tags(src_den._units)) & {"wet", "dry"}
        dst_tags = set(ureg._semantic_tags(dst_den._units)) & {"wet", "dry"}
        tag_from = next(iter(src_tags)) if len(src_tags) == 1 else None
        tag_to = next(iter(dst_tags)) if len(dst_tags) == 1 else None

        if from_moisture_basis is not None and tag_from is not None and from_moisture_basis != tag_from:
            raise ValueError(f"Price unit {price.units} ({tag_from} basis) contradicts from_moisture_basis={from_moisture_basis!r}.")
        if to_moisture_basis is not None and tag_to is not None and to_moisture_basis != tag_to:
            raise ValueError(f"Target unit {dst} ({tag_to} basis) contradicts to_moisture_basis={to_moisture_basis!r}.")

        eff_from_moist = from_moisture_basis or tag_from or "wet"
        eff_to_moist = to_moisture_basis or tag_to or eff_from_moist

        mag = _materials_to_decimal(price.magnitude, "price")

        # 1. Purity basis adjustment
        if from_purity_basis != to_purity_basis:
            if self.purity is None:
                raise ValueError("purity is required to convert price between 'gross' and 'pure' bases.")
            mag = mag / self.purity if (from_purity_basis == "gross" and to_purity_basis == "pure") else mag * self.purity

        # 2. Moisture basis adjustment
        if eff_from_moist != eff_to_moist:
            if self.moisture is None:
                raise ValueError("moisture is required to convert price between 'wet' and 'dry' bases.")
            rebased = rebase_moisture_intensive(
                mag,
                source_moisture=self.moisture if eff_from_moist == "wet" else Decimal(0),
                target_moisture=Decimal(0) if eff_to_moist == "dry" else self.moisture,
            )
            mag = rebased if isinstance(rebased, Decimal) else Decimal(str(rebased))

        working = ureg.Quantity(mag, price.units)
        if target_unit is None:
            return ureg.Quantity(_materials_preserve_numeric_type(working.magnitude, price), working.units)

        # 3. Optional FX conversion when ISO currencies differ
        src_iso = _numerator_currency_iso(working)
        dst_iso = _numerator_currency_iso(ureg.Quantity(1, dst_num))
        if src_iso != dst_iso:
            if rates is None:
                raise MissingFXRateError(f"FX rates are required to convert {src_iso} to {dst_iso}.")
            working = convert_currency(working, dst_num, rates=dict(rates), base_currency=base_currency)

        # Project semantic wet/dry tags on denominators now that moisture scaling has been applied
        if src_tags or dst_tags:
            src_factor, src_phys_units = ureg._project_mass_units(src_den._units)
            dst_factor, dst_phys_units = ureg._project_mass_units(dst_den._units)
            _, cur_num_after_fx = _split_ratio_units(working)[0], working
            work_num, _ = _split_ratio_units(working)
            phys_working = ureg.Quantity(
                _materials_to_decimal(working.magnitude, "price") / _materials_to_decimal(src_factor, "src_factor"),
                work_num / ureg.Unit(src_phys_units),
            )
            phys_dst = dst_num / ureg.Unit(dst_phys_units)
            ctx = self.as_context(use_bulk_density=use_bulk_density, active_only=False)
            with ureg.context(ctx):
                if ureg.Unit(src_phys_units).dimensionality != ureg.Unit(dst_phys_units).dimensionality:
                    phys_conv = _convert_price_quantity(phys_working, phys_dst, contexts=[ctx])
                else:
                    phys_conv = phys_working.to(phys_dst)
            final_mag = _materials_to_decimal(phys_conv.magnitude, "price") * _materials_to_decimal(dst_factor, "dst_factor")
            return ureg.Quantity(_materials_preserve_numeric_type(final_mag, price), dst)

        ctx = self.as_context(use_bulk_density=use_bulk_density, active_only=False)
        with ureg.context(ctx):
            if src_den.dimensionality != dst_den.dimensionality:
                converted = _convert_price_quantity(working, dst, contexts=[ctx])
            else:
                converted = working.to(dst)
        return ureg.Quantity(_materials_preserve_numeric_type(_materials_to_decimal(converted.magnitude, "price"), price), dst)

    def blend_with(
        self,
        other: MaterialProperties,
        *,
        self_mass_fraction: Decimal | float,
    ) -> MaterialProperties:
        """Blend two MaterialProperties records by mass fraction.

        - Densities: Harmonic mean by mass (additive ideal volumes: ``1/rho = w1/rho1 + w2/rho2``).
        - Molar mass: Harmonic mean by mass (conserving total moles: ``1/M = w1/M1 + w2/M2``).
        - Moisture & purity: Mass-weighted arithmetic mean.
        - Viscosity: Non-linear Refutas Viscosity Blending Number (VBN) in cSt when available.
        """
        from core.units.registry import ureg
        if not isinstance(other, MaterialProperties):
            raise TypeError(f"Expected MaterialProperties, got {type(other).__name__}.")
        w1 = _materials_validate_fraction(self_mass_fraction, "self_mass_fraction")
        w2 = Decimal(1) - w1

        def harmonic_qty(q1: pint.Quantity | None, q2: pint.Quantity | None, base_unit: str) -> pint.Quantity | None:
            from core.units.registry import ureg
            if q1 is None or q2 is None:
                return None
            v1 = q1.to(base_unit).magnitude
            v2 = q2.to(base_unit).magnitude
            blended = Decimal(1) / ((w1 / v1) + (w2 / v2))
            return ureg.Quantity(blended, base_unit).to(q1.units)

        b_dens = harmonic_qty(self.effective_true_density, other.effective_true_density, "kg / meter ** 3")
        b_bulk = harmonic_qty(self.effective_bulk_density, other.effective_bulk_density, "kg / meter ** 3")
        b_mm = harmonic_qty(self.molar_mass, other.molar_mass, "kg / mole")

        b_moist = (w1 * self.moisture + w2 * other.moisture) if (self.moisture is not None and other.moisture is not None) else None
        b_pur = (w1 * self.purity + w2 * other.purity) if (self.purity is not None and other.purity is not None) else None

        b_visc: pint.Quantity | None = None
        if self.viscosity is not None and other.viscosity is not None:
            try:
                nu1 = float(self.kinematic_viscosity("centistokes").magnitude)
                nu2 = float(other.kinematic_viscosity("centistokes").magnitude)
                vbn = float(w1) * _materials_refutas_vbn(nu1) + float(w2) * _materials_refutas_vbn(nu2)
                nu_blend = Decimal(str(_materials_cst_from_refutas_vbn(vbn)))
                if self.viscosity.dimensionality == _materials_DYNAMIC_VISC_DIM and b_dens is not None:
                    b_visc = (ureg.Quantity(nu_blend, "centistokes") * b_dens).to(self.viscosity.units)
                else:
                    b_visc = ureg.Quantity(nu_blend, "centistokes")
            except ValueError:
                b_visc = None

        b_por: Decimal | None = None
        if b_dens is None or b_bulk is None:
            if self.porosity is not None and other.porosity is not None:
                b_por = w1 * self.porosity + w2 * other.porosity

        return MaterialProperties(
            density=b_dens,
            bulk_density=b_bulk,
            molar_mass=b_mm,
            viscosity=b_visc,
            moisture=b_moist,
            purity=b_pur,
            porosity=b_por,
        )





_gas_data_DATA = json.loads(((Path(__file__).resolve().parents[0] / "reference_data") / "assays" / "gas_components.json").read_text(encoding="utf-8"))



GAS_COMPONENT_DATA_VERSION = _gas_data_DATA["version"]



GAS_MOLAR_MASSES = MappingProxyType({name: Decimal(item["molar_mass_g_mol"]) for name, item in _gas_data_DATA["components"].items()})





# Conventional abridged atomic weights; not isotope-specific measurements.
# CIAAW Atomic Weights 2021, including the 2024 revisions, retrieved 2026-11-27.
# https://ciaaw.org/abridged-atomic-weights.htm
ATOMIC_MASSES = MappingProxyType({k: Decimal(v) for k, v in {
    "H": "1.0080", "He": "4.0026", "Li": "6.94", "B": "10.81", "C": "12.011", "N": "14.007",
    "O": "15.999", "F": "18.998", "Na": "22.990", "Mg": "24.305", "Al": "26.982", "Si": "28.085",
    "P": "30.974", "S": "32.06", "Cl": "35.45", "Ar": "39.95", "K": "39.098", "Ca": "40.078",
    "Ti": "47.867", "V": "50.942", "Cr": "51.996", "Mn": "54.938", "Fe": "55.845", "Co": "58.933",
    "Ni": "58.693", "Cu": "63.546", "Zn": "65.38", "Ga": "69.723", "Ge": "72.630", "As": "74.922",
    "Sc": "44.956", "Sr": "87.62", "Y": "88.906", "Zr": "91.222", "Nb": "92.906", "Mo": "95.95", "Ag": "107.87",
    "Cd": "112.41", "Sn": "118.71", "Sb": "121.76", "Ba": "137.33", "La": "138.91", "Ce": "140.12",
    "Pr": "140.91", "Nd": "144.24", "Sm": "150.36", "Eu": "151.96", "Gd": "157.25", "Tb": "158.93",
    "Dy": "162.50", "Ho": "164.93", "Er": "167.26", "Tm": "168.93", "Yb": "173.05", "Lu": "174.97",
    "W": "183.84", "Pt": "195.08", "Au": "196.97", "Hg": "200.59", "Pb": "207.2", "Bi": "208.98",
}.items()})





@dataclass(frozen=True, slots=True)
class StoichiometricCompound:
    """Simple molecular formula with nested parentheses; no inferred speciation.

    Write hydrates as e.g. CaSO4(H2O)2. Isotopes, ionic charges and dot-hydrate
    notation require explicit external speciation rather than guessed parsing.
    """
    formula: str
    atoms: Mapping[str, int] = field(init=False)

    def __post_init__(self):
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        formula = _converters_petroleum_text(self.formula, "formula")
        tokens = re.findall(r"[A-Z][a-z]?|[0-9]+|[()]", formula)
        if "".join(tokens) != formula:
            raise ValueError("Unsupported formula syntax; use element symbols and parentheses.")
        stack = [{}]
        index = 0
        while index < len(tokens):
            token = tokens[index]
            if token == "(":
                stack.append({})
                index += 1
                continue
            if token == ")":
                if len(stack) == 1 or not stack[-1]:
                    raise ValueError("Unbalanced or empty formula group.")
                group = stack.pop()
            elif token in ATOMIC_MASSES:
                group = {token: 1}
            else:
                raise ValueError(f"Unknown element or misplaced count {token!r}.")
            index += 1
            count = 1
            if index < len(tokens) and tokens[index].isdigit():
                count = int(tokens[index])
                index += 1
                if count <= 0:
                    raise ValueError("Atom counts must be positive integers.")
            for element, number in group.items():
                stack[-1][element] = stack[-1].get(element, 0) + number * count
        if len(stack) != 1 or not stack[0]:
            raise ValueError("Unbalanced or empty formula.")
        object.__setattr__(self, "formula", formula)
        object.__setattr__(self, "atoms", MappingProxyType(stack[0]))

    @property
    def molar_mass(self):
        from core.units.registry import ureg
        return ureg.Quantity(sum((ATOMIC_MASSES[e] * n for e, n in self.atoms.items()), Decimal(0)), "g/mol")

    def element_fraction(self, element: str) -> Decimal:
        if element not in self.atoms:
            raise ValueError(f"{element!r} is absent from {self.formula}.")
        return ATOMIC_MASSES[element] * self.atoms[element] / self.molar_mass.magnitude

    def contained_element(self, quantity, *, element):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, self.element_fraction(element))

    def equivalent_mass(self, quantity, target_formula: str, *, element: str):
        from core.units.assays.base import scale_mass
        destination = StoichiometricCompound(target_formula)
        return scale_mass(quantity, self.element_fraction(element) / destination.element_fraction(element))





@dataclass(frozen=True, slots=True)
class ElementComposition:
    """Independent assayed mass fractions; unassayed balance is retained.

    Keys may be elements or explicitly reported species such as CaF2. Do not
    mix overlapping elemental and oxide analyses into the same mass balance.
    """
    mass_fractions: Mapping[str, Decimal]
    normalize: bool = False
    source_total: Decimal = field(init=False)

    def __post_init__(self):
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        if not isinstance(self.normalize, bool):
            raise TypeError("normalize must be boolean.")
        if not isinstance(self.mass_fractions, Mapping) or not self.mass_fractions:
            raise ValueError("Supply at least one assayed component.")
        raw = {_converters_petroleum_text(k, "component"): _converters_petroleum_number(v, "component fraction", minimum=0) for k, v in self.mass_fractions.items()}
        total = sum(raw.values(), Decimal(0))
        if self.normalize:
            if total <= 0:
                raise ValueError("Cannot normalize an empty composition.")
            raw = {k: v / total for k, v in raw.items()}
            largest = max(raw, key=raw.get)
            raw[largest] += 1 - sum(raw.values(), Decimal(0))
        elif total > 1:
            raise ValueError("Independent mass fractions cannot exceed one; normalize explicitly for proportional input.")
        object.__setattr__(self, "mass_fractions", MappingProxyType(raw))
        object.__setattr__(self, "source_total", total)

    @classmethod
    def from_grades(cls, grades: Mapping[str, pint.Quantity]):
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        return cls({k: _converters_petroleum_quantity_number(v, "dimensionless", "assay grade", minimum=0) for k, v in grades.items()})

    @property
    def unassayed_fraction(self):
        return 1 - sum(self.mass_fractions.values(), Decimal(0))

    def fraction(self, component):
        if component not in self.mass_fractions:
            raise ValueError(f"No assayed fraction for {component!r}.")
        return self.mass_fractions[component]





@dataclass(frozen=True, slots=True)
class GasComposition:
    mole_fractions: Mapping[str, Decimal]
    normalize: bool = False
    source_total: Decimal = field(init=False)

    def __post_init__(self) -> None:
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        if not isinstance(self.normalize, bool):
            raise TypeError("normalize must be a boolean.")
        if not isinstance(self.mole_fractions, Mapping) or not self.mole_fractions:
            raise ValueError("Supply a non-empty component-to-mole-fraction mapping.")
        unknown = set(self.mole_fractions) - set(GAS_MOLAR_MASSES)
        if unknown:
            raise ValueError(f"Unknown gas components: {sorted(unknown, key=str)}")
        values = {name: _converters_petroleum_number(value, f"mole fraction {name}", minimum=0) for name, value in self.mole_fractions.items()}
        total = sum(values.values(), Decimal(0))
        if total <= 0:
            raise ValueError("Gas composition must have positive total moles.")
        if total != 1 and not self.normalize:
            raise ValueError("Mole fractions must sum exactly to one; request normalize=True explicitly for proportional data.")
        values = {name: value / total if self.normalize else value for name, value in values.items() if value > 0}
        if self.normalize:
            # Decimal division can leave a last-place residual. Assign it to the
            # largest component; never normalize unknown or negative components.
            largest = max(values, key=values.get)
            values[largest] += 1 - sum(values.values(), Decimal(0))
        object.__setattr__(self, "mole_fractions", MappingProxyType(values))
        object.__setattr__(self, "source_total", total)

    @property
    def molar_mass(self) -> pint.Quantity:
        from core.units.registry import ureg
        return ureg.Quantity(sum((fraction * GAS_MOLAR_MASSES[name] for name, fraction in self.mole_fractions.items()), Decimal(0)), "g/mol")

    @property
    def ideal_relative_density(self) -> Decimal:
        from core.units.assays.energy.lng import _converters_gas_AIR_MOLAR_MASS
        return self.molar_mass.magnitude / _converters_gas_AIR_MOLAR_MASS

    @property
    def water_mole_fraction(self) -> Decimal:
        return self.mole_fractions.get("H2O", Decimal(0))

    def dry(self) -> GasComposition:
        remaining = {name: value for name, value in self.mole_fractions.items() if name != "H2O"}
        if not remaining:
            raise ValueError("Pure water has no dry-gas composition.")
        return GasComposition(remaining, normalize=True)

    def with_water(self, water_mole_fraction: Decimal) -> GasComposition:
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        water = _converters_petroleum_fraction(water_mole_fraction, "water mole fraction", allow_one=False)
        values = {name: value * (1 - water) for name, value in self.dry().mole_fractions.items()}
        values["H2O"] = water
        return GasComposition(values, normalize=True)

    def calorific_values(self) -> GasCalorificValues:
        from core.units.heat_content import CalorificMethod
        from core.units.heat_content import ISO6976_COMPONENTS
        from core.units.heat_content import ThermodynamicBasis
        from core.units.registry import ureg
        missing = set(self.mole_fractions) - set(ISO6976_COMPONENTS)
        if missing:
            raise ValueError(f"No bundled calorific data for components {sorted(missing)}; supply measured calorific values explicitly.")
        # Incoming water is spectator vapour, not a combustible constituent.
        # The reference table's water phase-change enthalpy is not gas HHV.
        # ASTM D3588 wet-gas convention: incoming water remains vapour; only
        # water formed by combustion contributes to the HHV/LHV difference.
        combustion = ((name, fraction) for name, fraction in self.mole_fractions.items() if name != "H2O")
        hhv = sum((fraction * Decimal(str(ISO6976_COMPONENTS[name]["molar_hhv_kj_mol_15c"])) for name, fraction in combustion), Decimal(0))
        lhv = sum((fraction * Decimal(str(ISO6976_COMPONENTS[name]["molar_lhv_kj_mol_15c"])) for name, fraction in self.mole_fractions.items() if name != "H2O"), Decimal(0))
        # kJ/mol divided by g/mol is MJ/kg.
        return GasCalorificValues(ureg.Quantity(hhv / self.molar_mass.magnitude, "MJ/kg"),
                                  ureg.Quantity(lhv / self.molar_mass.magnitude, "MJ/kg"),
                                  ureg.Quantity(15, "degC"), "energy_specs.json component sums; rounded 15 C combustion data; spectator water remains vapour",
                                  method=CalorificMethod.COMPONENT_SUM, thermodynamic_basis=ThermodynamicBasis.CONSTANT_PRESSURE)

    def blend(self, other: GasComposition, *, self_mole_fraction: Decimal) -> GasComposition:
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        if not isinstance(other, GasComposition):
            raise TypeError("Gas blending requires another GasComposition.")
        fraction = _converters_petroleum_fraction(self_mole_fraction, "self mole fraction")
        values = {name: self.mole_fractions.get(name, Decimal(0)) * fraction + other.mole_fractions.get(name, Decimal(0)) * (1 - fraction)
                  for name in self.mole_fractions.keys() | other.mole_fractions.keys()}
        return GasComposition(values, normalize=True)

    def quality(self, *, minimum_methane: Decimal, maximum_co2: Decimal, maximum_h2s: Decimal,
                maximum_water: Decimal) -> Mapping[str, bool]:
        """Evaluate explicit mole-fraction limits; no universal biomethane grade is assumed."""
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        limits = {name: _converters_petroleum_fraction(value, name) for name, value in
                  (("CH4", minimum_methane), ("CO2", maximum_co2), ("H2S", maximum_h2s), ("H2O", maximum_water))}
        return MappingProxyType({"methane": self.mole_fractions.get("CH4", Decimal(0)) >= limits["CH4"],
                                 "co2": self.mole_fractions.get("CO2", Decimal(0)) <= limits["CO2"],
                                 "h2s": self.mole_fractions.get("H2S", Decimal(0)) <= limits["H2S"],
                                 "water": self.water_mole_fraction <= limits["H2O"]})

    def to_gas_assay(self, *, name: str = "composition gas"):
        """Compatibility assay using explicit ideal reference-volume assumptions."""
        from core.units.reference_conditions import GasReferenceConditions
        from core.units.registry import ureg
        from core.units.thermodynamics import builtin_gas_backend
        from core.units.thermodynamics import gas_properties
        from core.units.assays.energy.natural_gas import GasAssay
        values = self.calorific_values()
        reference = GasReferenceConditions.from_standard("ferc", compressibility_factor=Decimal(1))
        backend = builtin_gas_backend("ideal", temperature_range=(ureg.Quantity(270, "K"), ureg.Quantity(330, "K")),
                                      pressure_range=(ureg.Quantity(90000, "Pa"), ureg.Quantity(120000, "Pa")))
        properties = gas_properties(self, reference, backend=backend, calorific_values=values)
        hhv = properties.hhv.to("btu_it/meter**3").magnitude * Decimal("0.3048") ** 3
        lhv = properties.lhv.to("btu_it/meter**3").magnitude * Decimal("0.3048") ** 3
        if hhv <= 0:
            raise ValueError("A noncombustible industrial gas has no positive-heating GasAssay.")
        return GasAssay(name=name, gcv_btu_per_scf=float(hhv), ncv_btu_per_scf=float(lhv),
                        relative_density=float(self.ideal_relative_density), molar_composition=self.mole_fractions)


from core.units.heat_content import CalorificMethod
from core.units.heat_content import ThermodynamicBasis
from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GasCalorificValues(_safety_QuantitySnapshot):
    hhv: pint.Quantity
    lhv: pint.Quantity
    combustion_temperature: pint.Quantity
    source: str
    method: CalorificMethod = CalorificMethod.MEASURED
    thermodynamic_basis: ThermodynamicBasis = ThermodynamicBasis.UNSPECIFIED

    def __post_init__(self) -> None:
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.heat_content import CalorificMethod
        from core.units.heat_content import ThermodynamicBasis
        from core.units.registry import ureg
        object.__setattr__(self, "method", CalorificMethod(self.method))
        object.__setattr__(self, "thermodynamic_basis", ThermodynamicBasis(self.thermodynamic_basis))
        hhv = _converters_petroleum_quantity_number(self.hhv, "MJ/kg", "HHV", minimum=0)
        lhv = _converters_petroleum_quantity_number(self.lhv, "MJ/kg", "LHV", minimum=0)
        if lhv > hhv:
            raise ValueError("LHV cannot exceed HHV.")
        object.__setattr__(self, "hhv", ureg.Quantity(hhv, "MJ/kg"))
        object.__setattr__(self, "lhv", ureg.Quantity(lhv, "MJ/kg"))
        object.__setattr__(self, "combustion_temperature", _converters_petroleum_temperature(self.combustion_temperature))
        object.__setattr__(self, "source", _converters_petroleum_text(self.source, "calorific source"))


__all__ = ('ATOMIC_MASSES', 'CompositionBasis', 'ElementComposition', 'GasCalorificValues', 'GasComposition', 'MaterialComposition', 'MaterialProperties', 'StoichiometricCompound')
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
    from core.units.registry import Quantity
    from core.units.registry import Unit
