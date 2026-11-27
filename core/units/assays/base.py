"""Commodex assays base implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from dataclasses import dataclass, field
from decimal import Decimal
from enum import Enum
from itertools import permutations
from types import MappingProxyType
from typing import Any
from typing import Any, Final, Literal
from typing import ClassVar
import inspect
import math
import pint
import re

from core.units.composition import MaterialProperties


# ---------------------------------------------------------------------------
# 4. Unified MaterialSpec & Non-Destructive Assay Adapters
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class MaterialSpec:
    """Composable commodity specification combining properties, composition, and processing."""

    name: str
    properties: MaterialProperties = field(default_factory=MaterialProperties)
    composition: MaterialComposition | None = None
    processing: ProcessingSpecification | None = None

    def __post_init__(self) -> None:
        from core.units.composition import MaterialComposition
        from core.units.composition import MaterialProperties
        from core.units.composition import _materials_GRAVIMETRIC_BASES
        from core.units.composition import _materials_clean_name
        from core.units.processing import ProcessingSpecification
        clean = _materials_clean_name(self.name, "name")
        if not isinstance(self.properties, MaterialProperties):
            raise TypeError(f"properties must be a MaterialProperties, got {type(self.properties).__name__}.")
        if self.composition is not None and not isinstance(self.composition, MaterialComposition):
            raise TypeError(f"composition must be a MaterialComposition or None, got {type(self.composition).__name__}.")
        if self.processing is not None and not isinstance(self.processing, ProcessingSpecification):
            raise TypeError(f"processing must be a ProcessingSpecification or None, got {type(self.processing).__name__}.")

        props = self.properties
        comp = self.composition
        if comp is not None:
            if props.moisture is not None and comp.moisture_fraction is not None:
                if abs(props.moisture - comp.moisture_fraction) > Decimal("1e-9"):
                    raise ValueError(
                        f"MaterialProperties.moisture ({props.moisture}) conflicts with MaterialComposition.moisture_fraction ({comp.moisture_fraction})."
                    )
            elif props.moisture is None and comp.moisture_fraction is not None:
                from dataclasses import replace
                props = replace(props, moisture=comp.moisture_fraction)
            elif comp.moisture_fraction is None and props.moisture is not None and comp.basis in _materials_GRAVIMETRIC_BASES:
                from dataclasses import replace
                comp = replace(comp, moisture_fraction=props.moisture)

        object.__setattr__(self, "name", clean)
        object.__setattr__(self, "properties", props)
        object.__setattr__(self, "composition", comp)

    def as_context(
        self,
        *,
        use_bulk_density: bool = False,
        active_only: bool = False,
        name: str | None = None,
    ) -> pint.Context:
        ctx_name = name or f"material_{self.name.lower().replace(' ', '_')}"
        return self.properties.as_context(
            use_bulk_density=use_bulk_density,
            active_only=active_only,
            name=ctx_name,
        )

    def volume_to_mass(self, volume: pint.Quantity, *, target_unit: str | Unit = "kilogram", use_bulk: bool = False) -> pint.Quantity:
        return self.properties.volume_to_mass(volume, target_unit=target_unit, use_bulk=use_bulk)

    def mass_to_volume(self, mass: pint.Quantity, *, target_unit: str | Unit = "meter ** 3", use_bulk: bool = False) -> pint.Quantity:
        return self.properties.mass_to_volume(mass, target_unit=target_unit, use_bulk=use_bulk)

    def gross_to_pure(self, gross_quantity: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        return self.properties.gross_to_pure(gross_quantity, target_unit=target_unit)

    def pure_to_gross(self, pure_quantity: pint.Quantity, *, target_unit: str | Unit | None = None) -> pint.Quantity:
        return self.properties.pure_to_gross(pure_quantity, target_unit=target_unit)

    def wet_to_dry(self, wet_mass: pint.Quantity, *, target_unit: str | Unit = "dmt") -> pint.Quantity:
        if self.composition is not None and self.composition.moisture_fraction is not None:
            return self.composition.wet_to_dry_mass(wet_mass, target_unit=target_unit)
        return self.properties.wet_to_dry(wet_mass, target_unit=target_unit)

    def dry_to_wet(self, dry_mass: pint.Quantity, *, target_unit: str | Unit = "wmt") -> pint.Quantity:
        if self.composition is not None and self.composition.moisture_fraction is not None:
            return self.composition.dry_to_wet_mass(dry_mass, target_unit=target_unit)
        return self.properties.dry_to_wet(dry_mass, target_unit=target_unit)

    def contained_mass(
        self,
        mass: pint.Quantity,
        component: str,
        *,
        input_basis: CompositionBasis | str | None = None,
        target_unit: str | Unit = "contained_metal_tonne",
    ) -> pint.Quantity:
        if self.composition is None:
            raise ValueError(f"{self.name!r} does not define a MaterialComposition.")
        return self.composition.contained_mass(mass, component, input_basis=input_basis, target_unit=target_unit)

    def required_mass_for_contained(
        self,
        contained: pint.Quantity,
        component: str,
        *,
        output_basis: CompositionBasis | str | None = None,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        if self.composition is None:
            raise ValueError(f"{self.name!r} does not define a MaterialComposition.")
        return self.composition.required_mass_for_contained(
            contained, component, output_basis=output_basis, target_unit=target_unit
        )

    def input_to_output(
        self,
        input_mass: pint.Quantity,
        output: str,
        *,
        apply_recovery: bool = True,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        if self.processing is None:
            raise ValueError(f"{self.name!r} does not define a ProcessingSpecification.")
        return self.processing.input_to_output(
            input_mass, output, apply_recovery=apply_recovery, target_unit=target_unit
        )

    def output_to_input(
        self,
        output_mass: pint.Quantity,
        output: str,
        *,
        apply_recovery: bool = True,
        target_unit: str | Unit | None = None,
    ) -> pint.Quantity:
        if self.processing is None:
            raise ValueError(f"{self.name!r} does not define a ProcessingSpecification.")
        return self.processing.output_to_input(
            output_mass, output, apply_recovery=apply_recovery, target_unit=target_unit
        )

    def normalize_price(self, price: pint.Quantity, target_unit: str | Unit | None = None, **kwargs: Any) -> pint.Quantity:
        return self.properties.normalize_price(price, target_unit=target_unit, **kwargs)

    def blend_with(
        self,
        other: MaterialSpec,
        *,
        self_mass_fraction: Decimal | float,
        blend_name: str | None = None,
    ) -> MaterialSpec:
        from core.units.composition import CompositionBasis
        from core.units.composition import _materials_validate_fraction
        if not isinstance(other, MaterialSpec):
            raise TypeError(f"Expected MaterialSpec, got {type(other).__name__}.")
        w1 = _materials_validate_fraction(self_mass_fraction, "self_mass_fraction")
        b_props = self.properties.blend_with(other.properties, self_mass_fraction=w1)
        b_comp = (
            self.composition.blend_with(
                other.composition,
                self_weight=w1,
                weight_basis=(
                    CompositionBasis.AS_RECEIVED
                    if (
                        self.composition.basis == CompositionBasis.DRY
                        and self.composition.moisture_fraction is not None
                        and other.composition.moisture_fraction is not None
                    )
                    else None
                ),
            )
            if (self.composition is not None and other.composition is not None)
            else None
        )
        b_proc = (
            self.processing.blend_with(other.processing, self_input_fraction=w1)
            if (self.processing is not None and other.processing is not None)
            else None
        )
        return MaterialSpec(
            name=blend_name or f"{self.name}/{other.name} Blend",
            properties=b_props,
            composition=b_comp,
            processing=b_proc,
        )

    @classmethod
    def from_assay(cls, assay: Any) -> MaterialSpec:
        """Extract composable MaterialSpec building blocks from any specialised assay instance."""
        return material_spec_from_assay(assay)





def material_spec_from_assay(assay: Any) -> MaterialSpec:
    """Adapt an existing specialised assay into ``MaterialSpec(properties, composition, processing)``."""
    from core.units.composition import CompositionBasis
    from core.units.composition import MaterialComposition
    from core.units.composition import MaterialProperties
    from core.units.composition import _materials_to_decimal
    from core.units.processing import ProcessingSpecification
    from core.units.registry import ureg
    from core.units import agricultural_specs as ag
    from core.units import conversions as conv

    if isinstance(assay, MaterialSpec):
        return assay
    if callable(getattr(assay, "to_material_spec", None)):
        converted_spec = assay.to_material_spec()
        if isinstance(converted_spec, MaterialSpec):
            return converted_spec
    if isinstance(assay, MaterialProperties):
        return MaterialSpec(name="material", properties=assay)
    if isinstance(assay, MaterialComposition):
        return MaterialSpec(
            name="material",
            properties=MaterialProperties(moisture=assay.moisture_fraction),
            composition=assay,
        )

    name = getattr(assay, "name", type(assay).__name__)

    if isinstance(assay, conv.CrudeAssay):
        imp = {"S": Decimal(str(assay.sulfur_percentage)) / Decimal(100)} if assay.sulfur_percentage is not None else {}
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(density=ureg.Quantity(Decimal(str(assay.density_kg_per_m3)), "kg / meter ** 3")),
            composition=MaterialComposition(basis=CompositionBasis.AS_RECEIVED, impurities=imp) if imp else None,
        )

    if isinstance(assay, conv.RefineryCrudeAssay):
        base = material_spec_from_assay(assay.bulk)
        raw_yields = {c.name: Decimal(str(c.yield_wt_pct)) / Decimal(100) for c in assay.cuts if c.yield_wt_pct > 0}
        total_y = sum(raw_yields.values(), Decimal(0))
        norm_yields = {k: v / total_y for k, v in raw_yields.items()} if total_y > Decimal(1) else raw_yields
        visc = ureg.Quantity(Decimal(str(assay.viscosity_cst_50c)), "centistokes") if assay.viscosity_cst_50c is not None else None
        from dataclasses import replace
        return MaterialSpec(
            name=name,
            properties=replace(base.properties, viscosity=visc),
            composition=base.composition,
            processing=ProcessingSpecification(input_material=assay.bulk.name, output_yields=norm_yields),
        )

    if isinstance(assay, conv.OreAssay):
        moist = Decimal(str(assay.moisture_pct)) / Decimal(100)
        grade = Decimal(str(assay.grade_pct)) / Decimal(100)
        pay_grade = Decimal(str(assay.payable_grade_pct)) / Decimal(100)
        rec = min(Decimal(1), pay_grade / grade) if grade > Decimal(0) else Decimal(0)
        bulk_rho = getattr(assay, "bulk_density_kg_per_m3", None)
        imp_raw = getattr(assay, "impurities_pct", None) or getattr(assay, "penalty_impurities_pct", None) or {}
        imp = {
            str(k): Decimal(str(v)) / Decimal(100)
            for k, v in imp_raw.items()
            if k != assay.primary_element and Decimal(0) <= Decimal(str(v)) <= Decimal(100)
        }
        if grade + sum(imp.values(), Decimal(0)) > Decimal(1):
            imp = {}
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                bulk_density=ureg.Quantity(Decimal(str(bulk_rho)), "kg / meter ** 3") if bulk_rho is not None else None,
                moisture=moist,
            ),
            composition=MaterialComposition(
                components={assay.primary_element: grade},
                basis=CompositionBasis.DRY,
                moisture_fraction=moist,
                impurities=imp,
            ),
            processing=ProcessingSpecification(
                input_material=name,
                output_yields={assay.primary_element: grade},
                recovery_rates={assay.primary_element: rec},
            ),
        )

    if isinstance(assay, conv.BulkOreAssay):
        moist = _materials_to_decimal(assay.moisture_fraction, "moisture_fraction") if assay.moisture_fraction is not None else None
        grades = {str(k): _materials_to_decimal(v, k) for k, v in assay.grades.items()}
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=moist),
            composition=MaterialComposition(components=grades, basis=CompositionBasis.DRY, moisture_fraction=moist),
        )

    if isinstance(assay, conv.BullionAssay):
        fin = Decimal(str(assay.fineness))
        pay = Decimal(str(assay.payable_factor))
        comps = {assay.primary_metal: fin}
        for k, v in (assay.byproduct_fineness or {}).items():
            comps[str(k)] = Decimal(str(v))
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(purity=fin),
            composition=MaterialComposition(components=comps, basis=CompositionBasis.AS_RECEIVED),
            processing=ProcessingSpecification(
                input_material=name,
                output_yields={assay.primary_metal: fin},
                recovery_rates={assay.primary_metal: pay},
            ),
        )

    if isinstance(assay, (conv.PetrochemAssay, conv.ChemicalAssay)):
        pur_pct = getattr(assay, "purity_pct", None)
        if pur_pct is None:
            pur_pct = getattr(assay, "concentration_pct", None)
        if pur_pct is None:
            pur_pct = getattr(assay, "concentration_wt_pct", None)
        pur = Decimal(str(pur_pct)) / Decimal(100) if pur_pct is not None else None
        moist = (Decimal(1) - pur) if isinstance(assay, conv.ChemicalAssay) and pur is not None and pur < Decimal(1) else None
        comp_key = getattr(assay, "formula", None) or name
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                density=ureg.Quantity(Decimal(str(assay.density_kg_per_m3)), "kg / meter ** 3"),
                molar_mass=ureg.Quantity(Decimal(str(assay.molecular_weight_g_per_mol)), "g / mole"),
                purity=pur,
                moisture=moist,
            ),
            composition=(
                MaterialComposition(
                    components={str(comp_key): pur},
                    basis=CompositionBasis.AS_RECEIVED,
                    moisture_fraction=moist,
                )
                if pur is not None
                else None
            ),
        )

    from core.units.assays.base import MaterialSpecification
    if isinstance(assay, MaterialSpecification):
        mm_attr = getattr(assay, "molar_mass", None)
        mm_qty = mm_attr if isinstance(mm_attr, pint.Quantity) else None
        pur_attr = getattr(assay, "active_mass_fraction", None)
        if pur_attr is None:
            pur_attr = getattr(assay, "purity_fraction", None)
        pur_dec = Decimal(str(pur_attr)) if pur_attr is not None and Decimal(str(pur_attr)) > 0 else None
        comp_obj: MaterialComposition | None = None
        if assay.composition is not None:
            comp_obj = MaterialComposition(
                components=dict(assay.composition.mass_fractions),
                basis=CompositionBasis.DRY,
                moisture_fraction=assay.moisture_fraction,
            )
        elif pur_dec is not None and getattr(assay, "solute_formula", None) is not None:
            comp_obj = MaterialComposition(
                components={str(getattr(assay, "solute_formula")): pur_dec},
                basis=CompositionBasis.AS_RECEIVED,
            )
        proc_obj: ProcessingSpecification | None = getattr(assay, "processing", None)
        if proc_obj is not None and not isinstance(proc_obj, ProcessingSpecification):
            raise TypeError("Declared processing must be ProcessingSpecification.")
        if proc_obj is None and assay.recovery_fractions and assay.composition is not None:
            proc_obj = ProcessingSpecification(
                input_material=name,
                output_yields={
                    k: assay.composition.mass_fractions[k]
                    for k in assay.recovery_fractions
                    if k in assay.composition.mass_fractions
                },
                recovery_rates=dict(assay.recovery_fractions),
            )
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                density=assay.density,
                molar_mass=mm_qty,
                moisture=assay.moisture_fraction,
                purity=pur_dec,
            ),
            composition=comp_obj,
            processing=proc_obj,
        )

    if isinstance(assay, conv.FertilizerAssay):
        nutrients = {
            k: Decimal(str(v)) / Decimal(100)
            for k, v in (("N", assay.n_pct), ("P2O5", assay.p2o5_pct), ("K2O", assay.k2o_pct))
            if v > 0
        }
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                density=ureg.Quantity(Decimal(str(assay.density_kg_per_m3)), "kg / meter ** 3")
                if assay.density_kg_per_m3 is not None
                else None
            ),
            composition=MaterialComposition(components=nutrients, basis=CompositionBasis.AS_RECEIVED),
        )

    if isinstance(assay, (conv.BiomassAssay, conv.PetcokeAssay)):
        moist = Decimal(str(assay.moisture_pct)) / Decimal(100)
        ash_db = Decimal(str(assay.ash_db_pct)) / Decimal(100)
        comps: dict[str, Decimal] = {}
        if hasattr(assay, "fixed_carbon_db_pct") and getattr(assay, "fixed_carbon_db_pct") is not None:
            comps["fixed_carbon"] = Decimal(str(assay.fixed_carbon_db_pct)) / Decimal(100)
        if hasattr(assay, "hydrogen_db_pct") and getattr(assay, "hydrogen_db_pct") is not None:
            comps["H"] = Decimal(str(assay.hydrogen_db_pct)) / Decimal(100)
        imps: dict[str, Decimal] = {}
        if hasattr(assay, "sulfur_db_pct") and getattr(assay, "sulfur_db_pct") is not None:
            imps["S"] = Decimal(str(assay.sulfur_db_pct)) / Decimal(100)
        if hasattr(assay, "chlorine_db_pct") and getattr(assay, "chlorine_db_pct") is not None:
            imps["Cl"] = Decimal(str(assay.chlorine_db_pct)) / Decimal(100)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                bulk_density=ureg.Quantity(Decimal(str(assay.bulk_density_kg_per_m3)), "kg / meter ** 3"),
                moisture=moist,
            ),
            composition=MaterialComposition(
                components=comps,
                basis=CompositionBasis.DRY,
                moisture_fraction=moist,
                ash_fraction=ash_db,
                impurities=imps,
            ),
        )

    if isinstance(assay, conv.CoalAssay):
        moist = Decimal(str(assay.total_moisture_pct)) / Decimal(100)
        ash_ar = Decimal(str(assay.ash_pct)) / Decimal(100) if assay.ash_pct is not None else None
        comps_coal = (
            {"H": Decimal(str(assay.hydrogen_ar_pct)) / Decimal(100)}
            if assay.hydrogen_ar_pct is not None
            else {}
        )
        imps = {"S": Decimal(str(assay.sulfur_pct)) / Decimal(100)} if assay.sulfur_pct is not None else {}
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(moisture=moist),
            composition=MaterialComposition(
                components=comps_coal,
                basis=CompositionBasis.AS_RECEIVED,
                moisture_fraction=moist,
                ash_fraction=ash_ar,
                impurities=imps,
            ),
        )

    if isinstance(assay, ag.OilseedCrushAssay):
        oil_y = Decimal(str(assay.oil_yield_pct)) / Decimal(100)
        meal_y = Decimal(str(assay.meal_yield_pct)) / Decimal(100)
        loss = max(Decimal(0), Decimal(1) - oil_y - meal_y)
        oil_rho = Decimal(str(assay.oil_density_kg_per_m3))
        seed_rho = Decimal(str(assay.bulk_seed_density_kg_per_m3))
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                density=ureg.Quantity(oil_rho, "kg / meter ** 3") if oil_rho >= seed_rho else None,
                bulk_density=ureg.Quantity(seed_rho, "kg / meter ** 3"),
            ),
            processing=(
                ProcessingSpecification(
                    input_material=assay.seed,
                    output_yields={"oil": oil_y, "meal": meal_y},
                    process_losses={"crush_shrink": loss} if loss > Decimal(0) else {},
                )
                if (oil_y + meal_y) <= Decimal(1)
                else None
            ),
        )

    if isinstance(assay, ag.SoftsAssay):
        p_y = Decimal(str(assay.processing_yield_pct)) / Decimal(100)
        s_y = Decimal(str(assay.secondary_yield_pct)) / Decimal(100)
        yields = {"primary": p_y}
        if s_y > Decimal(0):
            yields["secondary"] = s_y
        loss = max(Decimal(0), Decimal(1) - p_y - s_y)
        return MaterialSpec(
            name=name,
            properties=MaterialProperties(
                bulk_density=ureg.Quantity(Decimal(str(assay.bulk_density_kg_per_m3)), "kg / meter ** 3")
                if assay.bulk_density_kg_per_m3 is not None
                else None
            ),
            processing=(
                ProcessingSpecification(
                    input_material=assay.commodity,
                    output_yields=yields,
                    process_losses={"processing_loss": loss} if loss > Decimal(0) else {},
                )
                if (p_y + s_y) <= Decimal(1)
                else None
            ),
        )

    # Generic structural extraction fallback for remaining specialised assays
    dens_val = getattr(assay, "density_kg_per_m3", None) or getattr(assay, "liquid_density_kg_per_m3", None)
    bulk_val = getattr(assay, "bulk_density_kg_per_m3", None)
    mm_val = getattr(assay, "molecular_weight_g_per_mol", None)
    moist_pct = getattr(assay, "moisture_pct", None)
    if moist_pct is None:
        moist_pct = getattr(assay, "actual_moisture_pct", None)
    pur_pct = getattr(assay, "purity_pct", None)

    dens_q = (
        ureg.Quantity(Decimal(str(dens_val)), "kg / meter ** 3")
        if dens_val is not None and Decimal(str(dens_val)) > 0
        else None
    )
    bulk_q = (
        ureg.Quantity(Decimal(str(bulk_val)), "kg / meter ** 3")
        if bulk_val is not None and Decimal(str(bulk_val)) > 0
        else None
    )
    if dens_q is not None and bulk_q is not None and bulk_q.magnitude > dens_q.magnitude:
        dens_q = None

    moist_dec = Decimal(str(moist_pct)) / Decimal(100) if moist_pct is not None else None
    if moist_dec is not None and not (Decimal(0) <= moist_dec < Decimal(1)):
        moist_dec = None
    pur_dec_fb = Decimal(str(pur_pct)) / Decimal(100) if pur_pct is not None else None
    if pur_dec_fb is not None and not (Decimal(0) < pur_dec_fb <= Decimal(1)):
        pur_dec_fb = None

    props = MaterialProperties(
        density=dens_q,
        bulk_density=bulk_q,
        molar_mass=(
            ureg.Quantity(Decimal(str(mm_val)), "g / mole")
            if mm_val is not None and Decimal(str(mm_val)) > 0
            else None
        ),
        moisture=moist_dec,
        purity=pur_dec_fb,
    )
    comp_dict = getattr(assay, "molar_composition", None)
    comp = (
        MaterialComposition(
            components={str(k): Decimal(str(v)) for k, v in comp_dict.items()},
            basis=CompositionBasis.MOLAR,
        )
        if isinstance(comp_dict, Mapping) and comp_dict
        else None
    )
    return MaterialSpec(name=name, properties=props, composition=comp)





def beneficiated_mass(spec, quantity, *, component, target_fraction, recovery_fraction, mass_basis="dry"):
    from core.units.assays.energy.petroleum import _converters_petroleum_fraction
    grade = _converters_petroleum_fraction(target_fraction, "target grade")
    recovery = _converters_petroleum_fraction(recovery_fraction, "recovery")
    if grade == 0:
        raise ValueError("Target grade must be positive.")
    dry = spec._dry_feed(quantity, mass_basis)
    output = scale_mass(spec.contained_mass(dry, component=component), recovery / grade)
    if mass(output) > mass(dry):
        raise ValueError("Product exceeds dry feed; dilution/additional feed is not declared.")
    return output





_material_specs_RARE_EARTHS = frozenset({"La", "Ce", "Pr", "Nd", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Y", "Sc"})





def fractions(values: Mapping[str, Any], label: str) -> Mapping[str, Decimal]:
    from core.units.assays.energy.petroleum import _converters_petroleum_fraction
    from core.units.assays.energy.petroleum import _converters_petroleum_text
    if not isinstance(values, Mapping):
        raise TypeError(f"{label} must be a mapping.")
    result = {_converters_petroleum_text(k, label + " key"): _converters_petroleum_fraction(v, label + ": " + str(k)) for k, v in values.items()}
    return MappingProxyType(result)





class _spec_engine_QuantityMap(Mapping):
    """Immutable canonical quantities with defensive views and replace support."""
    __slots__ = ("_values", "_unit")

    def __init__(self, values, unit):
        object.__setattr__(self, "_values", MappingProxyType(dict(values)))
        object.__setattr__(self, "_unit", unit)

    def __setattr__(self, name, value):
        raise AttributeError("Quantity mappings are immutable.")

    def __getitem__(self, key):
        from core.units.registry import ureg
        return ureg.Quantity(self._values[key], self._unit)

    def __iter__(self):
        return iter(self._values)

    def __len__(self):
        return len(self._values)


from core.units.validation import authorized_conversion


@authorized_conversion
def mass(value: pint.Quantity) -> Decimal:
    from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
    return _converters_petroleum_quantity_number(value, "kg", "mass", minimum=0)


from core.units.validation import authorized_conversion


@authorized_conversion
def scale_mass(value: pint.Quantity, factor: Any) -> pint.Quantity:
    # A processed mass must not inherit the feed's WMT/DMT/contained label.
    from core.units.assays.energy.petroleum import _converters_petroleum_number
    from core.units.assays.energy.petroleum import _converters_petroleum_output
    return _converters_petroleum_output(mass(value) * _converters_petroleum_number(factor, "mass factor", minimum=0), "kg", value)





def require_mass_basis(quantity, basis):
    from core.units.exceptions import IncompatibleAnalyticalBasisError
    mass(quantity)
    from core.units.conversions import _ore_mass_factor
    _, _, labelled = _ore_mass_factor(quantity.units)
    if labelled is not None and labelled != basis:
        raise IncompatibleAnalyticalBasisError(f"Mass label {labelled!r} conflicts with the requested {basis!r} feed basis.")





def add_molar_context(context, molar_mass):
    """Bridge declared gross mass to the explicitly defined mole basis."""
    from core.units.assays.energy.petroleum import _converters_petroleum_output
    from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
    from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
    kg_per_mol = _converters_petroleum_positive_quantity(molar_mass, "kg/mol", "molar mass")

    def forward(_registry, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        return _converters_petroleum_output(mass(quantity) / kg_per_mol, "mol", quantity)

    def inverse(_registry, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        moles = _converters_petroleum_quantity_number(quantity, "mol", "moles", minimum=0)
        return _converters_petroleum_output(moles * kg_per_mol, "kg", quantity)

    context.add_transformation("[mass]", "[substance]", forward)
    context.add_transformation("[substance]", "[mass]", inverse)
    return context





@dataclass(frozen=True, slots=True)
class ElementRecoveryEngine:
    composition: ElementComposition
    recovery_fractions: Mapping[str, Decimal] = field(default_factory=dict)
    payable_fractions: Mapping[str, Decimal] = field(default_factory=dict)
    deduction_fractions: Mapping[str, Decimal] = field(default_factory=dict)

    def __post_init__(self):
        from core.units.composition import ElementComposition
        if not isinstance(self.composition, ElementComposition):
            raise TypeError("composition must be ElementComposition.")
        for name in ("recovery_fractions", "payable_fractions", "deduction_fractions"):
            values = fractions(getattr(self, name), name)
            if values.keys() - self.composition.mass_fractions.keys():
                raise ValueError(f"{name} contains an unassayed component.")
            object.__setattr__(self, name, values)

    def contained_mass(self, dry_mass, *, component):
        require_mass_basis(dry_mass, "dry")
        return scale_mass(dry_mass, self.composition.fraction(component))

    def recovered_mass(self, dry_mass, *, component):
        if component not in self.recovery_fractions:
            raise ValueError(f"Explicit recovery required for {component!r}.")
        return scale_mass(self.contained_mass(dry_mass, component=component), self.recovery_fractions[component])

    def payable_mass(self, dry_mass, *, component):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        if component not in self.payable_fractions:
            raise ValueError(f"Explicit payability required for {component!r}.")
        recovered = mass(self.recovered_mass(dry_mass, component=component))
        deduction = mass(dry_mass) * self.deduction_fractions.get(component, Decimal(0))
        return _converters_petroleum_output(max(recovered - deduction, Decimal(0)) * self.payable_fractions[component], "kg", dry_mass)

    def feed_required(self, output_mass, *, component, stage="recovered"):
        from core.units.registry import ureg
        if stage not in {"contained", "recovered", "payable"}:
            raise ValueError("stage must be contained, recovered or payable.")
        fraction = mass(getattr(self, stage + "_mass")(ureg.Quantity(Decimal(1), "kg"), component=component))
        if fraction == 0:
            raise ValueError("Cannot invert a zero component output.")
        return scale_mass(output_mass, 1 / fraction)

    def contained_masses(self, dry_mass):
        return MappingProxyType({k: self.contained_mass(dry_mass, component=k) for k in self.composition.mass_fractions})


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class ImpurityPenalty(_safety_QuantitySnapshot):
    component: str
    threshold_fraction: Decimal
    charge_per_dry_tonne_per_pct: pint.Quantity

    def __post_init__(self):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        object.__setattr__(self, "component", _converters_petroleum_text(self.component, "impurity component"))
        object.__setattr__(self, "threshold_fraction", _converters_petroleum_fraction(self.threshold_fraction, "impurity threshold"))
        validate_registry_quantity(self.charge_per_dry_tonne_per_pct)
        _converters_petroleum_number(self.charge_per_dry_tonne_per_pct.magnitude, "impurity charge", minimum=0)
        value = self.charge_per_dry_tonne_per_pct
        object.__setattr__(self, "charge_per_dry_tonne_per_pct", ureg.Quantity(Decimal(str(value.magnitude)), value.units))


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class ConcentrateSettlement(_safety_QuantitySnapshot):
    payable_kg: Mapping[str, Decimal]
    gross_value: pint.Quantity
    deductions: pint.Quantity
    net_value: pint.Quantity





def settle_concentrate(engine, dry_mass, prices, *, currency: str, treatment_charge=None,
                       refining_charges=None, impurity_penalties=()):
    """Explicit contract inputs; no metal-price, charge or FX fixings assumed.

    Absolute grade deductions are fractions of dry feed after recovery, floored
    at zero, followed by payability. Refining charges apply to payable metal.
    Impurity charges are per dry tonne per excess percentage point.
    """
    from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    target = parse_unit(currency)
    if len(target.dimensionality) != 1 or next(iter(target.dimensionality.values())) != 1 or not next(iter(target.dimensionality)).startswith("[currency_"):
        raise ValueError("Settlement currency must be a standalone currency unit.")
    if not prices:
        raise ValueError("Supply metal prices explicitly.")
    missing_prices = engine.payable_fractions.keys() - prices.keys()
    if missing_prices:
        raise ValueError(f"Missing prices for declared payable components: {sorted(missing_prices)}.")
    payable = {}
    gross = Decimal(0)
    deduction = Decimal(0)
    for component, price in prices.items():
        amount = engine.payable_mass(dry_mass, component=component)
        value = _converters_petroleum_quantity_number(price, target / ureg.kg, "metal price", minimum=0)
        payable[component] = mass(amount)
        gross += payable[component] * value
    for component, charge in (refining_charges or {}).items():
        if component not in payable:
            raise ValueError("Refining charge requires a priced payable component.")
        deduction += payable[component] * _converters_petroleum_quantity_number(charge, target / ureg.kg, "refining charge", minimum=0)
    if treatment_charge is not None:
        deduction += mass(dry_mass) * _converters_petroleum_quantity_number(treatment_charge, target / ureg.kg, "treatment charge", minimum=0)
    for penalty in impurity_penalties:
        if not isinstance(penalty, ImpurityPenalty):
            raise TypeError("Impurity penalties must be ImpurityPenalty records.")
        excess_pct = max(engine.composition.fraction(penalty.component) - penalty.threshold_fraction, Decimal(0)) * 100
        deduction += mass(dry_mass) / 1000 * excess_pct * _converters_petroleum_quantity_number(penalty.charge_per_dry_tonne_per_pct, target / ureg.metric_ton, "impurity penalty", minimum=0)
    return ConcentrateSettlement(MappingProxyType(payable), ureg.Quantity(gross, target),
                                ureg.Quantity(deduction, target), ureg.Quantity(gross - deduction, target))


from core.units.validation import _safety_QuantitySnapshot
from core.units.validation import authorized_conversion


@dataclass(frozen=True, slots=True, kw_only=True)
class MaterialSpecification(_safety_QuantitySnapshot):
    name: str
    density: pint.Quantity | None = None
    moisture_fraction: Decimal | None = None
    composition: ElementComposition | None = None
    recovery_fractions: Mapping[str, Decimal] = field(default_factory=dict)
    payable_fractions: Mapping[str, Decimal] = field(default_factory=dict)
    deduction_fractions: Mapping[str, Decimal] = field(default_factory=dict)
    heat_content: Any = None
    reference_temperature: pint.Quantity | None = None
    reference_pressure: pint.Quantity | None = None

    def __post_init__(self):
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.assays.energy.petroleum import _converters_petroleum_pressure
        from core.units.assays.energy.petroleum import _converters_petroleum_temperature
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.composition import ElementComposition
        from core.units.registry import ureg
        object.__setattr__(self, "name", _converters_petroleum_text(self.name, "specification name"))
        if self.density is not None:
            object.__setattr__(self, "density", ureg.Quantity(_converters_petroleum_positive_quantity(self.density, "kg/m^3", "density"), "kg/m^3"))
        if self.moisture_fraction is not None:
            object.__setattr__(self, "moisture_fraction", _converters_petroleum_fraction(self.moisture_fraction, "moisture", allow_one=False))
        if self.composition is not None and not isinstance(self.composition, ElementComposition):
            raise TypeError("composition must be ElementComposition; declare normalization explicitly there.")
        for name in ("recovery_fractions", "payable_fractions", "deduction_fractions"):
            object.__setattr__(self, name, fractions(getattr(self, name), name))
        if self.composition is None and any((self.recovery_fractions, self.payable_fractions, self.deduction_fractions)):
            raise ValueError("Recovery and settlement data require an assayed composition.")
        if self.composition is not None:
            self.recovery_engine  # Validate keys without requiring unspecified rates.
        if self.reference_temperature is not None:
            object.__setattr__(self, "reference_temperature", _converters_petroleum_temperature(self.reference_temperature))
        if self.reference_pressure is not None:
            object.__setattr__(self, "reference_pressure", _converters_petroleum_pressure(self.reference_pressure))
        if self.heat_content is not None:
            from core.units.heat_content import HeatContentAssay
            if not isinstance(self.heat_content, HeatContentAssay):
                raise TypeError("heat_content must be HeatContentAssay.")

    @property
    def recovery_engine(self):
        from core.units.exceptions import MissingAssayError
        if self.composition is None:
            raise MissingAssayError("An explicit composition is required.")
        return ElementRecoveryEngine(self.composition, self.recovery_fractions, self.payable_fractions, self.deduction_fractions)

    def wet_to_dry_mass(self, quantity):
        from core.units.analytical_basis import wet_to_dry_mass
        from core.units.exceptions import MissingAssayError
        if self.moisture_fraction is None:
            raise MissingAssayError("Explicit moisture is required.")
        require_mass_basis(quantity, "wet")
        return wet_to_dry_mass(quantity, self.moisture_fraction, target_unit="kg")

    def dry_to_wet_mass(self, quantity):
        from core.units.analytical_basis import dry_to_wet_mass
        from core.units.exceptions import MissingAssayError
        if self.moisture_fraction is None:
            raise MissingAssayError("Explicit moisture is required.")
        require_mass_basis(quantity, "dry")
        return dry_to_wet_mass(quantity, self.moisture_fraction, target_unit="kg")

    def mass_at_moisture(self, quantity, target_moisture_fraction):
        from core.units.analytical_basis import convert_moisture_mass
        from core.units.exceptions import MissingAssayError
        if self.moisture_fraction is None:
            raise MissingAssayError("Explicit moisture is required.")
        return convert_moisture_mass(
            quantity,
            source_moisture=self.moisture_fraction,
            target_moisture=target_moisture_fraction,
            target_unit="kg",
        )

    def _dry_feed(self, quantity, mass_basis):
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        if mass_basis == "dry":
            require_mass_basis(quantity, "dry")
            return quantity
        if mass_basis == "received":
            return self.wet_to_dry_mass(quantity)
        raise IncompatibleAnalyticalBasisError("mass_basis must be dry or received.")

    def contained_mass(self, quantity, *, component, mass_basis="dry"):
        return self.recovery_engine.contained_mass(self._dry_feed(quantity, mass_basis), component=component)

    def recovered_mass(self, quantity, *, component, mass_basis="dry"):
        return self.recovery_engine.recovered_mass(self._dry_feed(quantity, mass_basis), component=component)

    def payable_mass(self, quantity, *, component, mass_basis="dry"):
        return self.recovery_engine.payable_mass(self._dry_feed(quantity, mass_basis), component=component)

    def feed_required(self, quantity, *, component, stage="recovered", mass_basis="dry"):
        from core.units.exceptions import IncompatibleAnalyticalBasisError
        dry = self.recovery_engine.feed_required(quantity, component=component, stage=stage)
        if mass_basis == "received":
            return self.dry_to_wet_mass(dry)
        if mass_basis != "dry":
            raise IncompatibleAnalyticalBasisError("mass_basis must be dry or received.")
        return dry

    def volume_to_mass(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.exceptions import MissingAssayError
        if self.density is None:
            raise MissingAssayError("Explicit density at the declared state is required.")
        volume = _converters_petroleum_quantity_number(quantity, "m^3", "volume", minimum=0)
        return _converters_petroleum_output(volume * self.density.magnitude, "kg", quantity)

    def mass_to_volume(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.exceptions import MissingAssayError
        if self.density is None:
            raise MissingAssayError("Explicit density at the declared state is required.")
        return _converters_petroleum_output(mass(quantity) / self.density.magnitude, "m^3", quantity)

    def as_context(self, *, name=None):
        context = pint.Context(name or self.name)
        if self.density is not None:
            context.add_transformation("[volume]", "[mass]", lambda _r, q: self.volume_to_mass(q))
            context.add_transformation("[mass]", "[volume]", lambda _r, q: self.mass_to_volume(q))
        return context

    def as_heat_content(self, **overrides):
        from core.units.exceptions import MissingAssayError
        from core.units.heat_content import HeatContentAssay
        if self.heat_content is None:
            if not {"gcv", "ncv"} & overrides.keys():
                raise MissingAssayError("Supply measured calorific data as heat_content=HeatContentAssay(...).")
            metadata = dict(name=self.name, density=self.density, reference_temperature=self.reference_temperature,
                            reference_pressure=self.reference_pressure)
            metadata.update(overrides)
            return HeatContentAssay(**metadata)
        from dataclasses import replace
        updates = {}
        if self.heat_content.density is None and self.density is not None:
            updates["density"] = self.density
        updates.update(overrides)
        return replace(self.heat_content, **updates)

    @authorized_conversion
    def convert_quantity(self, quantity, target_unit, *, basis=None):
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.conversions import _is_thermal_endpoint
        from core.units.conversions import _convert_contextual_quantity
        from core.units.conversions import _ore_mass_factor
        from core.units.conversions import convert_bulk_ore
        target = parse_unit(target_unit)

        def mass_basis(unit):
            try:
                return _ore_mass_factor(unit)[2]
            except ValueError:
                return None

        source_basis, target_basis = mass_basis(quantity.units), mass_basis(target)
        cross_dimension = quantity.dimensionality != target.dimensionality
        if source_basis in {"wet", "dry"} and target_basis in {"wet", "dry"} and source_basis != target_basis:
            return convert_bulk_ore(quantity, target, moisture_fraction=self.moisture_fraction)
        if cross_dimension and "dry" in {source_basis, target_basis}:
            raise ValueError("Dry settlement mass needs an explicit dry_to_wet_mass()/wet_to_dry_mass() step before a gross-density or energy conversion.")
        if ({source_basis, target_basis} & {"contained", "payable"}) and (
            cross_dimension or source_basis is not None and target_basis is not None and source_basis != target_basis
        ):
            raise ValueError("Select contained_mass/recovered_mass/payable_mass and the component explicitly.")
        if _is_thermal_endpoint(quantity.units) or _is_thermal_endpoint(target):
            if quantity.dimensionality != target.dimensionality:
                from core.units.heat_content import convert_heat_content
                return convert_heat_content(quantity, target, assay=self.as_heat_content(), basis=basis or "lhv")
        context = self.as_context()
        with ureg.context(context):
            return _convert_contextual_quantity(quantity, target, [context])


from core.units.validation import authorized_conversion


@dataclass(frozen=True, slots=True, kw_only=True)
class _agricultural_families_AgriculturalFamilyAssay(MaterialSpecification):
    commodity: str
    package_mass: pint.Quantity | None = None
    foreign_matter_fraction: Decimal | None = None
    processing: ProcessingSpecification | None = None
    _commodities: ClassVar[frozenset[str] | None] = None

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.processing import ProcessingSpecification
        from core.units.registry import ureg
        MaterialSpecification.__post_init__(self)
        commodity = _converters_petroleum_text(self.commodity, "commodity").lower().replace(" ", "-")
        if self._commodities is not None and commodity not in self._commodities:
            raise ValueError(f"Unsupported {type(self).__name__} commodity {commodity!r}; choose from {sorted(self._commodities)}.")
        object.__setattr__(self, "commodity", commodity)
        if self.package_mass is not None:
            object.__setattr__(self, "package_mass", ureg.Quantity(_converters_petroleum_positive_quantity(self.package_mass, "kg", "package mass"), "kg"))
        object.__setattr__(self, "foreign_matter_fraction", _agricultural_families_yield(self.foreign_matter_fraction, "foreign matter"))
        if self.processing is not None and not isinstance(self.processing, ProcessingSpecification):
            raise TypeError("processing must be ProcessingSpecification.")

    def clean_mass(self, quantity):
        if self.foreign_matter_fraction is None:
            raise ValueError("Declare foreign matter, including explicit zero, before cleaning.")
        return scale_mass(quantity, 1 - self.foreign_matter_fraction)

    def moisture_standardised_mass(self, quantity, *, target_moisture_fraction, handling_loss_fraction=Decimal(0)):
        """Preserve dry matter, then apply declared proportional handling loss."""
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        target = _converters_petroleum_fraction(target_moisture_fraction, "target moisture", allow_one=False)
        loss = _converters_petroleum_fraction(handling_loss_fraction, "handling loss")
        dry = self.wet_to_dry_mass(quantity)
        return scale_mass(dry, (1 - loss) / (1 - target))

    @authorized_conversion
    def moisture_adjusted_price(self, quantity, *, target_moisture_fraction, handling_loss_fraction=Decimal(0)):
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        from core.units.conversions import _ore_mass_factor
        from core.units.conversions import _numerator_currency_iso
        validate_registry_quantity(quantity)
        _numerator_currency_iso(quantity)
        _, exponent, mass_basis = _ore_mass_factor(quantity.units)
        if exponent != -1 or mass_basis not in {None, "wet"}:
            raise ValueError("Moisture price must have a received-mass denominator.")
        magnitude = _converters_petroleum_number(quantity.magnitude, "price")  # signed quotes are valid
        factor = self.moisture_standardised_mass(ureg.Quantity(Decimal(1), "kg"),
                    target_moisture_fraction=target_moisture_fraction,
                    handling_loss_fraction=handling_loss_fraction).magnitude
        if factor == 0:
            raise ValueError("Cannot price zero standardised output.")
        return _converters_petroleum_output(magnitude / factor, quantity.units, quantity)

    def mass_to_packages(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        if self.package_mass is None:
            raise ValueError("Declare the package's net mass.")
        return _converters_petroleum_output(mass(quantity) / self.package_mass.magnitude, "dimensionless", quantity)

    def packages_to_mass(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_number
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.registry import ureg
        if self.package_mass is None:
            raise ValueError("Declare the package's net mass.")
        if isinstance(quantity, pint.Quantity):
            count = _converters_petroleum_quantity_number(quantity, "dimensionless", "package count", minimum=0)
            source = quantity
        else:
            count = _converters_petroleum_number(quantity, "package count", minimum=0)
            source = ureg.Quantity(quantity, "dimensionless")
        return _converters_petroleum_output(count * self.package_mass.magnitude, "kg", source)

    def product_mass(self, quantity, *, product):
        if self.processing is None:
            raise ValueError("Declare product yields in ProcessingSpecification.")
        return self.processing.input_to_output(quantity, product, target_unit="kg")

    def feedstock_required(self, quantity, *, product):
        if self.processing is None:
            raise ValueError("Declare product yields in ProcessingSpecification.")
        return self.processing.output_to_input(quantity, product, target_unit="kg")

    def all_product_masses(self, quantity):
        if self.processing is None:
            raise ValueError("Declare product yields in ProcessingSpecification.")
        return self.processing.all_outputs(quantity, target_unit="kg")





_spec_engine_MODES = {
    "volume_to_mass", "mass_to_volume", "wet_to_dry_mass", "dry_to_wet_mass", "contained_mass",
    "recovered_mass", "payable_mass", "feed_required", "fine_mass", "gross_mass_from_fine",
    "metallic_iron_mass", "oxide_bound_iron_mass", "beneficiated_product_mass", "treo_mass",
    "oxide_mass", "element_mass", "ndpr_element_mass", "calcined_mass", "uncalcined_mass",
    "cement_from_clinker", "clinker_from_cement", "calcination_co2_mass", "calcined_solid_mass",
    "compacted_volume", "loose_volume", "polished_mass", "rough_mass", "black_mass_from_pack",
    "clean_mass", "pellet_mass", "recycled_feed_required", "active_mass", "solution_mass",
    "active_to_moles", "moles_to_active", "solution_to_moles", "diluent_required", "diluted_mass",
    "pure_solvent_mass", "sulphur_mass", "acid_equivalent", "co_product_mass", "water_consumed",
    "ammonia_mass", "purified_mass", "dissolved_mass", "recoverable_mass", "brine_required",
    "gas_volume_to_mass", "mass_to_gas_volume", "upgraded_volume", "fuel_energy",
    "combustible_dry_mass", "co2_mass", "delivered_hydrogen", "carrier_required", "delivered_energy",
    "stored_energy", "output_energy", "input_required", "underlying_quantity", "certificate_count",
    "salt_mass", "dissolved_solids_mass", "solvent_water_mass",
    "piece_mass", "coating_mass", "unrecovered_mass", "reported_to_element", "fine_troy_ounces",
    "commercial_component_mass", "component_mass", "dilution_water_required", "equivalent_volume",
    "moisture_standardised_mass", "moisture_adjusted_price", "mass_to_packages", "packages_to_mass",
    "product_mass", "feedstock_required", "split_mass", "whole_pulse_required", "mass_to_bushels",
    "bushels_to_mass", "refined_mass", "crude_oil_required", "oil_mass", "dry_meal_mass", "meal_mass",
    "starch_mass", "recovered_starch_mass", "sugar_mass", "recovered_sugar_mass", "solids_mass",
    "reconstituted_mass", "reconstituted_volume", "milk_solids_mass", "milk_fat_mass", "protein_mass",
    "fat_equivalent_mass", "nutrient_mass", "saleable_mass", "raw_meat_required", "made_tea_mass",
    "fresh_leaf_required", "cured_mass", "green_leaf_required", "processed_spice_mass", "raw_spice_required",
    "kernel_mass", "shell_mass", "inshell_mass_required", "sawn_volume", "roundwood_required",
    "solid_volume", "stacked_volume", "area_to_mass", "mass_to_area", "roll_mass", "recovered_fibre_mass",
    "recovered_paper_required", "collectable_mass", "ash_mass", "fibre_feedstock_required", "degummed_mass",
}


from core.units.validation import authorized_conversion


@authorized_conversion
def convert_specification(
    quantity,
    target_unit,
    *,
    assay,
    mode="auto",
    return_metadata: bool = False,
    strict: bool = False,
    uncertainties=None,
    covariance=None,
    correlations=None,
    **kwargs,
):
    from core.units.exceptions import UnsupportedConversionError
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.precision import ensure_quantity_ieee754

    validate_registry_quantity(quantity)
    if mode == "auto":
        from core.units.conversions import convert
        options = dict(kwargs)
        if return_metadata:
            options["return_metadata"] = True
        if strict:
            options["strict"] = True
        for key, value in (("uncertainties", uncertainties), ("covariance", covariance), ("correlations", correlations)):
            if value is not None:
                options[key] = value
        return convert(
            quantity,
            target_unit,
            assay=assay,
            **options,
        )
    if mode not in _spec_engine_MODES or not callable(getattr(assay, mode, None)):
        raise UnsupportedConversionError(f"Unsupported conversion mode {mode!r} for {type(assay).__name__}.")
    method = getattr(assay, mode)
    call_kwargs = dict(kwargs)
    if "target_unit" in inspect.signature(method).parameters or mode == "upgraded_volume":
        # Tagged reference volumes need the declared model to construct the
        # destination; a subsequent plain Pint .to() cannot supply that state.
        call_kwargs.setdefault("target_unit", target_unit)
    result = method(quantity, **call_kwargs)
    if not isinstance(result, pint.Quantity):
        raise TypeError("Select a quantity-returning conversion method.")
    converted = ensure_quantity_ieee754(result.to(parse_unit(target_unit)))
    if not return_metadata:
        return converted

    from core.units.provenance import ConversionClassification
    from core.units.provenance import build_conversion_result

    yield_modes = {"recovered_mass", "payable_mass", "feed_required", "beneficiated_product_mass", "clean_mass", "pellet_mass", "recycled_feed_required", "product_mass", "feedstock_required", "split_mass", "whole_pulse_required", "refined_mass", "crude_oil_required", "oil_mass", "dry_meal_mass", "meal_mass", "recovered_starch_mass", "recovered_sugar_mass", "saleable_mass", "raw_meat_required", "made_tea_mass", "fresh_leaf_required", "cured_mass", "green_leaf_required", "processed_spice_mass", "raw_spice_required", "kernel_mass", "shell_mass", "inshell_mass_required", "sawn_volume", "roundwood_required", "recovered_fibre_mass", "recovered_paper_required", "collectable_mass", "fibre_feedstock_required", "degummed_mass"}
    cls_tag = (
        ConversionClassification.YIELD_DEPENDENT
        if mode in yield_modes
        else ConversionClassification.ASSAY_DEPENDENT
    )
    return build_conversion_result(
        converted,
        quantity,
        assay=assay,
        conversion_path=(str(quantity.units), mode, str(converted.units)),
        method=f"{type(assay).__name__}.{mode}",
        classification=cls_tag,
        uncertainties=uncertainties,
        covariance=covariance,
        correlations=correlations,
        parameter_evaluator=lambda p: convert_specification(
            ureg.Quantity(p["quantity"], quantity.units),
            target_unit,
            assay=assay,
            mode=mode,
            **kwargs,
        ),
        base_parameters={"quantity": Decimal(str(quantity.magnitude))},
    )





def install_api(namespace, classes):
    """Give each concrete specification its typed named context/converter pair."""
    exports = []
    for cls in classes:
        stem = re.sub(r"(?<!^)(?=[A-Z])", "_", cls.__name__.removesuffix("Assay").removesuffix("Specification")).lower()

        def make_factory(expected):
            def create(assay, *, name=None):
                if not isinstance(assay, expected):
                    raise TypeError(f"Expected {expected.__name__}.")
                return assay.as_context(name=name)
            return create

        def make_converter(expected):
            def convert(quantity, target_unit, *, assay, mode="auto", **kwargs):
                if not isinstance(assay, expected):
                    raise TypeError(f"Expected {expected.__name__}.")
                return convert_specification(quantity, target_unit, assay=assay, mode=mode, **kwargs)
            return convert

        for name, function in (("create_" + stem + "_context", make_factory(cls)), ("convert_" + stem, make_converter(cls))):
            function.__name__ = function.__qualname__ = name
            function.__module__ = namespace["__name__"]
            function.__doc__ = f"{'Create a physical context' if name.startswith('create_') else 'Convert quantities with explicit processing modes'} for {cls.__name__}."
            namespace[name] = function
            exports.append(name)
        exports.append(cls.__name__)
    return tuple(exports)


__all__ = ('ConcentrateSettlement', 'ElementRecoveryEngine', 'ImpurityPenalty', 'MaterialSpec', 'MaterialSpecification', 'convert_specification', 'material_spec_from_assay')
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
    from core.units.composition import CompositionBasis
    from core.units.composition import ElementComposition
    from core.units.composition import MaterialComposition
    from core.units.processing import ProcessingSpecification
    from core.units.registry import Unit
