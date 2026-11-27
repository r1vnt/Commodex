"""Commodex provenance implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field, fields, is_dataclass
from decimal import Context, Decimal, ROUND_HALF_EVEN
from enum import Enum
from functools import lru_cache
from hashlib import sha256
from types import MappingProxyType
from typing import Any, Final
import json
import math
import pint




_audit_HIGH_PREC_CTX: Final[Context] = Context(prec=50, rounding=ROUND_HALF_EVEN)





class ConversionClassification(str, Enum):
    """Reliability and epistemic classification of a unit/material conversion."""

    EXACT = "exact"
    MARKET_DEPENDENT = "market_dependent"
    ASSAY_DEPENDENT = "assay_dependent"
    COMPOSITION_DEPENDENT = "composition_dependent"
    YIELD_DEPENDENT = "estimated_or_measured_yield"
    STATE_AND_MODEL_DEPENDENT = "state_and_model_dependent"
    MOISTURE_AND_CALORIFIC_DEPENDENT = "moisture_and_calorific_dependent"





def _audit_dec(value: Any, label: str = "value") -> Decimal:
    if isinstance(value, bool) or not isinstance(value, (Decimal, int, float)):
        raise TypeError(f"{label} must be Decimal, int, or float, got {type(value).__name__}.")
    if isinstance(value, float) and not math.isfinite(value):
        raise ValueError(f"{label} must be finite, got {value!r}.")
    out = value if isinstance(value, Decimal) else Decimal(str(value))
    if not out.is_finite():
        raise ValueError(f"{label} must be finite, got {value!r}.")
    return out





def _audit_freeze_mapping(mapping: Mapping[str, Any] | None) -> Mapping[str, Any]:
    from core.units.registry import validate_registry_quantity
    from core.units.validation import _safety_copy_quantities
    if not mapping:
        return MappingProxyType({})
    frozen: dict[str, Any] = {}
    for key, val in mapping.items():
        if isinstance(val, pint.Quantity):
            validate_registry_quantity(val)
            frozen[str(key)] = _safety_copy_quantities(val)
        elif isinstance(val, Mapping):
            frozen[str(key)] = _audit_freeze_mapping(val)
        elif isinstance(val, list):
            frozen[str(key)] = tuple(_safety_copy_quantities(item) for item in val)
        else:
            frozen[str(key)] = _safety_copy_quantities(val)
    return MappingProxyType(frozen)





def _audit_serialize_audit_value(value: Any) -> Any:
    from core.units.precision import format_ieee754
    if isinstance(value, pint.Quantity):
        return {"magnitude": format_ieee754(value.magnitude), "unit": str(value.units)}
    if isinstance(value, Decimal):
        return format_ieee754(value)
    if isinstance(value, float):
        return format_ieee754(value)
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Mapping):
        return {str(k): _audit_serialize_audit_value(v) for k, v in sorted(value.items(), key=lambda kv: str(kv[0]))}
    if isinstance(value, (tuple, list)):
        return [_audit_serialize_audit_value(item) for item in value]
    if is_dataclass(value) and not isinstance(value, type):
        return {
            f.name: _audit_serialize_audit_value(getattr(value, f.name))
            for f in fields(value)
        }
    return value





def compute_assay_fingerprint(assay: Any) -> str:
    """Compute a deterministic SHA-256 digest of an assay's full parameter state."""
    payload = {
        "class": f"{type(assay).__module__}.{type(assay).__qualname__}",
        "state": _audit_serialize_audit_value(assay),
    }
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")
    return sha256(encoded).hexdigest()





@lru_cache(maxsize=1)
def _audit_reference_catalog_by_class() -> Mapping[str, tuple[tuple[Any, Any], ...]]:
    """Build a read-only lookup of pinned ReferenceRecords indexed by assay class name."""
    from core.units import agricultural_specs, conversions
    from core.units.standards import REFERENCE_VERSION
    from core.units.standards import reference_records

    categories = (
        "crude_oil",
        "refined_products",
        "natural_gas",
        "coal",
        "alternative_fuels",
        "metals",
        "chemicals",
        "agriculture",
    )
    class_lookup: dict[str, type[Any]] = {}
    for mod in (conversions, agricultural_specs):
        for name, obj in vars(mod).items():
            if isinstance(obj, type):
                class_lookup[name] = obj

    indexed: dict[str, list[tuple[Any, Any]]] = {}
    for category in categories:
        try:
            records = reference_records(category, version=REFERENCE_VERSION)
        except Exception:
            continue
        for rec in records:
            cls_name = rec.id.split(".", 1)[0]
            cls = class_lookup.get(cls_name)
            if cls is None:
                continue
            try:
                args, kwargs = rec.constructor_arguments()
                instance = cls(*args, **kwargs)
            except Exception:
                continue
            indexed.setdefault(cls_name, []).append((instance, rec))
    return MappingProxyType({k: tuple(v) for k, v in indexed.items()})


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class ConversionResult(_safety_QuantitySnapshot):
    """Complete, immutable audit trail for a physical or financial unit conversion."""

    quantity: Quantity

    original_quantity: Quantity
    commodity_id: str | None

    conversion_path: tuple[str, ...]
    assay_id: str | None
    assay_version: str | None

    reference_conditions: ReferenceConditions | None
    method: str
    standard: str | None
    standard_version: str | None

    registry_version: str
    registry_hash: str

    estimated: bool
    uncertainty: Decimal | None
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        from core.units.precision import pad_ieee754_decimal
        from core.units.registry import validate_registry_quantity
        from core.units.validation import _safety_copy_quantities
        if not isinstance(self.quantity, pint.Quantity):
            raise TypeError("ConversionResult.quantity must be a pint.Quantity.")
        if not isinstance(self.original_quantity, pint.Quantity):
            raise TypeError("ConversionResult.original_quantity must be a pint.Quantity.")
        validate_registry_quantity(self.quantity)
        validate_registry_quantity(self.original_quantity)
        object.__setattr__(self, "quantity", _safety_copy_quantities(self.quantity))
        object.__setattr__(self, "original_quantity", _safety_copy_quantities(self.original_quantity))
        object.__setattr__(self, "conversion_path", tuple(str(step) for step in self.conversion_path))
        object.__setattr__(self, "warnings", tuple(str(w) for w in self.warnings))
        if self.uncertainty is not None:
            u_dec = _audit_dec(self.uncertainty, "uncertainty")
            if u_dec < Decimal(0):
                raise ValueError(f"Conversion uncertainty must be non-negative, got {self.uncertainty!r}.")
            object.__setattr__(self, "uncertainty", pad_ieee754_decimal(u_dec))

    @property
    def classification(self) -> str:
        if self.reference_conditions is not None:
            return self.reference_conditions.classification
        return ConversionClassification.EXACT.value

    @property
    def uncertainty_quantity(self) -> Quantity | None:
        from core.units.registry import ureg
        if self.uncertainty is None:
            return None
        return ureg.Quantity(self.uncertainty, self.quantity.units)

    @property
    def relative_uncertainty(self) -> Decimal | None:
        from core.units.precision import pad_ieee754_decimal
        if self.uncertainty is None:
            return None
        mag = _audit_dec(self.quantity.magnitude, "quantity.magnitude")
        if mag.is_zero():
            return None
        return pad_ieee754_decimal(_audit_HIGH_PREC_CTX.divide(self.uncertainty, abs(mag)))

    def to_dict(self) -> dict[str, Any]:
        """Return a JSON-serializable audit dictionary reconstructing the conversion."""
        from core.units.precision import format_ieee754
        return {
            "quantity": _audit_serialize_audit_value(self.quantity),
            "original_quantity": _audit_serialize_audit_value(self.original_quantity),
            "commodity_id": self.commodity_id,
            "conversion_path": list(self.conversion_path),
            "assay_id": self.assay_id,
            "assay_version": self.assay_version,
            "reference_conditions": _audit_serialize_audit_value(self.reference_conditions),
            "classification": self.classification,
            "method": self.method,
            "standard": self.standard,
            "standard_version": self.standard_version,
            "registry_version": self.registry_version,
            "registry_hash": self.registry_hash,
            "estimated": self.estimated,
            "uncertainty": format_ieee754(self.uncertainty) if self.uncertainty is not None else None,
            "warnings": list(self.warnings),
        }





# ---------------------------------------------------------------------------
# Assay Provenance & Audit Builder
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AssayAuditInfo:
    """Resolved provenance, classification, and thermodynamic reference state for an assay."""

    commodity_id: str | None
    assay_id: str | None
    assay_version: str | None
    reference_conditions: ReferenceConditions
    standard: str | None
    standard_version: str | None
    estimated: bool
    warnings: tuple[str, ...] = ()





def inspect_assay_provenance(
    assay: Any,
    *,
    heating_basis: str | None = None,
    analytical_basis: str | None = None,
    extra_parameters: Mapping[str, Any] | None = None,
    classification_override: ConversionClassification | str | None = None,
) -> AssayAuditInfo:
    """Extract auditable identity, version, reference conditions, and reliability classification."""
    from core.units.reference_conditions import ReferenceConditions
    from core.units.registry import ureg
    if assay is None:
        cls_val = (
            classification_override.value
            if isinstance(classification_override, ConversionClassification)
            else (classification_override or ConversionClassification.EXACT.value)
        )
        return AssayAuditInfo(
            commodity_id=None,
            assay_id=None,
            assay_version=None,
            reference_conditions=ReferenceConditions(
                analytical_basis=analytical_basis,
                heating_basis=heating_basis,
                classification=cls_val,
                parameters=extra_parameters or {},
            ),
            standard=None,
            standard_version=None,
            estimated=False,
        )

    cls_name = type(assay).__name__
    fingerprint = compute_assay_fingerprint(assay)

    # Check pinned reference catalogue first
    matched_record = None
    for candidate, rec in _audit_reference_catalog_by_class().get(cls_name, ()):
        if candidate == assay:
            matched_record = rec
            break

    raw_name = getattr(assay, "name", None)
    if matched_record is not None:
        assay_id = matched_record.id
        assay_version = f"{matched_record.version}+sha256:{fingerprint[:16]}"
        estimated = matched_record.value_kind != "measured_value"
    else:
        assay_id = f"{cls_name}:{raw_name}" if raw_name else cls_name
        explicit_ver = getattr(assay, "assay_version", None) or getattr(assay, "version", None)
        assay_version = (
            f"{explicit_ver}+sha256:{fingerprint[:16]}"
            if explicit_ver
            else f"custom+sha256:{fingerprint[:16]}"
        )
        estimated = False

    temp: Quantity | None = None
    press: Quantity | None = None
    density: Quantity | None = None
    moisture: Decimal | None = None
    z_factor: Decimal | None = None
    a_basis = analytical_basis
    h_basis = heating_basis
    standard: str | None = None
    standard_version: str | None = None
    commodity_id: str | None = None
    classification = ConversionClassification.ASSAY_DEPENDENT
    params: dict[str, Any] = {"assay_class": cls_name, "assay_sha256": fingerprint}
    if matched_record is not None:
        params["reference_catalog_version"] = matched_record.version
        params["reference_value_kind"] = matched_record.value_kind
        params["reference_citation"] = matched_record.source.get("citation")
    warnings_list: list[str] = []

    # 1. CrudeAssay
    if cls_name == "CrudeAssay":
        grp = getattr(assay, "product_group", "crude")
        commodity_id = "crude_oil" if grp == "crude" else grp
        temp = ureg.Quantity(60.0, "degF")
        press = ureg.Quantity(101325.0, "Pa")
        density = ureg.Quantity(assay.density_kg_per_m3, "kg/m^3")
        standard = "API MPMS Ch. 11.1 / ASTM D1250"
        standard_version = "60F_101.325kPa"
        classification = ConversionClassification.ASSAY_DEPENDENT
        params.update({
            "name": assay.name,
            "api_gravity": assay.api_gravity,
            "specific_gravity_60F": assay.specific_gravity,
            "density_kg_per_m3": assay.density_kg_per_m3,
            "density_lb_per_gal": assay.density_lb_per_gal,
            "bbl_per_mt": assay.bbl_per_mt,
            "mt_per_bbl": assay.mt_per_bbl,
            "product_group": grp,
            "sulfur_percentage": assay.sulfur_percentage,
            "tan": assay.tan,
        })

    # 2. LiquidProductSpec
    elif cls_name == "LiquidProductSpec":
        commodity_id = str(getattr(assay.product, "value", assay.product))
        temp = assay.reference.temperature
        press = assay.reference.absolute_pressure
        density = assay.reference_density
        standard = "API MPMS Ch. 11.1"
        classification = ConversionClassification.ASSAY_DEPENDENT
        params["name"] = assay.name

    # 3. GasAssay & GasComposition
    elif cls_name == "GasAssay":
        commodity_id = "natural_gas"
        temp = ureg.Quantity(60.0, "degF")
        press = ureg.Quantity(14.73, "psia")
        classification = ConversionClassification.COMPOSITION_DEPENDENT
        standard = "ISO 6976 / GPA 2172"
        h_basis = h_basis or "gcv"
        params.update({
            "name": assay.name,
            "gcv_btu_per_scf": assay.gcv_btu_per_scf,
            "ncv_btu_per_scf": assay.ncv_btu_per_scf,
            "relative_density": assay.relative_density,
            "mmbtu_per_mt_lng": assay.mmbtu_per_mt_lng,
        })

    # 4. Solid / Liquid / Clean Fuels
    elif cls_name in {"BiomassAssay", "CoalAssay", "PetcokeAssay"}:
        commodity_id = cls_name.removesuffix("Assay").lower()
        classification = ConversionClassification.MOISTURE_AND_CALORIFIC_DEPENDENT
        moist_pct = getattr(assay, "moisture_pct", None) or getattr(assay, "moisture_percentage", None)
        if moist_pct is not None:
            moisture = _audit_dec(moist_pct, "moisture_pct") / Decimal(100)
        rho_val = getattr(assay, "bulk_density_kg_per_m3", None) or getattr(assay, "density_kg_per_m3", None)
        if rho_val is not None:
            density = ureg.Quantity(rho_val, "kg/m^3")
        a_basis = a_basis or "AR"
        h_basis = h_basis or "lhv"
        for field_info in fields(assay) if is_dataclass(assay) else ():
            params[field_info.name] = getattr(assay, field_info.name)

    elif cls_name in {"BiofuelAssay", "CleanFuelAssay", "NGLAssay", "HydrogenAssay"}:
        commodity_id = cls_name.removesuffix("Assay").lower()
        classification = (
            ConversionClassification.COMPOSITION_DEPENDENT
            if cls_name in {"CleanFuelAssay", "NGLAssay", "HydrogenAssay"}
            else ConversionClassification.ASSAY_DEPENDENT
        )
        rho_val = getattr(assay, "density_kg_per_m3", None)
        if rho_val is not None:
            density = ureg.Quantity(rho_val, "kg/m^3")
        h_basis = h_basis or ("hhv" if cls_name == "NGLAssay" else "lhv")
        for field_info in fields(assay) if is_dataclass(assay) else ():
            params[field_info.name] = getattr(assay, field_info.name)

    # 5. HeatContentAssay
    elif cls_name == "HeatContentAssay":
        commodity_id = getattr(assay, "name", None) or "fuel"
        a_basis = a_basis or getattr(assay, "basis", None)
        h_basis = h_basis or "lhv"
        temp = getattr(assay, "reference_temperature", None) or getattr(assay, "combustion_temperature", None)
        press = getattr(assay, "reference_pressure", None)
        density = getattr(assay, "density", None)
        moist_pct = getattr(assay, "moisture_pct", None)
        if moist_pct is not None:
            moisture = _audit_dec(moist_pct, "moisture_pct") / Decimal(100)
        if getattr(assay, "basis", None) in {"DRY_GAS", "WET_GAS"}:
            classification = ConversionClassification.COMPOSITION_DEPENDENT
            gas_state = assay.dry_gas_state if assay.basis == "DRY_GAS" else assay.wet_gas_state
            if gas_state is not None:
                temp = temp or gas_state.reference.temperature
                press = press or gas_state.reference.absolute_pressure
                z_factor = gas_state.reference.compressibility_factor
        else:
            classification = ConversionClassification.MOISTURE_AND_CALORIFIC_DEPENDENT
        standard = getattr(assay, "laboratory_standard", None) or getattr(assay, "method_source", None)
        cal_method = str(getattr(assay, "method", ""))
        if "latent_water" in cal_method or "explicit_ratio" in cal_method or "empirical" in cal_method:
            estimated = True
        params["calorific_method"] = cal_method
        if getattr(assay, "gcv", None) is not None:
            params["gcv"] = assay.gcv
        if getattr(assay, "ncv", None) is not None:
            params["ncv"] = assay.ncv

    # 6. Metals, Ores & Composable Material Building Blocks
    elif cls_name in {"OreAssay", "BulkOreAssay", "BatteryMaterialAssay", "BullionAssay", "UraniumAssay", "FerroalloyAssay", "MetalAlloyAssay"}:
        commodity_id = getattr(assay, "primary_element", None) or getattr(assay, "metal", None) or cls_name.removesuffix("Assay").lower()
        moist_frac = getattr(assay, "moisture_fraction", None)
        moist_pct = getattr(assay, "moisture_pct", None)
        if moist_frac is not None:
            moisture = _audit_dec(moist_frac, "moisture_fraction")
        elif moist_pct is not None:
            moisture = _audit_dec(moist_pct, "moisture_pct") / Decimal(100)
        rho_val = getattr(assay, "bulk_density_kg_per_m3", None) or getattr(assay, "density_kg_per_m3", None)
        if rho_val is not None:
            density = ureg.Quantity(rho_val, "kg/m^3")
        classification = ConversionClassification.ASSAY_DEPENDENT
        for field_info in fields(assay) if is_dataclass(assay) else ():
            params[field_info.name] = getattr(assay, field_info.name)

    elif cls_name in {"MaterialComposition", "MaterialProperties", "MaterialSpec"}:
        commodity_id = getattr(assay, "name", None) or cls_name.lower()
        if cls_name == "MaterialComposition":
            moisture = assay.moisture_fraction
            a_basis = a_basis or assay.basis.value
        elif cls_name == "MaterialProperties":
            density = assay.density or assay.bulk_density
            moisture = assay.moisture
        else:
            if assay.properties is not None:
                density = assay.properties.density or assay.properties.bulk_density
                moisture = assay.properties.moisture
            if assay.composition is not None:
                moisture = moisture if moisture is not None else assay.composition.moisture_fraction
                a_basis = a_basis or assay.composition.basis.value
        classification = ConversionClassification.ASSAY_DEPENDENT

    # 7. Agricultural & MaterialSpecification subclasses
    else:
        commodity_id = (
            getattr(assay, "commodity", None)
            or getattr(assay, "grain", None)
            or getattr(assay, "seed", None)
            or getattr(assay, "category", None)
            or getattr(assay, "species_group", None)
            or getattr(assay, "species", None)
            or getattr(assay, "variety", None)
            or getattr(assay, "wool_type", None)
            or getattr(assay, "feedstock", None)
            or getattr(assay, "product_type", None)
            or getattr(assay, "origin", None)
            or getattr(assay, "floral_source", None)
            or getattr(assay, "grade", None)
            or cls_name.removesuffix("Assay").lower()
        )
        if isinstance(commodity_id, str):
            commodity_id = commodity_id.lower()
        if isinstance(getattr(assay, "density", None), pint.Quantity):
            density = assay.density
        else:
            for rho_attr in ("density_kg_per_m3", "bulk_density_kg_per_m3", "oil_density_kg_per_m3", "milk_density_kg_per_m3", "ethanol_density_kg_per_m3", "seed_bulk_density_kg_per_m3"):
                rho_v = getattr(assay, rho_attr, None)
                if isinstance(rho_v, (int, float, Decimal)):
                    density = ureg.Quantity(rho_v, "kg/m^3")
                    break
        if isinstance(getattr(assay, "reference_temperature", None), pint.Quantity):
            temp = assay.reference_temperature
        elif isinstance(getattr(assay, "reference_temp_c", None), (int, float, Decimal)):
            temp = ureg.Quantity(getattr(assay, "reference_temp_c"), "degC")
        if isinstance(getattr(assay, "reference_pressure", None), pint.Quantity):
            press = assay.reference_pressure
        moist_f = getattr(assay, "moisture_fraction", None)
        moist_p = getattr(assay, "actual_moisture_pct", None)
        if moist_p is None:
            moist_p = getattr(assay, "moisture_pct", None)
        if moist_f is not None:
            moisture = _audit_dec(moist_f, "moisture_fraction")
        elif moist_p is not None:
            moisture = _audit_dec(moist_p, "moisture_pct") / Decimal(100)
        if is_dataclass(assay):
            for field_info in fields(assay):
                params[field_info.name] = getattr(assay, field_info.name)

    if extra_parameters:
        params.update(extra_parameters)

    final_cls = (
        classification_override.value
        if isinstance(classification_override, ConversionClassification)
        else (classification_override or classification.value)
    )
    ref_cond = ReferenceConditions(
        temperature=temp,
        pressure=press,
        density=density,
        moisture_fraction=moisture,
        compressibility_factor=z_factor,
        analytical_basis=a_basis,
        heating_basis=h_basis,
        classification=final_cls,
        parameters=params,
    )
    return AssayAuditInfo(
        commodity_id=commodity_id,
        assay_id=assay_id,
        assay_version=assay_version,
        reference_conditions=ref_cond,
        standard=standard,
        standard_version=standard_version,
        estimated=estimated,
        warnings=tuple(warnings_list),
    )





def build_conversion_result(
    quantity: pint.Quantity,
    original_quantity: pint.Quantity,
    *,
    assay: Any = None,
    commodity_id: str | None = None,
    conversion_path: Sequence[str] | None = None,
    method: str = "pint_unit_registry",
    standard: str | None = None,
    standard_version: str | None = None,
    classification: ConversionClassification | str | None = None,
    heating_basis: str | None = None,
    analytical_basis: str | None = None,
    reference_conditions: ReferenceConditions | None = None,
    extra_parameters: Mapping[str, Any] | None = None,
    estimated: bool | None = None,
    uncertainties: Mapping[str, Any] | Decimal | float | pint.Quantity | None = None,
    covariance: Mapping[tuple[str, str], Decimal | float | int] | Sequence[Sequence[Decimal | float | int]] | None = None,
    correlations: Mapping[tuple[str, str], Decimal | float | int] | Sequence[Sequence[Decimal | float | int]] | None = None,
    parameter_evaluator: Callable[[Mapping[str, Decimal]], Decimal | float | pint.Quantity] | None = None,
    base_parameters: Mapping[str, Decimal | float | int] | None = None,
    warnings: Sequence[str] = (),
) -> ConversionResult:
    """Assemble an immutable ``ConversionResult`` with optional uncertainty propagation."""
    from core.units.precision import pad_ieee754_decimal
    from core.units.registry import REGISTRY_VERSION
    from core.units.registry import get_registry_definition_hash
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    from core.units.uncertainty import propagate_callable_uncertainty
    info = inspect_assay_provenance(
        assay,
        heating_basis=heating_basis,
        analytical_basis=analytical_basis,
        extra_parameters=extra_parameters,
        classification_override=classification,
    )
    ref = reference_conditions if reference_conditions is not None else info.reference_conditions
    eff_cls = ref.classification
    eff_estimated = info.estimated if estimated is None else bool(estimated)
    if eff_cls == ConversionClassification.EXACT.value and estimated is None:
        eff_estimated = False

    path = (
        tuple(conversion_path)
        if conversion_path is not None
        else (str(original_quantity.units), str(quantity.units))
    )
    all_warnings = list(info.warnings) + list(warnings)

    # Propagate uncertainty
    propagated_u: Decimal | None = None
    if uncertainties is None:
        if eff_cls == ConversionClassification.EXACT.value:
            propagated_u = pad_ieee754_decimal(Decimal(0))
        else:
            propagated_u = None
    elif isinstance(uncertainties, (Decimal, int, float, pint.Quantity)) and not isinstance(uncertainties, bool):
        if isinstance(uncertainties, pint.Quantity):
            validate_registry_quantity(uncertainties)
            sigma_in = _audit_dec(uncertainties.to(original_quantity.units).magnitude, "uncertainty")
        else:
            sigma_in = _audit_dec(uncertainties, "uncertainty")
        if sigma_in < Decimal(0):
            raise ValueError("Uncertainty must be non-negative.")
        x_mag = _audit_dec(original_quantity.magnitude, "original_quantity")
        y_mag = _audit_dec(quantity.magnitude, "quantity")
        if not x_mag.is_zero():
            dydx = abs(_audit_HIGH_PREC_CTX.divide(y_mag, x_mag))
            propagated_u = pad_ieee754_decimal(_audit_HIGH_PREC_CTX.multiply(dydx, sigma_in))
        elif parameter_evaluator is not None and base_parameters and "quantity" in base_parameters:
            res = propagate_callable_uncertainty(
                parameter_evaluator,
                base_parameters,
                {"quantity": sigma_in},
            )
            propagated_u = res.standard_uncertainty
        else:
            propagated_u = pad_ieee754_decimal(Decimal(0))
    elif isinstance(uncertainties, Mapping):
        clean_sigmas: dict[str, Decimal] = {}
        for key, raw_u in uncertainties.items():
            if key in {"quantity", "input"} and isinstance(raw_u, pint.Quantity):
                validate_registry_quantity(raw_u)
                clean_sigmas["quantity"] = _audit_dec(raw_u.to(original_quantity.units).magnitude, f"uncertainties[{key}]")
            elif isinstance(raw_u, pint.Quantity):
                validate_registry_quantity(raw_u)
                clean_sigmas[str(key)] = _audit_dec(raw_u.magnitude, f"uncertainties[{key}]")
            else:
                norm_key = "quantity" if key == "input" else str(key)
                clean_sigmas[norm_key] = _audit_dec(raw_u, f"uncertainties[{key}]")

        if set(clean_sigmas) <= {"quantity", "relative"} and covariance is None and correlations is None:
            x_mag = _audit_dec(original_quantity.magnitude, "original_quantity")
            y_mag = _audit_dec(quantity.magnitude, "quantity")
            var_total = Decimal(0)
            if "quantity" in clean_sigmas:
                if clean_sigmas["quantity"] < Decimal(0):
                    raise ValueError("Standard uncertainty for 'quantity' must be non-negative.")
                dydx = _audit_HIGH_PREC_CTX.divide(y_mag, x_mag) if not x_mag.is_zero() else Decimal(0)
                term_q = _audit_HIGH_PREC_CTX.multiply(dydx * dydx, clean_sigmas["quantity"] ** 2)
                var_total = _audit_HIGH_PREC_CTX.add(var_total, term_q)
            if "relative" in clean_sigmas:
                if clean_sigmas["relative"] < Decimal(0):
                    raise ValueError("Relative uncertainty must be non-negative.")
                term_r = _audit_HIGH_PREC_CTX.multiply(y_mag * y_mag, clean_sigmas["relative"] ** 2)
                var_total = _audit_HIGH_PREC_CTX.add(var_total, term_r)
            propagated_u = pad_ieee754_decimal(_audit_HIGH_PREC_CTX.sqrt(var_total))
        elif parameter_evaluator is not None and base_parameters is not None:
            merged_base = dict(base_parameters)
            res = propagate_callable_uncertainty(
                parameter_evaluator,
                merged_base,
                clean_sigmas,
                covariance=covariance,
                correlations=correlations,
            )
            propagated_u = res.standard_uncertainty
        else:
            raise ValueError(
                f"Cannot propagate parameter uncertainties {sorted(clean_sigmas)} without a parameter evaluator."
            )
    else:
        raise TypeError(f"Unsupported uncertainties specification: {type(uncertainties).__name__}.")

    return ConversionResult(
        quantity=quantity,
        original_quantity=original_quantity,
        commodity_id=commodity_id if commodity_id is not None else info.commodity_id,
        conversion_path=path,
        assay_id=info.assay_id,
        assay_version=info.assay_version,
        reference_conditions=ref,
        method=method,
        standard=standard if standard is not None else info.standard,
        standard_version=standard_version if standard_version is not None else info.standard_version,
        registry_version=REGISTRY_VERSION,
        registry_hash=get_registry_definition_hash(ureg),
        estimated=eff_estimated,
        uncertainty=propagated_u,
        warnings=tuple(all_warnings),
    )


__all__ = ('AssayAuditInfo', 'ConversionClassification', 'ConversionResult', 'REFERENCE_DATA_HASHES', 'REFERENCE_VERSION', 'ReferenceRecord', 'build_conversion_result', 'compute_assay_fingerprint', 'get_registry_snapshot', 'inspect_assay_provenance', 'reference_record', 'reference_records')
_PUBLIC_ALIASES = {'reference_records': ('standards', 'reference_records'), 'REFERENCE_DATA_HASHES': ('registry', 'REFERENCE_DATA_HASHES'), 'get_registry_snapshot': ('registry', 'get_registry_snapshot'), 'ReferenceRecord': ('standards', 'ReferenceRecord'), 'REFERENCE_VERSION': ('standards', 'REFERENCE_VERSION'), 'reference_record': ('standards', 'reference_record')}

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
    from core.units.registry import REFERENCE_DATA_HASHES
    from core.units.standards import REFERENCE_VERSION
    from core.units.reference_conditions import ReferenceConditions
    from core.units.standards import ReferenceRecord
    from core.units.registry import get_registry_snapshot
    from core.units.standards import reference_record
    from core.units.standards import reference_records
