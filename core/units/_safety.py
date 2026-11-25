"""Internal, thread-local authorization for validated semantic conversions."""

from contextlib import contextmanager
from contextvars import ContextVar
from dataclasses import fields, is_dataclass
from decimal import Decimal
from functools import wraps
import math
from collections.abc import Mapping


_SEMANTIC_SCOPE = ContextVar("unit_semantic_scope", default=False)


@contextmanager
def semantic_conversion():
    token = _SEMANTIC_SCOPE.set(True)
    try:
        yield
    finally:
        _SEMANTIC_SCOPE.reset(token)


def authorized_conversion(function):
    @wraps(function)
    def call(*args, **kwargs):
        with semantic_conversion():
            return function(*args, **kwargs)
    return call


def _finite_tree(value, label):
    if isinstance(value, (float, Decimal)) and not math.isfinite(value):
        raise ValueError(f"{label} must be finite.")
    if isinstance(value, Mapping):
        for key, item in value.items():
            _finite_tree(item, f"{label}.{key}")
    elif isinstance(value, (tuple, list)):
        for item in value:
            _finite_tree(item, label)


def validated_assay(cls):
    """Reject nonfinite fields and keep explicit assay methods authorized."""
    initialize = cls.__init__

    @wraps(initialize)
    def checked_init(self, *args, **kwargs):
        for item in args:
            _finite_tree(item, cls.__name__)
        _finite_tree(kwargs, cls.__name__)
        initialize(self, *args, **kwargs)
        if is_dataclass(self):
            for field in fields(self):
                value = getattr(self, field.name)
                _finite_tree(value, field.name)
                if value is not None and isinstance(value, (int, float, Decimal)):
                    if field.name.endswith(("_pct", "_percentage")) and "_per_" not in field.name and not 0 <= value <= 100:
                        raise ValueError(f"{field.name} must be in [0, 100].")
        # Independent co-products cannot exceed their common feedstock mass.
        for names in (
            ("oil_yield_pct", "meal_yield_pct"),
            ("oil_yield_pct", "fishmeal_yield_pct"),
            ("processing_yield_pct", "secondary_yield_pct"),
            ("lint_turnout_pct", "cottonseed_turnout_pct"),
            ("fat_pct", "protein_pct"),
        ):
            values = [getattr(self, name, None) for name in names]
            if all(value is not None for value in values) and sum(values) > 100 + 1e-10:
                raise ValueError(f"Combined mass yields {names} cannot exceed 100%.")
        for child, parent in (("edible_yield_pct", "dressed_yield_pct"), ("boneless_meat_pct_of_live", "rtc_dressing_pct")):
            child_value, parent_value = getattr(self, child, None), getattr(self, parent, None)
            if child_value is not None and parent_value is not None and child_value > parent_value:
                raise ValueError(f"{child} cannot exceed {parent} on the same feedstock basis.")

    cls.__init__ = checked_init
    for name, member in tuple(vars(cls).items()):
        if name.startswith("_"):
            continue
        if callable(member) and not isinstance(member, (classmethod, staticmethod, property)):
            setattr(cls, name, _checked_method(member, name))
    return cls


def _checked_method(member, name):
    yield_methods = {
        "raw_to_processed_mass", "seed_mass_to_oil_mass", "seed_mass_to_meal_mass",
        "whole_to_dressed", "whole_to_edible", "live_to_carcass_mass",
        "live_to_rtc_mass", "live_to_boneless_mass", "gross_to_fine",
        "gross_to_payable_fine", "milk_to_solids_kg", "shell_to_liquid_egg_mass",
        "shell_to_dried_egg_mass", "wmt_to_dmt", "dmt_to_contained_metal",
    }

    @wraps(member)
    def checked(self, *args, **kwargs):
        _finite_tree(args, name)
        _finite_tree(kwargs, name)
        quantities = [value for value in args if hasattr(value, "magnitude") and hasattr(value, "dimensionality")]
        def is_mass(quantity):
            return len(quantity.dimensionality) == 1 and next(iter(quantity.dimensionality.values())) == 1 and next(iter(quantity.dimensionality)) in {"[mass]", "[wet_bulk_mass]", "[dry_bulk_mass]", "[contained_metal_mass]", "[payable_metal_mass]"}
        expected_basis = next((basis for prefix, basis in (("wmt_to_", "wet"), ("dmt_to_", "dry"), ("contained_metal_to_", "contained"), ("payable_metal_to_", "payable")) if name.startswith(prefix)), None)
        for quantity in quantities:
            _finite_tree(quantity.magnitude, name)
            if expected_basis is not None:
                tags = set(quantity._REGISTRY._semantic_tags(quantity.units._units))
                if tags and tags != {expected_basis}:
                    raise ValueError(f"{name} requires {expected_basis}-basis input, not {tags}.")
            if quantity.magnitude < 0 and is_mass(quantity):
                raise ValueError(f"{name} requires non-negative physical mass.")
        with semantic_conversion():
            result = member(self, *args, **kwargs)
            if hasattr(result, "magnitude"):
                _finite_tree(result.magnitude, name)
                if name in yield_methods and quantities and is_mass(result) and is_mass(quantities[0]):
                    registry = result._REGISTRY
                    value = registry._convert_mass_scale(result.magnitude, result.units._units, registry.kg._units)
                    source = registry._convert_mass_scale(quantities[0].magnitude, quantities[0].units._units, registry.kg._units)
                    if value < 0 or float(value) > float(source) + max(1e-10, abs(float(source)) * 1e-12):
                        raise ValueError(f"{name} violates feedstock mass balance.")
            else:
                _finite_tree(result, name)
            return result
    return checked
