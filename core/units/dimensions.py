"""Commodex dimensions implementation and public interfaces."""
from __future__ import annotations

from decimal import Decimal
from types import MappingProxyType
import pint

def energy_roles(unit):
    """Signed thermal/electrical exponents, including quote denominators."""
    return MappingProxyType(roles(unit))

_energy_roles_MARKERS = {"thermal": "thermal_energy_kind", "electrical": "electrical_energy_kind"}

def roles(unit):
    from core.units.registry import parse_unit
    dimensions = parse_unit(unit).dimensionality
    return {role: dimensions.get(f"[{marker}]", 0) for role, marker in _energy_roles_MARKERS.items()
            if dimensions.get(f"[{marker}]", 0)}

def physical_unit(unit, *, role):
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    parsed = parse_unit(unit)
    signature = roles(parsed)
    if signature and (set(signature) != {role} or abs(signature[role]) != 1):
        raise pint.DimensionalityError(parsed, f"{role} energy endpoint", extra_msg="; use an explicit generation, heat-pump or electrolysis process")
    exponent = signature.get(role, 0)
    return parsed / ureg.Unit(_energy_roles_MARKERS[role]) ** exponent

def physical_quantity(quantity, *, role):
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    if not isinstance(quantity, pint.Quantity):
        raise TypeError("Expected a Pint quantity.")
    validate_registry_quantity(quantity)
    return ureg.Quantity(quantity.magnitude, physical_unit(quantity.units, role=role))

def restore(quantity, target, *, role):
    from core.units.registry import parse_unit
    from core.units.registry import ureg
    target = parse_unit(target)
    physical_unit(target, role=role)  # Reject a conflicting output role.
    if roles(quantity.units):
        raise ValueError("Process output must be physical arithmetic before role attachment.")
    exponent = roles(target).get(role, 0)
    value = ureg.Quantity(quantity.magnitude, quantity.units * ureg.Unit(_energy_roles_MARKERS[role]) ** exponent)
    return value.to(target)

def typed_unit(quantity, role):
    from core.units.registry import ureg
    if quantity.dimensionality == ureg.joule.dimensionality:
        return ureg.Unit("MWh_th" if role == "thermal" else "MWh_e")
    if quantity.dimensionality == ureg.watt.dimensionality:
        return ureg.Unit("MW_th" if role == "thermal" else "MW_e")
    raise ValueError("Expected physical energy or power.")

def physical_heat_rate(quantity):
    """Accept a legacy physical ratio or typed thermal-input/electrical-output ratio."""
    from core.units.registry import ureg
    from core.units.registry import validate_registry_quantity
    validate_registry_quantity(quantity)
    signature = roles(quantity.units)
    if signature and signature != {"thermal": 1, "electrical": -1}:
        raise ValueError("Heat rate must be thermal input per electrical output.")
    unit = quantity.units
    if signature:
        unit = unit / ureg.thermal_energy_kind * ureg.electrical_energy_kind
    value = ureg.Quantity(Decimal(str(quantity.magnitude)), unit)
    return value.to(ureg.dimensionless)

__all__ = ('energy_roles',)
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
