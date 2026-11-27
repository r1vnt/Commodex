"""Commodex context implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from decimal import Decimal
from enum import Enum
from itertools import permutations
from types import MappingProxyType
from typing import Any, Final, Literal
import math
import pint




def create_material_context(properties: MaterialProperties, *, use_bulk_density: bool = False, active_only: bool = False, name: str | None = None,) -> pint.Context:
    
    """Build a Pint Context bridging [volume], [mass], and [substance] from ``MaterialProperties``."""
    from core.units.composition import MaterialProperties
    from core.units.registry import ureg
    if not isinstance(properties, MaterialProperties):
        raise TypeError(f"Expected MaterialProperties, got {type(properties).__name__}.")

    ctx = pint.Context(name or "material_properties")
    nodes: list[tuple[pint.Unit, Decimal]] = [(ureg.kilogram, Decimal(1))]

    try:
        rho = properties._select_density(use_bulk=use_bulk_density).to("kg / meter ** 3").magnitude
        nodes.append((ureg.meter ** 3, Decimal(1) / rho))
    except ValueError:
        pass

    if properties.molar_mass is not None:
        kg_per_mol = properties.molar_mass.to("kg / mole").magnitude
        if active_only:
            if properties.purity is None:
                raise ValueError("purity is required on MaterialProperties when active_only=True.")
            moles_per_kg = properties.purity / kg_per_mol
        else:
            moles_per_kg = Decimal(1) / kg_per_mol
        nodes.append((ureg.mole, moles_per_kg))

    def _make_edge(src_u: pint.Unit, dst_u: pint.Unit, ratio: Decimal):
        flt_ratio = float(ratio)

        def _apply(ureg_inst: pint.UnitRegistry, qty: pint.Quantity) -> pint.Quantity:
            q_src = qty.to(src_u)
            mult = ratio if isinstance(q_src.magnitude, Decimal) else flt_ratio
            return ureg_inst.Quantity(q_src.magnitude * mult, dst_u)

        return _apply

    # Each node stores (unit, units_of_node_per_1_kg); so 1 src_u = (1 / src_per_kg) * dst_per_kg dst_u
    for (src_unit, src_per_kg), (dst_unit, dst_per_kg) in permutations(nodes, 2):
        ctx.add_transformation(src_unit.dimensionality, dst_unit.dimensionality, _make_edge(src_unit, dst_unit, dst_per_kg / src_per_kg),)

    return ctx

__all__ = ('create_fx_context', 'create_gas_heating_context', 'create_heat_content_context', 'create_hydrocarbon_context', 'create_material_context')
_PUBLIC_ALIASES = {'create_heat_content_context': ('heat_content', 'create_heat_content_context'), 'create_gas_heating_context': ('conversions', 'create_gas_heating_context'), 'create_hydrocarbon_context': ('conversions', 'create_hydrocarbon_context'), 'create_fx_context': ('conversions', 'create_fx_context')}

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
    from core.units.composition import MaterialProperties
    from core.units.conversions import create_fx_context
    from core.units.conversions import create_gas_heating_context
    from core.units.heat_content import create_heat_content_context
    from core.units.conversions import create_hydrocarbon_context
