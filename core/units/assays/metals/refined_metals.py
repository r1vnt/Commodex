"""Commodex assays metals refined_metals implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class RefinedMetalAssay(MaterialSpecification):
    metal: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        MaterialSpecification.__post_init__(self)
        if self.metal not in {"Cu", "Zn", "Ni", "Pb", "Sn", "Al"}:
            raise ValueError("Refined metal must be Cu, Zn, Ni, Pb, Sn or Al.")
        if self.composition is None or self.composition.fraction(self.metal) <= 0:
            raise ValueError("Declare a positive assayed metal purity.")

    def fine_mass(self, quantity):
        return self.contained_mass(quantity, component=self.metal)

    def gross_mass_from_fine(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, 1 / self.composition.fraction(self.metal))


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (RefinedMetalAssay,))
__all__ = ('RefinedMetalAssay', 'convert_refined_metal', 'create_refined_metal_context')
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
    from core.units.assays.metals.refined_metals import convert_refined_metal
    from core.units.assays.metals.refined_metals import create_refined_metal_context
