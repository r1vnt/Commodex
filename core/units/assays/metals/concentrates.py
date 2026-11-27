"""Commodex assays metals concentrates implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class MultiMetalConcentrateAssay(MaterialSpecification):
    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        MaterialSpecification.__post_init__(self)
        if self.composition is None:
            raise ValueError("Concentrates require assayed components.")

    def settle(self, quantity, prices, *, currency, mass_basis="dry", **charges):
        from core.units.assays.base import settle_concentrate
        return settle_concentrate(self.recovery_engine, self._dry_feed(quantity, mass_basis), prices, currency=currency, **charges)

    def assayed_masses(self, quantity, *, mass_basis="dry"):
        return self.recovery_engine.contained_masses(self._dry_feed(quantity, mass_basis))





@dataclass(frozen=True, slots=True, kw_only=True)
class PreciousMetalConcentrateAssay(MultiMetalConcentrateAssay):
    def fine_troy_ounces(self, quantity, *, component, mass_basis="dry"):
        if component not in {"Au", "Ag", "Pt"}:
            raise ValueError("Select Au, Ag or Pt.")
        return self.contained_mass(quantity, component=component, mass_basis=mass_basis).to("troy_ounce")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (MultiMetalConcentrateAssay, PreciousMetalConcentrateAssay,))
__all__ = ('MultiMetalConcentrateAssay', 'PreciousMetalConcentrateAssay', 'convert_multi_metal_concentrate', 'convert_precious_metal_concentrate', 'create_multi_metal_concentrate_context', 'create_precious_metal_concentrate_context')
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
    from core.units.assays.metals.concentrates import convert_multi_metal_concentrate
    from core.units.assays.metals.concentrates import convert_precious_metal_concentrate
    from core.units.assays.metals.concentrates import create_multi_metal_concentrate_context
    from core.units.assays.metals.concentrates import create_precious_metal_concentrate_context
