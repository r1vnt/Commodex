"""Commodex assays environmental carbon_dioxide implementation and public interfaces."""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class CarbonDioxideAssay(MaterialSpecification):
    state: str
    purity_fraction: Decimal

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        MaterialSpecification.__post_init__(self)
        if self.state not in {"gas", "liquid", "dense-phase"}:
            raise ValueError("CO2 state must be gas, liquid or dense-phase.")
        object.__setattr__(self, "purity_fraction", _converters_petroleum_fraction(self.purity_fraction, "CO2 mass purity"))
        if self.density is None or self.reference_temperature is None or self.reference_pressure is None:
            raise ValueError("CO2 requires density and explicit temperature/absolute pressure for its state.")

    def co2_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, self.purity_fraction)

    def equivalent_volume(self, quantity, *, target):
        from core.units.exceptions import IncompatibleCommodityError
        if not isinstance(target, CarbonDioxideAssay):
            raise TypeError("Target must be CarbonDioxideAssay.")
        if target.purity_fraction != self.purity_fraction:
            raise IncompatibleCommodityError("State equivalence conserves the same composition; purity changes require a separation balance.")
        return target.mass_to_volume(self.volume_to_mass(quantity))


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (CarbonDioxideAssay,))
__all__ = ('CarbonDioxideAssay', 'convert_carbon_dioxide', 'create_carbon_dioxide_context')
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
    from core.units.assays.environmental.carbon_dioxide import convert_carbon_dioxide
    from core.units.assays.environmental.carbon_dioxide import create_carbon_dioxide_context
