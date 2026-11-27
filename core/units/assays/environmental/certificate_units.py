"""Commodex assays environmental certificate_units implementation and public interfaces."""
from __future__ import annotations

from dataclasses import dataclass
from decimal import Decimal
import pint

from core.units.validation import _safety_QuantitySnapshot
from core.units.validation import authorized_conversion


@dataclass(frozen=True, slots=True, kw_only=True)
class EnvironmentalUnitSpecification(_safety_QuantitySnapshot):
    scheme: str
    vintage: str
    certificate_type: str
    quantity_per_certificate: pint.Quantity
    issuer: str

    def __post_init__(self):
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        for name in ("scheme", "vintage", "certificate_type", "issuer"):
            object.__setattr__(self, name, _converters_petroleum_text(getattr(self, name), name))
        validate_registry_quantity(self.quantity_per_certificate)
        dims = self.quantity_per_certificate.dimensionality
        qualified_instrument = len(dims) == 1 and next(iter(dims.values())) == 1 and (
            next(iter(dims)).startswith("[instrument_") or next(iter(dims)) == "[renewable_identification_number]"
        )
        if dims not in {ureg.MWh.dimensionality, ureg.kg.dimensionality, (ureg.meter**3).dimensionality} and not qualified_instrument:
            raise ValueError("Certificate denomination must declare energy, mass, volume or a qualified environmental instrument.")
        _converters_petroleum_positive_quantity(self.quantity_per_certificate, self.quantity_per_certificate.units, "certificate denomination")
        value = self.quantity_per_certificate
        object.__setattr__(self, "quantity_per_certificate", ureg.Quantity(Decimal(str(value.magnitude)), value.units))

    @property
    def name(self):
        return f"{self.scheme}:{self.certificate_type}:{self.vintage}:{self.issuer}"

    def underlying_quantity(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.registry import ureg
        if isinstance(quantity, pint.Quantity):
            count = _converters_petroleum_quantity_number(quantity, "dimensionless", "certificate count", minimum=0)
            source = quantity
        else:
            from core.units.assays.energy.petroleum import _converters_petroleum_number as _number
            count = _number(quantity, "certificate count", minimum=0)
            source = ureg.Quantity(quantity, "dimensionless")
        value = count * Decimal(str(self.quantity_per_certificate.magnitude))
        return _converters_petroleum_output(value, self.quantity_per_certificate.units, source)

    def certificate_count(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        underlying = _converters_petroleum_quantity_number(quantity, self.quantity_per_certificate.units, "underlying quantity", minimum=0)
        return _converters_petroleum_output(underlying / Decimal(str(self.quantity_per_certificate.magnitude)), "dimensionless", quantity)

    def transfer_to(self, quantity, target):
        from core.units.exceptions import IncompatibleCommodityError
        if not isinstance(target, EnvironmentalUnitSpecification):
            raise TypeError("Target must be EnvironmentalUnitSpecification.")
        if any(getattr(self, k) != getattr(target, k) for k in ("scheme", "vintage", "certificate_type", "issuer")):
            raise IncompatibleCommodityError("Different certificate identities/vintages are not interchangeable.")
        return target.certificate_count(self.underlying_quantity(quantity))

    def as_context(self, *, name=None):
        context = pint.Context(name or self.name)
        context.add_transformation("[]", self.quantity_per_certificate.dimensionality, lambda _r, q: self.underlying_quantity(q))
        context.add_transformation(self.quantity_per_certificate.dimensionality, "[]", lambda _r, q: self.certificate_count(q))
        return context

    @authorized_conversion
    def convert_quantity(self, quantity, target_unit, *, basis=None):
        from core.units.registry import parse_unit
        from core.units.registry import ureg
        from core.units.conversions import _convert_contextual_quantity
        context = self.as_context()
        with ureg.context(context):
            return _convert_contextual_quantity(quantity, parse_unit(target_unit), [context])


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (EnvironmentalUnitSpecification,))
__all__ = ('EnvironmentalUnitSpecification', 'convert_environmental_unit', 'create_environmental_unit_context')
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
    from core.units.assays.environmental.certificate_units import convert_environmental_unit
    from core.units.assays.environmental.certificate_units import create_environmental_unit_context
