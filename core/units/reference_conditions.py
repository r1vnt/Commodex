"""Commodex reference_conditions implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from collections.abc import Mapping
from dataclasses import dataclass, field
from dataclasses import dataclass, field, fields, is_dataclass
from decimal import Context, Decimal, ROUND_HALF_EVEN
from decimal import Decimal, ROUND_HALF_EVEN, getcontext
from enum import Enum
from functools import lru_cache
from hashlib import sha256
from types import MappingProxyType
from typing import Any, Final
from typing import TYPE_CHECKING, Any, Final
import inspect
import json
import math
import pint
import warnings

from core.units.provenance import ConversionClassification
from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class ReferenceConditions(_safety_QuantitySnapshot):
    """Thermodynamic, analytical, and physical reference state of a conversion."""

    temperature: Quantity | None = None
    pressure: Quantity | None = None
    density: Quantity | None = None
    moisture_fraction: Decimal | None = None
    compressibility_factor: Decimal | None = None
    analytical_basis: str | None = None
    heating_basis: str | None = None
    classification: str = ConversionClassification.EXACT.value
    parameters: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        from core.units.precision import pad_ieee754_decimal
        from core.units.provenance import ConversionClassification
        from core.units.provenance import _audit_dec
        from core.units.provenance import _audit_freeze_mapping
        from core.units.registry import validate_registry_quantity
        from core.units.validation import _safety_copy_quantities
        for attr in ("temperature", "pressure", "density"):
            qty = object.__getattribute__(self, attr)
            if qty is not None:
                if not isinstance(qty, pint.Quantity):
                    raise TypeError(f"{attr} must be a pint.Quantity or None.")
                validate_registry_quantity(qty)
                object.__setattr__(self, attr, _safety_copy_quantities(qty))
        for dec_attr in ("moisture_fraction", "compressibility_factor"):
            val = object.__getattribute__(self, dec_attr)
            if val is not None:
                object.__setattr__(self, dec_attr, pad_ieee754_decimal(_audit_dec(val, dec_attr)))
        cls_val = (
            self.classification.value
            if isinstance(self.classification, ConversionClassification)
            else str(self.classification)
        )
        ConversionClassification(cls_val)  # validate membership
        object.__setattr__(self, "classification", cls_val)
        object.__setattr__(self, "parameters", _audit_freeze_mapping(self.parameters))

    def __getattribute__(self, name: str) -> Any:
        from core.units.validation import _safety_copy_quantities
        val = object.__getattribute__(self, name)
        if name == "parameters" and isinstance(val, Mapping):
            return MappingProxyType({k: _safety_copy_quantities(v) for k, v in val.items()})
        return _safety_copy_quantities(val)


from core.units.validation import _safety_QuantitySnapshot


@dataclass(frozen=True, slots=True)
class GasReferenceConditions(_safety_QuantitySnapshot):
    """Metering state and independent combustion metadata for a gas quantity.

    Neutral pressure units explicitly mean absolute here; gauge labels are
    rejected, never silently shifted. None means unknown, not ideal/dry.
    Z must describe the actual gas mixture at this state's temperature,
    pressure and water content. No equation of state is guessed.
    """
    temperature: Quantity
    absolute_pressure: Quantity
    compressibility_factor: Decimal | None = None
    water_mole_fraction: Decimal | None = None
    combustion_reference_temperature: Quantity | None = None
    standard_id: str | None = None

    def __post_init__(self) -> None:
        from core.units.conversions import _absolute_reference_temperature
        from core.units.conversions import _as_decimal
        from core.units.conversions import _classify_pressure_datum
        from core.units.registry import ureg
        from core.units.registry import validate_registry_quantity
        object.__setattr__(self, "temperature", _absolute_reference_temperature(self.temperature, "Metering temperature"))
        pressure = self.absolute_pressure
        if not isinstance(pressure, pint.Quantity):
            raise TypeError("absolute_pressure must be a pressure quantity.")
        validate_registry_quantity(pressure)
        if _classify_pressure_datum(pressure.units) == "gauge":
            raise ValueError("Gas reference pressure must be absolute, not gauge.")
        absolute = ureg.Quantity(_as_decimal(pressure.magnitude), ureg.Unit(pressure.units._units)).to("pascal")
        if not absolute.magnitude.is_finite() or absolute.magnitude <= 0:
            raise ValueError("Gas reference pressure must be finite and strictly positive.")
        object.__setattr__(self, "absolute_pressure", absolute)
        for name in ("compressibility_factor", "water_mole_fraction"):
            value = getattr(self, name)
            if value is None:
                continue
            if isinstance(value, bool):
                raise ValueError(f"{name} must be a measured number, not boolean.")
            number = _as_decimal(value)
            valid = number.is_finite() and (number > 0 if name == "compressibility_factor" else 0 <= number < 1)
            if not valid:
                raise ValueError(f"Invalid {name}: {value!r}")
            object.__setattr__(self, name, number)
        if self.combustion_reference_temperature is not None:
            object.__setattr__(self, "combustion_reference_temperature", _absolute_reference_temperature(self.combustion_reference_temperature, "Combustion temperature"))
        if self.standard_id is not None:
            if not isinstance(self.standard_id, str) or not self.standard_id.strip():
                raise ValueError("standard_id must be a non-empty string or None.")
            object.__setattr__(self, "standard_id", self.standard_id.strip())

    @classmethod
    def from_standard(cls, standard_id: str, *, compressibility_factor: Decimal | None = None,
                      water_mole_fraction: Decimal | None = Decimal(0),
                      combustion_reference_temperature: Quantity | None = None) -> GasReferenceConditions:
        """Construct a named dry reference; caller still supplies Z explicitly."""
        from core.units.conversions import _GAS_REFERENCE_PRESETS
        from core.units.registry import ureg
        key = standard_id.strip().lower()
        if key not in _GAS_REFERENCE_PRESETS:
            raise ValueError("Unknown P/T reference; Groningen requires a calorific model, not just gas-law scaling.")
        _, _, temperature, pressure = _GAS_REFERENCE_PRESETS[key]
        return cls(ureg.Quantity(Decimal(str(temperature)), "K"),
                   ureg.Quantity(Decimal(str(pressure)) * 1000, "Pa"),
                   compressibility_factor, water_mole_fraction,
                   combustion_reference_temperature, key)


__all__ = ('GasReferenceConditions', 'LiquidProductSpec', 'LiquidState', 'ReferenceConditions')
_PUBLIC_ALIASES = {'LiquidState': ('assays.energy.petroleum', 'LiquidState'), 'LiquidProductSpec': ('assays.energy.petroleum', 'LiquidProductSpec')}

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
    from core.units.assays.energy.petroleum import LiquidProductSpec
    from core.units.assays.energy.petroleum import LiquidState
    from core.units.registry import Quantity
