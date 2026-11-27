"""Commodex assays agriculture specialty implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from typing import ClassVar
import pint

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class SpecialtyGrainAssay(_agricultural_families_AgriculturalFamilyAssay):
    mass_per_bushel: pint.Quantity | None = None
    _commodities = frozenset({"millet", "sorghum", "rye", "triticale", "buckwheat", "quinoa"})

    def __post_init__(self):
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.registry import ureg
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        if self.mass_per_bushel is not None:
            object.__setattr__(self, "mass_per_bushel", ureg.Quantity(_converters_petroleum_positive_quantity(self.mass_per_bushel, "kg", "crop-weight bushel"), "kg"))

    def mass_to_bushels(self, quantity):
        """Crop-weight bushels; the registry's geometric bushel stays a volume."""
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        if self.mass_per_bushel is None:
            raise ValueError("Declare crop-weight mass_per_bushel.")
        return _converters_petroleum_output(mass(quantity) / self.mass_per_bushel.magnitude, "dimensionless", quantity)

    def bushels_to_mass(self, quantity):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if self.mass_per_bushel is None:
            raise ValueError("Declare crop-weight mass_per_bushel.")
        count = _converters_petroleum_quantity_number(quantity, "dimensionless", "crop bushel count", minimum=0)
        return _converters_petroleum_output(count * self.mass_per_bushel.magnitude, "kg", quantity)


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class TeaAssay(_agricultural_families_AgriculturalFamilyAssay):
    fresh_to_made_yield_fraction: Decimal | None = None
    _commodities = frozenset({"green", "black", "oolong", "white", "other"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "fresh_to_made_yield_fraction", _agricultural_families_yield(self.fresh_to_made_yield_fraction, "fresh-to-made tea yield"))

    def made_tea_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.fresh_to_made_yield_fraction, "fresh-to-made tea yield")

    def fresh_leaf_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.fresh_to_made_yield_fraction, "fresh-to-made tea yield")


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class TobaccoAssay(_agricultural_families_AgriculturalFamilyAssay):
    curing_yield_fraction: Decimal | None = None
    _commodities = frozenset({"leaf", "cured", "processed"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "curing_yield_fraction", _agricultural_families_yield(self.curing_yield_fraction, "curing yield"))

    def cured_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.curing_yield_fraction, "curing yield")

    def green_leaf_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.curing_yield_fraction, "curing yield")


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class SpiceAssay(_agricultural_families_AgriculturalFamilyAssay):
    processing_recovery_fraction: Decimal | None = None

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        object.__setattr__(self, "processing_recovery_fraction", _agricultural_families_yield(self.processing_recovery_fraction, "spice processing recovery"))

    def processed_spice_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.processing_recovery_fraction, "spice processing recovery")

    def raw_spice_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.processing_recovery_fraction, "spice processing recovery")


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class NutAssay(_agricultural_families_AgriculturalFamilyAssay):
    kernel_yield_fraction: Decimal | None = None
    shell_yield_fraction: Decimal | None = None
    _commodities = frozenset({"almonds", "cashews", "pistachios", "walnuts", "hazelnuts", "other"})

    def __post_init__(self):
        from core.units.analytical_basis import _agricultural_families_yield
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        for name in ("kernel_yield_fraction", "shell_yield_fraction"):
            object.__setattr__(self, name, _agricultural_families_yield(getattr(self, name), name))
        if self.kernel_yield_fraction is not None and self.shell_yield_fraction is not None and self.kernel_yield_fraction + self.shell_yield_fraction > 1:
            raise ValueError("Kernel plus shell yield cannot exceed inshell feed.")

    def kernel_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.kernel_yield_fraction, "kernel yield")

    def shell_mass(self, quantity):
        from core.units.analytical_basis import _agricultural_families_forward
        return _agricultural_families_forward(quantity, self.shell_yield_fraction, "shell yield")

    def inshell_mass_required(self, quantity):
        from core.units.analytical_basis import _agricultural_families_inverse
        return _agricultural_families_inverse(quantity, self.kernel_yield_fraction, "kernel yield")


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (SpecialtyGrainAssay, TeaAssay, TobaccoAssay, SpiceAssay, NutAssay,))
__all__ = ('NutAssay', 'SpecialtyGrainAssay', 'SpiceAssay', 'TeaAssay', 'TobaccoAssay', 'convert_nut', 'convert_specialty_grain', 'convert_spice', 'convert_tea', 'convert_tobacco', 'create_nut_context', 'create_specialty_grain_context', 'create_spice_context', 'create_tea_context', 'create_tobacco_context')
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
    from core.units.assays.agriculture.specialty import convert_nut
    from core.units.assays.agriculture.specialty import convert_specialty_grain
    from core.units.assays.agriculture.specialty import convert_spice
    from core.units.assays.agriculture.specialty import convert_tea
    from core.units.assays.agriculture.specialty import convert_tobacco
    from core.units.assays.agriculture.specialty import create_nut_context
    from core.units.assays.agriculture.specialty import create_specialty_grain_context
    from core.units.assays.agriculture.specialty import create_spice_context
    from core.units.assays.agriculture.specialty import create_tea_context
    from core.units.assays.agriculture.specialty import create_tobacco_context
