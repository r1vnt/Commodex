"""Commodex assays agriculture fruit implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from typing import ClassVar
import pint

from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class RootCropAssay(_agricultural_families_AgriculturalFamilyAssay):
    _commodities = frozenset({"potato", "cassava", "sweet-potato", "sugar-beet"})

    def starch_mass(self, quantity, *, mass_basis="received"):
        return self.contained_mass(quantity, component="starch", mass_basis=mass_basis)

    def sugar_mass(self, quantity, *, mass_basis="received"):
        return self.contained_mass(quantity, component="sugar", mass_basis=mass_basis)

    def recovered_sugar_mass(self, quantity, *, mass_basis="received"):
        return self.recovered_mass(quantity, component="sugar", mass_basis=mass_basis)


from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay


@dataclass(frozen=True, slots=True, kw_only=True)
class FruitConcentrateAssay(_agricultural_families_AgriculturalFamilyAssay):
    """Mass balances use declared soluble solids, independently of Brix readings.

    https://www.fao.org/4/v5030e/V5030E0i.htm describes sucrose-equivalent Brix.
    Fruit composition/calibration is needed to interpret refractometer readings.
    """
    soluble_solids_fraction: Decimal
    brix_reading: pint.Quantity | None = None
    solids_source: str = "declared soluble-solids mass fraction"

    def __post_init__(self):
        from core.units.assays.base import _agricultural_families_AgriculturalFamilyAssay
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        from core.units.assays.energy.petroleum import _converters_petroleum_text
        from core.units.registry import ureg
        _agricultural_families_AgriculturalFamilyAssay.__post_init__(self)
        solids = _converters_petroleum_fraction(self.soluble_solids_fraction, "soluble solids")
        if solids == 0:
            raise ValueError("Concentrate soluble solids must be positive.")
        if self.moisture_fraction is not None and solids + self.moisture_fraction > 1:
            raise ValueError("Soluble solids plus water cannot exceed total product mass.")
        object.__setattr__(self, "soluble_solids_fraction", solids)
        object.__setattr__(self, "solids_source", _converters_petroleum_text(self.solids_source, "solids source"))
        if self.brix_reading is not None:
            reading = _converters_petroleum_quantity_number(self.brix_reading, "degree_Brix", "Brix reading", minimum=0)
            if reading > 100:
                raise ValueError("Brix reading must lie in [0, 100].")
            object.__setattr__(self, "brix_reading", ureg.Quantity(reading, "degree_Brix"))

    @classmethod
    def from_brix(cls, *, brix, assume_brix_equals_solids=False, **metadata):
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if assume_brix_equals_solids is not True:
            raise ValueError("Declare calibrated soluble solids, or explicitly set assume_brix_equals_solids=True for the sucrose-equivalent approximation.")
        reading = _converters_petroleum_quantity_number(brix, "degree_Brix", "Brix reading", minimum=0)
        return cls(soluble_solids_fraction=reading / 100, brix_reading=brix,
                   solids_source="explicit Brix=solids sucrose-equivalent approximation", **metadata)

    def solids_mass(self, quantity):
        from core.units.assays.base import scale_mass
        return scale_mass(quantity, self.soluble_solids_fraction)

    def reconstituted_mass(self, quantity, *, target_solids_fraction):
        from core.units.assays.base import scale_mass
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        target = _converters_petroleum_fraction(target_solids_fraction, "target soluble solids")
        if target == 0 or target > self.soluble_solids_fraction:
            raise ValueError("Reconstitution target must be positive and no stronger than the concentrate.")
        return scale_mass(quantity, self.soluble_solids_fraction / target)

    def dilution_water_required(self, quantity, *, target_solids_fraction):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        output = self.reconstituted_mass(quantity, target_solids_fraction=target_solids_fraction)
        return _converters_petroleum_output(mass(output) - mass(quantity), "kg", quantity)

    def reconstituted_volume(self, quantity, *, target_solids_fraction, target_density):
        from core.units.assays.base import mass
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        output = self.reconstituted_mass(quantity, target_solids_fraction=target_solids_fraction)
        density = _converters_petroleum_positive_quantity(target_density, "kg/m^3", "reconstituted juice density")
        return _converters_petroleum_output(mass(output) / density, "m^3", quantity)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (RootCropAssay, FruitConcentrateAssay,))
__all__ = ('FruitConcentrateAssay', 'RootCropAssay', 'convert_fruit_concentrate', 'convert_root_crop', 'create_fruit_concentrate_context', 'create_root_crop_context')
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
    from core.units.assays.agriculture.fruit import convert_fruit_concentrate
    from core.units.assays.agriculture.fruit import convert_root_crop
    from core.units.assays.agriculture.fruit import create_fruit_concentrate_context
    from core.units.assays.agriculture.fruit import create_root_crop_context
