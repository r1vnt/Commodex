"""Commodex assays metals steel implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class SteelProductAssay(MaterialSpecification):
    product: str
    thickness: pint.Quantity | None = None
    width: pint.Quantity | None = None
    diameter: pint.Quantity | None = None
    coating_mass_per_area: pint.Quantity | None = None

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        from core.units.registry import ureg
        MaterialSpecification.__post_init__(self)
        if self.product not in {"billet", "slab", "HRC", "CRC", "rebar", "galvanised", "other"}:
            raise ValueError("Unknown steel product; declare product='other' for custom geometry.")
        for field in ("thickness", "width", "diameter"):
            if getattr(self, field) is not None:
                object.__setattr__(self, field, ureg.Quantity(_converters_petroleum_positive_quantity(getattr(self, field), "m", field), "m"))
        if self.coating_mass_per_area is not None:
            object.__setattr__(self, "coating_mass_per_area", ureg.Quantity(_converters_petroleum_positive_quantity(self.coating_mass_per_area, "kg/m^2", "coating"), "kg/m^2"))

    def piece_mass(self, length):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_positive_quantity
        size = _converters_petroleum_positive_quantity(length, "m", "piece length")
        if self.diameter is not None:
            if self.width is not None or self.thickness is not None:
                raise ValueError("Choose round diameter or rectangular dimensions.")
            area = Decimal("3.141592653589793238462643383") * self.diameter.magnitude ** 2 / 4
        elif self.width is not None and self.thickness is not None:
            area = self.width.magnitude * self.thickness.magnitude
        else:
            raise ValueError("Explicit round or rectangular dimensions are required.")
        return self.volume_to_mass(_converters_petroleum_output(area * size, "m^3", length))

    def coating_mass(self, coated_area):
        from core.units.assays.energy.petroleum import _converters_petroleum_output
        from core.units.assays.energy.petroleum import _converters_petroleum_quantity_number
        if self.coating_mass_per_area is None:
            raise ValueError("Explicit coating weight is required; area includes every coated face.")
        area = _converters_petroleum_quantity_number(coated_area, "m^2", "total coated area", minimum=0)
        return _converters_petroleum_output(area * self.coating_mass_per_area.magnitude, "kg", coated_area)


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class DirectReducedIronAssay(MaterialSpecification):
    total_iron_fraction: Decimal
    metallisation_fraction: Decimal

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.energy.petroleum import _converters_petroleum_fraction
        from core.units.composition import ElementComposition
        total = _converters_petroleum_fraction(self.total_iron_fraction, "total iron")
        metallisation = _converters_petroleum_fraction(self.metallisation_fraction, "metallisation")
        if self.composition is None:
            object.__setattr__(self, "composition", ElementComposition({"Fe": total}))
        elif self.composition.fraction("Fe") != total:
            raise ValueError("Total iron conflicts with composition.")
        object.__setattr__(self, "total_iron_fraction", total)
        object.__setattr__(self, "metallisation_fraction", metallisation)
        MaterialSpecification.__post_init__(self)

    def metallic_iron_mass(self, quantity, *, mass_basis="dry"):
        from core.units.assays.base import scale_mass
        return scale_mass(self._dry_feed(quantity, mass_basis), self.total_iron_fraction * self.metallisation_fraction)

    def oxide_bound_iron_mass(self, quantity, *, mass_basis="dry"):
        """Elemental Fe in oxides, not total oxide mass."""
        from core.units.assays.base import scale_mass
        return scale_mass(self._dry_feed(quantity, mass_basis), self.total_iron_fraction * (1 - self.metallisation_fraction))


from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class IronPelletAssay(MaterialSpecification):
    product: str

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        MaterialSpecification.__post_init__(self)
        if self.product not in {"pellet", "lump", "fines"}:
            raise ValueError("Iron product must be pellet, lump or fines.")
        if self.composition is None:
            raise ValueError("Explicit Fe grade is required.")
        self.composition.fraction("Fe")

    def beneficiated_product_mass(self, quantity, *, target_iron_fraction, iron_recovery_fraction, mass_basis="dry"):
        from core.units.assays.base import beneficiated_mass
        return beneficiated_mass(self, quantity, component="Fe", target_fraction=target_iron_fraction,
                                  recovery_fraction=iron_recovery_fraction, mass_basis=mass_basis)


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (SteelProductAssay, DirectReducedIronAssay, IronPelletAssay,))
__all__ = ('DirectReducedIronAssay', 'IronPelletAssay', 'SteelProductAssay', 'convert_direct_reduced_iron', 'convert_iron_pellet', 'convert_steel_product', 'create_direct_reduced_iron_context', 'create_iron_pellet_context', 'create_steel_product_context')
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
    from core.units.assays.metals.steel import convert_direct_reduced_iron
    from core.units.assays.metals.steel import convert_iron_pellet
    from core.units.assays.metals.steel import convert_steel_product
    from core.units.assays.metals.steel import create_direct_reduced_iron_context
    from core.units.assays.metals.steel import create_iron_pellet_context
    from core.units.assays.metals.steel import create_steel_product_context
