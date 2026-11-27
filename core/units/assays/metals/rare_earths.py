"""Commodex assays metals rare_earths implementation and public interfaces."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal
from types import MappingProxyType
import pint

from core.units.assays.base import MaterialSpecification


@dataclass(frozen=True, slots=True, kw_only=True)
class RareEarthAssay(MaterialSpecification):
    oxide_fractions: Mapping[str, Decimal]

    def __post_init__(self):
        from core.units.assays.base import MaterialSpecification
        from core.units.assays.base import _material_specs_RARE_EARTHS
        from core.units.composition import ElementComposition
        from core.units.composition import StoichiometricCompound
        MaterialSpecification.__post_init__(self)
        values = ElementComposition(self.oxide_fractions).mass_fractions
        for formula in values:
            atoms = StoichiometricCompound(formula).atoms
            if "O" not in atoms or len(atoms) != 2 or not (atoms.keys() & _material_specs_RARE_EARTHS):
                raise ValueError("Supply individual rare-earth oxide grades, not an undifferentiated TREO formula.")
        object.__setattr__(self, "oxide_fractions", values)

    def treo_mass(self, quantity, *, mass_basis="dry"):
        from core.units.assays.base import scale_mass
        return scale_mass(self._dry_feed(quantity, mass_basis), sum(self.oxide_fractions.values(), Decimal(0)))

    def oxide_mass(self, quantity, *, oxide, mass_basis="dry"):
        from core.units.assays.base import scale_mass
        if oxide not in self.oxide_fractions:
            raise ValueError("No assay for the selected oxide.")
        return scale_mass(self._dry_feed(quantity, mass_basis), self.oxide_fractions[oxide])

    def element_mass(self, quantity, *, element, mass_basis="dry"):
        from core.units.assays.base import _material_specs_RARE_EARTHS
        from core.units.assays.base import scale_mass
        from core.units.composition import StoichiometricCompound
        fraction = Decimal(0)
        for oxide, grade in self.oxide_fractions.items():
            compound = StoichiometricCompound(oxide)
            if element in compound.atoms:
                fraction += grade * compound.element_fraction(element)
        if element not in _material_specs_RARE_EARTHS:
            raise ValueError("Select a rare-earth element.")
        return scale_mass(self._dry_feed(quantity, mass_basis), fraction)

    def ndpr_element_mass(self, quantity, *, mass_basis="dry"):
        nd = self.element_mass(quantity, element="Nd", mass_basis=mass_basis)
        pr = self.element_mass(quantity, element="Pr", mass_basis=mass_basis)
        return nd + pr


from core.units.assays.base import install_api
_generated_exports = install_api(globals(), (RareEarthAssay,))
__all__ = ('RareEarthAssay', 'convert_rare_earth', 'create_rare_earth_context')
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
    from core.units.assays.metals.rare_earths import convert_rare_earth
    from core.units.assays.metals.rare_earths import create_rare_earth_context
