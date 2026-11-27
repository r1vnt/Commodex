"""Shared mathematical contracts; each leaf module selects its own specimens."""
from dataclasses import FrozenInstanceError, fields, is_dataclass
from decimal import Decimal as D
from functools import lru_cache
import pickle
import pint
import pytest
import core.units as u


@lru_cache(maxsize=1)
def specimens():
    from tests.core.units.test_sector_specifications import examples as sector_examples
    from tests.core.units.assays.agriculture.test_families import examples as agricultural_examples
    values = []
    for category in ("crude_oil","refined_products","natural_gas","coal","alternative_fuels","metals","chemicals","agriculture"):
        for record in u.reference_records(category):
            cls = getattr(u,record.id.split('.')[0])
            args,kwargs = record.constructor_arguments()
            values.append(cls(*args,**kwargs))
    values.extend(sector_examples())
    values.extend(agricultural_examples())
    values.append(u.RefineryCrudeAssay(u.CrudeAssay.brent(),(u.DistillationCut("whole",0,800,100,100,.84),)))
    values.append(u.LNGCargo("measured example",u.GasComposition({"CH4":D(1)}),u.LiquidState(u.Q_(112,"K"),u.Q_(150000,"Pa")),measured_density=u.Q_(D(450),"kg/m^3"),density_source="test density certificate"))
    unique = {}
    for value in values:
        unique.setdefault(type(value),value)
    return tuple(unique.values())


def cases(module):
    return [pytest.param(value,id=type(value).__name__) for value in specimens() if type(value).__module__ == module]


def assert_snapshot_and_ownership(value):
    if is_dataclass(value):
        first = fields(value)[0].name
        with pytest.raises((FrozenInstanceError,AttributeError)):
            setattr(value,first,getattr(value,first))
    foreign = pint.UnitRegistry().Quantity(1,"kg")
    with pytest.raises(u.ForeignRegistryError):
        u.convert(foreign,"g",assay=value)


def assert_physical_or_energy_contract(value):
    if isinstance(value,u.LNGCargo):
        assert value.mass(u.Q_(D(2),"m^3")).to("kg").magnitude == 900
        return
    if isinstance(value,u.RefineryCrudeAssay):
        # Refinery processing is a wrapper around its explicit bulk assay.
        return assert_physical_or_energy_contract(value.bulk)
    # Models with declared density have an independent, one-line mass oracle.
    density = getattr(value,"density",None)
    if isinstance(density,pint.Quantity) and density.is_compatible_with("kg/m^3"):
        amount = u.Q_(D(2),"m^3")
        converted = u.convert(amount,"kg",assay=value)
        expected = D(2)*D(str(density.to("kg/m^3").magnitude))
        assert float(converted.magnitude) == pytest.approx(float(expected),rel=2e-12)
        assert float(u.convert(converted,"m^3",assay=value).magnitude) == pytest.approx(2,rel=2e-12)
        return
    density = getattr(value,"density_kg_per_m3",None)
    if density is not None:
        assert float(u.convert(u.Q_(D(2),"m^3"),"kg",assay=value).magnitude) == pytest.approx(2*float(density),rel=2e-12)
        return
    # A missing physical specification must not prevent same-unit accounting,
    # and must not manufacture a commodity-to-currency exchange factor.
    assert u.convert(u.Q_(D(2),"kg"),"g",assay=value).magnitude == 2000
    with pytest.raises((ValueError,pint.PintError)):
        u.convert(u.Q_(D(2),"kg"),"USD",assay=value)
