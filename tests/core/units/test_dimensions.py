"""Physical, economic and process dimensions are independent contracts."""
from itertools import combinations
import pytest
import pint
import core.units as u
from core.units.dimensions import energy_roles


@pytest.mark.parametrize('unit,roles',[('MWh_th',{'thermal':1}),('USD/MWh_e',{'electrical':-1}),('MWh_th/MWh_e',{'thermal':1,'electrical':-1})])
def test_role_dimensions_survive_compound_quotes_and_base_unit_reduction(unit,roles):
    assert dict(energy_roles(unit)) == roles
    assert dict(energy_roles(u.Q_(1,unit).to_base_units().units)) == roles


@pytest.mark.parametrize('source,target',list(combinations(['USD','EUR','JPY','wmt','dmt','contained_metal_tonne','payable_metal_tonne','EUA','UKA','CCA','REC','GO','MWh_e','MWh_th'],2)))
def test_all_forbidden_cross_dimension_pairs_reject_plain_pint(source,target):
    with pytest.raises(pint.DimensionalityError): u.Q_(1,source).to(target)
