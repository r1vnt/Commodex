"""Solver-domain boundaries and explicit unsupported-model rejection."""
from decimal import Decimal as D
import pytest
import core.units as u


@pytest.mark.parametrize('pressure',[0.,.001,101.325,1000.,10000.])
def test_ideal_eos_is_pressure_independent_and_zero_pressure_is_finite(pressure):
    assert u.calculate_z_factor(u.Q_(pressure,'kPa'),u.Q_(300,'K'),method='ideal') == 1


@pytest.mark.parametrize('pressure,temperature',[(-1,300),(100,0),(100,-10),(float('nan'),300),(100,float('inf'))])
def test_invalid_eos_physical_states_are_rejected(pressure,temperature):
    with pytest.raises((ValueError,ArithmeticError)):
        u.calculate_z_factor(u.Q_(pressure,'kPa'),u.Q_(temperature,'K'),method='ideal')


def test_aga8_gross_is_not_silently_replaced_with_detail_or_heos():
    with pytest.raises(ValueError):
        u.calculate_z_factor(u.Q_(100,'kPa'),u.Q_(300,'K'),method='aga8_gross',molar_composition={'CH4':1})
