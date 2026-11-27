"""Independent processing balances and blend permutation/conservation contracts."""
from decimal import Decimal as D
import pytest
import core.units as u


def test_calorific_blend_is_mass_weighted_and_permutation_invariant():
    metadata={'thermodynamic_basis':u.ThermodynamicBasis.CONSTANT_PRESSURE,'combustion_temperature':u.Q_(25,'degC')}
    first=u.HeatContentAssay(gcv=u.Q_(D(20),'MJ/kg'),ncv=u.Q_(D(18),'MJ/kg'),**metadata)
    second=u.HeatContentAssay(gcv=u.Q_(D(40),'MJ/kg'),ncv=u.Q_(D(36),'MJ/kg'),**metadata)
    a=u.blend_heat_content([first,second],[D('.25'),D('.75')],weighting='mass')
    b=u.blend_heat_content([second,first],[D('.75'),D('.25')],weighting='mass')
    assert a.gcv.magnitude == b.gcv.magnitude == 35
    assert a.ncv.magnitude == b.ncv.magnitude == D('31.5')


def test_named_losses_are_multiplicative_and_close_energy_conservation():
    stages=[u.EnergyLoss('transport',D('.1')),u.EnergyLoss('storage',D('.2'))]
    balance=u.apply_energy_losses(u.Q_(D(100),'MWh_th'),stages)
    assert balance.delivered_energy.magnitude == 72
    lost=sum((quantity for _,quantity in balance.stage_losses),u.Q_(D(0),'MWh_th'))
    assert balance.delivered_energy+lost == balance.input_energy
    assert u.energy_before_losses(balance.delivered_energy,stages).magnitude == 100
