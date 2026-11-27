"""Physical-domain boundaries, parameter exclusivity and specific error contracts."""
from decimal import Decimal as D
import pytest
import pint
import core.units as u


@pytest.mark.parametrize('value',[float('nan'),float('inf'),float('-inf'),D('sNaN'),None,True,complex(1,2)])
def test_invalid_magnitude_input_is_rejected(value):
    with pytest.raises((TypeError,ValueError,ArithmeticError)):
        u.convert(u.Q_(value,'kg'),'g')


@pytest.mark.parametrize('constructor,args,kwargs',[
    (u.MaterialProperties,(),{'density':u.Q_(0,'kg/m^3')}),
    (u.EnergyLoss,('invalid',-1),{}),(u.EnergyLoss,('invalid',1.01),{}),
    (u.HeatPumpSpecification,(0,),{}),(u.GenerationProcess,(),{'efficiency':0}),
    (u.GenerationProcess,(),{'efficiency':1.01}),
    (u.GenerationProcess,(),{'efficiency':.5,'heat_rate':u.Q_(2,'MWh_th/MWh_e')}),
    (u.MaterialComposition,({'a':.8,'b':.8},),{}),
])
def test_invalid_density_efficiency_loss_composition_and_mutually_exclusive_parameters(constructor,args,kwargs):
    with pytest.raises((TypeError,ValueError,ArithmeticError)):
        constructor(*args,**kwargs)


def test_negative_physical_quantities_and_missing_assays_have_explicit_exception_types():
    with pytest.raises(u.InvalidPhysicalStateError): u.convert(u.Q_(-1,'kg'),'g')
    with pytest.raises((u.MissingAssayError,pint.DimensionalityError)): u.convert(u.Q_(1,'bbl'),'kg')
    with pytest.raises(pint.DimensionalityError): u.Q_(1,'EUA').to('UKA')
