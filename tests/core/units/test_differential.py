"""Independent Pint, Decimal, rational and analytical calculation oracles."""
from decimal import Decimal as D
from fractions import Fraction
import math
import pytest
import pint
import core.units as u


@pytest.mark.parametrize('source,target',[('kilogram','pound'),('meter','foot'),('joule','kilowatt_hour'),('pascal','bar'),('liter','meter**3')])
def test_independent_unmodified_pint_registry_agrees_for_physical_units(source,target):
    baseline=pint.UnitRegistry()
    assert float(u.Q_(D('12.5'),source).to(target).magnitude) == pytest.approx(float(baseline.Quantity(12.5,source).to(target).magnitude),rel=1e-14)


def test_fx_notional_and_moisture_against_independent_decimal_and_rational_arithmetic():
    rate=D('.8123456789')
    value=u.convert_currency(u.Q_(D(1250),'USD'),'EUR',rates={'EUR':rate})
    assert value.magnitude == D(1250)*rate
    ratio=Fraction(37,100)
    mass=D(100)*(1-D(ratio.numerator)/ratio.denominator)
    assert u.convert(u.Q_(D(100),'wmt'),'dmt',moisture_fraction=D('.37')).magnitude == mass
    notional=(u.Q_(D(50),'USD/bbl')*u.Q_(D(1000),'bbl')).to('USD')
    assert notional.magnitude == D(50000)


def test_analytical_heat_rate_and_numerical_generation_are_consistent():
    eta=D('.41')
    rate=u.heat_rate_from_efficiency(eta,unit='MWh_th/MWh_e')
    output=u.GenerationProcess(heat_rate=rate).convert(u.Q_(D(100),'MWh_th'),'MWh_e')
    assert abs(output.magnitude-D(100)*eta) < D('1e-12')


def test_ideal_gas_density_against_independent_equation_of_state():
    composition=u.GasComposition({'CH4':D(1)})
    state=u.GasReferenceConditions.from_standard('iso15',compressibility_factor=D(1))
    backend=u.builtin_gas_backend('ideal',temperature_range=(u.Q_(100,'K'),u.Q_(400,'K')),pressure_range=(u.Q_(1000,'Pa'),u.Q_(200000,'Pa')))
    actual=u.gas_density(composition,state,backend=backend).to('kg/m^3').magnitude
    expected=state.absolute_pressure.to('Pa').magnitude*composition.molar_mass.to('kg/mol').magnitude/(backend.gas_constant.to('J/mol/K').magnitude*state.temperature.to('K').magnitude)
    assert float(actual) == pytest.approx(float(expected),rel=2e-14)
