"""Independent dimensional factors, offsets, flow rates and commercial counts."""
from decimal import Decimal as D
from fractions import Fraction
import pytest
import pint
import core.units as u

FACTORS = [
    ("kg", "g", "1000"), ("lb", "kg", "0.45359237"), ("mt", "kg", "1000"),
    ("liter", "m^3", ".001"), ("bbl", "m^3", ".158987294928"),
    ("kWh", "J", "3600000"), ("MW", "W", "1000000"),
    ("kg/liter", "kg/m^3", "1000"), ("bar", "Pa", "100000"),
    ("kg/hour", "kg/second", str(D(1)/3600)), ("m^3/hour", "liter/hour", "1000"),
    ("kmol/hour", "mol/hour", "1000"), ("MWh/day", "kWh/day", "1000"),
    ("centipoise", "Pa*second", ".001"), ("mg/liter", "g/m^3", "1"),
    ("g/mol", "kg/mol", ".001"), ("hectare", "m^2", "10000"),
    ("mile", "meter", "1609.344"), ("ampere_hour", "coulomb", "3600"),
    ("percent", "dimensionless", ".01"), ("dozen", "item", "12"),
    ("FEU", "TEU", "2"), ("short_ton_mile", "kg*meter", "1459972.31821056"),
]


@pytest.mark.parametrize("source,target,factor", FACTORS)
@pytest.mark.parametrize("number", [2.5, D("2.5")])
def test_physical_conversion_factors_against_independent_decimal_arithmetic(source, target, factor, number):
    output = u.Q_(number, source).to(target)
    expected = D(str(number)) * D(factor)
    assert float(output.magnitude) == pytest.approx(float(expected), rel=2e-13)
    assert float(output.to(source).magnitude) == pytest.approx(float(number), rel=2e-13)
    assert isinstance(output.magnitude, type(number))


@pytest.mark.parametrize("unit,value,target,expected", [("degC",0,"kelvin",273.15), ("degF",32,"degC",0), ("delta_degF",9,"delta_degC",5)])
def test_absolute_temperature_and_temperature_differences(unit, value, target, expected):
    assert u.Q_(value, unit).to(target).magnitude == pytest.approx(expected)


def test_rational_scaling_and_compound_denominators():
    expected = Fraction(45359237, 100000000)
    assert D(str(u.Q_(D(1), "lb").to("kg").magnitude)) == D(expected.numerator) / expected.denominator
    amount = u.Q_(D("-3.25"), "USD/kg/hour")
    restored = amount.to("USD/mt/day").to(amount.units)
    # Pint's registry factors use binary64; Decimal magnitudes do not make
    # the underlying conversion constants arbitrary-precision rationals.
    assert abs(restored.magnitude - amount.magnitude) < D("1e-14")


@pytest.mark.parametrize("source,target", [("WS", "USD/mt"), ("price_tick", "USD"), ("contract", "kg"), ("TEU", "m^3"), ("gross_tonnage", "kg")])
def test_contract_and_shipping_indices_need_explicit_contract_factors(source, target):
    with pytest.raises(pint.DimensionalityError):
        u.Q_(1, source).to(target)
