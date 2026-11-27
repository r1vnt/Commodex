"""Explicit binary64 tolerances, Decimal contexts and non-finite boundaries."""
from decimal import Decimal as D, localcontext, getcontext, ROUND_HALF_EVEN, ROUND_HALF_UP
import math
import pytest
import core.units as u
from core.units.precision import decimal_precision


@pytest.mark.parametrize("magnitude", [D('1e-250'),D('1e250'),D('1e-25'),D('1e25'),D(0),D('-1e-25')])
def test_extreme_decimal_magnitudes_relative_tolerance_and_repeated_conversion_drift(magnitude):
    result = u.Q_(magnitude,"kg")
    for _ in range(50): result=result.to('lb').to('kg')
    error = abs(result.magnitude-magnitude)
    assert error <= abs(magnitude)*D('2e-13')+D('1e-270')


def test_decimal_context_isolation_rounding_and_cancellation():
    original = getcontext().copy()
    with decimal_precision(50):
        big = D('1e30')
        assert (big+D(1))-big == 1
        assert sum((D('.1') for _ in range(1000)),D(0)) == 100
    assert getcontext().prec == original.prec and getcontext().rounding == original.rounding
    with localcontext() as context:
        context.rounding=ROUND_HALF_EVEN
        assert D('1.005').quantize(D('.01')) == D('1.00')
        context.rounding=ROUND_HALF_UP
        assert D('1.005').quantize(D('.01')) == D('1.01')


@pytest.mark.parametrize("value", [float('nan'),float('inf'),float('-inf'),D('sNaN'),None,True,complex(1,2),'not-a-number'])
def test_invalid_numeric_types_and_nonfinite_values_fail(value):
    with pytest.raises((TypeError,ValueError,ArithmeticError)):
        u.parse_quantity(value,'kg')


def test_float_overflow_is_rejected_while_finite_decimal_extremes_remain_representable():
    with pytest.raises((ValueError,OverflowError)):
        u.parse_quote('1e999999 USD/kg',as_decimal=False)
    assert u.parse_quote('1e999999 USD/kg').magnitude.is_finite()
    assert u.parse_quantity(2,'kg',as_decimal=False).magnitude == 2.


def test_decimal_affine_temperature_compatibility_preserves_precision_and_inverse():
    from core.units.registry import _DecimalOffsetConverter
    converter = _DecimalOffsetConverter(1, 273.15)
    value = D('25.123456789012345678901234')
    assert converter.from_reference(converter.to_reference(value)) == value
    assert converter.to_reference(25.0) == pytest.approx(298.15)
    assert converter.from_reference(298.15) == pytest.approx(25.0)
    actual = u.Q_(value, 'degC').to('kelvin')
    assert actual.magnitude == value + D('273.15')
    assert actual.to('degC').magnitude == value
