"""Independent matrix algebra and seeded Monte Carlo checks of first-order uncertainty."""
from decimal import Decimal as D
import math
import random
import pytest
import core.units as u


def test_correlated_uncertainty_matches_seeded_monte_carlo_and_independent_formula():
    rng = random.Random(9173)
    sigma_x, sigma_y, rho = 2., 3., .4
    samples = []
    for _ in range(40000):
        z1, z2 = rng.gauss(0,1), rng.gauss(0,1)
        x, y = sigma_x*z1, sigma_y*(rho*z1+math.sqrt(1-rho*rho)*z2)
        samples.append(2*x-3*y)
    mean = sum(samples)/len(samples)
    observed = math.sqrt(sum((value-mean)**2 for value in samples)/(len(samples)-1))
    expected = math.sqrt((2*sigma_x)**2+(-3*sigma_y)**2+2*2*(-3)*rho*sigma_x*sigma_y)
    result = u.propagate_uncertainty({"x":2,"y":-3},{"x":sigma_x,"y":sigma_y},correlations={('x','y'):rho})
    assert float(result.standard_uncertainty) == pytest.approx(expected, rel=1e-14)
    assert observed == pytest.approx(expected, rel=.025)


@pytest.mark.parametrize("matrix", [[[1,2],[0,1]], [[1,2],[2,1]], [[-1,0],[0,1]], [[1]], [[1,float('nan')],[float('nan'),1]]])
def test_covariance_matrix_validation_rejects_asymmetry_negative_variance_and_non_psd(matrix):
    with pytest.raises((ValueError, ArithmeticError)):
        u.propagate_uncertainty({"x":1,"y":1},{"x":1,"y":1},covariance=matrix)


def test_absolute_relative_uncertainty_and_result_reproducibility():
    first = u.convert(u.Q_(D(10),"kg"),"g", return_metadata=True, uncertainties={"quantity":D(".1")})
    second = u.convert(u.Q_(D(10),"kg"),"g", return_metadata=True, uncertainties={"quantity":D(".1")})
    assert first.quantity == second.quantity
    assert first.uncertainty == D(100)
    assert first.uncertainty / first.quantity.magnitude == D(".01")


def test_indefinite_three_dimensional_matrix_is_rejected_even_for_positive_projection():
    matrix=[[1,-.75,-.75],[-.75,1,-.75],[-.75,-.75,1]]
    with pytest.raises(ValueError,match='positive semi-definite'):
        u.propagate_uncertainty({'a':1,'b':0,'c':0},{'a':1,'b':1,'c':1},correlations=matrix)


def test_singular_and_ill_scaled_valid_covariance_matrices_are_accepted():
    singular=u.propagate_uncertainty({'a':1,'b':-1},{'a':1,'b':1},correlations=[[1,1],[1,1]])
    assert singular.standard_uncertainty == 0
    scaled=u.propagate_uncertainty({'a':1,'b':1},{'a':D('1e-30'),'b':1},correlations=[[1,D('.1')],[D('.1'),1]])
    assert scaled.standard_uncertainty > 0
