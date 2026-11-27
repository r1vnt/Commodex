"""Generated algebraic and metamorphic tests, independent of implementation paths."""
import math
from decimal import Decimal as D
import pytest
import pint
pytest.importorskip("hypothesis", reason="Install commodex[dev] for generated property tests")
from hypothesis import given, strategies as st
import core.units as u

PAIRS = [("kg", "lb"), ("liter", "bbl"), ("J", "kWh"), ("bar", "Pa"), ("m", "ft"), ("kg/hour", "lb/day")]
numbers = st.floats(min_value=1e-8, max_value=1e8, allow_nan=False, allow_infinity=False)


@given(st.sampled_from(PAIRS), numbers, numbers, st.floats(min_value=.1, max_value=100, allow_nan=False))
def test_identity_inverse_roundtrip_linearity_additivity_homogeneity_and_monotonicity(pair, a, b, scale):
    source, target = pair
    qa, qb = u.Q_(a, source), u.Q_(b, source)
    ya, yb = qa.to(target), qb.to(target)
    assert qa.to(source).magnitude == a
    assert ya.to(source).magnitude == pytest.approx(a, rel=3e-13)
    assert (qa + qb).to(target).magnitude == pytest.approx((ya + yb).magnitude, rel=3e-13)
    assert (qa * scale).to(target).magnitude == pytest.approx((ya * scale).magnitude, rel=3e-13)
    assert (ya.magnitude <= yb.magnitude) if a <= b else (ya.magnitude >= yb.magnitude)
    assert qa.dimensionality == ya.dimensionality
    assert u.Q_(0, source).to(target).magnitude == 0
    assert (-qa).to(target).magnitude == pytest.approx(-ya.magnitude)


@given(numbers)
def test_conversion_composition_path_independence_and_scale_invariance(value):
    q = u.Q_(value, "kg")
    direct = q.to("lb")
    assert q.to("g").to("lb").magnitude == pytest.approx(direct.magnitude, rel=3e-13)
    assert q.to_base_units().to_base_units() == q.to_base_units()


@given(st.integers(min_value=0, max_value=10**12), st.integers(min_value=0, max_value=90))
def test_generated_moisture_roundtrip_conservation_and_boundaries(mass, moisture):
    wet = u.Q_(D(mass), "wmt")
    dry = u.convert(wet, "dmt", moisture_fraction=D(moisture)/100)
    assert dry.magnitude == D(mass)*(1-D(moisture)/100)
    assert abs(u.convert(dry, "wmt", moisture_fraction=D(moisture)/100).magnitude-D(mass)) < D("1e-12")


@given(st.sampled_from([("MWh_th","MWh_e"),("wmt","dmt"),("USD","EUR"),("REC","MWh"),("EUA","UKA")]), numbers)
def test_exhaustive_forbidden_families_never_gain_implicit_compatibility(pair, value):
    source, target = pair
    assert not u.Q_(value, source).is_compatible_with(target)
    with pytest.raises(pint.DimensionalityError):
        u.Q_(value, source).to_base_units().to(target)


@given(st.permutations([1,2,3,4,5]))
def test_permutation_invariance_of_linear_aggregation(order):
    expected = sum((u.Q_(value, "kg").to("g") for value in order), u.Q_(0,"g"))
    assert expected.magnitude == 15000


def test_exhaustive_compatible_pairs_in_the_versioned_unit_catalogue():
    import json
    from collections import defaultdict
    from core.units.standards import REFERENCE_DIRECTORY,REFERENCE_VERSION
    rows=json.loads((REFERENCE_DIRECTORY/'units'/REFERENCE_VERSION/'registry_units.json').read_text(encoding='utf-8'))['records']
    groups=defaultdict(list)
    for row in rows:
        unit=u.ureg.Unit(row['id'])
        groups[unit.dimensionality].append(unit)
    checked=0
    for units in groups.values():
        for source in units:
            for target in units:
                original=u.Q_(1.25,source)
                if u.are_convertible(source,target):
                    assert original.to(target).to(source).magnitude == pytest.approx(1.25,rel=2e-12)
                else:
                    with pytest.raises(pint.DimensionalityError): original.to(target)
                checked+=1
    assert checked > 2000
