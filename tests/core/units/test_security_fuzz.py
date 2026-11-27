"""Parser/configuration/payload abuse must fail without code execution or state changes."""
from decimal import Decimal
import json
import subprocess
import sys
import pytest
import pint
import core.units as u

ATTACKS = ["@import /etc/passwd", "__import__('os').system('x')", "meter.__class__", "meter ** 99999", "meter + second", "meter; print(1)", "(meter"*200, "x"*10000]


@pytest.mark.security
@pytest.mark.parametrize("text", ATTACKS)
def test_unit_parser_rejects_malicious_definitions_without_registry_mutation(text):
    before = u.get_registry_definition_hash()
    with pytest.raises((ValueError, pint.PintError, SyntaxError)):
        u.parse_unit(text)
    assert u.get_registry_definition_hash() == before


@pytest.mark.parametrize("text", ["kіlogram", "UЅD/kg", "ЕUR/MWh", "mеter"])
def test_unicode_homoglyphs_are_not_guessed_as_units_or_currencies(text):
    with pytest.raises((ValueError, pint.PintError)):
        u.parse_unit(text)


@pytest.mark.parametrize("version,category", [("../escape","crude_oil"),("2026-11-27","../crude_oil"),("C:\\windows","crude_oil"),(".","crude_oil")])
def test_reference_identifiers_cannot_traverse_paths(version, category):
    with pytest.raises(ValueError):
        u.reference_records(category, version=version)


@pytest.mark.parametrize("payload", [{}, {"m":"__import__('os')", "u":"kg"}, {"m":"NaN","u":"kg"}, {"m":"1","u":"meter.__class__"}, {"schema_version":999,"m":"1","u":"kg"}])
def test_unsafe_or_corrupted_quantity_payloads_are_rejected(payload):
    with pytest.raises((ValueError, TypeError, ArithmeticError, KeyError, pint.PintError)):
        u.quantity_from_dict(payload)


def test_recursive_configuration_fails_transactionally():
    payload = {"commodities": []}
    payload["commodities"].append(payload)
    before = u.get_registry_definition_hash()
    with pytest.raises(ValueError):
        u.load_desk_config(payload, freeze_after_load=False)
    assert u.get_registry_definition_hash() == before


def test_oversized_and_nested_parser_inputs_are_bounded_in_a_subprocess():
    program = """
import core.units as u
for text in ('x'*1000000, '('*10000+'kg'+')'*10000):
    try: u.parse_quote(text)
    except (ValueError, TypeError, ArithmeticError, Exception): pass
    else: raise AssertionError('oversized input accepted')
"""
    result = subprocess.run([sys.executable,"-B","-c",program], capture_output=True, text=True, timeout=15)
    assert result.returncode == 0, result.stderr


def test_generated_quote_and_configuration_fuzz_cases():
    pytest.importorskip("hypothesis")
    from hypothesis import given, settings, strategies as st
    @settings(max_examples=100, derandomize=True, deadline=None, database=None)
    @given(st.text(alphabet="@$;{}[]_/'\\\x000123αЖ", min_size=1,max_size=80))
    def check(text):
        before = u.get_registry_definition_hash()
        try:
            result = u.parse_quote(text)
        except (ValueError, TypeError, ArithmeticError, pint.PintError, SyntaxError):
            pass
        else:
            assert isinstance(result, pint.Quantity)
            assert Decimal(str(result.magnitude)).is_finite()
        assert u.get_registry_definition_hash() == before
    check()
