"""Deterministic properties and explicit gates for costly or external checks."""
import os
import random
from pathlib import Path
import pytest


def pytest_addoption(parser):
    for option in ("performance", "mutation", "network-audit", "release-install"):
        parser.addoption("--run-" + option, action="store_true", default=False)
    parser.addoption("--standards-fixtures", default=os.environ.get("COMMODEX_STANDARDS_FIXTURES"))
    parser.addoption("--test-order-seed",type=int,default=None)


def pytest_configure(config):
    (Path(__file__).resolve().parents[1] / 'artifacts').mkdir(exist_ok=True)
    for marker in ("property", "performance", "mutation", "network_audit", "release_install", "standards_reference", "security", "concurrency", "optional", "checklist(*ids)"):
        config.addinivalue_line("markers", marker + ": comprehensive unit-library test classification")
    try:
        from hypothesis import settings
    except ImportError:
        return
    settings.register_profile("commodex", max_examples=80, derandomize=True, deadline=None, database=None)
    settings.load_profile(os.environ.get("COMMODEX_HYPOTHESIS_PROFILE", "commodex"))


def pytest_collection_modifyitems(config, items):
    if config.getoption('--test-order-seed') is not None:
        random.Random(config.getoption('--test-order-seed')).shuffle(items)
    for item in items:
        for marker, option in (("performance", "performance"), ("mutation", "mutation"), ("network_audit", "network-audit"), ("release_install", "release-install")):
            if marker in item.keywords and not config.getoption("--run-" + option):
                item.add_marker(pytest.mark.skip(reason="Explicit gate: use --run-" + option))


@pytest.fixture(scope="session")
def repository_root():
    return Path(__file__).resolve().parents[1]


@pytest.fixture(autouse=True)
def restore_legacy_registry_extensions(request):
    """Keep the two legacy direct-definition tests independent of execution order."""
    mutations = {
        'tests/core/units/test_currencies_fx.py::FxConversionTests::test_settlement_float_multiplication_and_half_even_boundaries',
        'tests/core/units/test_pressure.py::PressureConversionTests::test_canonical_aliases_and_display_formats',
    }
    if request.node.nodeid not in mutations:
        yield
        return
    from core.units.registry import _registry_transaction

    class Rollback(Exception):
        pass

    try:
        with _registry_transaction():
            yield
            raise Rollback()
    except Rollback:
        pass
