"""Unsupported functionality is rejected rather than filled with guessed policies."""
from decimal import Decimal as D
import pytest
import core.units as u
from core.units.batch import convert_many


def test_bid_ask_objects_cannot_be_silently_used_as_single_fx_rates():
    with pytest.raises((ValueError,TypeError,ArithmeticError)):
        u.convert_fx(u.Q_(100,'USD'),'EUR',{'USD/EUR':{'bid':D('.8'),'ask':D('.9')}})


def test_gwp_policy_cannot_be_inferred_from_mass_unit_compatibility():
    with pytest.raises((u.UnsupportedConversionError,ValueError,TypeError)):
        u.convert(u.Q_(1,'kg'),'tonne_CO2e',gwp=25)


def test_portable_quantity_serialiser_rejects_assay_and_callback_objects():
    for value in (u.CrudeAssay.brent(),u.CrudeAssay.brent().as_context()):
        with pytest.raises(TypeError): u.quantity_to_dict(value)


def test_per_row_assay_dispatch_is_not_guessed_by_scalar_batch_api():
    with pytest.raises((u.UnsupportedConversionError,TypeError,ValueError)):
        convert_many([u.Q_(1,'bbl'),u.Q_(1,'bbl')],'kg',assays=[u.CrudeAssay.brent(),u.CrudeAssay.wti()])
