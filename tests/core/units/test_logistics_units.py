"""Contract multipliers are explicit; shipping indices retain their identities."""
from decimal import Decimal as D
import pytest
import pint
import core.units as u
from core.units.logistics_units import LOGISTICS_UNITS, logistics_unit


def test_contract_lot_tick_notional_and_capacity_quotation_contracts():
    contracts = u.Q_(D(5),'contract')
    multiplier = u.Q_(D(1000),'bbl/contract')
    price = u.Q_(D('-20'),'USD/bbl')
    assert (contracts*multiplier*price).to('USD').magnitude == -100000
    tick = u.Q_(D('.01'),'USD/bbl/price_tick')
    assert (u.Q_(3,'price_tick')*tick*contracts*multiplier).to('USD').magnitude == 150
    assert (u.Q_(2,'FEU')).to('TEU').magnitude == 4
    assert (u.Q_(30,'EUR/MW')*u.Q_(2,'MW')).to('EUR').magnitude == 60


def test_worldscale_and_packaging_require_their_declared_market_factors():
    flat_rate = u.Q_(D(12),'USD/mt')
    assert (u.Q_(D(75),'WS')/u.Q_(D(100),'WS')*flat_rate).to('USD/mt').magnitude == 9
    with pytest.raises(pint.DimensionalityError): u.Q_(75,'WS').to('USD/mt')
    assert all(logistics_unit(name) == u.parse_unit(name) for name in LOGISTICS_UNITS)
    with pytest.raises(ValueError): logistics_unit('made-up-package')
