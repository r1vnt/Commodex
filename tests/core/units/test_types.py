"""Unit contracts presented to future calendar, valuation, portfolio and risk modules."""
from datetime import datetime, timezone
from zoneinfo import ZoneInfo
from decimal import Decimal as D
import pytest
import pint
import core.units as u


@pytest.mark.parametrize('day,expected',[(datetime(2026,3,29),23),(datetime(2026,10,25),25)])
def test_calendar_delivery_duration_and_daylight_saving_unit_contract(day,expected):
    zone = ZoneInfo('Europe/London')
    start = day.replace(tzinfo=zone)
    end = day.replace(day=day.day+1,tzinfo=zone)
    duration = u.Q_((end.astimezone(timezone.utc)-start.astimezone(timezone.utc)).total_seconds(),'second').to('hour')
    assert duration.magnitude == expected
    assert (u.Q_(1,'MW_e')*duration).to('MWh_e').magnitude == expected


@pytest.mark.parametrize('name,source,target',[
    ('settlement','USD/MWh_e*MWh_e','USD'),('market_data','EUR/mt','EUR/kg'),
    ('forward_curve','USD/bbl','USD/liter'),('volatility','USD/MWh_e','USD/kWh_e'),
    ('fundamentals','kg/hour*hour','kg'),('power_valuation','USD/MWh_e*MW_e*hour','USD'),
    ('emissions_valuation','EUR/EUA*EUA','EUR'),('pricing_input','USD/mt*mt','USD'),
    ('spread','USD/MWh_e','USD/kWh_e'),('physical_asset','MW_e*hour','MWh_e'),
    ('trading','USD/bbl*bbl','USD'),('portfolio','mt','kg'),('greek','USD/(USD/MWh_e)','MWh_e'),
    ('risk_aggregation','USD','USD'),('instrument','contract*bbl/contract','bbl'),
    ('contract_multiplier','item*kg/item','kg')])
def test_future_module_quantity_contracts(name,source,target):
    assert u.Q_(D(1),source).to(target).magnitude > 0


def test_power_spread_contract_applies_generation_process_before_aggregation():
    fuel_cost = u.GenerationProcess(efficiency=D('.5')).convert(u.Q_(D(30),'EUR/MWh_th'),'EUR/MWh_e')
    spread = u.Q_(D(100),'EUR/MWh_e')-fuel_cost
    assert spread.magnitude == 40
    with pytest.raises(pint.DimensionalityError):
        u.Q_(100,'EUR/MWh_e')-u.Q_(30,'EUR/MWh_th')
