"""Economic energy roles survive every general-purpose conversion route."""
from decimal import Decimal
import pytest
import pint

from core.units import (Q_, ureg, convert, GenerationProcess, CrudeAssay,
    quantity_to_dict, quantity_from_dict, HeatPumpSpecification,
    CHPSpecification, thermal_to_output, output_to_thermal)
from core.units._safety import authorized_conversion


@pytest.mark.parametrize("source,target", [
    ("MWh_th", "MWh_e"), ("MWh_th", "MWh"), ("MW_e", "MW"),
    ("USD/MWh_th", "USD/MWh_e"), ("MWh_th/day", "MWh_e/day"),
    ("heat_rate", "dimensionless")])
def test_roles_cannot_be_erased(source, target):
    quantity = Q_(100, source)
    assert not quantity.is_compatible_with(target)
    for value in (quantity, quantity.to_base_units(), quantity.to_reduced_units(),
                  quantity_from_dict(quantity_to_dict(quantity))):
        with pytest.raises(pint.DimensionalityError):
            value.to(target)
    with pytest.raises((pint.DimensionalityError, ValueError)):
        convert(quantity, target, assay=CrudeAssay.brent())
    with pytest.raises(pint.DimensionalityError):
        authorized_conversion(lambda: quantity.to(target))()


def test_mixed_energy_cannot_be_added():
    with pytest.raises(pint.DimensionalityError):
        Q_(100, "MWh_th") + Q_(100, "MWh_e")


@pytest.mark.parametrize("number", [100., Decimal(100)])
def test_generation_amount_rate_and_quote_parity(number):
    process = GenerationProcess(efficiency=Decimal(".4"), output_basis="net", fuel_basis="lhv")
    for source, target in (("MWh_th", "MWh_e"), ("MW_th", "MW_e")):
        value = Q_(number, source)
        output = process.convert(value, target)
        assert output.magnitude == 40
        assert isinstance(output.magnitude, type(number))
        assert process.convert(output, source).magnitude == number
    quote = process.convert(Q_(number, "USD/MWh_th"), "USD/MWh_e")
    assert quote.magnitude == 250
    assert process.convert(quote, "USD/MWh_th").magnitude == number
    heat_rate = Q_(Decimal("2.5"), "MWh_th/MWh_e")
    assert GenerationProcess(heat_rate=heat_rate).thermal_efficiency == Decimal(".4")
    assert thermal_to_output(Q_(number, "MWh_th"), heat_rate=heat_rate).magnitude == 40
    assert output_to_thermal(Q_(number, "MWh_e"), heat_rate=heat_rate).magnitude == 250
    with pytest.raises(ValueError, match="FX"):
        process.convert(Q_(number, "USD/MWh_th"), "EUR/MWh_e")


def test_processes_reject_conflicting_roles_and_close_energy_balance():
    plant = CHPSpecification(.35, .45, "hhv", .1)
    balance = plant.outputs(Q_(100, "MWh_th"), fuel_basis="hhv")
    assert balance.energy_balance_residual.magnitude == 0
    assert balance.net_electricity.is_compatible_with("MWh_e")
    assert balance.useful_heat.is_compatible_with("MWh_th")
    with pytest.raises(pint.DimensionalityError):
        plant.outputs(Q_(100, "MWh_e"), fuel_basis="hhv")
    with pytest.raises(pint.DimensionalityError):
        HeatPumpSpecification(4).useful_heat(Q_(100, "MWh_th"))
    with pytest.raises(pint.DimensionalityError):
        thermal_to_output(Q_(100, "MWh_e"), efficiency=.4)
    with pytest.raises(ValueError):
        GenerationProcess(heat_rate=Q_(2.5, "MWh_e/MWh_th"))


def test_fuel_equivalents_accept_thermal_and_reject_electricity():
    from core.units.converters.energy import energy_to_gge, gge_to_energy
    quantity = Q_(Decimal(1), "MWh_lhv")
    equivalent = energy_to_gge(quantity, "EPA")
    assert abs(gge_to_energy(equivalent, target_unit="MWh_lhv").magnitude - 1) < Decimal("1e-12")
    with pytest.raises(pint.DimensionalityError):
        energy_to_gge(Q_(1, "MWh_e"), "EPA", basis="lhv")
