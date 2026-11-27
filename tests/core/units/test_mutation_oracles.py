"""Independent oracles used by the seeded mutation campaign and normal tests."""
from decimal import Decimal as D
import pytest
import core.units as u


def test_loss_arithmetic():
    assert u.apply_energy_losses(u.Q_(100,'MWh_th'),[u.EnergyLoss('loss',D('.1'))]).delivered_energy.magnitude == 90


def test_barrel_factor():
    # Independent US petroleum barrel definition: 42 gallons, each 231 in³.
    exact = D(42)*D(231)*D('.0254')**3
    assert D(str(u.BARREL_TO_M3)) == exact
    assay=u.CrudeAssay.brent()
    actual=u.convert(u.Q_(1,'bbl'),'kg',assay=assay).magnitude
    assert float(actual) == pytest.approx(.158987294928*assay.density_kg_per_m3,rel=1e-13)


def test_inverse_generation():
    assert u.GenerationProcess(efficiency=D('.5')).convert(u.Q_(D(50),'MWh_e'),'MWh_th').magnitude == 100


def test_valid_output_basis():
    assert u.GenerationProcess(efficiency=.5,output_basis='net').output_basis == 'net'


def test_efficiency_one_is_valid():
    assert u.GenerationProcess(efficiency=1).thermal_efficiency == 1


def test_generation_conditional():
    assert u.GenerationProcess(efficiency=D('.5')).convert(u.Q_(100,'MWh_th'),'MWh_e').magnitude == 50


def test_exactly_one_process_parameter():
    assert u.GenerationProcess(efficiency=.5).thermal_efficiency == D('.5')


def test_invalid_output_basis():
    with pytest.raises(ValueError): u.GenerationProcess(efficiency=.5,output_basis='unknown')


def test_zero_cop_rejected():
    with pytest.raises(ValueError): u.HeatPumpSpecification(0)


def test_currency_half_even_rounding():
    assert u.round_currency(u.Q_(D('1.005'),'USD')).magnitude == D('1.00')


def test_solver_residual():
    value=u.solve_z_factor(u.Q_(10000,'kPa'),u.Q_(300,'K'),method='dak')
    assert value.converged and value.residual < 1e-8


def test_fx_direction():
    assert u.convert_fx(u.Q_(D(100),'USD'),'EUR',{'USD/EUR':D('.8')}).magnitude == 80


def test_chp_conservation():
    balance=u.CHPSpecification(D('.3'),D('.4'),'hhv').outputs(u.Q_(D(100),'MWh_th'),fuel_basis='hhv')
    assert balance.rejected_heat.magnitude == 30
    assert balance.energy_balance_residual.magnitude == 0


def test_frozen_registry_is_reported():
    # Keep the normal suite's registry unfrozen; exercise freezing in a child.
    import subprocess,sys
    result=subprocess.run([sys.executable,'-B','-c','import core.units as u; u.freeze_registry(); assert u.is_registry_frozen()'],capture_output=True,text=True,timeout=25)
    assert result.returncode == 0,result.stderr
