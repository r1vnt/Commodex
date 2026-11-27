"""Explicit methods, thermodynamic reference states and closed energy balances."""
from dataclasses import replace
from decimal import Decimal as D
import importlib.util
import unittest

import pint
import core.units as u


Q = u.Q_
CP = u.ThermodynamicBasis.CONSTANT_PRESSURE
CV = u.ThermodynamicBasis.CONSTANT_VOLUME


def fuel(gross=20, net=18, **metadata):
    return u.HeatContentAssay(gcv=Q(D(gross), "MJ/kg"), ncv=Q(D(net), "MJ/kg"), basis="DB",
        thermodynamic_basis=CP, combustion_temperature=Q(D("298.15"), "K"), **metadata)


def ultimate():
    return u.UltimateAnalysis(carbon_pct=60, hydrogen_fuel_bound_pct=5, oxygen_pct=10,
        nitrogen_pct=2, sulphur_pct=3, ash_pct=10, moisture_pct=10)


class CalorificMethodTests(unittest.TestCase):
    def assertClose(self, actual, expected):
        self.assertAlmostEqual(float(actual.to(expected.units).magnitude), float(expected.magnitude), places=10)

    def test_no_implicit_missing_value_or_moisture_approximation(self):
        for metadata in ({}, {"hydrogen_pct": 5, "moisture_pct": 10}, {"ncv_to_gcv_ratio": D(".9")}):
            spec = u.HeatContentAssay(gcv=Q(20, "MJ/kg"), **metadata)
            with self.subTest(metadata=metadata), self.assertRaisesRegex(ValueError, "select method"):
                spec.calorific_value("ncv")
        spec = u.HeatContentAssay(ncv=Q(18, "MJ/kg"), moisture_pct=10)
        with self.assertRaisesRegex(ValueError, "moisture_correction_method"):
            spec.calorific_value("ncv", basis="DB")
        self.assertEqual(replace(spec, moisture_correction_method="latent_water").calorific_value("ncv", basis="DB").magnitude, D("20.27133333333333333333333333"))
        with self.assertRaisesRegex(ValueError, "method"):
            u.dry_to_as_received(Q(20, "MJ/kg"), 10, kind="ncv")

    def test_public_exports_and_temperature_dependent_opt_in_validation(self):
        from core.units import calorific_models, energy_systems
        for module in (calorific_models, energy_systems):
            self.assertEqual(len(module.__all__), len(set(module.__all__)))
            for name in module.__all__:
                self.assertIs(getattr(module, name), getattr(u, name))
        with self.assertRaisesRegex(ValueError, "Temperature-dependent"):
            u.HeatContentAssay(gcv=Q(20, "MJ/kg"), latent_heat_method="iapws95")
        with self.assertRaisesRegex(ValueError, "constant-volume"):
            u.HeatContentAssay(gcv=Q(20, "MJ/kg"), hydrogen_pct=5, moisture_pct=0, method="latent_water", thermodynamic_basis=CV).calorific_value()

    def test_selected_methods_have_distinct_results_and_provenance(self):
        measured = u.HeatContentAssay(gcv=Q(D(20), "MJ/kg"), ncv=Q(D(19), "MJ/kg"),
            hydrogen_pct=5, moisture_pct=10, method="latent_water", ncv_to_gcv_ratio=D(".8"))
        self.assertEqual(measured.calorific_value().magnitude, 19)
        self.assertEqual(measured.calorific_method(), u.CalorificMethod.MEASURED)
        latent = replace(measured, ncv=None)
        self.assertEqual(latent.calorific_value().magnitude, D("18.6647144"))
        self.assertEqual(latent.calorific_method(), u.CalorificMethod.LATENT_WATER_APPROXIMATION)
        ratio = replace(latent, method="explicit_ratio")
        self.assertEqual(ratio.calorific_value().magnitude, 16)
        resolved = latent.for_basis("DB")
        self.assertEqual(resolved.ncv_method, u.CalorificMethod.LATENT_WATER_APPROXIMATION)
        self.assertEqual(resolved.gcv_method, u.CalorificMethod.MEASURED)
        selected = replace(latent, basis="DB", moisture_pct=0).for_basis("DB", kind="ncv")
        self.assertIsNone(selected.gcv)
        self.assertEqual(selected.ncv_method, u.CalorificMethod.LATENT_WATER_APPROXIMATION)
        with self.assertRaises(ValueError):
            replace(latent, method="unknown")
        with self.assertRaises(ValueError):
            replace(latent, method="explicit_ratio", ncv_to_gcv_ratio=None)

    def test_gar_nar_enforce_as_received_in_density_context_and_router(self):
        spec = fuel(as_received_moisture_pct=20, moisture_correction_method="latent_water")
        self.assertEqual(spec.calorific_value("GAR").magnitude, 16)
        self.assertEqual(spec.calorific_value("NAR").magnitude, D("13.9116"))
        self.assertEqual(spec.energy_density("MJ/kg", kind="GAR").magnitude, 16)
        for kind in ("gar", "nar"):
            for number in (2., D(2)):
                amount = Q(number, "kg")
                expected = spec.calorific_value(kind) * D(2)
                actual = u.convert(amount, "MJ", assay=spec, basis=kind)
                self.assertClose(actual, expected * Q(1, "kg"))
                self.assertIsInstance(actual.magnitude, type(number))
        with self.assertRaisesRegex(ValueError, "AS RECEIVED"):
            spec.energy_density("MJ/kg", kind="gar", basis="DB")
        with self.assertRaisesRegex(ValueError, "Moisture"):
            fuel().calorific_value("GAR")
        gross_only = fuel(as_received_moisture_pct=20)
        self.assertEqual(u.convert(Q(2, "kg"), "MJ", assay=gross_only, basis="gar").magnitude, 32)
        self.assertEqual(u.convert(Q(50, "USD/GJ"), "USD/metric_ton", assay=gross_only, basis="gar").magnitude, 800)
        self.assertIsNone(gross_only.for_basis("AR", kind="gcv").ncv)
        with self.assertRaisesRegex(ValueError, "moisture_correction_method"):
            u.convert(Q(2, "kg"), "MJ", assay=gross_only, basis="nar")

    def test_combustion_temperature_is_independent_of_metering_temperature(self):
        ref = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1))
        spec = u.HeatContentAssay.from_gas_composition(u.GasComposition({"CH4": D(1)}),
            dry_reference=ref, wet_reference=replace(ref, water_mole_fraction=D(".02")))
        self.assertEqual(spec.method, u.CalorificMethod.COMPONENT_SUM)
        self.assertEqual(spec.thermodynamic_basis, CP)
        self.assertEqual(spec.combustion_temperature.to("degC").magnitude, 15)
        self.assertEqual(spec.reference_temperature.to("degC").magnitude, 15)
        warmer_meter = replace(ref, temperature=Q(D("303.15"), "K"))
        changed = replace(spec, dry_gas_state=replace(spec.dry_gas_state, reference=warmer_meter),
                          density=None, reference_temperature=None)
        self.assertEqual(changed.combustion_temperature, spec.combustion_temperature)
        self.assertNotEqual(changed.reference_temperature, spec.reference_temperature)

    def test_ultimate_and_proximate_analyses_validate_closure_and_hydrogen(self):
        analysis = ultimate()
        self.assertEqual(analysis.moisture_bound_hydrogen_pct, D(10) / D("8.936"))
        self.assertEqual(analysis.total_hydrogen_pct, D(5) + D(10) / D("8.936"))
        proximate = u.ProximateAnalysis(moisture_pct=10, ash_pct=10, volatile_matter_pct=30, fixed_carbon_pct=50)
        spec = u.HeatContentAssay(gcv=Q(25, "MJ/kg"), ultimate_analysis=analysis, proximate_analysis=proximate, method="latent_water")
        self.assertEqual(spec.hydrogen_pct, 5)
        self.assertEqual(spec.moisture_pct, 10)
        self.assertClose(spec.calorific_value(), Q(D(25) - D("2.442") * (D("8.936") * D(".05") + D(".1")), "MJ/kg"))
        for kwargs in ({"carbon_pct": 61}, {"moisture_pct": float("nan")}, {"hydrogen_fuel_bound_pct": True}, {"basis": "DB"}):
            with self.subTest(kwargs=kwargs), self.assertRaises((ValueError, TypeError)):
                replace(analysis, **kwargs)
        with self.assertRaises(ValueError):
            replace(spec, proximate_analysis=replace(proximate, moisture_pct=11, fixed_carbon_pct=49))
        with self.assertRaises(ValueError):
            replace(spec, hydrogen_pct=6)

    def test_empirical_correlation_requires_explicit_coefficients_and_ranges(self):
        correlation = u.EmpiricalHeatCorrelation("illustrative declared correlation", "test fixture; not a published standard",
            {"carbon_pct": D(30), "hydrogen_fuel_bound_pct": D(100)}, "AR", Q(25, "degC"), CP,
            {"carbon_pct": (D(50), D(70)), "hydrogen_fuel_bound_pct": (D(3), D(7))})
        spec = u.HeatContentAssay.from_ultimate_analysis(ultimate(), correlation=correlation)
        self.assertEqual(spec.gcv.magnitude, 23)
        self.assertEqual(spec.method, u.CalorificMethod.EMPIRICAL_CORRELATION)
        with self.assertRaisesRegex(ValueError, "select method"):
            spec.calorific_value("ncv")
        with self.assertRaisesRegex(ValueError, "applicability"):
            correlation.evaluate(replace(ultimate(), carbon_pct=40, other_pct=20))
        with self.assertRaises(ValueError):
            replace(correlation, valid_ranges_pct={})

    def test_laboratory_results_retain_method_report_and_thermodynamic_basis(self):
        report = u.LaboratoryCalorificResult("declared laboratory method/edition", "report-123", "laboratory source",
            "DB", CV, Q(25, "degC"), gcv=Q(D(25), "MJ/kg"))
        spec = u.HeatContentAssay.from_laboratory(report, name="reported sample")
        self.assertEqual(spec.method, u.CalorificMethod.STANDARD_LABORATORY)
        self.assertEqual(spec.laboratory_report_id, "report-123")
        self.assertEqual(spec.thermodynamic_basis, CV)
        with self.assertRaisesRegex(ValueError, "select method"):
            spec.calorific_value("ncv")
        with self.assertRaises(ValueError):
            u.HeatContentAssay.from_laboratory(report, gcv=Q(100, "MJ/kg"))
        with self.assertRaises(ValueError):
            replace(spec, laboratory_standard=None)

    def test_constant_pressure_volume_identity_sign_units_and_roundtrip(self):
        # CH4(g)+2 O2(g) -> CO2(g)+2 H2O(l), delta gas moles = -2.
        correction = u.ideal_gas_combustion_correction(Q(D(-2), "mol/mol"), Q(25, "degC"))
        for magnitude in (885., D(885)):
            internal = Q(magnitude, "kJ/mol")
            enthalpy = u.convert_combustion_heat(internal, source=CV, target=CP, pressure_volume_correction=correction)
            expected = D(885) + 2 * D("8.314462618") * D("298.15") / 1000
            self.assertClose(enthalpy, Q(expected, "kJ/mol"))
            restored = u.convert_combustion_heat(enthalpy, source=CP, target=CV, pressure_volume_correction=correction)
            self.assertClose(restored, internal)
            self.assertIsInstance(enthalpy.magnitude, type(magnitude))
        spec = replace(fuel(), thermodynamic_basis=CV)
        moved = spec.for_thermodynamic_basis(CP, gcv_pv_correction=Q(-1, "MJ/kg"), ncv_pv_correction=Q(0, "MJ/kg"))
        self.assertEqual(moved.gcv.magnitude, 21)
        self.assertEqual(moved.ncv.magnitude, 18)
        self.assertEqual(moved.thermodynamic_basis, CP)
        self.assertTrue(moved.derivations)
        with self.assertRaises(ValueError):
            spec.for_thermodynamic_basis(CP, gcv_pv_correction=Q(-1, "MJ/kg"))
        with self.assertRaises(ValueError):
            u.convert_combustion_heat(Q(10, "MJ/kg"), source="unspecified", target=CP, pressure_volume_correction=Q(1, "MJ/kg"))

    def test_iso_component_table_path_uses_declared_reference_data_and_spectator_water(self):
        table = u.ISO6976Table("2016/test-fixture", "explicit illustrative values", Q(25, "degC"),
            {"CH4": (D(890), D(802)), "N2": (D(0), D(0))})
        composition = u.GasComposition({"CH4": D(".9"), "N2": D(".1")})
        values = table.calculate(composition)
        self.assertEqual(values.method, u.CalorificMethod.COMPONENT_SUM)
        self.assertClose(values.hhv * composition.molar_mass, Q(D(801), "kJ/mol"))
        wet = table.calculate(composition.with_water(D(".1")))
        self.assertClose(wet.hhv * composition.with_water(D(".1")).molar_mass, Q(D("720.9"), "kJ/mol"))
        ref = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1))
        spec = u.HeatContentAssay.from_gas_composition(composition, dry_reference=ref,
            wet_reference=replace(ref, water_mole_fraction=D(".1")), component_table=table)
        self.assertEqual(spec.combustion_temperature.to("degC").magnitude, 25)
        self.assertEqual(spec.reference_temperature.to("degC").magnitude, 15)
        self.assertClose(u.convert(Q(D(1), "mol"), "kJ", assay=spec, basis="hhv"), Q(D(801), "kJ"))
        with self.assertRaises(ValueError):
            table.calculate(u.GasComposition({"C2H6": D(1)}))
        with self.assertRaises(ValueError):
            replace(table, component_values_kj_per_mol={"H2O": (D(44), D(0))})

    @unittest.skipUnless(importlib.util.find_spec("CoolProp"), "optional CoolProp not installed")
    def test_temperature_dependent_latent_heat_actual_iapws_backends(self):
        for method in ("iapws95", "if97"):
            at25 = u.water_latent_heat(Q(25, "degC"), method=method)
            at100 = u.water_latent_heat(Q(100, "degC"), method=method)
            self.assertAlmostEqual(float(at25.magnitude), 2.4417, places=3)
            self.assertAlmostEqual(float(at100.magnitude), 2.2564, places=3)
            spec = u.HeatContentAssay(gcv=Q(20, "MJ/kg"), hydrogen_pct=5, moisture_pct=10,
                method="latent_water", latent_heat_method=method, combustion_temperature=Q(100, "degC"))
            self.assertEqual(spec.latent_heat_mj_per_kg_water, at100.magnitude)
            self.assertClose(spec.calorific_value(), Q(D(20) - at100.magnitude * D(".5468"), "MJ/kg"))
        with self.assertRaises(ValueError):
            u.water_latent_heat(Q(700, "K"), method="iapws95")


class EnergySystemTests(unittest.TestCase):
    def assertClose(self, actual, expected):
        self.assertAlmostEqual(float(actual.to(expected.units).magnitude), float(expected.magnitude), places=10)

    def test_efficiency_hhv_lhv_switch_and_gross_net_generation(self):
        spec = fuel(50, 45)
        eta = u.convert_thermal_efficiency(D(".4"), assay=spec, source_basis="hhv", target_basis="lhv")
        self.assertEqual(eta, D(".4") * 50 / 45)
        self.assertAlmostEqual(u.convert_thermal_efficiency(eta, assay=spec, source_basis="lhv", target_basis="hhv"), D(".4"), places=25)
        self.assertGreater(u.convert_thermal_efficiency(D(".99"), assay=spec, source_basis="hhv", target_basis="lhv"), 1)
        volume_only = u.HeatContentAssay(gcv=Q(D(1000), "btu_it/scf"), ncv=Q(D(900), "btu_it/scf"))
        self.assertIsNone(volume_only.standard_density)
        self.assertEqual(u.convert_thermal_efficiency(D(".4"), assay=volume_only, source_basis="hhv", target_basis="lhv"), D(".4") * 1000 / 900)
        for numeric in (100., D(100)):
            for unit in ("MW_e", "MWh_e"):
                gross = Q(numeric, unit)
                net = u.net_generation(gross, parasitic_fraction=D(".05"))
                self.assertEqual(net.magnitude, 95)
                self.assertEqual(u.gross_generation(net, parasitic_fraction=D(".05")).magnitude, 100)
                self.assertIsInstance(net.magnitude, type(numeric))
        self.assertClose(u.net_generation(Q(100, "MWh"), parasitic_load=Q(5000, "kWh")), Q(95, "MWh_e"))
        with self.assertRaises(ValueError):
            u.net_generation(Q(100, "MW"), parasitic_load=Q(101, "MW"))
        with self.assertRaises(ValueError):
            u.gross_generation(Q(100, "MW"), parasitic_fraction=1)

    def test_chp_energy_closure_and_basis_agreement(self):
        plant = u.CHPSpecification(D(".35"), D(".45"), "hhv", parasitic_fraction=D(".1"))
        for numeric in (100., D(100)):
            for unit in ("MWh", "MW"):
                balance = plant.outputs(Q(numeric, unit), fuel_basis="hhv")
                self.assertEqual(balance.gross_electricity.magnitude, 35)
                self.assertEqual(balance.net_electricity.magnitude, 31.5)
                self.assertEqual(balance.useful_heat.magnitude, 45)
                self.assertEqual(balance.rejected_heat.magnitude, 20)
                self.assertClose(balance.energy_balance_residual, Q(0, "MWh_th" if unit == "MWh" else "MW_th"))
                self.assertFalse(balance.net_electricity.is_compatible_with(balance.useful_heat))
                self.assertAlmostEqual(float(balance.total_useful_efficiency), .765)
                self.assertIsInstance(balance.useful_heat.magnitude, type(numeric))
        with self.assertRaises(ValueError):
            plant.outputs(Q(100, "MWh"), fuel_basis="lhv")
        with self.assertRaises(ValueError):
            u.CHPSpecification(.7, .5, "hhv")

    def test_heat_pump_cop_can_exceed_one_and_power_roundtrip(self):
        pump = u.HeatPumpSpecification(D(4))
        for value in (Q(10., "MWh_e"), Q(D(10), "MW_e")):
            unit = "MWh_th" if value.is_compatible_with("MWh_e") else "MW_th"
            heat = pump.useful_heat(value, target_unit=unit)
            self.assertEqual(heat.magnitude, 40)
            self.assertClose(pump.electrical_input(heat, target_unit=value.units), value)
        for value in (0, -1, float("nan"), True):
            with self.subTest(value=value), self.assertRaises((ValueError, TypeError)):
                u.HeatPumpSpecification(value)

    def test_electrolysis_energy_mass_and_flow_roundtrips_with_explicit_scope(self):
        h2 = fuel(144, 120)
        electrolyser = u.ElectrolysisSpecification(Q(D(50), "kWh/kg"), "system AC input including auxiliaries", Q(0, "kWh/kg"))
        self.assertAlmostEqual(electrolyser.efficiency(hydrogen_heat_content=h2, basis="hhv"), D(".8"), places=14)
        self.assertClose(u.ElectrolysisSpecification.from_efficiency(efficiency=D(".8"), hydrogen_heat_content=h2, basis="hhv", source="declared design").electricity_per_kg_h2, Q(50, "kWh_e/kg"))
        for numeric in (1., D(1)):
            for source, target in (("MWh", "kg"), ("MW", "kg/hour")):
                electricity = Q(numeric, source)
                hydrogen = electrolyser.hydrogen_output(electricity, target_unit=target)
                self.assertEqual(hydrogen.magnitude, 20)
                self.assertClose(electrolyser.electricity_required(hydrogen, target_unit=source), electricity)
                self.assertIsInstance(hydrogen.magnitude, type(numeric))
        with self.assertRaises(ValueError):
            u.ElectrolysisSpecification(Q(0, "kWh/kg"), "bad input")

    def test_carrier_combustion_energy_is_separate_from_hydrogen_recovery_and_cracking_load(self):
        cracker = u.HydrogenCrackingSpecification(D(".18"), D(".9"), fuel(144, 120), Q(D(3), "MJ/kg"), fuel(22, 19))
        for numeric in (100., D(100)):
            balance = cracker.balance(Q(numeric, "kg"), basis="lhv")
            self.assertClose(balance.delivered_hydrogen_mass, Q(D("16.2"), "kg"))
            self.assertClose(balance.delivered_hydrogen_energy, Q(D(1944), "MJ_th"))
            self.assertClose(balance.carrier_combustion_energy, Q(D(1900), "MJ_th"))
            self.assertClose(balance.cracking_energy_input, Q(D(300), "MJ_th"))
        self.assertIsNone(replace(cracker, carrier_heat_content=None).balance(Q(100, "kg"), basis="lhv").carrier_combustion_energy)
        carrier = u.HydrogenCarrierAssay(name="declared carrier", carrier="LOHC", h2_mass_fraction=D(".18"),
            cracking_recovery_fraction=D(".9"), hydrogen_heat_content=fuel(144, 120), heat_content=fuel(22, 19))
        self.assertEqual(carrier.as_heat_content().energy_role, "recovered_hydrogen")
        self.assertClose(carrier.combustion_energy(Q(100, "kg"), target_unit="MJ"), Q(1900, "MJ"))
        self.assertClose(carrier.delivered_energy(Q(100, "kg"), target_unit="MJ"), Q(1944, "MJ"))
        self.assertClose(carrier.cracking_balance(Q(100, "kg"), cracking_energy_per_kg_carrier=Q(3, "MJ/kg")).cracking_energy_input, Q(300, "MJ_th"))
        with self.assertRaises(ValueError):
            u.convert_thermal_efficiency(.5, assay=carrier, source_basis="hhv", target_basis="lhv")

    def test_sequential_named_losses_close_and_invert_without_adding_fractions(self):
        stages = (u.EnergyLoss("transport", D(".1")), u.EnergyLoss("conversion", D(".2")), u.EnergyLoss("storage", D(".05")))
        for numeric in (100., D(100)):
            balance = u.apply_energy_losses(Q(numeric, "MWh"), stages)
            self.assertAlmostEqual(float(balance.delivered_energy.magnitude), 68.4, places=12)
            self.assertClose(balance.delivered_energy + sum((value for _, value in balance.stage_losses), Q(0, "MWh")), balance.input_energy)
            self.assertClose(u.energy_before_losses(balance.delivered_energy, stages), balance.input_energy)
        with self.assertRaises(ValueError):
            u.energy_before_losses(Q(1, "MJ"), (u.EnergyLoss("total", 1),))
        with self.assertRaises(ValueError):
            u.apply_energy_losses(Q(1, "MJ"), stages + (stages[0],))

    def test_mass_and_molar_blending_agree_when_weights_represent_same_amounts(self):
        first = fuel(20, 18, molar_mass=Q(D(".02"), "kg/mol"))
        second = fuel(40, 36, molar_mass=Q(D(".04"), "kg/mol"))
        by_mass = u.blend_heat_content((first, second), (Q(1, "kg"), Q(3, "kg")), weighting="mass")
        by_moles = u.blend_heat_content((first, second), (Q(50, "mol"), Q(75, "mol")), weighting="molar")
        self.assertEqual(by_mass.gcv.magnitude, 35)
        self.assertEqual(by_mass.ncv.magnitude, D("31.5"))
        self.assertClose(by_mass.energy_density("MJ/kg"), by_moles.energy_density("MJ/kg"))
        self.assertEqual(by_mass.molar_mass.magnitude, D(".032"))
        self.assertClose(u.convert(Q(4, "kg"), "MJ", assay=by_mass, basis="hhv"), Q(140, "MJ"))
        self.assertEqual(by_mass.method, u.CalorificMethod.COMPONENT_SUM)
        with self.assertRaises(ValueError):
            u.blend_heat_content((first, replace(second, thermodynamic_basis=CV)), (.25, .75), weighting="mass")
        with self.assertRaises(ValueError):
            u.blend_heat_content((first, second), (.5, .6), weighting="mass")

    def test_native_molar_blend_and_measured_ar_blend_do_not_infer_unknown_metadata(self):
        first = u.HeatContentAssay(gcv=Q(100, "kJ/mol"), basis="DB", thermodynamic_basis=CP, combustion_temperature=Q(25, "degC"))
        second = replace(first, gcv=Q(200, "kJ/mol"))
        blend = u.blend_heat_content((first, second), (.5, .5), weighting="molar")
        self.assertIsNone(blend.molar_mass)
        self.assertClose(blend.energy_density("kJ/mol", kind="hhv"), Q(150, "kJ/mol"))
        ar = [u.HeatContentAssay(gcv=Q(value, "MJ/kg"), thermodynamic_basis=CP, combustion_temperature=Q(25, "degC")) for value in (20, 30)]
        blend = u.blend_heat_content(ar, (.5, .5), weighting="mass")
        self.assertEqual(blend.gcv.magnitude, 25)
        self.assertIsNone(blend.moisture_pct)

    def test_wet_gas_blending_requires_new_z_state_and_respects_water_composition(self):
        ref = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1))
        specs = [u.HeatContentAssay.from_gas_composition(u.GasComposition(recipe), dry_reference=ref,
            wet_reference=replace(ref, water_mole_fraction=water), basis="WET_GAS")
            for recipe, water in (({"CH4": D(1)}, D(".1")), ({"CH4": D(".5"), "N2": D(".5")}, D(".2")))]
        molar = sum((spec.molar_mass for spec in specs), Q(D(0), "kg/mol")) / 2
        water = D(".15")
        # Use the registry's pinned molecular mass rather than another water table.
        from core.units._gas_data import GAS_MOLAR_MASSES
        dry_molar = (molar - Q(water * GAS_MOLAR_MASSES["H2O"] / 1000, "kg/mol")) / (1-water)
        state = u.GasMoistureState(replace(ref, water_mole_fraction=water, compressibility_factor=D(".98")), dry_molar)
        blend = u.blend_heat_content(specs, (.5, .5), weighting="molar", gas_state=state)
        self.assertClose(blend.energy_density("MJ/mol", kind="hhv"), sum((spec.energy_density("MJ/mol", kind="hhv") for spec in specs), Q(D(0), "MJ/mol")) / 2)
        with self.assertRaisesRegex(ValueError, "explicit blend"):
            u.blend_heat_content(specs, (.5, .5), weighting="molar")


if __name__ == "__main__":
    unittest.main()
