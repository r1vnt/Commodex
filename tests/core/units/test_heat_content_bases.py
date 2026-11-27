"""Analytical bases, explicit fuel states and wet/dry gas calorific parity."""

from dataclasses import FrozenInstanceError, replace
from decimal import Decimal as D
import unittest

import pint
import core.units as u


SOLID_BASES = ("AR", "AD", "DB", "DAF", "DMMF", "AS_FIRED", "AS_DELIVERED", "REFERENCE_MOISTURE")


def solid_spec():
    return u.HeatContentAssay(name="declared solid fuel", gcv=u.Q_(D(20), "MJ/kg"), ncv=u.Q_(D(18), "MJ/kg"),
        moisture_pct=D(20), air_dried_moisture_pct=D(10), ash_db_pct=D(10), mineral_matter_db_pct=D(20), method="latent_water",
        mineral_matter_gcv_correction=u.Q_(D(1), "MJ/kg"), mineral_matter_ncv_correction=u.Q_(D(".5"), "MJ/kg"),
        density=u.Q_(D(650), "kg/m^3"),
        as_fired_state=u.FuelCondition(D(5), temperature=u.Q_(D(40), "degC"), density=u.Q_(D(700), "kg/m^3"), label="burner feed"),
        as_delivered_state=u.FuelCondition(D(25), density=u.Q_(D(600), "kg/m^3"), label="delivery survey"),
        reference_moisture_state=u.FuelCondition(D(12), label="contract moisture 12%"))


def gas_spec(*, dry_z=D(1), wet_z=D(1), water=D(".1")):
    dry_reference = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=dry_z)
    wet_reference = replace(dry_reference, compressibility_factor=wet_z, water_mole_fraction=water)
    composition = u.GasComposition({"CH4": D(".9"), "N2": D(".1")})
    return u.HeatContentAssay.from_gas_composition(composition, dry_reference=dry_reference, wet_reference=wet_reference), composition


class SolidCalorificBasisTests(unittest.TestCase):
    def assertQuantityClose(self, actual, expected):
        self.assertAlmostEqual(float(actual.to(expected.units).magnitude), float(expected.magnitude), places=10)

    def test_dmmf_uses_mineral_matter_not_ash_and_supports_explicit_calorific_corrections(self):
        plain = u.HeatContentAssay(gcv=u.Q_(D(20), "MJ/kg"), ncv=u.Q_(D(18), "MJ/kg"), basis="DB", ash_db_pct=10, mineral_matter_db_pct=20)
        self.assertQuantityClose(plain.calorific_value("gcv", basis="DAF"), u.Q_(D(20) / D(".9"), "MJ/kg"))
        self.assertEqual(plain.calorific_value("gcv", basis="DMMF").magnitude, 25)
        corrected = replace(plain, mineral_matter_gcv_correction=u.Q_(D(1), "MJ/kg"), mineral_matter_ncv_correction=u.Q_(D(".5"), "MJ/kg"))
        self.assertEqual(corrected.calorific_value("gcv", basis="DMMF").magnitude, D("23.75"))
        self.assertEqual(corrected.calorific_value("ncv", basis="DMMF").magnitude, D("21.875"))
        round_trip = corrected.for_basis("DMMF").for_basis("DB")
        self.assertEqual(round_trip.gcv.magnitude, 20)
        self.assertEqual(round_trip.ncv.magnitude, 18)
        with self.assertRaisesRegex(ValueError, "ash is not mineral matter"):
            replace(plain, mineral_matter_db_pct=None).calorific_value("gcv", basis="DMMF")
        with self.assertRaisesRegex(ValueError, "ash_db_pct"):
            replace(plain, ash_db_pct=None).calorific_value("gcv", basis="DAF")

    def test_maf_and_adb_are_aliases_while_named_fuel_states_remain_distinct(self):
        assay = solid_spec()
        self.assertEqual(assay.calorific_value("gcv", basis="MAF"), assay.calorific_value("gcv", basis="DAF"))
        self.assertEqual(assay.for_basis("MAF").basis, "DAF")
        self.assertEqual(assay.calorific_value("ncv", basis="ADB"), assay.calorific_value("ncv", basis="AD"))
        self.assertEqual(assay.calorific_value("gcv", basis="as fired").magnitude, D("23.75"))
        self.assertEqual(assay.calorific_value("gcv", basis="as-delivered").magnitude, D("18.75"))
        self.assertAlmostEqual(assay.calorific_value("gcv", basis="reference moisture").magnitude, D(22), places=20)
        self.assertEqual(assay.calorific_value("ncv", basis="AF").magnitude, D("21.832875"))
        self.assertEqual(assay.calorific_value("ncv", basis="AS_DELIVERED").magnitude, D("16.722375"))
        self.assertEqual(assay.calorific_value("ncv", basis="RM").magnitude, D("20.0442"))

    def test_all_solid_basis_pairs_round_trip_both_heating_values(self):
        original = solid_spec()
        for source in SOLID_BASES:
            projected = original.for_basis(source)
            for target in SOLID_BASES:
                for kind in ("gcv", "ncv"):
                    with self.subTest(source=source, target=target, kind=kind):
                        actual = projected.calorific_value(kind, basis=target)
                        self.assertQuantityClose(actual, original.calorific_value(kind, basis=target))
                        self.assertIsInstance(actual.magnitude, D)
            self.assertQuantityClose(projected.for_basis("AR").gcv, original.gcv)
        self.assertEqual(original.moisture_pct, 20)
        self.assertEqual(original.density.magnitude, 650)

    def test_firing_density_and_contexts_use_target_condition_without_reusing_source_density(self):
        original = solid_spec()
        fired = original.for_basis("AS_FIRED")
        delivered = original.for_basis("AS_DELIVERED")
        self.assertEqual(fired.density.magnitude, 700)
        self.assertEqual(delivered.density.magnitude, 600)
        self.assertEqual(fired.reference_temperature.to("degC").magnitude, 40)
        self.assertIsNone(original.for_basis("DB").density)
        self.assertQuantityClose(original.energy_density("MJ/m^3", kind="ncv", basis="AS_FIRED"), u.Q_(D("21.832875") * 700, "MJ/m^3"))
        for numeric in (1000., D(1000)):
            for projected in (fired, delivered, original.for_basis("DMMF")):
                with self.subTest(basis=projected.basis, numeric=type(numeric)):
                    amount = u.Q_(numeric, "kg")
                    energy = u.convert(amount, "GJ", assay=projected)
                    restored = u.convert(energy, "kg", assay=projected)
                    self.assertQuantityClose(restored, amount)
                    self.assertIsInstance(energy.magnitude, type(numeric))
            volume = u.Q_(numeric, "m^3")
            energy = u.convert(volume, "GJ", assay=fired)
            self.assertQuantityClose(energy, u.Q_(D(str(numeric)) * 700 * D("21.832875") / 1000, "GJ"))
        quote = u.convert(u.Q_(D(50), "USD/GJ"), "EUR/metric_ton", assay=fired, rates={"USD/EUR": D(".9")})
        self.assertQuantityClose(quote, u.Q_(D(50) * D(".9") * D("21.832875"), "EUR/metric_ton"))
        with self.assertRaisesRegex(ValueError, "new density"):
            original.energy_density("MJ/m^3", basis="DB")

    def test_master_router_analytical_density_conversion_and_condition_projection(self):
        assay = solid_spec()
        for numeric in (20., D(20)):
            original = u.Q_(numeric, "MJ/kg")
            result = u.convert(original, "MJ/kg", assay=assay, basis="gcv", analytical_basis="DMMF")
            self.assertEqual(result.magnitude, 30)
            self.assertIsInstance(result.magnitude, type(numeric))
            fired = u.convert(original, "MJ/m^3", assay=assay, basis="gcv", analytical_basis="AS_FIRED")
            self.assertEqual(fired.magnitude, D("23.75") * 700)
        with self.assertRaisesRegex(ValueError, "for_basis"):
            u.convert(u.Q_(1, "kg"), "GJ", assay=assay, analytical_basis="AS_FIRED")

    def test_missing_counterpart_is_resolved_before_rebasing_when_estimation_data_exists(self):
        assay = u.HeatContentAssay(gcv=u.Q_(D(20), "MJ/kg"), moisture_pct=10, hydrogen_pct=5,
                                  as_fired_state=u.FuelCondition(5), method="latent_water")
        fired = assay.for_basis("AS_FIRED")
        self.assertQuantityClose(fired.ncv, assay.calorific_value("ncv", basis="AS_FIRED"))
        self.assertIsNone(fired.hydrogen_pct)
        only_gross = replace(assay, hydrogen_pct=None).for_basis("AS_FIRED")
        self.assertIsNone(only_gross.ncv)

    def test_conditions_are_explicit_immutable_and_validate_metadata(self):
        with self.assertRaisesRegex(ValueError, "explicit"):
            u.HeatContentAssay(gcv=u.Q_(20, "MJ/kg"), basis="AS_FIRED", moisture_pct=5)
        with self.assertRaisesRegex(ValueError, "conflicts"):
            u.HeatContentAssay(gcv=u.Q_(20, "MJ/kg"), basis="AS_DELIVERED", moisture_pct=5, as_delivered_state=u.FuelCondition(10))
        density = u.Q_(D(700), "kg/m^3")
        condition = u.FuelCondition(D(10), density=density)
        density.ito("g/liter")
        condition.density.ito("g/liter")
        self.assertEqual(condition.density.units, u.ureg.kg / u.ureg.meter ** 3)
        with self.assertRaises(FrozenInstanceError):
            condition.moisture_pct = D(20)
        with self.assertRaises(ValueError):
            u.FuelCondition(10, absolute_pressure=u.Q_(1, "psig"))
        for moisture in (-1, 100, float("nan"), True):
            with self.subTest(moisture=moisture), self.assertRaises((ValueError, TypeError)):
                u.FuelCondition(moisture)
        for kwargs in ({"mineral_matter_db_pct": 100}, {"mineral_matter_db_pct": True},
                       {"mineral_matter_gcv_correction": u.Q_(1, "kg")}, {"mineral_matter_gcv_correction": u.Q_(float("nan"), "MJ/kg")}):
            with self.subTest(kwargs=kwargs), self.assertRaises((ValueError, TypeError, pint.PintError)):
                u.HeatContentAssay(gcv=u.Q_(20, "MJ/kg"), **kwargs)


class GasCalorificBasisTests(unittest.TestCase):
    def assertRelativeClose(self, actual, expected):
        self.assertAlmostEqual(float(actual.to(expected.units).magnitude / expected.magnitude), 1., places=11)

    def test_gross_and_net_molar_and_volume_dilution_without_spectator_latent_deduction(self):
        assay, composition = gas_spec()
        for kind in ("gcv", "ncv"):
            dry_molar = assay.energy_density("MJ/mol", kind=kind)
            wet_molar = assay.energy_density("MJ/mol", kind=kind, basis="WET_GAS")
            self.assertRelativeClose(wet_molar, dry_molar * D(".9"))
            dry_volume = assay.energy_density("MJ/Sm3", kind=kind)
            wet_volume = assay.energy_density("MJ/m^3", kind=kind, basis="WET_GAS")
            self.assertAlmostEqual(float(wet_volume.magnitude / dry_volume.magnitude), .9, places=12)
            wet_composition = composition.with_water(D(".1"))
            expected = wet_composition.calorific_values().hhv if kind == "gcv" else wet_composition.calorific_values().lhv
            self.assertRelativeClose(assay.calorific_value(kind, basis="WET_GAS"), expected)
        dry_ratio = assay.gcv.magnitude / assay.ncv.magnitude
        wet = assay.for_basis("WET_GAS")
        self.assertAlmostEqual(wet.gcv.magnitude / wet.ncv.magnitude, dry_ratio, places=25)
        self.assertEqual(u.GasComposition({"H2O": 1}).calorific_values().hhv.magnitude, 0)

    def test_compressibility_temperature_pressure_and_molar_mass_are_all_respected(self):
        assay, _ = gas_spec(dry_z=D(".99"), wet_z=D(".95"))
        wet_state = replace(assay.wet_gas_state, reference=replace(assay.wet_gas_state.reference, temperature=u.Q_(D(310), "K"), absolute_pressure=u.Q_(D(200000), "Pa")))
        assay = replace(assay, wet_gas_state=wet_state)
        expected_ratio = D(".9") * D(200000) / D(101325) * D("288.15") / 310 * D(".99") / D(".95")
        for kind in ("gcv", "ncv"):
            dry = assay.energy_density("MJ/m^3", kind=kind)
            wet = assay.energy_density("MJ/m^3", kind=kind, basis="WET_GAS")
            self.assertRelativeClose(wet, dry * expected_ratio)
        original = assay.for_basis("WET_GAS").for_basis("DRY_GAS")
        self.assertRelativeClose(original.gcv, assay.gcv)
        self.assertRelativeClose(original.ncv, assay.ncv)

    def test_contexts_amounts_flows_quotes_float_decimal_and_original_gas_engine_agree(self):
        dry, composition = gas_spec()
        wet = dry.for_basis("WET_GAS")
        backend = u.builtin_gas_backend("ideal", temperature_range=(u.Q_(270, "K"), u.Q_(330, "K")),
                                        pressure_range=(u.Q_(50000, "Pa"), u.Q_(200000, "Pa")))
        for assay, volume_unit, gas_composition, reference in (
            (dry, "Sm3", composition, dry.dry_gas_state.reference),
            (wet, "m^3", composition.with_water(D(".1")), wet.wet_gas_state.reference),
        ):
            for kind in ("hhv", "lhv"):
                for numeric in (1000., D(1000)):
                    with self.subTest(state=assay.basis, kind=kind, numeric=type(numeric)):
                        quantity = u.Q_(numeric, volume_unit)
                        actual = u.convert(quantity, "MWh", assay=assay, basis=kind)
                        expected = u.composition_gas_energy(quantity, gas_composition, reference, backend=backend, basis=kind, target_unit="MWh")
                        self.assertRelativeClose(actual, expected)
                        self.assertRelativeClose(u.convert(actual, volume_unit, assay=assay, basis=kind), quantity)
                        self.assertIsInstance(actual.magnitude, type(numeric))
                        flow = u.convert(u.Q_(numeric, volume_unit + "/day"), "MW", assay=assay, basis=kind)
                        self.assertRelativeClose(flow, u.Q_(actual.magnitude / 24, "MW"))
                        quote = u.convert(u.Q_(D(50), "EUR/MWh"), "EUR/" + volume_unit, assay=assay, basis=kind)
                        self.assertRelativeClose(quote, u.Q_(D(50) * D(str(expected.magnitude)) / D(str(numeric)), "EUR/" + volume_unit))
                        for amount_unit in ("kg", "mol"):
                            converted = u.convert(quantity, amount_unit, assay=assay, basis=kind)
                            self.assertRelativeClose(u.convert(converted, volume_unit, assay=assay, basis=kind), quantity)

    def test_molar_only_and_volume_only_gas_values_do_not_require_unrelated_metadata(self):
        dry_ref = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=None)
        wet_ref = replace(dry_ref, water_mole_fraction=D(".1"))
        dry_state, wet_state = u.GasMoistureState(dry_ref), u.GasMoistureState(wet_ref)
        molar = u.HeatContentAssay(gcv=u.Q_(D(".9"), "MJ/mol"), ncv=u.Q_(D(".8"), "MJ/mol"), basis="DRY_GAS", dry_gas_state=dry_state, wet_gas_state=wet_state)
        self.assertEqual(molar.calorific_value("gcv", basis="WET_GAS").magnitude, D(".81"))
        self.assertEqual(u.convert(u.Q_(D(2), "mol"), "MJ", assay=molar).magnitude, D("1.6"))
        with self.assertRaisesRegex(ValueError, "dry_gas_molar_mass"):
            molar.energy_density("MJ/kg")
        with self.assertRaisesRegex(ValueError, "compressibility"):
            molar.energy_density("MJ/m^3")
        volumetric = replace(molar, gcv=u.Q_(D(40), "MJ/Sm3"), ncv=u.Q_(D(36), "MJ/Sm3"))
        self.assertEqual(volumetric.energy_density("MJ/m^3").magnitude, 36)
        self.assertEqual(u.convert(u.Q_(D(2), "m^3"), "MJ", assay=volumetric).magnitude, 72)
        with self.assertRaisesRegex(ValueError, "compressibility"):
            volumetric.calorific_value("ncv", basis="WET_GAS")

    def test_gas_specific_energy_master_conversion_and_tag_validation(self):
        assay, _ = gas_spec()
        original = u.Q_(D(40), "MJ/Sm3")
        for numeric in (40., D(40)):
            incoming = u.Q_(numeric, "MJ/Sm3")
            wet = u.convert(incoming, "MJ/m^3", assay=assay, basis="gcv", analytical_basis="WET_GAS")
            self.assertAlmostEqual(float(wet.magnitude), 36., places=11)
            self.assertIsInstance(wet.magnitude, type(numeric))
        wet_assay = assay.for_basis("WET_GAS")
        with self.assertRaisesRegex(ValueError, "contradicts"):
            wet_assay.energy_density("MJ/Sm3")
        with self.assertRaisesRegex(ValueError, "contradicts"):
            replace(wet_assay, gcv=original, ncv=u.Q_(D(36), "MJ/Sm3"))
        with self.assertRaisesRegex(ValueError, "separate"):
            assay.calorific_value("gcv", basis="DB")
        with self.assertRaisesRegex(ValueError, "separate"):
            u.HeatContentAssay(gcv=u.Q_(D(20), "MJ/kg"), moisture_pct=10, dry_gas_state=assay.dry_gas_state).calorific_value("gcv", basis="DRY_GAS")

    def test_wet_gas_hhv_lhv_estimation_excludes_incoming_spectator_water(self):
        assay, _ = gas_spec()
        wet = assay.for_basis("WET_GAS")
        known_hydrogen_pct = D(10)
        gross = replace(wet, gcv=u.Q_(D(50), "MJ/kg"), ncv=None, hydrogen_pct=known_hydrogen_pct, method="latent_water")
        expected = D(50) - D("2.442") * D("8.936") * D(".1")
        self.assertEqual(gross.calorific_value("ncv").magnitude, expected)

    def test_missing_gas_states_unknown_water_and_inconsistent_composition_are_rejected(self):
        ref = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1))
        with self.assertRaisesRegex(ValueError, "water_mole_fraction"):
            u.GasMoistureState(replace(ref, water_mole_fraction=None))
        with self.assertRaisesRegex(ValueError, "explicit"):
            u.HeatContentAssay(gcv=u.Q_(20, "MJ/kg"), basis="WET_GAS")
        with self.assertRaisesRegex(ValueError, "zero water"):
            u.HeatContentAssay(gcv=u.Q_(20, "MJ/kg"), basis="DRY_GAS", dry_gas_state=u.GasMoistureState(replace(ref, water_mole_fraction=D(".1"))))
        dry, _ = gas_spec()
        with self.assertRaisesRegex(ValueError, "same dry-gas"):
            replace(dry, wet_gas_state=replace(dry.wet_gas_state, dry_gas_molar_mass=u.Q_(D(40), "g/mol")))
        with self.assertRaisesRegex(ValueError, "mole fractions"):
            replace(dry, moisture_pct=10)
        with self.assertRaisesRegex(ValueError, "conflicts"):
            u.GasMoistureState.from_composition(u.GasComposition({"CH4": D(".9"), "H2O": D(".1")}), ref)

    def test_original_seven_gas_reference_routes_survive_explicit_moisture_profiles(self):
        composition = u.GasComposition({"CH4": D(".9"), "N2": D(".1")})
        gas = composition.to_gas_assay()
        reference = u.GasReferenceConditions.from_standard("ferc", compressibility_factor=D(1))
        dry_state = u.GasMoistureState.from_composition(composition, reference)
        wet_state = u.GasMoistureState.from_composition(composition.with_water(D(".05")), replace(reference, water_mole_fraction=D(".05")))
        explicit = u.heat_content_from_assay(gas, basis="DRY_GAS", dry_gas_state=dry_state, wet_gas_state=wet_state)
        for projected in (explicit, explicit.for_basis("WET_GAS").for_basis("DRY_GAS")):
            for source in ("Nm3", "Sm3", "scf", "scf_tx", "scf_60F_14696psi_dry", "Sm3_20C", "m3_Groningen"):
                for kind in ("gcv", "ncv"):
                    with self.subTest(source=source, kind=kind):
                        quantity = u.Q_(D(1000), source)
                        actual = u.convert(quantity, "MWh", assay=projected, basis=kind)
                        expected = u.convert_gas_energy(quantity, "MWh", assay=gas, basis=kind)
                        self.assertRelativeClose(actual, expected)
                        self.assertRelativeClose(u.convert(actual, source, assay=projected, basis=kind), quantity)


if __name__ == "__main__":
    unittest.main()
