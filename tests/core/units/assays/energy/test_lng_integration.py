"""Gas composition, thermodynamic/cargo accounting and thermal-route regressions."""

from dataclasses import replace
from decimal import Decimal
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import pint

import core.units as u
from core.units import conversions as c
from core.units import thermodynamics
from core.units.converters import gas as model


D = Decimal


class GasCompositionTests(unittest.TestCase):
    def setUp(self):
        self.composition = u.GasComposition({"CH4": D(".9"), "C2H6": D(".05"), "N2": D(".05")})
        self.reference = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1))
        self.backend = u.builtin_gas_backend("ideal", temperature_range=(u.Q_(100, "K"), u.Q_(400, "K")),
                                            pressure_range=(u.Q_(10000, "Pa"), u.Q_(5000000, "Pa")))

    def test_normalization_is_explicit_and_original_total_is_retained(self):
        for values in ({"CH4": 90, "N2": 10}, {"CH4": D(".999999")}, {"CH4": -1}, {"CH4": True},
                       {"CH4": float("nan")}, {"UNKNOWN": 0, "CH4": 1}, {"CH4": 0}):
            with self.subTest(values=values), self.assertRaises((TypeError, ValueError)):
                u.GasComposition(values)
        composition = u.GasComposition({"CH4": 90, "N2": 10}, normalize=True)
        self.assertEqual(composition.source_total, 100)
        self.assertEqual(composition.mole_fractions, {"CH4": D(".9"), "N2": D(".1")})
        with self.assertRaises(TypeError):
            composition.mole_fractions["CH4"] = 1
        with self.assertRaises(ValueError):
            u.GasAssay(gcv_btu_per_scf=1000, molar_composition={"CH4": 99, "N2": 1})
        assay = u.GasAssay(gcv_btu_per_scf=1000, molar_composition={"CH4": 99, "N2": 1}, normalize_composition=True)
        self.assertAlmostEqual(assay.molar_composition["CH4"], .99)
        with self.assertRaises(ValueError):
            u.calculate_gas_calorific_value({"CH4": D(".9999")})
        self.assertGreater(u.calculate_gas_calorific_value({"CH4": 99, "N2": 1}, normalize_composition=True).magnitude, 0)

    def test_composition_mass_heating_values_and_wobbe_are_recalculated(self):
        methane = u.GasComposition({"CH4": 1})
        self.assertEqual(methane.molar_mass.magnitude, D("16.0425"))
        props = u.gas_properties(methane, self.reference, backend=self.backend)
        expected_density = D(101325) * D(".0160425") / (D(str(c.UNIVERSAL_GAS_CONSTANT_J_MOL_K)) * D("288.15"))
        self.assertAlmostEqual(float(props.density.magnitude), float(expected_density), places=12)
        self.assertGreater(props.hhv, props.lhv)
        self.assertEqual(props.gross_wobbe, props.hhv / props.relative_density.sqrt())
        self.assertEqual(props.net_wobbe, props.lhv / props.relative_density.sqrt())
        self.assertEqual(props.combustion_temperature.magnitude, D("288.15"))
        assay = methane.to_gas_assay()
        self.assertEqual(u.GasAssay.from_composition({"CH4": 1}).gcv_btu_per_scf, assay.gcv_btu_per_scf)
        self.assertAlmostEqual(assay.molar_mass_g_per_mol, 16.0425)
        self.assertLess(assay.net_wobbe_index_mj_per_sm3, assay.wobbe_index_mj_per_sm3)
        mass_energy = u.convert(u.Q_(D(1000), "kg"), "MJ", assay=assay, basis="hhv")
        composition_mass_energy = methane.calorific_values().hhv.magnitude * 1000
        self.assertAlmostEqual(float(mass_energy.magnitude), float(composition_mass_energy), places=7)
        with self.assertRaises(ValueError):
            u.calculate_gas_calorific_value({"CH4": 1}, reference_conditions=replace(self.reference, compressibility_factor=None))
        unmeasured = u.GasAssay(gcv_btu_per_scf=1000)
        with self.assertRaisesRegex(ValueError, "no default gas"):
            unmeasured.z_factor(u.Q_(1, "MPa"), u.Q_(300, "K"), method="dak")

    def test_industrial_gases_are_valid_noncombustible_compositions(self):
        for component in ("O2", "N2", "Ar", "He", "CO2"):
            with self.subTest(component=component):
                composition = u.GasComposition({component: 1})
                properties = u.gas_properties(composition, self.reference, backend=self.backend)
                self.assertGreater(properties.density.magnitude, 0)
                self.assertEqual(properties.hhv.magnitude, 0)
                self.assertEqual(properties.gross_wobbe.magnitude, 0)
                with self.assertRaises(ValueError):
                    composition.to_gas_assay()

    def test_hydrogen_and_biomethane_blending_uses_moles_and_explicit_quality_limits(self):
        methane, hydrogen = u.GasComposition({"CH4": 1}), u.GasComposition({"H2": 1})
        blend = methane.blend(hydrogen, self_mole_fraction=D(".8"))
        self.assertEqual(blend.mole_fractions["H2"], D(".2"))
        pure_properties = u.gas_properties(methane, self.reference, backend=self.backend)
        blend_properties = u.gas_properties(blend, self.reference, backend=self.backend)
        self.assertLess(blend_properties.hhv, pure_properties.hhv)
        self.assertLess(blend_properties.molar_mass, pure_properties.molar_mass)
        quality = blend.quality(minimum_methane=D(".95"), maximum_co2=D(".02"), maximum_h2s=D(".00001"), maximum_water=D(".0001"))
        self.assertFalse(quality["methane"])
        self.assertTrue(quality["co2"])
        dirty = u.GasComposition({"CH4": D(".96"), "CO2": D(".03"), "H2S": D(".01")})
        self.assertFalse(dirty.quality(minimum_methane=D(".95"), maximum_co2=D(".02"), maximum_h2s=D(".00001"), maximum_water=D(".0001"))["h2s"])
        self.assertGreater(dirty.molar_mass.magnitude, methane.molar_mass.magnitude)

    def test_missing_heavy_component_calorific_data_requires_measurements(self):
        composition = u.GasComposition({"CH4": D(".99"), "C7H16": D(".01")})
        self.assertGreater(composition.molar_mass.magnitude, 16)
        with self.assertRaisesRegex(ValueError, "No bundled calorific data"):
            composition.calorific_values()
        values = u.GasCalorificValues(u.Q_(55, "MJ/kg"), u.Q_(50, "MJ/kg"), u.Q_(15, "degC"), "measured example")
        self.assertGreater(u.gas_properties(composition, self.reference, backend=self.backend, calorific_values=values).hhv.magnitude, 0)

    def test_real_gas_density_wobbe_air_and_combustion_states_are_explicit(self):
        calls = []
        def z(composition, state):
            calls.append(state)
            return D(".95")
        backend = replace(self.backend, method_id="specified alternative", ideal=False, z_calculator=z)
        reference = replace(self.reference, compressibility_factor=None)
        with self.assertRaisesRegex(ValueError, "air compressibility"):
            u.gas_properties(self.composition, reference, backend=backend)
        calls.clear()
        properties = u.gas_properties(self.composition, reference, backend=backend, air_compressibility=D(".99"))
        self.assertEqual(len(calls), 1)
        self.assertEqual(properties.relative_density, self.composition.ideal_relative_density * D(".99") / D(".95"))
        self.assertEqual(properties.metering_reference.compressibility_factor, D(".95"))
        self.assertEqual(properties.metering_reference.water_mole_fraction, 0)
        with self.assertRaisesRegex(ValueError, "metering Z conflicts"):
            u.gas_density(self.composition, self.reference, backend=backend)
        with self.assertRaisesRegex(ValueError, "Combustion reference"):
            u.gas_properties(self.composition, replace(self.reference, combustion_reference_temperature=u.Q_(25, "degC")), backend=self.backend)
        with self.assertRaisesRegex(ValueError, "outside backend"):
            u.gas_density(self.composition, replace(self.reference, temperature=u.Q_(450, "K")), backend=self.backend)

    def test_explicit_eos_selection_preserves_model_identity_and_sour_fractions(self):
        with self.assertRaises(ValueError):
            u.builtin_gas_backend("auto", temperature_range=self.backend.temperature_range, pressure_range=self.backend.pressure_range)
        sour = u.GasComposition({"CH4": D(".93"), "CO2": D(".05"), "H2S": D(".02")})
        backend = u.builtin_gas_backend("dak", temperature_range=self.backend.temperature_range, pressure_range=self.backend.pressure_range)
        reference = replace(self.reference, absolute_pressure=u.Q_(1, "MPa"), temperature=u.Q_(300, "K"), compressibility_factor=None)
        z = backend.z_factor(sour, reference)
        expected = u.calculate_z_factor(reference.absolute_pressure, reference.temperature, float(sour.ideal_relative_density), method="dak",
                                         molar_composition=sour.mole_fractions)
        self.assertAlmostEqual(float(z), expected)
        with self.assertRaisesRegex(ValueError, "conflicts"):
            u.calculate_z_factor(reference.absolute_pressure, reference.temperature, method="dak", molar_composition=sour.mole_fractions, co2_mol_frac=0)
        with self.assertRaisesRegex(ValueError, "hydrocarbon approximation"):
            backend.z_factor(u.GasComposition({"H2": 1}), reference)
        for method, r_value in (("aga8", "8.31451"), ("gerg2008", "8.314472"), ("heos", str(c.UNIVERSAL_GAS_CONSTANT_J_MOL_K))):
            selected = u.builtin_gas_backend(method, temperature_range=self.backend.temperature_range, pressure_range=self.backend.pressure_range)
            self.assertEqual(selected.gas_constant.magnitude, D(r_value))
            with patch.object(thermodynamics, "calculate_z_factor", side_effect=ImportError("backend unavailable")):
                # The factory's closure binds the function at construction.
                selected = u.builtin_gas_backend(method, temperature_range=self.backend.temperature_range, pressure_range=self.backend.pressure_range)
                with self.assertRaisesRegex(ImportError, "backend unavailable"):
                    selected.z_factor(self.composition, reference)
        with patch.object(thermodynamics, "calculate_z_factor", return_value=.95), patch.dict("sys.modules", {
            "CoolProp": SimpleNamespace(), "CoolProp.CoolProp": SimpleNamespace(PhaseSI=lambda *args: "liquid")
        }):
            selected = u.builtin_gas_backend("heos", temperature_range=self.backend.temperature_range, pressure_range=self.backend.pressure_range)
            with self.assertRaisesRegex(ValueError, "not a single gas"):
                selected.z_factor(self.composition, reference)

    def test_composition_energy_matches_volume_and_mass_routes_at_matching_states(self):
        for unit, key in (("Nm3", "din0"), ("Sm3", "iso15"), ("scf", "ferc")):
            reference = u.GasReferenceConditions.from_standard(key, compressibility_factor=D(1))
            volume = u.Q_(D(1000), unit)
            actual = u.composition_gas_energy(volume, self.composition, reference, backend=self.backend, basis="HHV")
            geometric = model._geometric_volume(volume, reference)
            mass = geometric * u.gas_density(self.composition, reference, backend=self.backend).magnitude
            expected = mass * self.composition.calorific_values().hhv.magnitude
            self.assertAlmostEqual(float(actual.magnitude), float(expected), places=9)
        with self.assertRaisesRegex(ValueError, "contradicts"):
            u.composition_gas_energy(u.Q_(1, "Nm3"), self.composition, self.reference, backend=self.backend, basis="HHV")


class WaterAndDewPointTests(unittest.TestCase):
    def setUp(self):
        self.composition = u.GasComposition({"CH4": D(".95"), "N2": D(".05")})
        self.reference = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1))
        self.water = u.ideal_water_saturation_model()
        self.backend = u.builtin_gas_backend("ideal", temperature_range=(u.Q_(270, "K"), u.Q_(330, "K")),
                                            pressure_range=(u.Q_(10000, "Pa"), u.Q_(200000, "Pa")))

    def test_dry_wet_saturated_gas_and_dew_point_round_trip(self):
        saturated = self.water.saturate(self.composition, self.reference)
        self.assertIsNone(saturated.reference.compressibility_factor)
        self.assertGreater(saturated.composition.water_mole_fraction, 0)
        self.assertAlmostEqual(float(saturated.composition.dry().mole_fractions["CH4"]), .95)
        dew = self.water.dew_point(saturated.composition, saturated.reference.absolute_pressure)
        self.assertAlmostEqual(float(dew.temperature.magnitude), 288.15, places=9)
        drier = saturated.composition.with_water(saturated.composition.water_mole_fraction / 2)
        self.assertLess(self.water.dew_point(drier, self.reference.absolute_pressure).temperature, dew.temperature)
        content = u.gas_water_content(saturated.composition, saturated.reference, backend=self.backend)
        self.assertGreater(content.magnitude, 0)
        with self.assertRaisesRegex(ValueError, "water fraction conflicts"):
            u.gas_density(saturated.composition, self.reference, backend=self.backend)

    def test_water_equilibrium_does_not_extrapolate_to_ice_or_pipeline_pressure(self):
        for reference in (replace(self.reference, temperature=u.Q_(260, "K")), replace(self.reference, absolute_pressure=u.Q_(1, "MPa"))):
            with self.subTest(reference=reference), self.assertRaisesRegex(ValueError, "applicability"):
                self.water.saturate(self.composition, reference)
        with self.assertRaises(ValueError):
            self.water.dew_point(self.composition, self.reference.absolute_pressure)
        with self.assertRaises(ValueError):
            self.water.dew_point(self.composition.with_water(D("1e-6")), self.reference.absolute_pressure)

    def test_phase_backends_distinguish_water_and_hydrocarbon_dew_points(self):
        calls = []
        def dew(composition, pressure, kind):
            calls.append(kind)
            return u.Q_(280 if kind == "water" else 290, "K")
        backend = replace(self.backend, dew_point_calculator=dew, method_id="specified equilibrium")
        wet = self.composition.with_water(D(".01"))
        self.assertEqual(u.gas_dew_point(wet, u.Q_(101325, "Pa"), kind="water", backend=backend).temperature.magnitude, 280)
        self.assertEqual(u.gas_dew_point(self.composition, u.Q_(101325, "Pa"), kind="hydrocarbon", backend=backend).temperature.magnitude, 290)
        self.assertEqual(calls, ["water", "hydrocarbon"])
        with self.assertRaisesRegex(ValueError, "does not provide"):
            u.gas_dew_point(self.composition, u.Q_(101325, "Pa"), kind="hydrocarbon", backend=self.backend)
        invalid = replace(backend, dew_point_calculator=lambda composition, pressure, kind: u.Q_(500, "K"))
        with self.assertRaisesRegex(ValueError, "applicability"):
            u.gas_dew_point(self.composition, u.Q_(101325, "Pa"), kind="hydrocarbon", backend=invalid)


class LNGCargoTests(unittest.TestCase):
    def setUp(self):
        self.composition = u.GasComposition({"CH4": D(".95"), "N2": D(".05")})
        self.state = u.LiquidState(u.Q_(112, "K"), u.Q_(150000, "Pa"))
        self.cargo = u.LNGCargo("example cargo", self.composition, self.state, measured_density=u.Q_(D(450), "kg/meter**3"), density_source="example density certificate")
        self.reference = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1))
        self.backend = u.builtin_gas_backend("ideal", temperature_range=(u.Q_(90, "K"), u.Q_(350, "K")),
                                            pressure_range=(u.Q_(10000, "Pa"), u.Q_(5000000, "Pa")))

    def test_cargo_density_is_measured_or_composition_temperature_dependent(self):
        with self.assertRaisesRegex(ValueError, "exactly one"):
            u.LNGCargo("missing", self.composition, self.state)
        calls = []
        def density(composition, state):
            calls.append((composition, state))
            return u.Q_(D(450) + (D(112) - state.temperature.magnitude) + D(100) * composition.mole_fractions.get("N2", 0), "kg/meter**3")
        backend = replace(self.backend, liquid_density_calculator=density, method_id="synthetic composition/T density model")
        cargo = u.LNGCargo("calculated", self.composition, self.state, density_backend=backend)
        self.assertEqual(cargo.liquid_density.magnitude, 455)
        self.assertLess(replace(cargo, liquid_state=replace(self.state, temperature=u.Q_(115, "K"))).liquid_density.magnitude, 455)
        self.assertNotEqual(replace(cargo, composition=u.GasComposition({"CH4": 1})).liquid_density.magnitude, 455)
        self.assertTrue(calls)

    def test_liquid_gas_mass_and_energy_routes_agree_without_a_fixed_expansion_ratio(self):
        volume = u.Q_(D(1000), "meter**3")
        gas_volume = self.cargo.gas_equivalent(volume, self.reference, backend=self.backend, target_unit="Sm3")
        restored = self.cargo.liquid_equivalent(gas_volume, self.reference, backend=self.backend)
        self.assertAlmostEqual(float(restored.magnitude), 1000, places=9)
        gas_energy = u.composition_gas_energy(gas_volume, self.composition, self.reference, backend=self.backend, basis="HHV", target_unit="MMBtu")
        self.assertAlmostEqual(float(gas_energy.magnitude), float(self.cargo.energy(volume, basis="HHV").magnitude), places=8)
        self.assertNotAlmostEqual(float(self.cargo.expansion_ratio(self.reference, backend=self.backend)), 600, places=4)
        self.assertNotEqual(self.cargo.expansion_ratio(self.reference, backend=self.backend), replace(self.cargo, measured_density=u.Q_(460, "kg/meter**3")).expansion_ratio(self.reference, backend=self.backend))

    def test_loading_unloading_and_adjustments_have_explicit_directions_and_bases(self):
        opening, closing = u.LNGInventory(self.cargo, u.Q_(D(100), "meter**3")), u.LNGInventory(self.cargo, u.Q_(D(1000), "meter**3"))
        adjustment = u.LNGTransferAdjustment("vapor return", u.Q_(D(100), "MJ"), "HHV", "debit", "contract example", u.Q_(15, "degC"))
        result = u.calculate_lng_transfer(opening, closing, direction="loading", basis="HHV", adjustments=(adjustment,), target_unit="MJ")
        self.assertEqual(result.net_transferred_energy.magnitude, result.gross_transferred_energy.magnitude - 100)
        reverse = u.calculate_lng_transfer(closing, opening, direction="unloading", basis="HHV", adjustments=(), target_unit="MJ")
        self.assertEqual(reverse.gross_transferred_energy, result.gross_transferred_energy)
        self.assertEqual(result.opening_density.magnitude, 450)
        for invalid in (replace(adjustment, basis="LHV"), replace(adjustment, combustion_temperature=u.Q_(25, "degC"))):
            with self.subTest(invalid=invalid), self.assertRaises(ValueError):
                u.calculate_lng_transfer(opening, closing, direction="loading", basis="HHV", adjustments=(invalid,))
        with self.assertRaisesRegex(ValueError, "Duplicate"):
            u.calculate_lng_transfer(opening, closing, direction="loading", basis="HHV", adjustments=(adjustment, adjustment))
        with self.assertRaisesRegex(ValueError, "direction"):
            u.calculate_lng_transfer(opening, closing, direction="unloading", basis="HHV", adjustments=())

    def test_temperature_density_inventory_changes_do_not_create_energy(self):
        other_cargo = replace(self.cargo, liquid_state=replace(self.state, temperature=u.Q_(110, "K")), measured_density=u.Q_(D(500), "kg/meter**3"))
        opening = u.LNGInventory(self.cargo, u.Q_(D(1000), "meter**3"))
        closing = u.LNGInventory(other_cargo, u.Q_(D(900), "meter**3"))
        result = u.calculate_lng_transfer(opening, closing, direction="unloading", basis="HHV", adjustments=(), target_unit="MJ")
        self.assertEqual(result.net_transferred_energy.magnitude, 0)

    def test_wet_recipes_and_unidentified_density_or_energy_bases_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "water-free"):
            replace(self.cargo, composition=self.composition.with_water(D(".001")))
        with self.assertRaises(ValueError):
            replace(self.cargo, density_source=None)
        with self.assertRaises(ValueError):
            self.cargo.energy(u.Q_(D(1), "meter**3"), basis="unknown")
        with self.assertRaises(ValueError):
            replace(self.cargo, measured_density=u.Q_(float("nan"), "kg/meter**3"))
        gas = self.composition.to_gas_assay()
        for converter, quantity in ((gas.lng_liquid_to_energy, u.Q_(D(1), "meter**3")),
                                    (gas.energy_to_lng_liquid, u.Q_(D(1), "MMBtu")),
                                    (gas.lng_liquid_to_gas_volume, u.Q_(D(1), "meter**3"))):
            with self.subTest(converter=converter.__name__), self.assertRaisesRegex(ValueError, "cargo-specific"):
                converter(quantity)
        original = u.Q_(D(1000), "meter**3")
        energy = gas.lng_liquid_to_energy(original, liquid_density_kg_m3=D(450))
        restored = gas.energy_to_lng_liquid(energy, liquid_density_kg_m3=D(450))
        self.assertAlmostEqual(float(restored.magnitude), 1000, places=9)
        self.assertGreater(gas.lng_liquid_to_gas_volume(original, liquid_density_kg_m3=D(450)).magnitude, 0)


class GasThermalDispatcherRegressionTests(unittest.TestCase):
    def test_all_requested_amount_flow_and_price_routes_match_original_gas_context(self):
        gas = u.GasAssay(gcv_btu_per_scf=1000, relative_density=.6)
        for magnitude in (D(1000), 1000.):
            for source, target in (("Nm3", "MWh"), ("Sm3", "MMBtu"), ("scf", "MWh")):
                original = u.Q_(magnitude, source)
                for basis in ("gcv", "ncv"):
                    with self.subTest(source=source, target=target, kind=type(magnitude), basis=basis):
                        expected = c.convert_gas_energy(original, target, assay=gas, basis=basis)
                        actual = u.convert(original, target, assay=gas, basis=basis)
                        self.assertIsInstance(actual.magnitude, type(magnitude))
                        self.assertAlmostEqual(float(actual.magnitude), float(expected.magnitude), places=9)
                        self.assertAlmostEqual(float(u.convert(actual, source, assay=gas, basis=basis).magnitude), float(magnitude), places=8)
        amount = u.convert(u.Q_(D(1000), "Nm3"), "MWh", assay=gas)
        flow = u.convert(u.Q_(D(1000), "Nm3/day"), "MW", assay=gas)
        self.assertAlmostEqual(float(flow.magnitude), float(amount.magnitude / 24), places=10)
        price = u.Q_(D(50), "EUR/MWh")
        actual = u.convert(price, "EUR/Sm3", assay=gas)
        expected = c.convert_commodity_quote(price, "EUR/Sm3", assay=gas)
        self.assertAlmostEqual(float(actual.magnitude), float(expected.magnitude), places=12)
        self.assertAlmostEqual(float(actual.magnitude), 50 * gas.gcv_mj_per_sm3 / 3600, places=12)
        self.assertAlmostEqual(float(u.convert(actual, "EUR/MWh", assay=gas).magnitude), 50, places=10)

    def test_shared_heat_adapter_retains_all_seven_reference_families(self):
        gas = u.GasAssay(gcv_btu_per_scf=1030, ncv_btu_per_scf=930, relative_density=.62)
        heat = u.heat_content_from_assay(gas)
        for unit in ("Nm3", "Sm3", "scf", "scf_tx", "scf_60F_14696psi_dry", "Sm3_20C", "m3_Groningen"):
            for basis in ("hhv", "lhv"):
                with self.subTest(unit=unit, basis=basis):
                    original = u.Q_(D(1000), unit)
                    expected = c.convert_gas_energy(original, "MJ", assay=gas, basis=basis)
                    actual = u.convert(original, "MJ", assay=gas, heat_content=heat, basis=basis)
                    self.assertAlmostEqual(float(actual.magnitude), float(expected.magnitude), places=8)


if __name__ == "__main__":
    unittest.main()
