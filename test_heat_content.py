"""Shared calorific data against the real registry and commodity assays."""

from decimal import Decimal
import inspect
import unittest

import pint

import core.units as package
from core.units import agricultural_specs as agricultural
from core.units import conversions as c
from core.units import heat_content as h
from core.units.registry import ureg
from test_registry_compatibility import assay_presets
from test_agricultural_integration import agricultural_presets


def Q(value, unit):
    return ureg.Quantity(value, unit)


class HeatContentTests(unittest.TestCase):
    def assertQuantityClose(self, actual, expected):
        self.assertAlmostEqual(float(actual.to(expected.units).magnitude), float(expected.magnitude), places=8)

    def test_exports(self):
        self.assertEqual(len(h.__all__), len(set(h.__all__)))
        for name in h.__all__:
            self.assertIs(getattr(h, name), getattr(package, name))

    def test_measured_pair_and_explicit_estimates(self):
        gross = Q(Decimal("20"), "MJ/kg")
        net = h.gcv_to_ncv(gross, moisture_pct=10, hydrogen_pct=5)
        self.assertEqual(net.magnitude, Decimal("18.6647144"))
        self.assertEqual(h.ncv_to_gcv(net, moisture_pct=10, hydrogen_pct=5).magnitude, gross.magnitude)
        # Measured pairs take priority over either estimation mechanism.
        assay = h.HeatContentAssay(gcv=gross, ncv=Q(19, "MJ/kg"),
                                  moisture_pct=10, hydrogen_pct=5, ncv_to_gcv_ratio=.8)
        self.assertEqual(assay.calorific_value("LHV").magnitude, 19)
        self.assertEqual(assay.calorific_value("HHV").magnitude, 20)
        self.assertEqual(h.gcv_to_ncv(gross, ncv_to_gcv_ratio=Decimal("0.9")).magnitude, 18)
        with self.assertRaisesRegex(ValueError, "hydrogen"):
            h.gcv_to_ncv(gross)
        with self.assertRaisesRegex(ValueError, "Moisture"):
            h.gcv_to_ncv(gross, hydrogen_pct=5)
        for magnitude in (20., Decimal(20)):
            gross = Q(magnitude, "MJ/kg")
            net = h.gcv_to_ncv(gross, moisture_pct=10, hydrogen_pct=5)
            restored = h.ncv_to_gcv(net, moisture_pct=10, hydrogen_pct=5)
            dry = h.as_received_to_dry(net, 10, kind="ncv")
            wet = h.dry_to_as_received(dry, 10, kind="ncv")
            for quantity in (net, restored, dry, wet):
                self.assertIsInstance(quantity.magnitude, type(magnitude))
            self.assertQuantityClose(restored, gross)
            self.assertQuantityClose(wet, net)

    def test_analytical_bases_with_water_loss_and_ash(self):
        assay = h.HeatContentAssay(gcv=Q(20, "MJ/kg"), ncv=Q(18, "MJ/kg"),
                                  moisture_pct=20, air_dried_moisture_pct=10, ash_db_pct=10)
        self.assertEqual(assay.calorific_value("gcv", basis="DB").magnitude, Decimal(25))
        self.assertEqual(assay.calorific_value("ncv", basis="DB").magnitude, Decimal("23.1105"))
        self.assertEqual(assay.calorific_value("gcv", basis="AD").magnitude, Decimal("22.5"))
        self.assertEqual(assay.calorific_value("ncv", basis="AD").magnitude, Decimal("20.55525"))
        self.assertQuantityClose(assay.calorific_value("gcv", basis="DAF"), Q(20 / .72, "MJ/kg"))
        for kind in ("ncv", "gcv"):
            original = Q(Decimal("18"), "MJ/kg")
            dry = h.as_received_to_dry(original, 20, kind=kind)
            restored = h.dry_to_as_received(dry, 20, kind=kind)
            self.assertEqual(restored.magnitude, original.magnitude)
        daf = h.HeatContentAssay(ncv=Q(30, "MJ/kg"), basis="DAF", ash_db_pct=10,
                                air_dried_moisture_pct=10)
        self.assertEqual(daf.calorific_value("ncv", basis="DB").magnitude, Decimal(27))
        self.assertEqual(daf.calorific_value("ncv", basis="AD").magnitude, Decimal("24.0558"))
        with self.assertRaisesRegex(ValueError, "ash_db_pct"):
            h.HeatContentAssay(gcv=Q(20, "MJ/kg"), moisture_pct=0).calorific_value("gcv", basis="DAF")
        # A dry assay can carry AR/AD moisture metadata for all reverse paths.
        for source_basis in ("AR", "AD", "DB", "DAF"):
            for kind in ("gcv", "ncv"):
                source = assay.calorific_value(kind, basis=source_basis)
                other = h.HeatContentAssay(**{kind: source}, basis=source_basis,
                                          as_received_moisture_pct=20, air_dried_moisture_pct=10, ash_db_pct=10)
                for target_basis in ("AR", "AD", "DB", "DAF"):
                    with self.subTest(source=source_basis, target=target_basis, kind=kind):
                        self.assertQuantityClose(other.calorific_value(kind, basis=target_basis), assay.calorific_value(kind, basis=target_basis))

    def test_mass_volume_molar_energy_densities_and_contexts(self):
        assay = h.HeatContentAssay(gcv=Q(50, "MJ/kg"), ncv=Q(45, "MJ/kg"),
                                  density=Q(800, "kg/m^3"), molar_mass=Q(16, "g/mol"),
                                  standard_density=Q(Decimal("0.8"), "kg/Nm3"))
        self.assertEqual(assay.energy_density("MJ/m^3").magnitude, 36000)
        self.assertEqual(assay.energy_density("MJ/Nm3").magnitude, Decimal("36.0"))
        self.assertEqual(assay.energy_density("kJ/mol").magnitude, 720)
        self.assertQuantityClose(assay.energy_density("kWh/kg"), Q(12.5, "kWh/kg"))
        pairs = (("kg", "MWh"), ("m^3", "GJ"), ("Nm3", "MMBtu_it"),
                 ("mol", "kWh"), ("kg", "mol"), ("m^3", "Nm3"),
                 ("kg/hour", "MW"), ("Nm3/day", "kg/hour"),
                 ("USD/metric_ton", "USD/MMBtu_it"), ("USD/Nm3", "USD/kWh"))
        for source, destination in pairs:
            for magnitude in (1.25, Decimal("1.25")):
                with self.subTest(source=source, destination=destination, magnitude=magnitude):
                    original = Q(magnitude, source)
                    converted = c.convert(original, destination, assay=assay)
                    restored = c.convert(converted, source, assay=assay)
                    self.assertQuantityClose(restored, original)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(converted.magnitude, Decimal)
        with ureg.context(h.create_heat_content_context(assay, basis="HHV")):
            self.assertEqual(Q(Decimal(1), "kg").to("MJ").magnitude, 50)
        for source in (Q(36, "MJ/Nm3"), Q(720, "kJ/mol"), Q(36000, "MJ/m^3")):
            adapted = h.HeatContentAssay(ncv=source, density=assay.density,
                                        standard_density=assay.standard_density, molar_mass=assay.molar_mass)
            self.assertEqual(adapted.calorific_value().magnitude, Decimal(45))

    def test_calorific_density_dispatch_preserves_input_value(self):
        assay = h.HeatContentAssay(ncv=Q(20, "MJ/kg"), moisture_pct=20, density=Q(650, "kg/m^3"))
        # Caller quantity is 18, while assay benchmark is 20.
        result = h.convert_heat_content(Q(Decimal(18), "MJ/kg"), "MJ/kg", assay=assay, analytical_basis="DB")
        self.assertEqual(result.magnitude, Decimal("23.1105"))
        self.assertEqual(c.convert(Q(Decimal(18), "MJ/kg"), "MJ/m^3", assay=assay).magnitude, 11700)
        with self.assertRaisesRegex(ValueError, "new density"):
            h.convert_heat_content(Q(18, "MJ/kg"), "MJ/m^3", assay=assay, analytical_basis="DB")

    def test_volume_only_gas_and_reference_state_separation(self):
        original = c.GasAssay(gcv_btu_per_scf=1000)
        assay = h.heat_content_from_assay(original)
        self.assertQuantityClose(h.convert_heat_content(Q(1000000, "scf"), "MMBtu_it", assay=assay, basis="hhv"), Q(1000, "MMBtu_it"))
        with self.assertRaises(pint.DimensionalityError):
            h.convert_heat_content(Q(1, "scf"), "kg", assay=assay)
        with self.assertRaises(pint.DimensionalityError):
            h.convert_heat_content(Q(1, "scf"), "Nm3", assay=assay)
        molar = h.HeatContentAssay(ncv=Q(800, "kJ/mol"))
        self.assertQuantityClose(h.convert_heat_content(Q(1, "mol"), "kJ", assay=molar), Q(800, "kJ"))
        for unit in ("scf", "scf_tx", "scf_60F_14696psi_dry", "Sm3", "Nm3", "Sm3_20C", "m3_Groningen"):
            with self.subTest(reference=unit):
                spec = h.HeatContentAssay(ncv=Q(Decimal(36), ureg.MJ / ureg.Unit(unit)))
                self.assertQuantityClose(h.convert_heat_content(Q(Decimal(2), unit), "MJ", assay=spec), Q(72, "MJ"))
        gas = c.GasAssay(gcv_btu_per_scf=1000, relative_density=.6)
        mass_heat = h.convert_heat_content(Q(1, "metric_ton"), "MMBtu_it", assay=gas, basis="hhv")
        # GasAssay stores this cargo benchmark rounded to four decimals.
        self.assertEqual(round(float(mass_heat.magnitude), 4), gas.mmbtu_per_mt_lng)

    def test_efficiency_heat_rate_and_power(self):
        rate = h.heat_rate_from_efficiency(Q(40, "percent"))
        self.assertAlmostEqual(float(rate.magnitude), 8530.35408281985, places=7)
        self.assertAlmostEqual(float(h.efficiency_from_heat_rate(rate)), .4, places=14)
        for magnitude in (100., Decimal(100)):
            energy = Q(magnitude, "MWh")
            output = h.thermal_to_output(energy, efficiency=Decimal("0.4"))
            self.assertEqual(output.magnitude, 40)
            restored = h.output_to_thermal(output, heat_rate=rate)
            self.assertQuantityClose(restored, energy)
            if isinstance(magnitude, Decimal):
                self.assertIsInstance(restored.magnitude, Decimal)
        self.assertEqual(h.thermal_to_output(Q(100, "MW"), efficiency=.4, target_unit="MW").magnitude, 40)
        self.assertEqual(h.thermal_to_output(Q(100, "MWh_th"), efficiency=.4, target_unit="MWh_e").magnitude, 40)
        self.assertEqual(h.output_to_thermal(Q(40, "MW_e"), efficiency=.4, target_unit="MW_th").magnitude, 100)
        for rate_unit in ("heat_rate_us", "heat_rate_macro", "heat_rate_metric", "heat_rate_gas_turbine"):
            self.assertAlmostEqual(float(h.efficiency_from_heat_rate(rate.to(rate_unit))), .4, places=14)
        for kwargs in ({}, {"efficiency": .4, "heat_rate": rate}, {"efficiency": 40}, {"efficiency": 0}, {"efficiency": float("nan")}):
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                h.thermal_to_output(Q(1, "MWh"), **kwargs)
        with self.assertRaises(ValueError):
            h.efficiency_from_heat_rate(Q(1000, "btu_it/kWh"))

    def test_all_fuel_presets_match_existing_energy_semantics(self):
        fuel_types = (c.GasAssay, c.CoalAssay, c.CleanFuelAssay, c.BiofuelAssay,
                      c.NGLAssay, c.BiomassAssay, c.PetcokeAssay, c.HydrogenAssay)
        checked = 0
        for name, original in assay_presets():
            if not isinstance(original, fuel_types):
                continue
            for basis in ("hhv", "lhv"):
                for magnitude in (1.25, Decimal("1.25")):
                    with self.subTest(preset=name, basis=basis, magnitude=magnitude):
                        spec = h.heat_content_from_assay(original)
                        source = Q(magnitude, "scf" if isinstance(original, c.GasAssay) else "kg")
                        expected = c.convert(source, "MJ", assay=original, basis=basis)
                        actual = h.convert_heat_content(source, "MJ", assay=original, basis=basis)
                        self.assertQuantityClose(actual, expected)
                        restored = c.convert(actual, source.units, assay=spec, basis=basis)
                        self.assertQuantityClose(restored, source)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(actual.magnitude, Decimal)
                        checked += 1
        self.assertGreater(checked, 150)

    def test_explicit_heat_data_for_every_other_assay(self):
        fuels = (c.GasAssay, c.CoalAssay, c.CleanFuelAssay, c.BiofuelAssay,
                 c.NGLAssay, c.BiomassAssay, c.PetcokeAssay, c.HydrogenAssay)
        seen = set()
        refinery = c.RefineryCrudeAssay(c.CrudeAssay.brent(),
                                      (c.DistillationCut("whole", 0, 800, 100, 100, .84),))
        examples = list(assay_presets()) + list(agricultural_presets()) + [
            ("RefineryCrudeAssay", refinery), ("BulkOreAssay", c.BulkOreAssay()),
        ]
        for name, original in examples:
            if isinstance(original, fuels) or type(original) in seen:
                continue
            seen.add(type(original))
            with self.subTest(preset=name):
                with self.assertRaisesRegex(ValueError, "Supply at least one"):
                    h.heat_content_from_assay(original)
                spec = h.HeatContentAssay.from_assay(original, ncv=Q(Decimal(18), "MJ/kg"))
                self.assertEqual(c.convert(Q(Decimal(2), "kg"), "MJ", assay=original, heat_content=spec).magnitude, 36)
                self.assertEqual(h.convert_heat_content(Q(Decimal(2), "kg"), "MJ", assay=original, ncv=Q(18, "MJ/kg")).magnitude, 36)
        expected = {value for module in (c, agricultural) for value in vars(module).values()
                    if inspect.isclass(value) and value.__module__ == module.__name__
                    and value.__name__.endswith("Assay") and value not in fuels}
        self.assertEqual(seen, expected)

    def test_fx_quote_router_with_heat_layer(self):
        spec = h.HeatContentAssay(ncv=Q(Decimal(20), "MJ/kg"))
        converted = c.convert(Q(Decimal(100), "USD/metric_ton"), "EUR/GJ", heat_content=spec,
                              rates={"USD/EUR": Decimal("0.9")})
        self.assertEqual(converted.magnitude, Decimal("4.5"))

    def test_invalid_metadata(self):
        invalid = ({"gcv": Q(0, "MJ/kg")}, {"gcv": Q(float("nan"), "MJ/kg")},
                   {"gcv": Q(10, "MJ/kg"), "ncv": Q(11, "MJ/kg")},
                   {"ncv": Q(10, "MJ/kg"), "moisture_pct": 100},
                   {"ncv": Q(10, "MJ/kg"), "hydrogen_pct": 101},
                   {"ncv": Q(10, "MJ/kg"), "basis": "DB", "moisture_pct": 10},
                   {"ncv": Q(10, "MJ/kg"), "density": Q(1, "kg")},
                   {"ncv": Q(10, "MJ/kg"), "molar_mass": Q(-1, "g/mol")},
                   {"ncv": Q(10, "MJ/kg"), "reference_pressure": Q(1, "psig")},
                   {"ncv": Q(10, "MJ/kg"), "reference_temperature": Q(0, "kelvin")},
                   {"ncv": Q(10, "MJ/kg"), "ncv_to_gcv_ratio": 1.1},
                   {"ncv": Q(10, "MJ/kg"), "basis": "unknown"},
                   {"ncv": Q(10, "MJ/kg"), "density": Q(1, "kg/m^3"), "standard_density": Q(2, "kg/m^3")})
        for kwargs in invalid:
            with self.subTest(kwargs=kwargs), self.assertRaises(ValueError):
                h.HeatContentAssay(**kwargs)


if __name__ == "__main__":
    unittest.main()
