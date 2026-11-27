"""Regression coverage for commodity assay Pint contexts."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import conversions as c
from core.units.registry import ureg


class AssayContextTests(unittest.TestCase):
    def test_physical_bridges_and_round_trips(self):
        cases = (
            (c.BiofuelAssay.ethanol(), ("m^3", "kg", "MJ")),
            (c.NGLAssay.propane(), ("m^3", "kg", "MJ")),
            (c.OreAssay.iron_ore_62fe(), ("m^3", "kg")),
            (c.UraniumAssay.leu_4_5pct(), ("kg", "MWh")),
            (c.BiomassAssay.industrial_pellets_i2(), ("m^3", "kg", "MJ")),
            (c.PetrochemAssay.benzene(), ("m^3", "kg", "mol")),
            (c.FertilizerAssay.uan_32(), ("m^3", "kg")),
            (c.BitumenAssay.pen_60_70(), ("m^3", "kg")),
            (c.PetcokeAssay.usgc_fuel_grade_65s(), ("m^3", "kg", "MJ")),
            (c.ChemicalAssay.caustic_soda_50pct(), ("m^3", "kg", "mol")),
        )
        for assay, units in cases:
            for magnitude in (2.5, Decimal("2.5")):
                for source in units:
                    for target in units:
                        if source == target:
                            continue
                        with self.subTest(assay=assay.name, source=source, target=target, magnitude=magnitude):
                            quantity = ureg.Quantity(magnitude, source)
                            with self.assertRaises(pint.DimensionalityError):
                                quantity.to(target)
                            with ureg.context(assay.as_context()):
                                result = quantity.to(target)
                                self.assertGreater(result.magnitude, 0)
                                self.assertAlmostEqual(float(result.to(source).magnitude), 2.5)
                                if isinstance(magnitude, Decimal):
                                    self.assertIsInstance(result.magnitude, Decimal)

    def test_physical_reference_values(self):
        ethanol = c.BiofuelAssay.ethanol()
        self.assertAlmostEqual(float(ureg.Quantity(1, "m^3").to("kg", ethanol.as_context()).magnitude), ethanol.density_kg_per_m3)
        uranium = c.UraniumAssay.leu_4_5pct()
        thermal = ureg.Quantity(1, "kg").to("MWh", uranium.as_context())
        self.assertAlmostEqual(float(thermal.magnitude), 1080.0)
        electrical = c.create_uranium_context(uranium, burnup_gwd_per_mtu=50, thermal_efficiency=0.33)
        self.assertAlmostEqual(float(ureg.Quantity(1, "kg").to("MWh", electrical).magnitude), 396.0)

    def test_named_mass_only_contexts(self):
        for assay in (c.BatteryMaterialAssay.spodumene_sc6(), c.BullionAssay.gold_dore()):
            with self.subTest(assay=assay.name):
                context = assay.as_context(context_name="custom_assay")
                self.assertEqual(context.name, "custom_assay")
                self.assertFalse(context.funcs)
                quantity = ureg.Quantity(1, "kg")
                self.assertEqual(quantity.to("g", context).magnitude, 1000)

    def test_missing_optional_density(self):
        assays = (
            replace(c.OreAssay.iron_ore_62fe(), bulk_density_wmt_per_m3=None),
            replace(c.FertilizerAssay.uan_32(), density_kg_per_m3=None),
        )
        for assay in assays:
            with self.subTest(assay=assay.name):
                context = assay.as_context()
                self.assertFalse(context.funcs)
                with self.assertRaises(pint.DimensionalityError):
                    ureg.Quantity(1, "m^3").to("kg", context)

    def test_options_and_invalid_inputs(self):
        fuels = (c.BiofuelAssay.ethanol(), c.NGLAssay.propane(), c.BiomassAssay.industrial_pellets_i2(), c.PetcokeAssay.usgc_fuel_grade_65s())
        for assay in fuels:
            with self.subTest(assay=assay.name):
                with self.assertRaises(ValueError):
                    assay.as_context(basis="invalid")
                mass = ureg.Quantity(1, "kg")
                self.assertGreater(mass.to("MJ", assay.as_context(basis="hhv")).magnitude, mass.to("MJ", assay.as_context(basis="lhv")).magnitude)
        uranium = c.UraniumAssay.leu_4_5pct()
        for value in (0, -1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                c.create_uranium_context(uranium, burnup_gwd_per_mtu=value)
        for value in (0, -1, 1.1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                c.create_uranium_context(uranium, thermal_efficiency=value)
        assay = replace(c.PetrochemAssay.benzene(), purity_pct=90)
        context = c.create_petrochem_context(assay, active_only=True)
        mass = ureg.Quantity(Decimal("1000"), "kg")
        result = mass.to("mol", context)
        self.assertEqual(result, assay.mass_to_moles(mass, active_only=True))
        self.assertAlmostEqual(float(result.to("kg", context).magnitude), 1000)

    def test_all_factories_are_exported(self):
        names = ("biofuel", "ngl", "ore", "battery_material", "uranium", "biomass", "petrochem", "bullion", "fertilizer", "bitumen", "petcoke", "chemical")
        for name in names:
            self.assertIn(f"create_{name}_context", c.__all__)


class ChemicalAssayTests(unittest.TestCase):
    def test_presets_and_contexts(self):
        self.assertIn("ChemicalAssay", c.__all__)
        presets = ("caustic_soda_50pct", "caustic_soda_flakes", "sulphuric_acid_98pct", "sulphuric_acid_66be", "hydrochloric_acid_32pct", "methanol", "soda_ash_dense", "molten_sulphur", "hydrogen_peroxide_50pct", "liquid_chlorine")
        for preset in presets:
            assay = getattr(c.ChemicalAssay, preset)()
            for active_only in (True, False):
                with self.subTest(preset=preset, active_only=active_only):
                    context = assay.as_context(active_only=active_only, context_name="chemical_test")
                    self.assertEqual(context.name, "chemical_test")
                    volume = ureg.Quantity(Decimal("1"), "m^3")
                    mass = volume.to("kg", context)
                    self.assertAlmostEqual(float(mass.magnitude), assay.density_kg_per_m3)
                    expected_kmol = assay.density_kg_per_m3 / assay.molecular_weight_g_per_mol
                    if active_only:
                        expected_kmol *= assay.active_mass_fraction
                    amount = volume.to("kmol", context)
                    self.assertAlmostEqual(float(amount.magnitude), expected_kmol)
                    self.assertIsInstance(amount.magnitude, Decimal)
                    self.assertAlmostEqual(float(amount.to("m^3", context).magnitude), 1)
                    if assay.is_liquid:
                        self.assertIsNotNone(assay.baume_degrees)
                    else:
                        self.assertIsNone(assay.baume_degrees)

    def test_active_mass_and_molar_basis(self):
        assay = c.ChemicalAssay.caustic_soda_50pct()
        mass = ureg.Quantity(Decimal("10"), "metric_ton")
        active = assay.solution_to_active_mass(mass)
        self.assertEqual(active.magnitude, 5)
        self.assertEqual(assay.active_to_solution_mass(active), mass)
        volume = ureg.Quantity(Decimal("1"), "m^3")
        active = assay.volume_to_active_mass(volume)
        self.assertEqual(active.to("kg").magnitude, Decimal("762.5"))
        self.assertEqual(assay.active_mass_to_volume(active), volume)
        self.assertAlmostEqual(float(assay.mass_to_moles(ureg.Quantity(39.997, "kg"), active_only=False).magnitude), 1)
        self.assertEqual(mass.to("kg", assay.as_context()).magnitude, 10000)

    def test_dilution_and_blending_mass_balance(self):
        assay = c.ChemicalAssay.caustic_soda_50pct()
        mass = ureg.Quantity(Decimal("10"), "metric_ton")
        water = assay.water_required_to_dilute(mass, 25)
        self.assertEqual(water.magnitude, 10)
        diluted = assay.dilute_to(25)
        self.assertEqual(diluted.solution_to_active_mass(mass + water), assay.solution_to_active_mass(mass))
        self.assertEqual(assay.dilute_to(25, new_density_kg_per_m3=1275).density_kg_per_m3, 1275)
        blend = assay.blend_with(diluted, self_mass_fraction=0.4)
        self.assertEqual(blend.concentration_pct, 35)
        self.assertEqual(blend.solution_to_active_mass(mass).magnitude, Decimal("3.5"))
        self.assertGreater(assay.dilute_to(0.00001).concentration_pct, 0)

    def test_price_parity_and_properties(self):
        assay = c.ChemicalAssay.caustic_soda_50pct()
        self.assertEqual(assay.solution_price_to_active_price(250), 500)
        self.assertEqual(assay.active_price_to_solution_price(500), 250)
        self.assertAlmostEqual(assay.mt_price_to_gal_price(assay.gal_price_to_mt_price(1.25)), 1.25, places=4)
        self.assertAlmostEqual(assay.mt_price_to_st_price(assay.st_price_to_mt_price(300)), 300, places=4)
        self.assertEqual(assay.active_kg_per_m3, 762.5)
        self.assertAlmostEqual(assay.specific_gravity, 1525 / c.WATER_DENSITY_15C_KG_M3, places=6)
        self.assertAlmostEqual(assay.lb_per_gal, float(assay.volume_to_mass(ureg.Quantity(1, "gallon")).to("lb").magnitude), places=4)

    def test_invalid_specs_and_operations(self):
        assay = c.ChemicalAssay.caustic_soda_50pct()
        invalid = {"formula": ("", " "), "molecular_weight_g_per_mol": (0, -1, float("nan")), "concentration_pct": (0, 101, float("inf")), "density_kg_per_m3": (299, 3001, float("nan")), "reference_temp_c": (-101, 201, float("inf"))}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for target in (0, 50, 75, float("nan")):
            with self.assertRaises(ValueError):
                assay.dilute_to(target)
            with self.assertRaises(ValueError):
                assay.water_required_to_dilute(ureg.Quantity(1, "metric_ton"), target)
        for other in (c.ChemicalAssay.methanol(), c.ChemicalAssay.caustic_soda_flakes()):
            with self.assertRaises(ValueError):
                assay.blend_with(other, self_mass_fraction=0.5)
        with self.assertRaises(ValueError):
            assay.blend_with(assay, self_mass_fraction=1.1)
        with self.assertRaises(TypeError):
            assay.volume_to_mass(1)
        with self.assertRaises(pint.DimensionalityError):
            assay.volume_to_mass(ureg.Quantity(1, "kg"))
        for method in (assay.solution_price_to_active_price, assay.active_price_to_solution_price, assay.gal_price_to_mt_price, assay.mt_price_to_gal_price, assay.st_price_to_mt_price, assay.mt_price_to_st_price):
            with self.assertRaises(ValueError):
                method(float("nan"))


if __name__ == "__main__":
    unittest.main()
