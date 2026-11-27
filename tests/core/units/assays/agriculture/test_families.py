"""Physical/financial routing and process balances for 18 agricultural families."""

from dataclasses import FrozenInstanceError, replace
from decimal import Decimal as D
import inspect
import re
import subprocess
import sys
import unittest

import pint

import core.units as u
from core.units import agricultural_specs as original, agricultural_families as families


def examples():
    common = dict(name="declared sample", density=u.Q_(D(800), "kg/m^3"), moisture_fraction=D(".12"),
                  package_mass=u.Q_(D(25), "kg"), foreign_matter_fraction=D(".02"))
    yield u.PulseAssay(**common, commodity="lentils", split_yield_fraction=D(".8"))
    yield u.SpecialtyGrainAssay(**common, commodity="sorghum", mass_per_bushel=u.Q_(D(56), "lb"))
    yield u.VegetableOilAssay(**{**common, "moisture_fraction": D(".001")}, commodity="palm", refining_recovery_fraction=D(".97"))
    yield u.OilseedAssay(**common, commodity="sesame", composition=u.ElementComposition({"oil": D(".4")}), oil_recovery_fraction=D(".9"), dry_meal_yield_fraction=D(".6"))
    yield u.StarchAssay(**common, commodity="corn", composition=u.ElementComposition({"starch": D(".98")}), recovery_fractions={"starch": D(".95")})
    yield u.RootCropAssay(**{**common, "moisture_fraction": D(".75")}, commodity="sugar-beet", composition=u.ElementComposition({"sugar": D(".6"), "starch": D(".1")}), recovery_fractions={"sugar": D(".9")})
    yield u.FruitConcentrateAssay(**common, commodity="orange", soluble_solids_fraction=D(".65"))
    yield u.DairyProductAssay(**common, commodity="WMP", composition=u.ElementComposition({"fat": D(".3"), "protein": D(".3"), "lactose": D(".35")}))
    yield u.AnimalFeedAssay(**common, commodity="protein-meal", composition=u.ElementComposition({"protein": D(".5"), "fat": D(".1")}))
    yield u.MeatProductAssay(**common, commodity="beef", saleable_yield_fraction=D(".7"))
    yield u.TeaAssay(**common, commodity="black", fresh_to_made_yield_fraction=D(".22"))
    yield u.TobaccoAssay(**common, commodity="leaf", curing_yield_fraction=D(".2"))
    yield u.SpiceAssay(**common, commodity="pepper", processing_recovery_fraction=D(".9"))
    yield u.NutAssay(**common, commodity="almonds", kernel_yield_fraction=D(".6"), shell_yield_fraction=D(".35"))
    yield u.TimberAssay(**{**common, "moisture_fraction": None}, commodity="softwood", moisture_content_dry_basis=D("1.25"), lumber_recovery_fraction=D(".5"))
    yield u.PaperAssay(**common, commodity="OCC", basis_weight=u.Q_(D(200), "g/m^2"), pulping_recovery_fraction=D(".8"))
    yield u.AgriculturalResidueAssay(**common, commodity="bagasse", collectable_fraction=D(".8"), composition=u.ElementComposition({"ash": D(".05")}))
    yield u.SilkAndNaturalFibreAssay(**common, commodity="silk", primary_recovery_fraction=D(".2"), degumming_recovery_fraction=D(".75"))


class AdditionalAgriculturalFamilyTests(unittest.TestCase):
    def test_all_18_families_exports_and_existing_import_paths(self):
        samples = list(examples())
        defined = {value for name, value in vars(families).items() if not name.startswith("_")
                   and inspect.isclass(value) and name in families.__all__ and name.endswith("Assay")}
        self.assertEqual(len(defined), 18)
        self.assertEqual({type(sample) for sample in samples}, defined)
        self.assertEqual(len(original.__all__), 99)
        self.assertEqual(len(original.__all__), len(set(original.__all__)))
        namespace = {}
        exec("from core.units.agricultural_specs import *", namespace)
        for name in families.__all__:
            self.assertIs(namespace[name], getattr(families, name))
            self.assertIs(getattr(u, name), getattr(families, name))
        legacy = {value for value in vars(original).values() if inspect.isclass(value)
                  and value.__name__ in original._LEGACY_EXPORTS and value.__name__.endswith("Assay")}
        self.assertEqual(len(legacy | defined), 33)

    def test_all_families_float_decimal_amount_flow_and_fx_quote_round_trips(self):
        for sample in examples():
            stem = re.sub(r"(?<!^)(?=[A-Z])", "_", type(sample).__name__.removesuffix("Assay")).lower()
            converter = getattr(u, "convert_" + stem)
            context = getattr(u, "create_" + stem + "_context")(sample)
            self.assertIsInstance(context, pint.Context)
            for numeric in (1.25, D("1.25")):
                for source, target in (("m^3", "kg"), ("m^3/day", "kg/hour"), ("USD/metric_ton", "USD/m^3")):
                    with self.subTest(family=type(sample).__name__, source=source, numeric=type(numeric)):
                        original_q = u.Q_(numeric, source)
                        output = converter(original_q, target, assay=sample)
                        restored = u.convert(output, source, assay=sample)
                        self.assertAlmostEqual(float(restored.magnitude), float(numeric), places=9)
                        self.assertIsInstance(output.magnitude, type(numeric))
                        self.assertIsInstance(restored.magnitude, type(numeric))
                fx_quote = u.convert(u.Q_(numeric, "USD/metric_ton"), "EUR/m^3", assay=sample, rates={"USD/EUR": D(".9")})
                self.assertAlmostEqual(float(fx_quote.magnitude), float(numeric) * .8 * .9, places=10)

    def test_every_listed_commodity_label_and_new_family_immutability(self):
        count = 0
        for sample in examples():
            commodities = sample._commodities or ({"orange", "apple", "grape", "pineapple"} if isinstance(sample, u.FruitConcentrateAssay) else {"pepper", "cardamom", "vanilla", "ginger", "turmeric"})
            for commodity in commodities:
                with self.subTest(family=type(sample).__name__, commodity=commodity):
                    changed = replace(sample, commodity=commodity)
                    self.assertEqual(changed.volume_to_mass(u.Q_(D(1), "m^3")).magnitude, 800)
                    count += 1
            with self.assertRaises(FrozenInstanceError):
                sample.commodity = "mutated"
            sample.density.ito("g/liter")
            self.assertEqual(sample.density.units, u.ureg.kg / u.ureg.meter ** 3)
        self.assertGreater(count, 65)

    def test_moisture_cleaning_packaging_and_price_invariants_all_families(self):
        for sample in examples():
            for numeric in (1000., D(1000)):
                with self.subTest(family=type(sample).__name__, numeric=type(numeric)):
                    quantity = u.Q_(numeric, "kg")
                    target = D(".1")
                    output = sample.moisture_standardised_mass(quantity, target_moisture_fraction=target)
                    expected_dry = D(str(numeric)) * (1 - sample.moisture_fraction)
                    self.assertAlmostEqual(float(output.magnitude) * .9, float(expected_dry), places=8)
                    price = sample.moisture_adjusted_price(u.Q_(numeric, "USD/metric_ton"), target_moisture_fraction=target)
                    self.assertAlmostEqual(float(price.magnitude) * float(output.magnitude), float(numeric) ** 2, places=5)
                    self.assertAlmostEqual(float(sample.clean_mass(quantity).magnitude), float(numeric) * .98, places=8)
                    packages = sample.mass_to_packages(quantity)
                    restored = sample.packages_to_mass(packages)
                    self.assertAlmostEqual(float(restored.magnitude), float(numeric), places=8)
                    self.assertIsInstance(output.magnitude, type(numeric))
                    self.assertIsInstance(packages.magnitude, type(numeric))
                    self.assertIsInstance(price.magnitude, type(numeric))
        pulse = next(examples())
        self.assertEqual(pulse.packages_to_mass(D("1.5")).magnitude, D("37.5"))
        with self.assertRaises(ValueError):
            pulse.moisture_adjusted_price(u.Q_(D(1), "USD/dmt"), target_moisture_fraction=D(".1"))

    def test_pulse_specialty_bushels_and_oilseed_mass_balances(self):
        samples = list(examples())
        pulse, grain, _, seeds = samples[:4]
        quantity = u.Q_(D(1000), "kg")
        self.assertEqual(pulse.split_mass(quantity).magnitude, 800)
        self.assertEqual(pulse.whole_pulse_required(pulse.split_mass(quantity)).magnitude, 1000)
        count = grain.mass_to_bushels(quantity)
        self.assertAlmostEqual(grain.bushels_to_mass(count).magnitude, 1000, places=20)
        self.assertEqual(u.Q_(1, "bushel").dimensionality, (u.ureg.meter ** 3).dimensionality)
        oil, dry_meal = seeds.oil_mass(quantity), seeds.dry_meal_mass(quantity)
        self.assertEqual(oil.magnitude, D("316.8"))
        self.assertEqual(dry_meal.magnitude, 528)
        self.assertLess(oil.magnitude + dry_meal.magnitude, 880)
        self.assertAlmostEqual(seeds.meal_mass(quantity, meal_moisture_fraction=D(".12")).magnitude, 600, places=20)
        for numeric in (1000., D(1000)):
            result = u.convert_oilseed(u.Q_(numeric, "kg"), "kg", assay=seeds, mode="oil_mass")
            self.assertAlmostEqual(float(result.magnitude), 316.8)
            self.assertIsInstance(result.magnitude, type(numeric))

    def test_vegetable_oil_refining_and_observed_temperature_routes(self):
        oil = u.VegetableOilAssay(name="palm", commodity="palm", density=u.Q_(D(900), "kg/m^3"),
              reference_temperature=u.Q_(D(20), "degC"), thermal_expansion_per_kelvin=D(".0008"),
              temperature_range=(u.Q_(D(10), "degC"), u.Q_(D(80), "degC")), refining_recovery_fraction=D(".97"))
        for numeric in (10., D(10)):
            quantity = u.Q_(numeric, "m^3")
            mass = u.convert(quantity, "kg", assay=oil, observed_temp=u.Q_(50, "degC"))
            self.assertAlmostEqual(float(mass.magnitude), float(numeric) * 900 / 1.024, places=8)
            restored = u.convert(mass, "m^3", assay=oil, observed_temp=50)
            self.assertAlmostEqual(float(restored.magnitude), float(numeric), places=9)
            self.assertIsInstance(mass.magnitude, type(numeric))
        input_mass = u.Q_(D(1000), "kg")
        self.assertEqual(oil.refined_mass(input_mass).magnitude, 970)
        self.assertAlmostEqual(oil.crude_oil_required(oil.refined_mass(input_mass)).magnitude, 1000, places=20)
        with self.assertRaises(ValueError):
            oil.volume_correction_factor(100)

    def test_starch_roots_dairy_and_animal_feed_blends(self):
        samples = list(examples())
        starch, root, dairy, feed = samples[4], samples[5], samples[7], samples[8]
        quantity = u.Q_(D(1000), "kg")
        self.assertEqual(starch.starch_mass(quantity).magnitude, D("862.4"))
        self.assertEqual(starch.recovered_starch_mass(quantity).magnitude, D("819.28"))
        self.assertEqual(root.sugar_mass(quantity).magnitude, 150)
        self.assertEqual(root.recovered_sugar_mass(quantity).magnitude, 135)
        self.assertEqual(dairy.milk_solids_mass(quantity).magnitude, 880)
        self.assertEqual(dairy.milk_fat_mass(quantity).magnitude, 264)
        self.assertEqual(dairy.fat_equivalent_mass(quantity, target_fat_fraction=D(".8"), fat_recovery_fraction=D(".95")).magnitude, D("313.5"))
        other = replace(feed, moisture_fraction=D(".5"), composition=u.ElementComposition({"protein": D(".2"), "fat": D(".05")}))
        blend = feed.blend(other, quantity, quantity, name="balanced feed")
        self.assertEqual(blend.moisture_fraction, D(".31"))
        protein = blend.nutrient_mass(u.Q_(D(2000), "kg"), nutrient="protein")
        original_protein = feed.nutrient_mass(quantity, nutrient="protein") + other.nutrient_mass(quantity, nutrient="protein")
        self.assertAlmostEqual(protein.magnitude, original_protein.magnitude, places=20)
        self.assertIsNone(blend.density)
        complete_a = replace(feed, composition=u.ElementComposition({"protein": D(".3"), "fat": D(".7")}))
        complete_b = replace(other, composition=u.ElementComposition({"protein": D(".4"), "fat": D(".6")}))
        complete_blend = complete_a.blend(complete_b, u.Q_(D(1), "kg"), u.Q_(D(2), "kg"), name="complete")
        self.assertEqual(complete_blend.composition.unassayed_fraction, 0)

    def test_fruit_reconstitution_solids_and_explicit_brix_assumption(self):
        juice = u.FruitConcentrateAssay(name="orange concentrate", commodity="orange", soluble_solids_fraction=D(".65"))
        quantity = u.Q_(D(1000), "kg")
        output = juice.reconstituted_mass(quantity, target_solids_fraction=D(".13"))
        self.assertEqual(output.magnitude, 5000)
        self.assertEqual(juice.dilution_water_required(quantity, target_solids_fraction=D(".13")).magnitude, 4000)
        self.assertEqual(juice.solids_mass(quantity).magnitude, output.magnitude * D(".13"))
        volume = juice.reconstituted_volume(quantity, target_solids_fraction=D(".13"), target_density=u.Q_(D(1040), "kg/m^3"))
        self.assertAlmostEqual(volume.magnitude * 1040, 5000, places=20)
        with self.assertRaisesRegex(ValueError, "explicitly"):
            u.FruitConcentrateAssay.from_brix(name="orange", commodity="orange", brix=u.Q_(D(65), "degree_Brix"))
        approximate = u.FruitConcentrateAssay.from_brix(name="orange", commodity="orange", brix=u.Q_(D(65), "degree_Brix"), assume_brix_equals_solids=True)
        self.assertEqual(approximate.soluble_solids_fraction, D(".65"))
        self.assertIn("approximation", approximate.solids_source)
        with self.assertRaises(ValueError):
            juice.reconstituted_mass(quantity, target_solids_fraction=D(".7"))

    def test_meat_tea_tobacco_spice_nuts_and_fibre_process_inverses(self):
        samples = list(examples())
        quantity = u.Q_(D(1000), "kg")
        for model, forward, inverse in (
            (samples[9], "saleable_mass", "raw_meat_required"),
            (samples[10], "made_tea_mass", "fresh_leaf_required"),
            (samples[11], "cured_mass", "green_leaf_required"),
            (samples[12], "processed_spice_mass", "raw_spice_required"),
            (samples[13], "kernel_mass", "inshell_mass_required"),
            (samples[17], "recovered_fibre_mass", "fibre_feedstock_required"),
        ):
            with self.subTest(family=type(model).__name__):
                output = getattr(model, forward)(quantity)
                self.assertAlmostEqual(getattr(model, inverse)(output).magnitude, 1000, places=20)
                self.assertLessEqual(output.magnitude, 1000)
        nut = samples[13]
        self.assertEqual((nut.kernel_mass(quantity) + nut.shell_mass(quantity)).magnitude, 950)
        fibre = samples[17]
        self.assertEqual(fibre.degummed_mass(fibre.recovered_fibre_mass(quantity)).magnitude, 150)

    def test_timber_dry_basis_moisture_and_volume_basis_inverses(self):
        timber = list(examples())[14]
        quantity = u.Q_(D(2250), "kg")
        self.assertAlmostEqual(timber.wet_to_dry_mass(quantity).magnitude, 1000, places=20)
        self.assertAlmostEqual(timber.dry_to_wet_mass(u.Q_(D(1000), "kg")).magnitude, 2250, places=20)
        self.assertGreater(timber.moisture_content_dry_basis, 1)
        volume = u.Q_(D(100), "m^3")
        self.assertEqual(timber.roundwood_required(timber.sawn_volume(volume)).magnitude, 100)
        self.assertEqual(timber.stacked_volume(timber.solid_volume(volume, solid_volume_fraction=D(".65")), solid_volume_fraction=D(".65")).magnitude, 100)
        with self.assertRaises(ValueError):
            replace(timber, moisture_fraction=D(".1"))

    def test_paper_area_rolls_pulping_and_cross_dimensional_quotes(self):
        paper = list(examples())[15]
        for numeric in (1000., D(1000)):
            area = u.Q_(numeric, "m^2")
            output = u.convert(area, "kg", assay=paper)
            self.assertAlmostEqual(float(output.magnitude), float(numeric) * .2, places=8)
            self.assertAlmostEqual(float(u.convert(output, "m^2", assay=paper).magnitude), float(numeric), places=8)
            self.assertIsInstance(output.magnitude, type(numeric))
        self.assertEqual(paper.roll_mass(u.Q_(D(1000), "m"), width=u.Q_(D(2), "m")).magnitude, 400)
        quantity = u.Q_(D(1000), "kg")
        self.assertEqual(paper.recovered_fibre_mass(quantity).magnitude, 800)
        self.assertEqual(paper.recovered_paper_required(paper.recovered_fibre_mass(quantity)).magnitude, 1000)
        quote = u.convert(u.Q_(D(100), "EUR/metric_ton"), "EUR/m^2", assay=paper)
        self.assertAlmostEqual(quote.magnitude, D(".02"), places=12)

    def test_residue_heat_and_reusable_material_processing_projections(self):
        residue = list(examples())[16]
        quantity = u.Q_(D(1000), "kg")
        self.assertEqual(residue.collectable_mass(quantity).magnitude, 800)
        self.assertEqual(residue.ash_mass(quantity).magnitude, 44)
        heat = u.HeatContentAssay(ncv=u.Q_(D(8), "MJ/kg"))
        fuel = replace(residue, heat_content=heat)
        self.assertEqual(u.convert(quantity, "MJ", assay=fuel).magnitude, 8000)
        process = u.ProcessingSpecification(input_material="root crop", output_yields={"starch": D(".2"), "pulp": D(".15")}, recovery_rates={"starch": D(".9")})
        root = u.RootCropAssay(name="cassava", commodity="cassava", processing=process)
        self.assertEqual(root.product_mass(quantity, product="starch").magnitude, 180)
        self.assertAlmostEqual(root.feedstock_required(root.product_mass(quantity, product="starch"), product="starch").magnitude, 1000, places=20)
        self.assertIs(u.material_spec_from_assay(root).processing, process)
        for sample in examples():
            projection = u.material_spec_from_assay(sample)
            self.assertEqual(projection.properties.density.magnitude, 800)

    def test_missing_inputs_bad_grades_and_unsupported_modes_are_rejected(self):
        with self.assertRaises(ValueError):
            u.PulseAssay(name="unknown", commodity="wheat")
        pulse = u.PulseAssay(name="peas", commodity="peas")
        for method in (pulse.split_mass, pulse.clean_mass, pulse.mass_to_packages, pulse.wet_to_dry_mass):
            with self.assertRaises(ValueError):
                method(u.Q_(D(1), "kg"))
        with self.assertRaises(ValueError):
            u.OilseedAssay(name="seeds", commodity="sesame", composition=u.ElementComposition({"oil": D(".5")}), oil_recovery_fraction=D(1), dry_meal_yield_fraction=D(".6"))
        with self.assertRaises(ValueError):
            u.NutAssay(name="nuts", commodity="almonds", kernel_yield_fraction=D(".7"), shell_yield_fraction=D(".4"))
        with self.assertRaises(ValueError):
            u.TeaAssay(name="tea", commodity="black", fresh_to_made_yield_fraction=float("nan"))
        with self.assertRaises(TypeError):
            u.StarchAssay(name="starch", commodity="corn", moisture_fraction=True)
        with self.assertRaises(ValueError):
            u.convert_pulse(u.Q_(D(1), "kg"), "kg", assay=pulse, mode="__init__")
        with self.assertRaises(TypeError):
            u.convert_pulse(u.Q_(D(1), "kg"), "kg", assay=list(examples())[1])

    def test_fresh_import_orders_do_not_create_cycles(self):
        for module in ("agricultural_specs", "agricultural_families", "_spec_engine"):
            script = f"from core.units import {module}; from core.units.agricultural_specs import PulseAssay; from core.units import Q_,convert; a=PulseAssay(name='peas',commodity='peas',density=Q_(800,'kg/m^3')); assert convert(Q_(1,'m^3'),'kg',assay=a).magnitude==800"
            result = subprocess.run([sys.executable, "-B", "-c", script], capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
