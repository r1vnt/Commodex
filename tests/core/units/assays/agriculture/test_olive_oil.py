"""Olive-oil acidity screens, thermal volume, refining and price regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class OliveOilAssayTests(unittest.TestCase):
    def test_exports_presets_and_acidity_only_screen(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("OliveOilAssay", "create_olive_oil_context", "convert_olive_oil"):
            self.assertIn(name, a.__all__)
        evoo = a.OliveOilAssay.extra_virgin_evoo()
        for grade, threshold in (("EVOO", 0.8), ("VIRGIN", 2.0), ("REFINED", 0.3)):
            self.assertTrue(replace(evoo, grade=grade, free_fatty_acid_pct=threshold).meets_grade_acidity_limit)
            self.assertFalse(replace(evoo, grade=grade, free_fatty_acid_pct=threshold + 0.0001).meets_grade_acidity_limit)
        self.assertIsNone(a.OliveOilAssay.lampante_crude().meets_grade_acidity_limit)
        self.assertIsNone(replace(evoo, grade="POMACE").meets_grade_acidity_limit)
        self.assertTrue(a.OliveOilAssay.refined_olive_oil().meets_grade_acidity_limit)

    def test_density_and_temperature_round_trips(self):
        assay = a.OliveOilAssay.extra_virgin_evoo()
        for magnitude in (2.5, Decimal("2.5")):
            for temperature in (None, -20, 20, Decimal("40"), 100):
                with self.subTest(magnitude=magnitude, temperature=temperature):
                    volume = ureg.Quantity(magnitude, "m^3")
                    mass = assay.volume_to_mass(volume, temp_c=temperature)
                    restored = assay.mass_to_volume(mass, temp_c=temperature)
                    self.assertAlmostEqual(float(restored.magnitude), 2.5)
                    vcf = 1 if temperature is None else 1 - 0.00072 * (float(temperature) - 20)
                    self.assertAlmostEqual(float(mass.to("kg").magnitude), 2.5 * 915 * vcf)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(mass.magnitude, Decimal)
        density = assay.volume_to_mass(ureg.Quantity(Decimal("1"), "m^3"), temp_c=40).to("kg")
        self.assertAlmostEqual(float(density.magnitude), 901.824)
        fahrenheit = assay.volume_to_mass(ureg.Quantity(1, "m^3"), temp_c=ureg.Quantity(68, "degF"))
        self.assertAlmostEqual(float(fahrenheit.to("kg").magnitude), 915)

    def test_converter_temperature_does_not_misroute_same_dimension(self):
        assay = a.OliveOilAssay.extra_virgin_evoo()
        for magnitude in (2.5, Decimal("2.5")):
            for source, target in (("kg", "lb"), ("m^3", "liter")):
                quantity = ureg.Quantity(magnitude, source)
                result = a.convert_olive_oil(quantity, target, assay=assay, temp_c=40)
                self.assertAlmostEqual(float(result.magnitude), float(quantity.to(target).magnitude))
            volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_olive_oil(volume, "kg", assay=assay, temp_c=40)
            expected = assay.volume_to_mass(volume, temp_c=40).to("kg")
            self.assertAlmostEqual(float(mass.magnitude), float(expected.magnitude))
            back = a.convert_olive_oil(mass, "m^3", assay=assay, temp_c=40)
            self.assertAlmostEqual(float(back.magnitude), 2.5)
            rate = a.convert_olive_oil(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay)
            self.assertAlmostEqual(float(rate.magnitude), float(magnitude) * 915 / 24)
            hot_rate = a.convert_olive_oil(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay, temp_c=40)
            self.assertAlmostEqual(float(hot_rate.magnitude), float(magnitude) * 901.824 / 24)
            restored = a.convert_olive_oil(hot_rate, "m^3/day", assay=assay, temp_c=40)
            self.assertAlmostEqual(float(restored.magnitude), float(magnitude))

    def test_refining_yield_and_explicit_modes(self):
        assay = a.OliveOilAssay.lampante_crude(refining_yield_pct=95)
        for magnitude in (100.0, Decimal("100")):
            crude = ureg.Quantity(magnitude, "mt")
            refined = a.convert_olive_oil(crude, "mt", assay=assay, mode="crude_to_refined")
            self.assertAlmostEqual(float(refined.magnitude), 95)
            restored = a.convert_olive_oil(refined, "mt", assay=assay, mode="refined_to_crude")
            self.assertAlmostEqual(float(restored.magnitude), 100)
            if isinstance(magnitude, Decimal):
                self.assertIsInstance(refined.magnitude, Decimal)
        self.assertEqual(a.convert_olive_oil(ureg.Quantity(100, "mt"), "mt", assay=assay).magnitude, 100)

    def test_price_helpers_and_combined_fx_thermal_quotes(self):
        assay = a.OliveOilAssay.extra_virgin_evoo()
        self.assertEqual(assay.price_per_100kg_to_price_per_mt(Decimal("500")), 5000)
        self.assertEqual(assay.price_per_mt_to_price_per_100kg(Decimal("5000")), 500)
        for price in (5.0, Decimal("5")):
            gallon_price = assay.price_per_kg_to_price_per_us_gallon(price, temp_c=40)
            self.assertAlmostEqual(float(gallon_price), 5 * 901.824 * 0.003785411784, places=4)
            self.assertAlmostEqual(float(assay.price_per_us_gallon_to_price_per_kg(gallon_price, temp_c=40)), 5, places=4)
            if isinstance(price, Decimal):
                self.assertIsInstance(gallon_price, Decimal)
        quote = c.convert(ureg.Quantity(Decimal("2"), "EUR/kg"), "USD/gallon", assay=assay, observed_temp=40, rates={"EUR/USD": Decimal("1.1")})
        self.assertAlmostEqual(float(quote.magnitude), 2.2 * 901.824 * 0.003785411784)
        quintal = c.convert(ureg.Quantity(Decimal("100"), "EUR/quintal"), "USD/mt", assay=assay, rates={"EUR/USD": Decimal("1.1")})
        self.assertEqual(quintal.magnitude, 1100)

    def test_invalid_values_specs_and_no_silent_temperature_clamping(self):
        assay = a.OliveOilAssay.extra_virgin_evoo()
        invalid = {"grade": ("unknown", None), "free_fatty_acid_pct": (-1, 31, Decimal("NaN")), "refining_yield_pct": (0, 101), "density_kg_per_m3": (849, 981), "reference_temp_c": (-21, 101)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for temperature in (-21, 101, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.as_context(temp_c=temperature)
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.volume_to_mass(ureg.Quantity(bad, "m^3"))
            with self.assertRaises(ValueError):
                assay.crude_to_refined_mass(ureg.Quantity(bad, "kg"))
        with self.assertRaises(ValueError):
            assay.price_per_kg_to_price_per_us_gallon(Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_olive_oil(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_olive_oil(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")
        with self.assertRaises(pint.DimensionalityError):
            ureg.Quantity(1, "m^3").to("kg")


if __name__ == "__main__":
    unittest.main()
