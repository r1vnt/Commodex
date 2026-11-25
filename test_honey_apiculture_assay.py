"""Honey density, constant-solids basis and nominal package regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class HoneyApicultureAssayTests(unittest.TestCase):
    def test_exports_presets_and_volume_round_trips(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("HoneyApicultureAssay", "create_honey_apiculture_context", "convert_honey_apiculture"):
            self.assertIn(name, a.__all__)
        for assay in (a.HoneyApicultureAssay.bulk_polyfloral_export(), a.HoneyApicultureAssay.us_grade_a_clover()):
            for magnitude in (2.5, Decimal("2.5")):
                for unit in ("m^3", "liter", "gallon"):
                    with self.subTest(floral_source=assay.floral_source, unit=unit, magnitude=magnitude):
                        volume = ureg.Quantity(magnitude, unit)
                        mass = assay.volume_to_mass(volume)
                        restored = assay.mass_to_volume(mass).to(unit)
                        self.assertAlmostEqual(float(restored.magnitude), 2.5)
                        expected = volume.to("m^3").magnitude * (Decimal(str(assay.density_kg_per_m3)) if isinstance(magnitude, Decimal) else assay.density_kg_per_m3)
                        self.assertAlmostEqual(float(mass.to("kg").magnitude), float(expected))
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(mass.magnitude, Decimal)

    def test_moisture_solids_conservation_and_reverse(self):
        assay = a.HoneyApicultureAssay.bulk_polyfloral_export()
        for magnitude in (100.0, Decimal("100")):
            for actual in (Decimal("17.5"), Decimal("18"), Decimal("20.123456789")):
                mass = ureg.Quantity(magnitude, "mt")
                normalised = a.convert_honey_apiculture(mass, "mt", assay=assay, mode="normalise_moisture", actual_moisture_pct=actual)
                expected = float(magnitude) * float(100 - actual) / 82
                self.assertAlmostEqual(float(normalised.magnitude), expected)
                self.assertAlmostEqual(float(normalised.magnitude) * 0.82, float(magnitude) * float(100 - actual) / 100)
                restored = a.convert_honey_apiculture(normalised, "mt", assay=assay, mode="reverse_moisture", actual_moisture_pct=actual)
                self.assertAlmostEqual(float(restored.magnitude), float(magnitude))
                if isinstance(magnitude, Decimal):
                    self.assertIsInstance(normalised.magnitude, Decimal)
        self.assertGreater(assay.solids_normalisation_factor, 1)
        self.assertEqual(assay.density_kg_per_m3, 1420)

    def test_approximate_brix_and_beeswax(self):
        assay = a.HoneyApicultureAssay.bulk_polyfloral_export()
        self.assertEqual(assay.brix_solids_pct, 82.5)
        for moisture in (17.5, Decimal("17.123456789")):
            brix = assay.moisture_to_brix(moisture)
            self.assertAlmostEqual(float(assay.brix_to_moisture(brix)), float(moisture))
            if isinstance(moisture, Decimal):
                self.assertIsInstance(brix, Decimal)
        mass = ureg.Quantity(Decimal("100"), "kg")
        wax = assay.honey_to_beeswax_mass(mass)
        self.assertEqual(wax.magnitude, Decimal("1.8"))
        self.assertEqual(a.convert_honey_apiculture(mass, "kg", assay=assay, mode="honey_to_beeswax"), wax)
        self.assertEqual(replace(assay, beeswax_ratio_pct_of_honey=0).honey_to_beeswax_mass(mass).magnitude, 0)

    def test_drum_tote_and_pail_counts(self):
        assay = a.HoneyApicultureAssay.bulk_polyfloral_export()
        cases = ((assay.drums_to_mass, assay.mass_to_drums, 300), (assay.ibc_totes_to_mass, assay.mass_to_ibc_totes, 1400), (assay.pails_to_mass, assay.mass_to_pails, 27.2155422))
        for count in (10.0, Decimal("10.123456789")):
            for forward, reverse, unit_mass in cases:
                mass = forward(count)
                self.assertAlmostEqual(float(mass.to("kg").magnitude), float(count) * unit_mass)
                recovered = reverse(mass)
                self.assertAlmostEqual(float(recovered), float(count))
                if isinstance(count, Decimal):
                    self.assertIsInstance(mass.magnitude, Decimal)
                    self.assertIsInstance(recovered, Decimal)
        self.assertAlmostEqual(float(assay.pails_to_mass(Decimal("1")).to("lb").magnitude), 60)
        custom = replace(assay, drum_mass_kg=290)
        self.assertEqual(custom.drums_to_mass(Decimal("10")).magnitude, Decimal("2.9"))

    def test_prices_rates_and_quote_router(self):
        assay = a.HoneyApicultureAssay.bulk_polyfloral_export()
        for price in (2.0, Decimal("2")):
            converted = assay.usd_per_lb_to_usd_per_mt(price)
            expected = ureg.Quantity(price, "USD/lb").to("USD/mt")
            self.assertAlmostEqual(float(converted), float(expected.magnitude), places=4)
            self.assertAlmostEqual(float(assay.usd_per_mt_to_usd_per_lb(converted)), 2, places=4)
            if isinstance(price, Decimal):
                self.assertIsInstance(converted, Decimal)
        volume = ureg.Quantity(Decimal("24"), "m^3/day")
        rate = a.convert_honey_apiculture(volume, "kg/hour", assay=assay)
        self.assertAlmostEqual(float(rate.magnitude), 1420)
        quote = c.convert_commodity_quote(ureg.Quantity(Decimal("1000"), "USD/m^3"), "EUR/kg", assay=assay, rates={"EUR": Decimal("0.9")})
        self.assertAlmostEqual(float(quote.magnitude), 900 / 1420)
        with self.assertRaises(pint.DimensionalityError):
            ureg.Quantity(1, "m^3").to("kg")

    def test_invalid_specs_values_and_modes(self):
        assay = a.HoneyApicultureAssay.bulk_polyfloral_export()
        invalid = {"floral_source": ("", None), "moisture_pct": (9, 31, Decimal("NaN")), "base_moisture_pct": (9, 31), "density_kg_per_m3": (1299, 1501), "drum_mass_kg": (0, -1), "beeswax_ratio_pct_of_honey": (-1, 26)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.volume_to_mass(ureg.Quantity(bad, "m^3"))
            with self.assertRaises(ValueError):
                assay.honey_to_beeswax_mass(ureg.Quantity(bad, "kg"))
            with self.assertRaises(ValueError):
                assay.drums_to_mass(bad)
        for bad in (9, 31, Decimal("NaN")):
            with self.assertRaises(ValueError):
                assay.normalise_to_base_moisture(ureg.Quantity(1, "mt"), actual_moisture_pct=bad)
        with self.assertRaises(ValueError):
            assay.brix_to_moisture(60)
        with self.assertRaises(ValueError):
            assay.usd_per_lb_to_usd_per_mt(Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_honey_apiculture(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_honey_apiculture(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()
