"""Composable thermal routing through the master conversion dispatcher."""

from dataclasses import asdict
from decimal import Decimal
import subprocess
import sys
import unittest
from unittest.mock import patch

from core.units import agricultural_specs as agricultural
from core.units import conversions as conversions
from core.units import heat_content as thermal
from core.units.registry import ureg


class HeatContentDispatchTests(unittest.TestCase):
    def setUp(self):
        self.grain = agricultural.GrainAssay.corn_56lb()
        self.fuel = thermal.HeatContentAssay(
            gcv=ureg.Quantity(Decimal("20"), "MJ/kg"),
            ncv=ureg.Quantity(Decimal("18"), "MJ/kg"),
            density=ureg.Quantity(Decimal("800"), "kg/m^3"),
            molar_mass=ureg.Quantity(Decimal("16"), "g/mol"),
        )

    def test_mass_volume_and_moles_round_trip_in_both_numeric_types(self):
        for source in ("kg", "m^3", "mol"):
            for target in ("MJ", "GJ", "MMBtu", "MWh", "therm"):
                for magnitude in (1.25, Decimal("1.25")):
                    with self.subTest(source=source, target=target, magnitude=magnitude):
                        original = ureg.Quantity(magnitude, source)
                        result = conversions.convert(original, target, assay=self.grain, heat_content=self.fuel)
                        expected = thermal.convert_heat_content(original, target, assay=self.fuel)
                        restored = conversions.convert(result, source, assay=self.grain, heat_content=self.fuel)
                        self.assertAlmostEqual(float(result.magnitude), float(expected.magnitude), places=8)
                        self.assertAlmostEqual(float(restored.magnitude), float(magnitude), places=8)
                        self.assertIsInstance(result.magnitude, type(magnitude))
                        self.assertIsInstance(restored.magnitude, type(magnitude))

    def test_explicit_spec_preserves_all_metadata_and_inputs(self):
        original_grain = asdict(self.grain)
        original_fuel = asdict(self.fuel)
        for assay in (self.grain, agricultural.EthanolFeedstockAssay.us_corn_ethanol()):
            result = conversions.convert(ureg.Quantity(Decimal(2), "m^3"), "MJ", assay=assay, heat_content=self.fuel)
            self.assertEqual(result.magnitude, Decimal(28800))
        self.assertEqual(asdict(self.grain), original_grain)
        self.assertEqual(asdict(self.fuel), original_fuel)

    def test_overrides_and_explicit_heat_basis(self):
        amount = ureg.Quantity(Decimal(100), "kg")
        self.assertEqual(conversions.convert(amount, "MJ", assay=self.grain, ncv=ureg.Quantity(17, "MJ/kg")).magnitude, 1700)
        self.assertEqual(conversions.convert(amount, "MJ", assay=self.grain, heat_content=self.fuel, basis="hhv").magnitude, 2000)
        self.assertEqual(conversions.convert(amount, "MJ", assay=self.fuel, ncv=ureg.Quantity(15, "MJ/kg")).magnitude, 1500)
        self.assertEqual(self.fuel.ncv.magnitude, 18)
        density = conversions.convert(
            ureg.Quantity(Decimal(18), "MJ/kg"), "MJ/kg", assay=self.fuel,
            moisture_pct=20, analytical_basis="DB",
        )
        self.assertEqual(density.magnitude, Decimal("23.1105"))

    def test_nonfuel_requires_calorific_data_but_plain_energy_scaling_does_not(self):
        with self.assertRaisesRegex(ValueError, "Supply at least one"):
            conversions.convert(ureg.Quantity(100, "kg"), "MJ", assay=self.grain)
        with self.assertRaises(TypeError):
            conversions.convert(ureg.Quantity(100, "kg"), "MJ", heat_content=object())
        result = conversions.convert(ureg.Quantity(Decimal(1000), "MJ"), "GJ", assay=self.grain)
        self.assertEqual(result.magnitude, 1)
        with self.assertRaises(TypeError):
            conversions.convert(ureg.Quantity(100, "kg"), "MJ", assay=self.fuel, misspelled_heat_field=1)

    def test_thermal_spec_does_not_replace_nonthermal_commodity_density(self):
        original = ureg.Quantity(Decimal(1), "m^3")
        expected = conversions.convert(original, "kg", assay=self.grain)
        with patch.object(thermal, "convert_heat_content", side_effect=AssertionError("unexpected thermal dispatch")):
            result = conversions.convert(original, "kg", assay=self.grain, heat_content=self.fuel)
            self.assertEqual(result, expected)
            pressure = conversions.convert(ureg.Quantity(0, "psig"), "psia", heat_content=self.fuel)
            self.assertGreater(pressure.magnitude, 0)
            fx = conversions.convert(ureg.Quantity(100, "USD"), "EUR", heat_content=self.fuel, rates={"USD/EUR": .9})
            self.assertAlmostEqual(fx.magnitude, 90)

    def test_energy_rate_and_fx_price_composition(self):
        rate = conversions.convert(ureg.Quantity(Decimal(1), "kg/hour"), "kW", assay=self.grain, heat_content=self.fuel)
        self.assertAlmostEqual(rate.magnitude, Decimal(5), places=12)
        quote = conversions.convert(
            ureg.Quantity(Decimal(100), "USD/metric_ton"), "EUR/GJ",
            assay=self.grain, heat_content=self.fuel, rates={"USD/EUR": Decimal("0.9")},
            quantize_settlement=True,
        )
        self.assertEqual(quote.magnitude, Decimal("5.00"))

    def test_nuclear_burnup_retains_its_noncalorific_path(self):
        assay = conversions.UraniumAssay.leu_4_5pct()
        original = ureg.Quantity(Decimal(10), "kg")
        with patch.object(thermal, "convert_heat_content", side_effect=AssertionError("unexpected calorific dispatch")):
            energy = conversions.convert(original, "MWh", assay=assay)
            restored = conversions.convert(energy, "kg", assay=assay)
        self.assertGreater(energy.magnitude, 0)
        self.assertAlmostEqual(restored.magnitude, original.magnitude, places=12)

    def test_scale_factors_preserve_decimal_precision_and_float_type(self):
        factor = Decimal("1.123456789012345678901234567")
        self.assertEqual(conversions._scale_factor(Decimal(1), factor), factor)
        self.assertIsInstance(conversions._scale_factor(1., factor), float)

    def test_import_orders_and_thermal_execution_in_fresh_processes(self):
        for first, second in (("conversions", "heat_content"), ("heat_content", "conversions")):
            script = (
                f"from core.units import {first}, {second}; "
                "from core.units.conversions import convert; "
                "from core.units.heat_content import HeatContentAssay; "
                "from core.units.registry import Q_; "
                "result = convert(Q_(2., 'kg'), 'MJ', assay=HeatContentAssay(ncv=Q_(18, 'MJ/kg'))); "
                "assert result.magnitude == 36. and isinstance(result.magnitude, float)"
            )
            result = subprocess.run([sys.executable, "-c", script], capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
