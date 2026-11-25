"""Pressure datum, vacuum, barometric and differential regressions."""

from decimal import Decimal
import unittest
from unittest.mock import patch

import pint

from core.units import conversions as c
from core.units.registry import ureg


class PressureConversionTests(unittest.TestCase):
    def test_decimal_source_survives_float_returning_conversion(self):
        source = ureg.Quantity(Decimal("0"), "psig")
        original_to = type(source).to

        def float_source_conversion(quantity, *args, **kwargs):
            result = original_to(quantity, *args, **kwargs)
            if quantity is source:
                return ureg.Quantity(float(result.magnitude), result.units)
            return result

        atmosphere = Decimal("101325.123456789123456789")
        with patch.object(type(source), "to", float_source_conversion):
            result = c.convert_pressure(source, "Pa", p_atm=atmosphere)
            same_datum = c.convert_pressure(source, "barg", p_atm=atmosphere)
        self.assertEqual(result.magnitude, atmosphere)
        self.assertIsInstance(same_datum.magnitude, Decimal)

    def test_unqualified_barometric_units(self):
        for unit in ("torr", "atm", "mmHg", "inHg"):
            with self.subTest(unit=unit):
                pressure = ureg.Quantity(1, unit)
                result = c.convert_pressure(pressure, "psig")
                expected = pressure.to("Pa").magnitude - 101325
                self.assertAlmostEqual(float(c.convert_pressure(result, "Pa", p_atm=101325).magnitude) - 101325, float(expected))
                with self.assertRaises(ValueError):
                    c.convert_pressure(pressure, "psig", neutral_is_absolute=False)

    def test_canonical_aliases_and_display_formats(self):
        ureg.define("@alias psig = section5_gauge_alias")
        previous = ureg.formatter.default_format
        try:
            for display in ("", "~P", "L"):
                ureg.formatter.default_format = display
                for alias in ("psig", "PSIG", "lb_per_sq_in_gauge", "section5_gauge_alias"):
                    with self.subTest(display=display, alias=alias):
                        result = c.convert_pressure(ureg.Quantity(0, alias), "lb_per_sq_in_abs")
                        self.assertAlmostEqual(float(result.magnitude), c.STANDARD_ATMOSPHERE_PSI)
                for alias in ("barg", "bar_gauge", "BARG"):
                    self.assertAlmostEqual(float(c.convert_pressure(ureg.Quantity(0, alias), "bar_absolute").magnitude), 1.01325)
        finally:
            ureg.formatter.default_format = previous

    def test_neutral_units_and_strict_mode(self):
        self.assertAlmostEqual(float(c.convert_pressure(ureg.Quantity(100, "psig"), "psi").magnitude), 114.6959487755134)
        self.assertAlmostEqual(float(c.convert_pressure(ureg.Quantity(114.6959487755134, "psi"), "psig").magnitude), 100)
        for source, target in (("psig", "psi"), ("psi", "psig"), ("Pa", "kPa")):
            with self.subTest(source=source, target=target), self.assertRaises(ValueError):
                c.convert_pressure(ureg.Quantity(1, source), target, neutral_is_absolute=False)
        result = c.convert_pressure(ureg.Quantity(1, "barg"), "psig", neutral_is_absolute=False)
        self.assertAlmostEqual(float(result.magnitude), float(ureg.Quantity(1, "bar").to("psi").magnitude))

    def test_vacuum_and_same_datum_validation(self):
        vacuum = c.convert_pressure(ureg.Quantity(-10, "psig"), "psia")
        self.assertAlmostEqual(float(vacuum.magnitude), 4.6959487755134)
        self.assertEqual(c.convert_pressure(ureg.Quantity(-1.01325, "barg"), "bara").magnitude, 0)
        self.assertAlmostEqual(float(c.convert_pressure(ureg.Quantity(0, "bara"), "barg").magnitude), -1.01325)
        for source, magnitude, target in (("psia", -5, "bara"), ("psig", -20, "barg"), ("psig", -20, "psia"), ("Pa", -1, "kPa")):
            with self.subTest(source=source, target=target), self.assertRaises(ValueError):
                c.convert_pressure(ureg.Quantity(magnitude, source), target)
            c.convert_pressure(ureg.Quantity(magnitude, source), target, allow_negative_absolute=True)

    def test_atmosphere_precision_and_legacy_alias(self):
        atmosphere = Decimal("101325.123456789123456789")
        zero = ureg.Quantity(Decimal("0"), "psig")
        result = c.convert_pressure(zero, "Pa", p_atm=atmosphere)
        self.assertEqual(result.magnitude, atmosphere)
        self.assertIsInstance(result.magnitude, Decimal)
        self.assertEqual(c.convert_pressure(zero, "Pa", atmospheric_pressure=atmosphere), result)
        self.assertEqual(c.convert_pressure_datum(zero, "Pa", atmospheric_pressure=atmosphere), result)
        self.assertEqual(c.convert_pressure(zero, "Pa", p_atm=atmosphere, atmospheric_pressure=90000), result)
        self.assertEqual(c.convert_pressure_datum(zero, "Pa", p_atm=atmosphere, atmospheric_pressure=90000), result)
        self.assertEqual(c.convert_pressure(zero, "Pa", p_atm=ureg.Quantity(atmosphere, "Pa")), result)
        for value in (0, -1, float("nan"), float("inf"), Decimal("NaN"), ureg.Quantity(float("inf"), "Pa"), ureg.Quantity(1, "barg")):
            with self.subTest(value=value), self.assertRaises(ValueError):
                c.convert_pressure(zero, "Pa", p_atm=value)

    def test_elevation_reference_values_and_units(self):
        self.assertEqual(c.atmospheric_pressure_at_elevation(0).magnitude, 101325)
        high = c.atmospheric_pressure_at_elevation(2000)
        self.assertAlmostEqual(float(high.magnitude), 79495, delta=2)
        feet = ureg.Quantity(2000, "m").to("ft")
        self.assertAlmostEqual(c.atmospheric_pressure_at_elevation(feet).magnitude, high.magnitude)
        self.assertAlmostEqual(c.atmospheric_pressure_at_elevation(feet.magnitude, elevation_unit="ft").magnitude, high.magnitude)
        self.assertAlmostEqual(c.atmospheric_pressure_at_elevation(feet.magnitude, unit="ft").magnitude, high.magnitude)
        self.assertAlmostEqual(c.atmospheric_pressure_at_elevation(11000).magnitude, 22632, delta=2)
        self.assertGreater(c.atmospheric_pressure_at_elevation(-500).magnitude, 101325)
        self.assertEqual(c.atmospheric_pressure_at_elevation(0, sea_level_pa=100000).magnitude, 100000)
        self.assertAlmostEqual(c.convert_pressure(ureg.Quantity(0, "barg"), "Pa", p_atm=high).magnitude, high.magnitude)
        for height in (-501, 11001, float("nan"), ureg.Quantity(float("inf"), "m")):
            with self.assertRaises(ValueError):
                c.atmospheric_pressure_at_elevation(height)
        for baseline in (0, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                c.atmospheric_pressure_at_elevation(0, sea_level_pa=baseline)

    def test_differential_pressure_cancels_shared_atmosphere(self):
        p1 = ureg.Quantity(Decimal("2"), "barg")
        p2 = ureg.Quantity(Decimal("100"), "kPa")
        mixed = c.differential_pressure(p1, p2, "kPa")
        self.assertEqual(mixed.magnitude, Decimal("201.325"))
        self.assertIsInstance(mixed.magnitude, Decimal)
        for atmosphere in (90000, Decimal("101325.123456789")):
            result = c.differential_pressure(p1, ureg.Quantity(1, "barg"), "bar", p_atm=atmosphere)
            self.assertEqual(result.magnitude, 1)
            self.assertIsInstance(result.magnitude, Decimal)
        self.assertEqual(c.differential_pressure(ureg.Quantity(1, "bara"), ureg.Quantity(2, "bara"), "bar").magnitude, -1)

    def test_type_dimension_and_finite_checks(self):
        with self.assertRaises(TypeError):
            c.convert_pressure(1, "Pa")
        for source, target in (("kg", "kg"), ("Pa", "kg")):
            with self.assertRaises(pint.DimensionalityError):
                c.convert_pressure(ureg.Quantity(1, source), target)
        with self.assertRaises(pint.DimensionalityError):
            c.atmospheric_pressure_at_elevation(ureg.Quantity(1, "kg"))
        with self.assertRaises(TypeError):
            c.differential_pressure(1, ureg.Quantity(1, "Pa"))
        for value in (float("nan"), float("inf"), Decimal("NaN")):
            with self.assertRaises(ValueError):
                c.convert_pressure(ureg.Quantity(value, "psia"), "bara")
        with self.assertRaises(ValueError):
            c.convert_pressure(ureg.Quantity(1e308, "barg"), "bara", p_atm=1e308)
        with self.assertRaises(ValueError):
            c.convert_pressure(ureg.Quantity(0.0, "barg"), "bara", p_atm=Decimal("1e999"))
        self.assertIn("atmospheric_pressure_at_elevation", c.__all__)
        self.assertIn("differential_pressure", c.__all__)


if __name__ == "__main__":
    unittest.main()
