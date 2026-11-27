"""Ore moisture, grade, prices, rates and scalar helper regressions."""

from decimal import Decimal
import unittest

import pint

from core.units import conversions as c
from core.units.registry import ureg


class BulkOreConversionTests(unittest.TestCase):
    def test_moisture_and_reverse_extraction(self):
        wet = ureg.Quantity(Decimal("100000"), "wmt")
        dry = c.convert_wet_to_dry_mass(wet, Decimal("8"))
        self.assertEqual(dry.magnitude, 92000)
        self.assertEqual(c.convert_dry_to_wet_mass(dry, Decimal("8")), wet)
        for target in ("mt", "dmtu", "dltu", "lb"):
            metal = c.extract_contained_metal(dry, Decimal("62"), target_unit=target)
            restored = c.extract_dry_mass_from_contained(metal, Decimal("62"))
            self.assertAlmostEqual(float(restored.magnitude), 92000)
            self.assertIsInstance(metal.magnitude, Decimal)
        self.assertEqual(c.extract_contained_metal(ureg.Quantity(100000, "dmt"), 62, target_unit="dmtu").magnitude, 6200000)
        self.assertEqual(c.extract_contained_metal(dry, 0).magnitude, 0)
        with self.assertRaises(ValueError):
            c.extract_dry_mass_from_contained(ureg.Quantity(1, "mt"), 0)
        self.assertIn("extract_dry_mass_from_contained", c.__all__)

    def test_prices_preserve_actual_currency_and_mass_scale(self):
        for currency in ("USD", "GBP_pence", "US_cent"):
            price = ureg.Quantity(Decimal("110"), currency + " / dmt")
            result = c.convert_ore_moisture_price(price, "wmt", Decimal("8"))
            self.assertEqual(result.units, ureg.Unit(currency + " / wmt"))
            self.assertEqual(result.magnitude, Decimal("101.2"))
            self.assertEqual(c.convert_ore_moisture_price(result, "dmt", 8), price)
        generic = ureg.Quantity(Decimal("110"), "USD / metric_ton")
        result = c.convert_ore_moisture_price(generic, "wmt", 8)
        self.assertEqual(result.units, ureg.USD / ureg.wmt)
        self.assertEqual(result.magnitude, 110)
        result = c.convert_ore_moisture_price(generic, "kg", 8)
        self.assertEqual(result.units, ureg.USD / ureg.kg)
        self.assertEqual(result.magnitude, Decimal("0.11"))

    def test_structural_basis_in_amounts_rates_and_quotes(self):
        assay = c.BulkOreAssay(moisture_fraction=Decimal("0.08"), grades={"Fe": Decimal("0.62")})
        cases = (("wmt", "dmt", "92"), ("wmt/day", "dmt/hour", str(Decimal("92") / 24)), ("USD/dmt", "USD/wmt", "92"), ("GBP_pence/dmt", "GBP_pence/wmt", "92"), ("wmt", "dmtu", "5704"), ("dmtu", "dmt", str(Decimal("1") / Decimal("0.62"))))
        previous = ureg.formatter.default_format
        try:
            for display in ("", "~P", "L"):
                ureg.formatter.default_format = display
                for source, target, expected in cases:
                    with self.subTest(display=display, source=source, target=target):
                        quantity = ureg.Quantity(Decimal("100"), source)
                        result = c.convert_bulk_ore(quantity, target, assay=assay)
                        self.assertAlmostEqual(float(result.magnitude), float(expected))
                        self.assertIsInstance(result.magnitude, Decimal)
                        back = c.convert_bulk_ore(result, source, assay=assay)
                        self.assertAlmostEqual(float(back.magnitude), 100)
        finally:
            ureg.formatter.default_format = previous
        contained_quote = c.convert_bulk_ore(ureg.Quantity(Decimal("62"), "USD/dmt"), "USD/dmtu", assay=assay)
        self.assertEqual(contained_quote.magnitude, 1)
        scaled = c.convert_bulk_ore(ureg.Quantity(100, "wmt"), "kg", assay=assay)
        self.assertEqual(scaled.magnitude, 100000)

    def test_scalar_helper_freezes_grades_and_bridges_ore_assay(self):
        grades = {"Fe": Decimal("0.62")}
        assay = c.BulkOreAssay("iron", Decimal("0.08"), grades)
        grades["Fe"] = 0.1
        self.assertEqual(assay.grades["Fe"], Decimal("0.62"))
        with self.assertRaises(TypeError):
            assay.grades["Fe"] = 0.1
        self.assertEqual(hash(assay), hash(c.BulkOreAssay("iron", Decimal("0.08"), {"Fe": Decimal("0.62")})))
        self.assertEqual(assay.wmt_to_dmt(Decimal("100")), 92)
        self.assertEqual(assay.dmt_to_wmt(Decimal("92")), 100)
        self.assertEqual(assay.contained_metal_dmt(Decimal("100"), "Fe"), 62)
        self.assertEqual(assay.contained_dmtu(Decimal("100")), 6200)
        ore = assay.to_ore_assay()
        self.assertEqual(ore.grade_pct, 62)
        self.assertEqual(ore.moisture_pct, 8)
        quantity = ureg.Quantity(Decimal("100"), "wmt")
        self.assertEqual(c.convert_bulk_ore(quantity, "dmtu", assay=ore), c.convert_bulk_ore(quantity, "dmtu", assay=assay))
        with self.assertRaises(ValueError):
            c.BulkOreAssay(moisture_fraction=0.6, grades={"Fe": 0.62}).to_ore_assay()

    def test_decimal_percentages_retain_precision(self):
        percentage = Decimal("8.123456789123456789")
        wet = ureg.Quantity(Decimal("100"), "wmt")
        self.assertEqual(c.convert_wet_to_dry_mass(wet, percentage).magnitude, 100 - percentage)
        self.assertEqual(c.extract_contained_metal(ureg.Quantity(Decimal(100), "dmt"), percentage).magnitude, percentage)
        with self.assertRaises(pint.DimensionalityError):
            c.extract_contained_metal(wet, percentage)

    def test_invalid_percentages_fractions_dimensions_and_grades(self):
        for value in (float("nan"), float("inf"), Decimal("NaN"), Decimal("Infinity"), -1, 101):
            with self.subTest(value=value):
                for operation in (
                    lambda: c.convert_wet_to_dry_mass(ureg.Quantity(1, "wmt"), value),
                    lambda: c.convert_dry_to_wet_mass(ureg.Quantity(1, "dmt"), value),
                    lambda: c.convert_ore_moisture_price(ureg.Quantity(1, "USD/dmt"), "wmt", value),
                    lambda: c.extract_contained_metal(ureg.Quantity(1, "dmt"), value),
                    lambda: c.extract_dry_mass_from_contained(ureg.Quantity(1, "mt"), value),
                ):
                    with self.assertRaises(ValueError):
                        operation()
        for value in (Decimal("NaN"), float("inf"), -0.1, 1.1):
            with self.assertRaises(ValueError):
                c.BulkOreAssay(moisture_fraction=value)
            with self.assertRaises(ValueError):
                c.BulkOreAssay(grades={"Fe": value})
        with self.assertRaises(ValueError):
            c.convert_bulk_ore(ureg.Quantity(1, "wmt"), "dmtu")
        self.assertEqual(c.convert_bulk_ore(ureg.Quantity(1, "wmt"), "dmtu", grade=0, moisture_fraction=0).magnitude, 0)
        with self.assertRaises(ValueError):
            c.convert_bulk_ore(ureg.Quantity(1, "dmtu"), "dmt", grade=0)
        with self.assertRaises(ValueError):
            c.convert_bulk_ore(ureg.Quantity(1, "USD/dmt"), "USD/dmtu", grade=0)
        with self.assertRaises(ValueError):
            c.convert_ore_moisture_price(ureg.Quantity(1, "USD/dmt"), "meter", 8)
        with self.assertRaises(ValueError):
            c.convert_ore_moisture_price(ureg.Quantity(1, "dmt"), "wmt", 8)
        with self.assertRaises(pint.DimensionalityError):
            c.convert_bulk_ore(ureg.Quantity(1, "wmt/day"), "dmt", moisture_fraction=0.08)
        with self.assertRaises(TypeError):
            c.convert_bulk_ore(1, "dmt")


if __name__ == "__main__":
    unittest.main()

# Mirrored family contract coverage
import pytest as _contract_pytest
from tests.core.units._family_contracts import cases as _contract_cases
from tests.core.units._family_contracts import assert_snapshot_and_ownership as _snapshot_contract
from tests.core.units._family_contracts import assert_physical_or_energy_contract as _physical_contract

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.metals.ores'))
def test_ores_immutable_snapshot_and_foreign_registry_contract(specimen):
    _snapshot_contract(specimen)

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.metals.ores'))
def test_ores_independent_density_or_missing_factor_contract(specimen):
    _physical_contract(specimen)
