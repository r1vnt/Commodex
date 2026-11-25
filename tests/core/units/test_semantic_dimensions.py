"""Semantic dimensions are an invariant, not merely a .to() runtime check."""

from decimal import Decimal
import unittest

import pint

import core.units as units
from core.units import conversions as c


class SemanticDimensionTests(unittest.TestCase):
    def test_instruments_and_mass_bases_are_distinct_dimensions(self):
        groups = (
            ("EUA", "UKA", "CCA", "VCU", "kg"),
            ("REC", "REGO", "MWh"),
            ("wmt", "dmt", "contained_metal_tonne", "payable_metal_tonne", "kg"),
        )
        for group in groups:
            dimensions = [units.parse_unit(name).dimensionality for name in group]
            self.assertEqual(len(set(dimensions)), len(group))
            for source in group:
                for target in group:
                    if source != target:
                        with self.subTest(source=source, target=target), self.assertRaises(pint.DimensionalityError):
                            units.Q_(100, source).to(target)

    def test_aliases_and_contained_subdivisions_remain_consistent(self):
        self.assertEqual(units.Q_(1, "WMT").to("wet_metric_tonne").magnitude, 1)
        self.assertEqual(units.Q_(100, "DMTU").to("contained_metal_tonne").magnitude, 1)

    def test_wet_dry_contained_and_payable_round_trip(self):
        assay = c.OreAssay.copper_concentrate()
        original = units.Q_(Decimal(100), "wmt")
        for target in ("dmt", "contained_metal_tonne", "payable_metal_tonne"):
            converted = units.convert(original, target, assay=assay)
            restored = units.convert(converted, original.units, assay=assay)
            self.assertIsInstance(converted.magnitude, Decimal)
            self.assertAlmostEqual(float(restored.magnitude), 100, places=10)
        dry = assay.wmt_to_dmt(original)
        self.assertEqual(dry.units, units.ureg.dmt)
        self.assertEqual(assay.dmt_to_contained_metal(dry).units, units.ureg.contained_metal_tonne)
        self.assertEqual(assay.dmt_to_payable_metal(dry).units, units.ureg.payable_metal_tonne)

    def test_moisture_grade_and_payability_are_required(self):
        for target in ("dmt", "contained_metal_tonne", "payable_metal_tonne"):
            with self.assertRaises((ValueError, pint.DimensionalityError)):
                units.convert(units.Q_(1, "wmt"), target)
        with self.assertRaisesRegex(ValueError, "Payable"):
            c.convert_bulk_ore(units.Q_(1, "dmt"), "payable_metal_tonne", grade=.25)
        assay = c.OreAssay.copper_concentrate()
        with self.assertRaisesRegex(ValueError, "dry-basis"):
            assay.dmt_to_contained_metal(units.Q_(1, "wmt"))
        with self.assertRaisesRegex(ValueError, "Contained-metal"):
            c.extract_dry_mass_from_contained(units.Q_(1, "dmt"), 25)
        with self.assertRaisesRegex(ValueError, "output unit"):
            c.extract_contained_metal(units.Q_(1, "dmt"), 25, target_unit="wmt")
        with self.assertRaisesRegex(ValueError, "output unit"):
            c.extract_dry_mass_from_contained(units.Q_(1, "contained_metal_tonne"), 25, target_unit="payable_metal_tonne")

    def test_semantic_rates_and_fx_quotes(self):
        assay = c.OreAssay.copper_concentrate()
        for source, target in (("wmt/day", "dmt/hour"), ("USD/wmt", "EUR/contained_metal_tonne")):
            original = units.Q_(Decimal(24), source)
            kwargs = {"rates": {"EUR": Decimal(".9")}} if "USD" in source else {}
            converted = units.convert(original, target, assay=assay, **kwargs)
            restored = units.convert(converted, source, assay=assay, **kwargs)
            self.assertAlmostEqual(float(restored.magnitude), 24, places=10)

    def test_thermal_rate_retains_semantic_mass_on_reverse_conversion(self):
        assay = c.PetcokeAssay.usgc_fuel_grade_65s()
        original = units.Q_(Decimal(24), "wmt/day")
        energy = units.convert(original, "MW", assay=assay)
        restored = units.convert(energy, original.units, assay=assay)
        self.assertIsInstance(restored.magnitude, Decimal)
        self.assertAlmostEqual(float(restored.magnitude), 24, places=9)
