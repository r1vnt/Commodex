"""Natural rubber DRC, nominal bales, bulk and price parity regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class RubberAssayTests(unittest.TestCase):
    def test_exports_presets_and_drc_round_trips(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("RubberAssay", "create_rubber_context", "convert_rubber"):
            self.assertIn(name, a.__all__)
        presets = (("tsr20_standard", 100), ("rss3_bale_standard", 100), ("field_latex_60drc", 60), ("cup_lump_55drc", 55))
        for preset, expected in presets:
            assay = getattr(a.RubberAssay, preset)()
            for magnitude in (100.0, Decimal("100")):
                with self.subTest(preset=preset, magnitude=magnitude):
                    wet = ureg.Quantity(magnitude, "mt")
                    dry = assay.wet_to_dry_mass(wet)
                    self.assertAlmostEqual(float(dry.magnitude), expected)
                    self.assertAlmostEqual(float(assay.dry_to_wet_mass(dry).magnitude), 100)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(dry.magnitude, Decimal)

    def test_lot_drc_override_precision_and_modes(self):
        assay = a.RubberAssay.field_latex_60drc()
        wet = ureg.Quantity(Decimal("100"), "mt")
        override = Decimal("55.123456789123456789")
        dry = assay.wet_to_dry_mass(wet, drc_pct=override)
        self.assertEqual(dry.magnitude, override)
        self.assertAlmostEqual(float(assay.dry_to_wet_mass(dry, drc_pct=override).magnitude), 100)
        result = a.convert_rubber(wet, "kg", assay=assay, mode="wet_to_dry", drc_pct=50)
        self.assertEqual(result.magnitude, 50000)
        restored = a.convert_rubber(result, "mt", assay=assay, mode="dry_to_wet", drc_pct=50)
        self.assertEqual(restored.magnitude, 100)
        self.assertEqual(a.convert_rubber(wet, "kg", assay=assay).magnitude, 100000)

    def test_bale_counts_preserve_precision_and_custom_sizes(self):
        for assay in (a.RubberAssay.tsr20_standard(), a.RubberAssay.rss3_bale_standard(), replace(a.RubberAssay.tsr20_standard(), bale_mass_kg=35)):
            for count in (30.0, Decimal("30.123456789")):
                with self.subTest(grade=assay.grade, bale=assay.bale_mass_kg, count=count):
                    mass = assay.bales_to_mass(count)
                    self.assertAlmostEqual(float(mass.to("kg").magnitude), float(count) * assay.bale_mass_kg)
                    self.assertAlmostEqual(float(assay.mass_to_bales(mass)), float(count))
                    if isinstance(count, Decimal):
                        self.assertIsInstance(mass.magnitude, Decimal)
                        self.assertIsInstance(assay.mass_to_bales(mass), Decimal)
        self.assertAlmostEqual(a.RubberAssay.tsr20_standard().bales_to_mass(30).magnitude, 1)
        self.assertEqual(replace(a.RubberAssay.tsr20_standard(), bale_mass_kg=35).bales_to_mass(Decimal("100")).magnitude, Decimal("3.5"))

    def test_quote_and_netback_parity_matches_registry(self):
        assay = a.RubberAssay.field_latex_60drc()
        for price in (150.0, Decimal("150")):
            converted = assay.cents_per_kg_to_usd_per_mt(price)
            self.assertEqual(converted, 1500)
            self.assertEqual(assay.usd_per_mt_to_cents_per_kg(converted), 150)
            expected = ureg.Quantity(price, "US_cent/kg").to("USD/metric_ton")
            self.assertAlmostEqual(float(converted), float(expected.magnitude))
            if isinstance(price, Decimal):
                self.assertIsInstance(converted, Decimal)
        self.assertEqual(assay.wet_price_to_dry_price_per_mt(Decimal("600")), 1000)
        self.assertEqual(assay.wet_price_to_dry_price_per_mt(Decimal("600"), drc_pct=50), 1200)

    def test_bulk_context_rates_master_and_missing_density(self):
        assay = a.RubberAssay.field_latex_60drc()
        for magnitude in (2.5, Decimal("2.5")):
            source = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_rubber(source, "kg", assay=assay)
            self.assertEqual(mass.magnitude, 2375)
            self.assertAlmostEqual(float(a.convert_rubber(mass, "m^3", assay=assay).magnitude), 2.5)
            rate = a.convert_rubber(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay)
            self.assertAlmostEqual(float(rate.magnitude), 2375 / 24)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/m^3"), "EUR/kg", assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 * 0.9 / 950)
        missing = replace(assay, bulk_density_kg_per_m3=None)
        self.assertFalse(missing.as_context().funcs)
        with self.assertRaises(ValueError):
            missing.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(pint.DimensionalityError):
            a.convert_rubber(ureg.Quantity(1, "m^3"), "kg", assay=missing)

    def test_invalid_specs_overrides_and_values(self):
        assay = a.RubberAssay.tsr20_standard()
        invalid = {"grade": ("", " ", None), "dry_rubber_content_pct": (0, 101, Decimal("NaN")), "bale_mass_kg": (0, -1, float("inf")), "bulk_density_kg_per_m3": (99, 2001)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for override in (0, -1, 101, Decimal("NaN"), Decimal("Infinity")):
            for method in (assay.wet_to_dry_mass, assay.dry_to_wet_mass):
                with self.assertRaises(ValueError):
                    method(ureg.Quantity(1, "mt"), drc_pct=override)
            with self.assertRaises(ValueError):
                assay.wet_price_to_dry_price_per_mt(100, drc_pct=override)
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.wet_to_dry_mass(ureg.Quantity(bad, "kg"))
            with self.assertRaises(ValueError):
                assay.bales_to_mass(bad)
        for method in (assay.cents_per_kg_to_usd_per_mt, assay.usd_per_mt_to_cents_per_kg, assay.wet_price_to_dry_price_per_mt):
            with self.assertRaises(ValueError):
                method(Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_rubber(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_rubber(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()

# Mirrored family contract coverage
import pytest as _contract_pytest
from tests.core.units._family_contracts import cases as _contract_cases
from tests.core.units._family_contracts import assert_snapshot_and_ownership as _snapshot_contract
from tests.core.units._family_contracts import assert_physical_or_energy_contract as _physical_contract

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.rubber'))
def test_rubber_immutable_snapshot_and_foreign_registry_contract(specimen):
    _snapshot_contract(specimen)

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.rubber'))
def test_rubber_independent_density_or_missing_factor_contract(specimen):
    _physical_contract(specimen)
