"""Crop-weight bushels, contractual shrink and drying charge regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class GrainAssayTests(unittest.TestCase):
    def test_exports_and_existing_agricultural_imports(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("SoftsAssay", "OilseedCrushAssay", "LivestockDairyAssay", "GrainAssay", "create_grain_context", "convert_grain"):
            self.assertIn(name, a.__all__)
        namespace = {}
        exec("from core.units.agricultural_specs import *", namespace)
        self.assertIs(namespace["GrainAssay"], a.GrainAssay)

    def test_crop_bushel_weights_match_registry(self):
        cases = (("corn_56lb", "bushel_corn", 56), ("wheat_60lb", "bushel_wheat", 60), ("soybeans_60lb", "bushel_soybeans", 60), ("barley_48lb", "bushel_barley", 48), ("oats_32lb", "bushel_oats", 32))
        for preset, bushel_unit, weight in cases:
            assay = getattr(a.GrainAssay, preset)()
            for count in (100.0, Decimal("100.123456789")):
                with self.subTest(preset=preset, count=count):
                    mass = assay.bushels_to_mass(count)
                    expected = ureg.Quantity(count, bushel_unit).to("metric_ton")
                    self.assertAlmostEqual(float(mass.magnitude), float(expected.magnitude))
                    recovered = assay.mass_to_bushels(mass)
                    self.assertAlmostEqual(float(recovered), float(count))
                    self.assertAlmostEqual(float(assay.mass_to_bushels(ureg.Quantity(weight, "lb"))), 1)
                    if isinstance(count, Decimal):
                        self.assertIsInstance(mass.magnitude, Decimal)
                        self.assertIsInstance(recovered, Decimal)
        self.assertNotEqual(ureg.bushel.dimensionality, ureg.bushel_corn.dimensionality)

    def test_shrink_handling_loss_and_reverse(self):
        assay = replace(a.GrainAssay.corn_56lb(actual_moisture_pct=20), handling_loss_pct=0.5)
        factor = (Decimal("80") / Decimal("84.5")) - Decimal("0.005")
        self.assertAlmostEqual(assay.moisture_adjustment_factor(), float(factor))
        for magnitude in (100.0, Decimal("100")):
            wet = ureg.Quantity(magnitude, "metric_ton")
            adjusted = assay.apply_moisture_deduction(wet)
            restored = assay.reverse_moisture_deduction(adjusted)
            self.assertAlmostEqual(float(adjusted.magnitude), 100 * float(factor))
            self.assertAlmostEqual(float(restored.magnitude), 100)
            if isinstance(magnitude, Decimal):
                self.assertIsInstance(adjusted.magnitude, Decimal)
            self.assertAlmostEqual(float(a.convert_grain(wet, "kg", assay=assay, mode="apply_moisture").magnitude), 100000 * float(factor))
            self.assertAlmostEqual(float(a.convert_grain(adjusted, "mt", assay=assay, mode="reverse_moisture").magnitude), 100)

    def test_base_dry_grain_and_zero_factor_rules(self):
        assay = replace(a.GrainAssay.corn_56lb(), handling_loss_pct=0.5)
        self.assertEqual(assay.moisture_adjustment_factor(), 1)
        self.assertEqual(assay.moisture_adjustment_factor(12), 1)
        enabled = replace(assay, allow_moisture_gain=True)
        self.assertGreater(enabled.moisture_adjustment_factor(12), 1)
        self.assertEqual(enabled.moisture_adjustment_factor(15.5), 1)
        self.assertEqual(assay.moisture_adjustment_factor(99.9), 0)
        with self.assertRaises(ValueError):
            assay.reverse_moisture_deduction(ureg.Quantity(1, "mt"), actual_moisture_pct=99.9)

    def test_drying_fee_price_and_decimal_moisture(self):
        assay = replace(a.GrainAssay.corn_56lb(actual_moisture_pct=20), drying_cost_per_mt_per_pct=2)
        for magnitude in (100.0, Decimal("100")):
            mass = ureg.Quantity(magnitude, "mt")
            charge = assay.drying_charge_usd(mass)
            self.assertEqual(charge, 900)
            self.assertEqual(assay.drying_charge_usd(mass, charge_per_mt_per_pct=Decimal("3")), 1350)
            self.assertEqual(assay.drying_charge_usd(mass, actual_moisture_pct=12), 0)
            if isinstance(magnitude, Decimal):
                self.assertIsInstance(charge, Decimal)
        self.assertAlmostEqual(
            assay.adjust_price_for_moisture(Decimal("200")),
            Decimal("200") * Decimal("80") / Decimal("84.5"),
            places=16,
        )
        precise = Decimal("20.123456789123456789")
        adjusted = assay.apply_moisture_deduction(ureg.Quantity(Decimal("84.5"), "mt"), actual_moisture_pct=precise)
        self.assertAlmostEqual(float(adjusted.magnitude), float(100 - precise))

    def test_bulk_context_rates_and_master_router(self):
        assay = a.GrainAssay.corn_56lb()
        for magnitude in (2.5, Decimal("2.5")):
            volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_grain(volume, "kg", assay=assay)
            self.assertAlmostEqual(float(mass.magnitude), 1800)
            self.assertAlmostEqual(float(a.convert_grain(mass, "m^3", assay=assay).magnitude), 2.5)
            rate = a.convert_grain(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay)
            self.assertAlmostEqual(float(rate.magnitude), 75)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/m^3"), "EUR/kg", assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 * 0.9 / 720)
        missing = replace(assay, bulk_density_kg_per_m3=None)
        self.assertFalse(missing.as_context().funcs)
        with self.assertRaises(ValueError):
            missing.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(pint.DimensionalityError):
            a.convert_grain(ureg.Quantity(1, "m^3"), "kg", assay=missing)

    def test_invalid_inputs_and_specs(self):
        assay = a.GrainAssay.corn_56lb()
        invalid = {"grain": ("", " "), "bushel_weight_lb": (0, -1, float("nan")), "base_moisture_pct": (-1, 100, float("inf")), "actual_moisture_pct": (-1, 100, Decimal("NaN")), "handling_loss_pct": (-1, 50), "drying_cost_per_mt_per_pct": (-1, float("nan")), "bulk_density_kg_per_m3": (99, 2001)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for value in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.bushels_to_mass(value)
            with self.assertRaises(ValueError):
                assay.mass_to_bushels(ureg.Quantity(value, "kg"))
        with self.assertRaises(ValueError):
            assay.drying_charge_usd(ureg.Quantity(1, "mt"), charge_per_mt_per_pct=Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_grain(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_grain(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()

# Mirrored family contract coverage
import pytest as _contract_pytest
from tests.core.units._family_contracts import cases as _contract_cases
from tests.core.units._family_contracts import assert_snapshot_and_ownership as _snapshot_contract
from tests.core.units._family_contracts import assert_physical_or_energy_contract as _physical_contract

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.grains'))
def test_grains_immutable_snapshot_and_foreign_registry_contract(specimen):
    _snapshot_contract(specimen)

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.grains'))
def test_grains_independent_density_or_missing_factor_contract(specimen):
    _physical_contract(specimen)
