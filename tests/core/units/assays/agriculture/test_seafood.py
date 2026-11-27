"""Seafood processing, glaze, reduction and registry integration regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class SeafoodAssayTests(unittest.TestCase):
    def test_exports_and_all_presets(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("SeafoodAssay", "create_seafood_context", "convert_seafood"):
            self.assertIn(name, a.__all__)
        for preset in ("atlantic_salmon", "american_lobster", "snow_crab", "vannamei_shrimp", "anchoveta_reduction"):
            assay = getattr(a.SeafoodAssay, preset)()
            for magnitude in (100.0, Decimal("100")):
                whole = ureg.Quantity(magnitude, "kg")
                for forward, inverse, yield_pct in ((assay.whole_to_dressed, assay.dressed_to_whole, assay.dressed_yield_pct), (assay.whole_to_edible, assay.edible_to_whole, assay.edible_yield_pct), (assay.net_to_glazed_mass, assay.glazed_to_net_mass, None)):
                    with self.subTest(preset=preset, method=forward.__name__, magnitude=magnitude):
                        result = forward(whole)
                        self.assertAlmostEqual(float(inverse(result).magnitude), 100)
                        if yield_pct is not None:
                            self.assertAlmostEqual(float(result.magnitude), yield_pct)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(result.magnitude, Decimal)

    def test_glaze_basis_and_explicit_modes(self):
        assay = a.SeafoodAssay.vannamei_shrimp(glaze_pct=10)
        self.assertEqual(assay.glazed_to_net_mass(ureg.Quantity(Decimal("100"), "kg")).magnitude, 90)
        self.assertEqual(assay.net_to_glazed_mass(ureg.Quantity(Decimal("90"), "kg")).magnitude, 100)
        # Edible yield is a percentage of whole mass, not of dressed mass.
        whole = ureg.Quantity(Decimal("100"), "kg")
        dressed = a.convert_seafood(whole, "kg", assay=assay, mode="whole_to_dressed")
        edible = a.convert_seafood(whole, "kg", assay=assay, mode="whole_to_edible")
        self.assertEqual(dressed.magnitude, 65)
        self.assertEqual(edible.magnitude, 50)
        for forward, reverse in (("whole_to_dressed", "dressed_to_whole"), ("whole_to_edible", "edible_to_whole"), ("net_to_glazed", "glazed_to_net")):
            changed = a.convert_seafood(whole, "lb", assay=assay, mode=forward)
            back = a.convert_seafood(changed, "kg", assay=assay, mode=reverse)
            self.assertAlmostEqual(float(back.magnitude), 100)
        self.assertEqual(a.convert_seafood(whole, "kg", assay=assay).magnitude, 100)

    def test_reduction_mass_balance(self):
        assay = a.SeafoodAssay.anchoveta_reduction()
        for magnitude in (100.0, Decimal("100")):
            whole = ureg.Quantity(magnitude, "mt")
            yields = assay.whole_to_reduction_yields(whole)
            self.assertEqual(yields["meal"].magnitude, 22.5)
            self.assertEqual(yields["oil"].magnitude, 4.5)
            self.assertLessEqual(sum(q.magnitude for q in yields.values()), magnitude)
            self.assertAlmostEqual(float(a.convert_seafood(whole, "kg", assay=assay, mode="whole_to_meal").magnitude), 22500)
            self.assertAlmostEqual(float(a.convert_seafood(whole, "kg", assay=assay, mode="whole_to_oil").magnitude), 4500)
            if isinstance(magnitude, Decimal):
                self.assertIsInstance(yields["oil"].magnitude, Decimal)
        with self.assertRaises(ValueError):
            replace(assay, dressed_yield_pct=90, oil_yield_pct=20)
        with self.assertRaises(ValueError):
            a.SeafoodAssay.atlantic_salmon().whole_to_reduction_yields(ureg.Quantity(1, "mt"))

    def test_piece_weight_and_price_parity(self):
        shrimp = a.SeafoodAssay.vannamei_shrimp(count_per_lb=20)
        self.assertAlmostEqual(shrimp.piece_mass_grams(), 453.59237 / 20, places=4)
        self.assertEqual(shrimp.whole_price_to_edible_price(Decimal("5")), 10)
        self.assertIsInstance(shrimp.whole_price_to_edible_price(Decimal("5")), Decimal)
        self.assertEqual(shrimp.whole_price_to_edible_price(5.0), 10)
        with self.assertRaises(ValueError):
            a.SeafoodAssay.atlantic_salmon().piece_mass_grams()

    def test_bulk_context_rates_quotes_and_missing_density(self):
        assay = a.SeafoodAssay.anchoveta_reduction()
        for magnitude in (2.5, Decimal("2.5")):
            volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_seafood(volume, "kg", assay=assay)
            self.assertEqual(mass.magnitude, 1600)
            self.assertAlmostEqual(float(a.convert_seafood(mass, "m^3", assay=assay).magnitude), 2.5)
            rate = a.convert_seafood(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=assay)
            self.assertAlmostEqual(float(rate.magnitude), 1600 / 24)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/m^3"), "EUR/kg", assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 * 0.9 / 640)
        missing = a.SeafoodAssay.atlantic_salmon()
        self.assertFalse(missing.as_context().funcs)
        with self.assertRaises(ValueError):
            missing.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(pint.DimensionalityError):
            a.convert_seafood(ureg.Quantity(1, "m^3"), "kg", assay=missing)

    def test_invalid_specs_masses_and_modes(self):
        assay = a.SeafoodAssay.atlantic_salmon()
        invalid = {"species_group": ("unknown", None), "dressed_yield_pct": (0, 101, float("nan")), "edible_yield_pct": (0, 101), "oil_yield_pct": (-1, 101), "glaze_pct": (-1, 100, Decimal("NaN")), "count_per_lb": (0, float("inf")), "bulk_density_kg_per_m3": (99, 2001)}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for value in (-1, float("nan"), Decimal("Infinity")):
            for method in (assay.whole_to_dressed, assay.dressed_to_whole, assay.whole_to_edible, assay.edible_to_whole, assay.net_to_glazed_mass, assay.glazed_to_net_mass):
                with self.assertRaises(ValueError):
                    method(ureg.Quantity(value, "kg"))
        with self.assertRaises(ValueError):
            assay.whole_price_to_edible_price(Decimal("NaN"))
        with self.assertRaises(TypeError):
            a.convert_seafood(1, "kg", assay=assay)
        with self.assertRaises(ValueError):
            a.convert_seafood(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()

# Mirrored family contract coverage
import pytest as _contract_pytest
from tests.core.units._family_contracts import cases as _contract_cases
from tests.core.units._family_contracts import assert_snapshot_and_ownership as _snapshot_contract
from tests.core.units._family_contracts import assert_physical_or_energy_contract as _physical_contract

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.seafood'))
def test_seafood_immutable_snapshot_and_foreign_registry_contract(specimen):
    _snapshot_contract(specimen)

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.seafood'))
def test_seafood_independent_density_or_missing_factor_contract(specimen):
    _physical_contract(specimen)
