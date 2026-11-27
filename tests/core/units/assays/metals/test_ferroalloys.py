"""Ferroalloy composition, settlement, stowage and dosing regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import conversions as c
from core.units.registry import ureg


class FerroalloyAssayTests(unittest.TestCase):
    def test_presets_exports_and_stowage_contexts(self):
        self.assertIn("FerroalloyAssay", c.__all__)
        self.assertIn("create_ferroalloy_context", c.__all__)
        presets = ("hc_ferrochrome", "lc_ferrochrome", "hc_ferromanganese", "silicomanganese", "ferrosilicon_75", "ferromolybdenum_60", "ferrovanadium_80", "ferroniobium_65")
        for name in presets:
            assay = getattr(c.FerroalloyAssay, name)()
            context = assay.as_context(context_name="test_ferroalloy")
            self.assertEqual(context.name, "test_ferroalloy")
            for magnitude in (2.5, Decimal("2.5")):
                with self.subTest(preset=name, magnitude=magnitude):
                    mass = ureg.Quantity(magnitude, "metric_ton")
                    contained = assay.alloy_to_contained_metal(mass)
                    self.assertAlmostEqual(float(contained.magnitude), 2.5 * assay.primary_element_pct / 100)
                    self.assertAlmostEqual(float(assay.contained_metal_to_alloy(contained).magnitude), 2.5)
                    self.assertEqual(mass.to("kg", context).magnitude, 2500)
                    if assay.bulk_density_kg_per_m3 is None:
                        self.assertFalse(context.funcs)
                        continue
                    volume = mass.to("m^3", context)
                    self.assertAlmostEqual(float(volume.magnitude), 2500 / assay.bulk_density_kg_per_m3)
                    self.assertAlmostEqual(float(volume.to("metric_ton", context).magnitude), 2.5)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(volume.magnitude, Decimal)
                        self.assertIsInstance(contained.magnitude, Decimal)
                    with self.assertRaises(pint.DimensionalityError):
                        mass.to("m^3")

    def test_secondary_elements_immutability_and_hash(self):
        source = {" Si ": 16.0}
        assay = c.FerroalloyAssay("SiMn", " Mn ", 65, carbon_pct=2, secondary_elements=source)
        source[" Si "] = 99
        self.assertEqual(assay.primary_element, "Mn")
        self.assertEqual(assay.secondary_elements["Si"], 16)
        with self.assertRaises(TypeError):
            assay.secondary_elements["Si"] = 20
        same = c.FerroalloyAssay("SiMn", "Mn", 65, carbon_pct=2, secondary_elements={"Si": 16})
        self.assertEqual(assay, same)
        self.assertEqual(hash(assay), hash(same))
        self.assertEqual(assay.implied_iron_balance_pct, 17)
        mass = ureg.Quantity(Decimal("10"), "metric_ton")
        silicon = assay.alloy_to_contained_metal(mass, element="Si")
        self.assertEqual(silicon.magnitude, Decimal("1.6"))
        self.assertEqual(assay.contained_metal_to_alloy(silicon, element="Si"), mass)
        with self.assertRaises(KeyError):
            assay.alloy_to_contained_metal(mass, element="Mo")
        zero = replace(assay, secondary_elements={"Si": 0})
        with self.assertRaises(ValueError):
            zero.contained_metal_to_alloy(mass, element="Si")

    def test_price_parity(self):
        assay = c.FerroalloyAssay.ferrovanadium_80()
        self.assertEqual(assay.gross_mt_price_to_contained_kg_price(24000), 30)
        self.assertEqual(assay.contained_kg_price_to_gross_mt_price(30), 24000)
        self.assertAlmostEqual(assay.gross_mt_price_to_contained_lb_price(24000), 30 * 0.45359237, places=4)
        self.assertAlmostEqual(assay.contained_lb_price_to_gross_mt_price(10), 10 * 1000 / 0.45359237 * 0.8, places=4)
        self.assertAlmostEqual(assay.contained_lb_price_to_contained_kg_price(10), 10 / 0.45359237, places=4)
        self.assertAlmostEqual(assay.contained_kg_price_to_contained_lb_price(30), 30 * 0.45359237, places=4)

    def test_pro_rata_settlement(self):
        assay = replace(c.FerroalloyAssay.ferrosilicon_75(), primary_element_pct=72)
        self.assertEqual(assay.pro_rata_price(1500), 1440)
        self.assertEqual(assay.pro_rata_price(1500, basis_grade_pct=80), 1350)
        with self.assertRaises(ValueError):
            replace(assay, basis_grade_pct=None).pro_rata_price(1500)
        for grade in (0, -1, 101, float("nan")):
            with self.assertRaises(ValueError):
                assay.pro_rata_price(1500, basis_grade_pct=grade)

    def test_ladle_dosing(self):
        assay = c.FerroalloyAssay.silicomanganese()
        steel = ureg.Quantity(Decimal("100"), "metric_ton")
        manganese_dose = assay.ladle_addition_required(steel, 0.5, recovery_rate=0.85)
        self.assertIsInstance(manganese_dose.magnitude, Decimal)
        self.assertAlmostEqual(float(manganese_dose.magnitude), 500 / (0.65 * 0.85))
        silicon_dose = assay.ladle_addition_required(steel, 0.1, recovery_rate=0.9, element="Si")
        self.assertAlmostEqual(float(silicon_dose.magnitude), 100 / (0.16 * 0.9))
        for pickup in (0, -1, 101, float("nan")):
            with self.assertRaises(ValueError):
                assay.ladle_addition_required(steel, pickup)
        for recovery in (0, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                assay.ladle_addition_required(steel, 0.5, recovery_rate=recovery)

    def test_blending_preserves_mass_balance(self):
        alloy = c.FerroalloyAssay.silicomanganese()
        other = replace(alloy, primary_element_pct=70, carbon_pct=1, secondary_elements={"Si": 12, "P": 0.1}, bulk_density_kg_per_m3=3500)
        blend = alloy.blend_with(other, self_mass_fraction=0.4)
        self.assertEqual(blend.primary_element_pct, 68)
        self.assertAlmostEqual(blend.secondary_elements["Si"], 13.6)
        self.assertAlmostEqual(blend.secondary_elements["P"], 0.06)
        self.assertAlmostEqual(blend.carbon_pct, 1.4)
        self.assertAlmostEqual(blend.bulk_density_kg_per_m3, 1 / (0.4 / 3300 + 0.6 / 3500), places=2)
        # Intermediate rounding used to make this valid composition exceed 100%.
        thirds = c.FerroalloyAssay("rounding edge", "Mn", 33.33336, carbon_pct=0, secondary_elements={"Si": 33.33336, "Cr": 33.33328})
        self.assertAlmostEqual(thirds.blend_with(thirds, self_mass_fraction=0.5).implied_iron_balance_pct, 0)
        with self.assertRaises(ValueError):
            alloy.blend_with(c.FerroalloyAssay.ferrosilicon_75(), self_mass_fraction=0.5)
        with self.assertRaises(ValueError):
            alloy.blend_with(other, self_mass_fraction=1.1)

    def test_missing_density_and_invalid_specs(self):
        assay = c.FerroalloyAssay.ferroniobium_65()
        self.assertIsNone(assay.bulk_density_kg_per_m3)
        self.assertIsNone(assay.carbon_pct)
        self.assertIsNone(assay.stowage_factor_m3_per_mt)
        with self.assertRaises(ValueError):
            assay.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(ValueError):
            assay.mass_to_bulk_volume(ureg.Quantity(1, "kg"))
        invalid = {"primary_element": ("", " ", None), "primary_element_pct": (0, 101, float("nan")), "basis_grade_pct": (0, 101, float("inf")), "carbon_pct": (-1, 16, float("nan")), "bulk_density_kg_per_m3": (499, 7001, float("inf")), "secondary_elements": ({"Nb": 1}, {" Nb ": 1}, {"Si": -1}, {"Si": 40}, {"": 1}, {"Si": 1, " Si ": 1})}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        with self.assertRaises(ValueError):
            replace(assay, carbon_pct=1, secondary_elements={"C": 1})
        with self.assertRaises(TypeError):
            assay.alloy_to_contained_metal(1)
        with self.assertRaises(pint.DimensionalityError):
            assay.alloy_to_contained_metal(ureg.Quantity(1, "m^3"))
        for method in (assay.gross_mt_price_to_contained_kg_price, assay.contained_kg_price_to_gross_mt_price, assay.gross_mt_price_to_contained_lb_price, assay.contained_lb_price_to_gross_mt_price, assay.contained_lb_price_to_contained_kg_price, assay.contained_kg_price_to_contained_lb_price, assay.pro_rata_price):
            with self.assertRaises(ValueError):
                method(float("nan"))


if __name__ == "__main__":
    unittest.main()

# Mirrored family contract coverage
import pytest as _contract_pytest
from tests.core.units._family_contracts import cases as _contract_cases
from tests.core.units._family_contracts import assert_snapshot_and_ownership as _snapshot_contract
from tests.core.units._family_contracts import assert_physical_or_energy_contract as _physical_contract

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.metals.ferroalloys'))
def test_ferroalloys_immutable_snapshot_and_foreign_registry_contract(specimen):
    _snapshot_contract(specimen)

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.metals.ferroalloys'))
def test_ferroalloys_independent_density_or_missing_factor_contract(specimen):
    _physical_contract(specimen)
