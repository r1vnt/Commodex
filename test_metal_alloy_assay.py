"""General alloy storage, composition, remelt and valuation regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import conversions as c
from core.units.registry import ureg


class MetalAlloyAssayTests(unittest.TestCase):
    def test_presets_exports_and_contexts(self):
        self.assertIn("MetalAlloyAssay", c.__all__)
        self.assertIn("create_metal_alloy_context", c.__all__)
        presets = ("cartridge_brass_c260", "aluminium_6061", "aluminium_7075", "monel_400", "inconel_625", "titanium_6al_4v")
        for name in presets:
            assay = getattr(c.MetalAlloyAssay, name)()
            for magnitude in (2.5, Decimal("2.5")):
                with self.subTest(preset=name, magnitude=magnitude):
                    volume = ureg.Quantity(magnitude, "m^3")
                    context = assay.as_context(context_name="test_metal_alloy")
                    self.assertEqual(context.name, "test_metal_alloy")
                    result = volume.to("kg", context)
                    self.assertAlmostEqual(float(result.magnitude), 2.5 * assay.solid_density_kg_per_m3)
                    self.assertAlmostEqual(float(result.to("m^3", context).magnitude), 2.5)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(result.magnitude, Decimal)
                    with self.assertRaises(pint.DimensionalityError):
                        volume.to("kg")

    def test_composition_immutability_and_contained_mass(self):
        composition = {" Cu ": 70, "Zn": 30}
        assay = c.MetalAlloyAssay("brass", composition, 8530)
        composition["Zn"] = 40
        self.assertEqual(assay.composition_pct["Zn"], 30)
        self.assertEqual(assay.base_metal, "Cu")
        self.assertEqual(assay.base_metal_pct, 70)
        self.assertEqual(assay.balance_other_pct, 0)
        same = c.MetalAlloyAssay("brass", {"Zn": 30, "Cu": 70}, 8530)
        self.assertEqual(assay, same)
        self.assertEqual(hash(assay), hash(same))
        with self.assertRaises(TypeError):
            assay.composition_pct["Cu"] = 80
        mass = ureg.Quantity(Decimal("10"), "metric_ton")
        for element, expected in (("Cu", 7), ("Zn", 3)):
            contained = assay.alloy_to_contained_element(mass, element)
            self.assertEqual(contained.magnitude, expected)
            self.assertIsInstance(contained.magnitude, Decimal)
            self.assertEqual(assay.contained_element_to_alloy(contained, element), mass)
        self.assertEqual(mass.to("kg", assay.as_context()).magnitude, 10000)
        self.assertEqual(c.MetalAlloyAssay.aluminium_6061().balance_other_pct, 0.5)
        with self.assertRaises(KeyError):
            assay.element_fraction("Ni")
        with self.assertRaises(ValueError):
            replace(assay, composition_pct={"Cu": 70, "Zn": 0}).element_fraction("Zn")

    def test_remelt_recovery_and_piece_counts(self):
        assay = replace(c.MetalAlloyAssay.aluminium_6061(), remelt_yield_pct=92)
        mass = ureg.Quantity(Decimal("10"), "metric_ton")
        recovered = assay.gross_to_recovered_melt_mass(mass)
        self.assertEqual(recovered.magnitude, Decimal("9.2"))
        self.assertEqual(assay.recovered_melt_to_gross_mass(recovered), mass)
        self.assertIsInstance(recovered.magnitude, Decimal)
        self.assertAlmostEqual(float(assay.piece_count_to_mass(10).magnitude), 226.8)
        self.assertAlmostEqual(assay.mass_to_piece_count(ureg.Quantity(226.8, "kg")), 10)
        self.assertEqual(assay.piece_count_to_mass(0).magnitude, 0)
        for count in (-1, float("nan"), float("inf")):
            with self.assertRaises(ValueError):
                assay.piece_count_to_mass(count)
        no_piece = c.MetalAlloyAssay.monel_400()
        with self.assertRaises(ValueError):
            no_piece.mass_to_piece_count(mass)
        with self.assertRaises(ValueError):
            no_piece.piece_count_to_mass(1)

    def test_melt_valuation_and_price_parity(self):
        assay = replace(c.MetalAlloyAssay.cartridge_brass_c260(), remelt_yield_pct=90)
        prices = {"Cu": 10000, "Zn": 3000}
        self.assertEqual(assay.gross_melt_value_per_mt(prices), 7900)
        self.assertEqual(assay.gross_melt_value_per_mt(prices, apply_remelt_yield=True, treatment_deduction_usd_per_mt=100), 7010)
        self.assertEqual(assay.gross_melt_value_per_mt({"Cu": 10000}), 7000)
        self.assertEqual(assay.gross_melt_value_per_mt({}), 0)
        self.assertAlmostEqual(assay.lb_price_to_mt_price(1), 1000 / 0.45359237, places=4)
        self.assertAlmostEqual(assay.mt_price_to_lb_price(1000), 0.45359237, places=6)
        with self.assertRaises(ValueError):
            assay.gross_melt_value_per_mt({"Cu": float("nan")})
        with self.assertRaises(ValueError):
            assay.gross_melt_value_per_mt(prices, treatment_deduction_usd_per_mt=float("inf"))
        for method in (assay.lb_price_to_mt_price, assay.mt_price_to_lb_price):
            with self.assertRaises(ValueError):
                method(float("nan"))

    def test_blend_mass_balance_and_rounding_edges(self):
        brass = c.MetalAlloyAssay.cartridge_brass_c260()
        other = replace(brass, composition_pct={"Cu": 60, "Zn": 35, "Sn": 5}, solid_density_kg_per_m3=8700, remelt_yield_pct=80)
        blend = brass.blend_with(other, self_mass_fraction=0.4)
        self.assertEqual(dict(blend.composition_pct), {"Cu": 64, "Sn": 3, "Zn": 33})
        self.assertEqual(blend.remelt_yield_pct, 88)
        self.assertAlmostEqual(blend.solid_density_kg_per_m3, 1 / (0.4 / 8530 + 0.6 / 8700), places=2)
        thirds = c.MetalAlloyAssay("rounding edge", {"Cu": 33.33336, "Zn": 33.33336, "Sn": 33.33328}, 8530, remelt_yield_pct=0.00001)
        result = thirds.blend_with(thirds, self_mass_fraction=0.5)
        self.assertAlmostEqual(sum(result.composition_pct.values()), 100)
        self.assertGreater(result.remelt_yield_pct, 0)
        with self.assertRaises(ValueError):
            brass.blend_with(other, self_mass_fraction=1.1)
        with self.assertRaises(TypeError):
            brass.blend_with(None, self_mass_fraction=0.5)

    def test_invalid_specs_types_and_units(self):
        assay = c.MetalAlloyAssay.cartridge_brass_c260()
        invalid = {"composition_pct": ({}, {"Cu": 0}, {"Cu": 101}, {"Cu": float("nan")}, {"Cu": 70, "Zn": 40}, {"": 70}, {None: 70}, {"Cu": 60, " Cu ": 10}), "solid_density_kg_per_m3": (999, 22001, float("nan")), "base_metal": ("Ni", "", 42), "remelt_yield_pct": (0, 101, float("inf")), "nominal_piece_mass_kg": (0, -1, float("nan"))}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        with self.assertRaises(ValueError):
            replace(assay, composition_pct={"Cu": 0, "Zn": 30}, base_metal="Cu")
        with self.assertRaises(TypeError):
            assay.solid_volume_to_mass(1)
        with self.assertRaises(TypeError):
            assay.gross_to_recovered_melt_mass(1)
        with self.assertRaises(pint.DimensionalityError):
            assay.solid_volume_to_mass(ureg.Quantity(1, "kg"))
        with self.assertRaises(pint.DimensionalityError):
            assay.alloy_to_contained_element(ureg.Quantity(1, "m^3"), "Cu")


if __name__ == "__main__":
    unittest.main()
