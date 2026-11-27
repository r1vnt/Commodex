"""Public export integrity and high-level converter entry points."""

from decimal import Decimal
import unittest

from core.units import conversions as c
from core.units.registry import ureg


class ConversionExportTests(unittest.TestCase):
    def test_exports_exist_and_are_unique(self):
        self.assertEqual(len(c.__all__), len(set(c.__all__)))
        missing = [name for name in c.__all__ if name not in c.__dict__]
        self.assertEqual(missing, [])

    def test_star_import_matches_public_exports(self):
        namespace = {}
        exec("from core.units.conversions import *", namespace)
        self.assertEqual(set(namespace) - {"__builtins__"}, set(c.__all__))
        for name in c.__all__:
            self.assertIs(namespace[name], c.__dict__[name])

    def test_domain_converter_entry_points(self):
        cases = (
            ("biofuel", c.BiofuelAssay.ethanol(), "kg", "MJ"),
            ("ngl", c.NGLAssay.propane(), "kg", "MJ"),
            ("ore", c.OreAssay.iron_ore_62fe(), "wmt", "dmt"),
            ("battery_material", c.BatteryMaterialAssay.spodumene_sc6(), "kg", "g"),
            ("uranium", c.UraniumAssay.leu_4_5pct(), "kg", "MWh"),
            ("biomass", c.BiomassAssay.industrial_pellets_i2(), "kg", "MJ"),
            ("petrochem", c.PetrochemAssay.benzene(), "m^3", "kg"),
            ("bullion", c.BullionAssay.gold_dore(), "kg", "g"),
            ("fertilizer", c.FertilizerAssay.uan_32(), "m^3", "kg"),
            ("bitumen", c.BitumenAssay.pen_60_70(), "m^3", "kg"),
            ("petcoke", c.PetcokeAssay.usgc_fuel_grade_65s(), "kg", "MJ"),
            ("chemical", c.ChemicalAssay.caustic_soda_50pct(), "kg", "mol"),
            ("ferroalloy", c.FerroalloyAssay.silicomanganese(), "m^3", "kg"),
            ("metal_alloy", c.MetalAlloyAssay.cartridge_brass_c260(), "m^3", "kg"),
        )
        for name, assay, source, target in cases:
            with self.subTest(converter=name):
                converter = getattr(c, "convert_" + name)
                self.assertIn("convert_" + name, c.__all__)
                quantity = ureg.Quantity(Decimal("10"), source)
                result = converter(quantity, target, assay=assay)
                self.assertGreater(result.magnitude, 0)
                self.assertIsInstance(result.magnitude, Decimal)
                restored = converter(result, source, assay=assay)
                self.assertAlmostEqual(float(restored.magnitude), 10)
                with self.assertRaises(TypeError):
                    converter(quantity, target, assay=object())


if __name__ == "__main__":
    unittest.main()
