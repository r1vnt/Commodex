"""Master router regressions across pressure, currencies and physical assays."""

from decimal import Decimal
import unittest

from core.units import conversions as c
from core.units.registry import ureg


class MasterConversionTests(unittest.TestCase):
    def test_pressure_alias_neutral_and_same_datum_validation(self):
        zero = ureg.Quantity(Decimal("0"), "lb_per_sq_in_gauge")
        result = c.convert(zero, "Pa", p_atm=Decimal("90000.123"), atmospheric_pressure=101325)
        self.assertEqual(result.magnitude, Decimal("90000.123"))
        self.assertAlmostEqual(c.convert(ureg.Quantity(1, "bar"), "barg").magnitude, -0.01325)
        with self.assertRaises(ValueError):
            c.convert(ureg.Quantity(-1, "psia"), "bara")

    def test_sparse_and_anchor_fx_with_physical_quote_conversion(self):
        assay = c.CrudeAssay.brent()
        source = ureg.Quantity(Decimal("75.5"), "USD/oil_barrel")
        expected = c.convert_commodity_quote(source, "EUR/mt", assay=assay, rates={"EUR": Decimal("0.9")})
        maps = ({"EUR": Decimal("0.9")}, {"USD/EUR": Decimal("0.9")}, {("USD", "GBP"): Decimal("0.75"), ("GBP", "EUR"): Decimal("1.2")})
        for rates in maps:
            with self.subTest(rates=rates):
                result = c.convert(source, "EUR/mt", assay=assay, rates=rates)
                self.assertEqual(result.units, expected.units)
                self.assertAlmostEqual(float(result.magnitude), float(expected.magnitude))
                self.assertIsInstance(result.magnitude, Decimal)
        raw = c.convert(source, "EUR/mt", api_gravity=assay.api_gravity, rates={"USD/EUR": 0.9})
        self.assertAlmostEqual(float(raw.magnitude), float(expected.magnitude), places=3)

    def test_fx_currency_targets_overrides_and_settlement(self):
        source = ureg.Quantity(Decimal("100"), "USD/oil_barrel")
        result = c.convert(source, "GBP_pence", rates={"USD/GBP": Decimal("0.8")})
        self.assertEqual(result.units, ureg.GBP_pence / ureg.oil_barrel)
        self.assertEqual(result.magnitude, 8000)
        result = c.convert(source, "EUR/oil_barrel", target_currency="GBP_pence", rates={"GBP": Decimal("0.8")})
        self.assertEqual(result.magnitude, 8000)
        self.assertEqual(result.units, ureg.GBP_pence / ureg.oil_barrel)
        result = c.convert(ureg.Quantity(1.245, "GBP/therm_it"), "GBP_pence/therm_it", quantize_settlement=True)
        self.assertEqual(result.magnitude, 124)
        result = c.convert(ureg.Quantity(0.07, "USD"), "EUR", rates={"USD/EUR": 7.5}, quantize_settlement=True)
        self.assertEqual(result.magnitude, Decimal("0.52"))

    def test_broader_assays_use_contexts_and_sensible_defaults(self):
        for assay in (c.BiofuelAssay.ethanol(), c.NGLAssay.propane(), c.BiomassAssay.industrial_pellets_i2(), c.PetcokeAssay.usgc_fuel_grade_65s(), c.HydrogenAssay.liquid_h2()):
            with self.subTest(assay=assay.name):
                mass = ureg.Quantity(Decimal("1"), "kg")
                result = c.convert(mass, "MJ", assay=assay)
                expected = mass.to("MJ", assay.as_context())
                self.assertAlmostEqual(float(result.magnitude), float(expected.magnitude))
        chemical = c.ChemicalAssay.caustic_soda_50pct()
        active = c.convert(ureg.Quantity(1, "kg"), "mol", assay=chemical)
        pure = c.convert(ureg.Quantity(1, "kg"), "mol", assay=chemical, active_only=False)
        self.assertAlmostEqual(pure.magnitude, 2 * active.magnitude)
        gas = c.convert(ureg.Quantity(1, "scf_60F_14_73psia"), "btu_it", btu_per_scf=1000)
        self.assertAlmostEqual(gas.magnitude, 1000)

    def test_ore_price_and_mass_routes_include_fx(self):
        result = c.convert(ureg.Quantity(Decimal("110"), "USD/dmt"), "EUR/wmt", moisture_fraction=Decimal("0.08"), rates={"USD/EUR": Decimal("0.9")})
        self.assertEqual(result.magnitude, Decimal("91.08"))
        self.assertEqual(result.units, ureg.EUR / ureg.wmt)
        ore = c.OreAssay.copper_concentrate()
        result = c.convert(ureg.Quantity(100, "wmt"), "dmtu", assay=ore)
        self.assertAlmostEqual(result.magnitude, 100 * ore.dry_mass_fraction * ore.grade_pct)

    def test_observed_temperatures_for_amounts_and_quotes(self):
        for assay, temp in ((c.PetrochemAssay.benzene(), 40), (c.BitumenAssay.pen_60_70(), 150)):
            with self.subTest(assay=assay.name):
                volume = ureg.Quantity(Decimal("10"), "m^3")
                mass = c.convert(volume, "kg", assay=assay, observed_temp=temp)
                expected = assay.observed_volume_to_mass(volume, temp).to("kg")
                self.assertAlmostEqual(float(mass.magnitude), float(expected.magnitude))
                restored = c.convert(mass, "m^3", assay=assay, observed_temp=temp)
                self.assertAlmostEqual(float(restored.magnitude), 10)
                price = ureg.Quantity(Decimal("10"), "USD/m^3")
                result = c.convert(price, "EUR/kg", assay=assay, observed_temp=temp, rates={"USD/EUR": Decimal("0.9")})
                self.assertAlmostEqual(float(result.magnitude), 90 / float(mass.magnitude))
        brent = c.CrudeAssay.brent()
        volume = ureg.Quantity(10, "oil_barrel")
        result = c.convert(volume, "mt", assay=brent, observed_temp=100)
        standard = brent.correct_volume_to_60f(volume, 100)
        self.assertAlmostEqual(result.magnitude, c.convert_hydrocarbon(standard, "mt", assay=brent).magnitude)

    def test_plain_scaling_and_clear_failures(self):
        self.assertEqual(c.convert(ureg.Quantity(1, "meter"), "cm").magnitude, 100)
        with self.assertRaisesRegex(ValueError, "FX rates are required"):
            c.convert(ureg.Quantity(1, "USD/kg"), "EUR/kg")
        with self.assertRaises(ValueError):
            c.convert(ureg.Quantity(1, "USD"), "EUR", source_currency="GBP", rates={"USD/EUR": 0.9})
        with self.assertRaisesRegex(ValueError, "volume-correction"):
            c.convert(ureg.Quantity(1, "kg"), "MJ", assay=c.HydrogenAssay.liquid_h2(), observed_temp=20)
        with self.assertRaises(TypeError):
            c.convert(ureg.Quantity(1, "kg"), "kg", assay=object())
        with self.assertRaises(TypeError):
            c.convert(1, "kg")


if __name__ == "__main__":
    unittest.main()
