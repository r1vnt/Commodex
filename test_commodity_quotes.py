"""Combined FX and physical-denominator quote conversion regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import conversions as c
from core.units.registry import ureg


class CommodityQuoteTests(unittest.TestCase):
    def test_crude_fx_and_round_trip(self):
        assay = c.CrudeAssay.brent()
        quote = ureg.Quantity(Decimal("75.50"), "USD/oil_barrel")
        result = c.convert_commodity_quote(quote, "EUR/metric_ton", crude_assay=assay, rates={"EUR": Decimal("0.9")})
        expected = Decimal("75.5") * Decimal("0.9") / (Decimal(str(assay.density_kg_per_m3)) * Decimal("0.158987294928") / 1000)
        self.assertAlmostEqual(float(result.magnitude), float(expected))
        self.assertIsInstance(result.magnitude, Decimal)
        restored = c.convert_commodity_quote(result, quote.units, assay=assay, rates={"EUR": Decimal("0.9")})
        self.assertAlmostEqual(float(restored.magnitude), 75.5)

    def test_minor_currencies_energy_scaling_and_settlement(self):
        source = ureg.Quantity(Decimal("85"), "GBP_pence/therm_it")
        result = c.convert_commodity_quote(source, "USD/MMBtu_it", rates={"GBP": Decimal("0.8")})
        self.assertAlmostEqual(float(result.magnitude), 10.625)
        self.assertEqual(result.units, ureg.USD / ureg.MMBtu_it)
        for target, expected in (("GBP/therm_it", Decimal("1.24")), ("GBP_pence/therm_it", Decimal("124"))):
            result = c.convert_commodity_quote(ureg.Quantity(1.245, "GBP/therm_it"), target, quantize_settlement=True)
            self.assertEqual(result.magnitude, expected)
        result = c.convert_commodity_quote(ureg.Quantity(0.07, "USD/MWh"), "EUR/MWh", rates={"EUR": 7.5}, quantize_settlement=True)
        self.assertEqual(result.magnitude, Decimal("0.52"))

    def test_broader_assay_physical_legs(self):
        cases = (
            (c.BiofuelAssay.ethanol(), "MJ", lambda a, q: a.volume_to_energy(q)),
            (c.NGLAssay.propane(), "MJ", lambda a, q: a.volume_to_energy(q)),
            (c.BiomassAssay.industrial_pellets_i2(), "MJ", lambda a, q: a.bulk_volume_to_energy(q)),
            (c.PetcokeAssay.usgc_fuel_grade_65s(), "MJ", lambda a, q: a.mass_to_energy(a.bulk_volume_to_wmt(q))),
            (c.BitumenAssay.pen_60_70(), "kg", lambda a, q: a.volume_to_mass(q)),
            (c.HydrogenAssay.liquid_h2(), "MJ", lambda a, q: a.volume_to_energy(q)),
            (c.PetrochemAssay.benzene(), "mol", lambda a, q: a.mass_to_moles(a.volume_to_mass(q))),
            (c.ChemicalAssay.caustic_soda_50pct(), "mol", lambda a, q: a.mass_to_moles(a.volume_to_mass(q))),
        )
        for assay, target_den, physical in cases:
            with self.subTest(assay=assay.name):
                source = ureg.Quantity(Decimal("10"), "USD/m^3")
                result = c.convert_commodity_quote(source, "EUR/" + target_den, assay=assay, rates={"EUR": Decimal("0.9")})
                factor = physical(assay, ureg.Quantity(Decimal("1"), "m^3")).to(target_den).magnitude
                self.assertAlmostEqual(float(result.magnitude), float(Decimal("9") / factor))
                restored = c.convert_commodity_quote(result, source.units, assay=assay, rates={"EUR": Decimal("0.9")})
                self.assertAlmostEqual(float(restored.magnitude), 10)

    def test_chemical_and_petrochemical_active_options(self):
        chemical = c.ChemicalAssay.caustic_soda_50pct()
        source = ureg.Quantity(Decimal("10"), "USD/kg")
        active = c.convert_commodity_quote(source, "USD/kmol", assay=chemical)
        pure = c.convert_commodity_quote(source, "USD/kmol", assay=chemical, active_only=False)
        self.assertAlmostEqual(float(active.magnitude), 799.94)
        self.assertAlmostEqual(float(pure.magnitude), 399.97)
        petrochem = replace(c.PetrochemAssay.benzene(), purity_pct=90)
        full = c.convert_commodity_quote(source, "USD/mol", assay=petrochem)
        active = c.convert_commodity_quote(source, "USD/mol", assay=petrochem, active_only=True)
        self.assertAlmostEqual(float(active.magnitude), float(full.magnitude) / 0.9)

    def test_gas_context_and_energy_basis_selection(self):
        gas = c.GasAssay(gcv_btu_per_scf=1000, name="test_gas")
        result = c.convert_commodity_quote(ureg.Quantity(Decimal("1"), "USD/scf_60F_14_73psia"), "USD/MMBtu_it", gas_assay=gas)
        self.assertAlmostEqual(float(result.magnitude), 1000)
        assay = c.BiofuelAssay.ethanol()
        source = ureg.Quantity(1, "USD/kg")
        lhv = c.convert_commodity_quote(source, "USD/MJ", assay=assay, energy_basis="lhv")
        hhv = c.convert_commodity_quote(source, "USD/MJ", assay=assay, energy_basis="hhv")
        self.assertGreater(lhv.magnitude, hhv.magnitude)

    def test_context_composition_and_zero_negative_quotes(self):
        class MolarAssay:
            def as_context(self):
                ctx = pint.Context("test_molar_bridge")
                ctx.add_transformation("[mass]", "[substance]", lambda _u, q: ureg.Quantity(q.to("kg").magnitude * 2, "mol"))
                ctx.add_transformation("[substance]", "[mass]", lambda _u, q: ureg.Quantity(q.to("mol").magnitude / 2, "kg"))
                return ctx

        crude = c.CrudeAssay.brent()
        for price in (Decimal("0"), Decimal("-10")):
            result = c.convert_commodity_quote(ureg.Quantity(price, "USD/m^3"), "USD/mol", assay=MolarAssay(), crude_assay=crude)
            self.assertAlmostEqual(float(result.magnitude), float(price) / (2 * crude.density_kg_per_m3))
        with self.assertRaises(pint.DimensionalityError):
            ureg.Quantity(1, "m^3").to("mol")

    def test_clear_errors(self):
        with self.assertRaisesRegex(ValueError, "FX rates are required"):
            c.convert_commodity_quote(ureg.Quantity(1, "USD/MWh"), "EUR/MWh")
        with self.assertRaisesRegex(ValueError, "requires an assay"):
            c.convert_commodity_quote(ureg.Quantity(1, "USD/m^3"), "USD/kg")
        with self.assertRaises(ValueError):
            c.convert_commodity_quote(ureg.Quantity(1, "kg"), "USD/kg")
        with self.assertRaises(TypeError):
            c.convert_commodity_quote(1, "USD/kg")
        with self.assertRaises(TypeError):
            c.convert_commodity_quote(ureg.Quantity(1, "USD/kg"), "USD/kg", assay=object())
        for amount in (float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                c.convert_commodity_quote(ureg.Quantity(amount, "USD/kg"), "USD/kg")


if __name__ == "__main__":
    unittest.main()
