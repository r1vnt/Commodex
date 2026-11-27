"""Grain ethanol yield, DDGS, density and gross crush economics regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class EthanolFeedstockAssayTests(unittest.TestCase):
    def test_exports_presets_and_bushel_weights(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("EthanolFeedstockAssay", "create_ethanol_feedstock_context", "convert_ethanol_feedstock"):
            self.assertIn(name, a.__all__)
        for assay, bushel_unit in ((a.EthanolFeedstockAssay.us_corn_ethanol(), "bushel_corn"), (a.EthanolFeedstockAssay.wheat_ethanol(), "bushel_wheat")):
            for count in (100.0, Decimal("100.123456789")):
                with self.subTest(feedstock=assay.feedstock, count=count):
                    mass = assay.bushels_to_feedstock_mass(count)
                    expected = ureg.Quantity(count, bushel_unit).to("mt")
                    self.assertAlmostEqual(float(mass.magnitude), float(expected.magnitude))
                    back = assay.feedstock_mass_to_bushels(mass)
                    self.assertAlmostEqual(float(back), float(count))
                    if isinstance(count, Decimal):
                        self.assertIsInstance(mass.magnitude, Decimal)
                        self.assertIsInstance(back, Decimal)

    def test_ethanol_and_ddgs_reference_yields_and_inverse(self):
        for assay in (a.EthanolFeedstockAssay.us_corn_ethanol(), a.EthanolFeedstockAssay.wheat_ethanol()):
            for magnitude in (1.0, Decimal("1")):
                feed = ureg.Quantity(magnitude * (Decimal(str(assay.bushel_weight_lb)) if isinstance(magnitude, Decimal) else assay.bushel_weight_lb), "lb")
                ethanol = assay.feedstock_mass_to_ethanol_volume(feed)
                ddgs = assay.feedstock_mass_to_ddgs_mass(feed)
                self.assertAlmostEqual(float(ethanol.magnitude), assay.ethanol_gal_per_bushel)
                self.assertAlmostEqual(float(ddgs.to("lb").magnitude), assay.ddgs_lb_per_bushel)
                for output, reverse in ((ethanol, assay.ethanol_volume_to_feedstock_mass), (ddgs, assay.ddgs_mass_to_feedstock_mass)):
                    self.assertAlmostEqual(float(reverse(output).to("lb").magnitude), assay.bushel_weight_lb)
                    if isinstance(magnitude, Decimal):
                        self.assertIsInstance(output.magnitude, Decimal)
                one_tonne = ureg.Quantity(magnitude, "mt")
                litres = assay.feedstock_mass_to_ethanol_volume(one_tonne).to("liter").magnitude
                self.assertAlmostEqual(float(litres), assay.ethanol_litres_per_mt_feedstock, places=4)
        no_ddgs = replace(a.EthanolFeedstockAssay.us_corn_ethanol(), ddgs_lb_per_bushel=0)
        self.assertEqual(no_ddgs.feedstock_mass_to_ddgs_mass(ureg.Quantity(1, "mt")).magnitude, 0)
        with self.assertRaises(ValueError):
            no_ddgs.ddgs_mass_to_feedstock_mass(ureg.Quantity(1, "mt"))

    def test_ethanol_density_is_separate_from_feedstock_yield(self):
        assay = a.EthanolFeedstockAssay.us_corn_ethanol()
        for magnitude in (1.0, Decimal("1")):
            ethanol = ureg.Quantity(magnitude, "m^3")
            mass = assay.ethanol_volume_to_mass(ethanol)
            self.assertAlmostEqual(float(mass.to("kg").magnitude), 789)
            self.assertAlmostEqual(float(assay.ethanol_mass_to_volume(mass).magnitude), 1)
            if isinstance(magnitude, Decimal):
                self.assertIsInstance(mass.magnitude, Decimal)
            feed = a.convert_ethanol_feedstock(ethanol, "mt", assay=assay)
            self.assertNotAlmostEqual(float(feed.magnitude), float(mass.magnitude))
            direct = a.convert_ethanol_feedstock(ethanol, "kg", assay=assay, mode="ethanol_vol_to_mass")
            self.assertAlmostEqual(float(direct.magnitude), 789)

    def test_crush_margin_and_breakeven_prices(self):
        assay = a.EthanolFeedstockAssay.us_corn_ethanol()
        for grain, ethanol, ddgs in ((4.5, 2.0, 150.0), (Decimal("4.5"), Decimal("2"), Decimal("150"))):
            margin = assay.gross_crush_margin_usd_per_bushel(grain, ethanol, ddgs)
            self.assertAlmostEqual(float(margin), 2.475)
            per_tonne = assay.gross_crush_margin_usd_per_mt(grain, ethanol, ddgs)
            bushels = ureg.Quantity(1, "mt").to("lb").magnitude / 56
            self.assertAlmostEqual(float(per_tonne), 2.475 * bushels, places=4)
            self.assertAlmostEqual(float(assay.breakeven_grain_usd_per_bushel(ethanol, ddgs)), 6.975)
            self.assertAlmostEqual(float(assay.breakeven_ethanol_usd_per_gallon(grain, ddgs)), 3.225 / 2.85, places=12)
            if isinstance(grain, Decimal):
                self.assertIsInstance(margin, Decimal)
                self.assertIsInstance(per_tonne, Decimal)
        self.assertLess(assay.gross_crush_margin_usd_per_bushel(10, 1, 0), 0)

    def test_modes_context_rates_and_quote_router(self):
        assay = a.EthanolFeedstockAssay.us_corn_ethanol()
        for magnitude in (100.0, Decimal("100")):
            feed = ureg.Quantity(magnitude, "mt")
            for forward, reverse, target in (("feedstock_to_ethanol", "ethanol_to_feedstock", "gallon"), ("feedstock_to_ddgs", "ddgs_to_feedstock", "mt")):
                changed = a.convert_ethanol_feedstock(feed, target, assay=assay, mode=forward)
                restored = a.convert_ethanol_feedstock(changed, "mt", assay=assay, mode=reverse)
                self.assertAlmostEqual(float(restored.magnitude), 100)
            rate = a.convert_ethanol_feedstock(ureg.Quantity(magnitude, "mt/day"), "gallon/hour", assay=assay)
            back = a.convert_ethanol_feedstock(rate, "mt/day", assay=assay)
            self.assertAlmostEqual(float(back.magnitude), 100)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/mt"), "EUR/gallon", assay=assay, rates={"EUR": Decimal("0.9")})
            inverse = c.convert_commodity_quote(quote, "USD/mt", assay=assay, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(inverse.magnitude), 100)
        self.assertEqual(a.convert_ethanol_feedstock(ureg.Quantity(100, "kg"), "kg", assay=assay).magnitude, 100)

    def test_validation(self):
        assay = a.EthanolFeedstockAssay.us_corn_ethanol()
        invalid = {"feedstock": ("", " ", None), "bushel_weight_lb": (0, float("nan")), "ethanol_gal_per_bushel": (0, -1), "ddgs_lb_per_bushel": (-1, 57), "ethanol_density_kg_per_m3": (649, 951, float("inf"))}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.bushels_to_feedstock_mass(bad)
            with self.assertRaises(ValueError):
                assay.feedstock_mass_to_ethanol_volume(ureg.Quantity(bad, "mt"))
            with self.assertRaises(ValueError):
                assay.ethanol_volume_to_mass(ureg.Quantity(bad, "m^3"))
        with self.assertRaises(ValueError):
            assay.gross_crush_margin_usd_per_bushel(Decimal("NaN"), 2, 150)
        with self.assertRaises(ValueError):
            assay.breakeven_ethanol_usd_per_gallon(4.5, Decimal("NaN"))
        with self.assertRaises(pint.DimensionalityError):
            assay.feedstock_mass_to_ethanol_volume(ureg.Quantity(1, "m^3"))
        with self.assertRaises(ValueError):
            a.convert_ethanol_feedstock(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()

# Mirrored family contract coverage
import pytest as _contract_pytest
from tests.core.units._family_contracts import cases as _contract_cases
from tests.core.units._family_contracts import assert_snapshot_and_ownership as _snapshot_contract
from tests.core.units._family_contracts import assert_physical_or_energy_contract as _physical_contract

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.feedstocks'))
def test_feedstocks_immutable_snapshot_and_foreign_registry_contract(specimen):
    _snapshot_contract(specimen)

@_contract_pytest.mark.parametrize("specimen", _contract_cases('core.units.assays.agriculture.feedstocks'))
def test_feedstocks_independent_density_or_missing_factor_contract(specimen):
    _physical_contract(specimen)
