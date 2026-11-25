"""Pulp fibre-basis and lumber nominal-volume registry integration tests."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units import agricultural_specs as a
from core.units import conversions as c
from core.units.registry import ureg


class PulpLumberAssayTests(unittest.TestCase):
    def test_exports_and_presets(self):
        self.assertEqual(len(a.__all__), len(set(a.__all__)))
        self.assertFalse(set(a.__all__) - vars(a).keys())
        for name in ("PulpLumberAssay", "create_pulp_lumber_context", "convert_pulp_lumber"):
            self.assertIn(name, a.__all__)
        for preset in ("nbsk_softwood_pulp", "bek_hardwood_pulp", "spf_dimensional_lumber"):
            assay = getattr(a.PulpLumberAssay, preset)()
            self.assertTrue(assay.as_context().funcs)

    def test_exact_board_foot_and_mbf_registry_scales(self):
        self.assertAlmostEqual(ureg.Quantity(1, "board_foot").to("inch^3").magnitude, 144)
        self.assertAlmostEqual(ureg.Quantity(1, "MBF").to("foot^3").magnitude, 1000 / 12)
        self.assertAlmostEqual(ureg.Quantity(1, "MBF").to("m^3").magnitude, 2.359737216)
        self.assertAlmostEqual(ureg.Quantity(1, "MMBF").to("MBF").magnitude, 1000)
        lumber = a.PulpLumberAssay.spf_dimensional_lumber()
        for count in (100.0, Decimal("100.123456789")):
            for dressed in (False, True):
                with self.subTest(count=count, dressed=dressed):
                    volume = lumber.mbf_to_cubic_metres(count, apply_dressed_recovery=dressed)
                    recovered = lumber.cubic_metres_to_mbf(volume, apply_dressed_recovery=dressed)
                    self.assertAlmostEqual(float(recovered), float(count))
                    if isinstance(count, Decimal):
                        self.assertIsInstance(volume.magnitude, Decimal)
                        self.assertIsInstance(recovered, Decimal)
        self.assertEqual(lumber.mbf_to_cubic_metres(Decimal("1")).magnitude, Decimal("2.359737216"))
        self.assertEqual(lumber.mbf_to_cubic_metres(Decimal("1"), apply_dressed_recovery=True).magnitude, Decimal("2.359737216") * Decimal("0.92"))

    def test_pulp_moisture_bases_and_explicit_modes(self):
        for moisture in (0, 10, 55):
            assay = replace(a.PulpLumberAssay.nbsk_softwood_pulp(), moisture_pct=moisture)
            for magnitude in (100.0, Decimal("100")):
                with self.subTest(moisture=moisture, magnitude=magnitude):
                    actual = ureg.Quantity(magnitude, "mt")
                    dry = assay.mass_to_bdmt(actual)
                    air = assay.mass_to_admt(actual)
                    self.assertEqual(dry.units, ureg.bdmt)
                    self.assertEqual(air.units, ureg.admt)
                    self.assertAlmostEqual(float(dry.magnitude), 100 - moisture)
                    self.assertAlmostEqual(float(air.magnitude), (100 - moisture) / 0.9)
                    self.assertAlmostEqual(float(assay.bdmt_to_mass(dry).magnitude), 100)
                    self.assertAlmostEqual(float(assay.admt_to_mass(air).magnitude), 100)
                    for forward, reverse, target in (("to_bdmt", "from_bdmt", "BDMT"), ("to_admt", "from_admt", "ADMT")):
                        changed = a.convert_pulp_lumber(actual, target, assay=assay, mode=forward)
                        restored = a.convert_pulp_lumber(changed, "mt", assay=assay, mode=reverse)
                        self.assertAlmostEqual(float(restored.magnitude), 100)
                        if isinstance(magnitude, Decimal):
                            self.assertIsInstance(changed.magnitude, Decimal)
        # Labels alone do not apply the 90% fibre adjustment.
        with self.assertRaises(pint.DimensionalityError):
            ureg.Quantity(1, "ADMT").to("BDMT")

    def test_packages_and_prices(self):
        pulp = a.PulpLumberAssay.nbsk_softwood_pulp()
        for count in (10.0, Decimal("10.123456789")):
            mass = pulp.packages_to_mass(count)
            self.assertAlmostEqual(float(mass.to("kg").magnitude), float(count) * 250)
            self.assertAlmostEqual(float(pulp.mass_to_packages(mass)), float(count))
            if isinstance(count, Decimal):
                self.assertIsInstance(mass.magnitude, Decimal)
                self.assertIsInstance(pulp.mass_to_packages(mass), Decimal)
        self.assertEqual(pulp.admt_price_to_bdmt_price(Decimal("900")), 1000)
        self.assertEqual(pulp.bdmt_price_to_admt_price(Decimal("1000")), 900)
        mbf_price = pulp.m3_price_to_mbf_price(Decimal("100"))
        self.assertEqual(mbf_price, Decimal("235.9737"))
        self.assertAlmostEqual(float(pulp.mbf_price_to_m3_price(mbf_price)), 100, places=4)
        with self.assertRaises(ValueError):
            a.PulpLumberAssay.spf_dimensional_lumber().mass_to_packages(ureg.Quantity(1, "mt"))

    def test_bulk_volume_context_master_and_quotes(self):
        pulp = a.PulpLumberAssay.nbsk_softwood_pulp()
        for magnitude in (2.5, Decimal("2.5")):
            volume = ureg.Quantity(magnitude, "m^3")
            mass = a.convert_pulp_lumber(volume, "kg", assay=pulp)
            self.assertEqual(mass.magnitude, 1875)
            self.assertAlmostEqual(float(a.convert_pulp_lumber(mass, "m^3", assay=pulp).magnitude), 2.5)
            rate = a.convert_pulp_lumber(ureg.Quantity(magnitude, "m^3/day"), "kg/hour", assay=pulp)
            self.assertAlmostEqual(float(rate.magnitude), 1875 / 24)
            quote = c.convert_commodity_quote(ureg.Quantity(magnitude, "USD/m^3"), "EUR/kg", assay=pulp, rates={"EUR": Decimal("0.9")})
            self.assertAlmostEqual(float(quote.magnitude), 2.5 * 0.9 / 750)
        lumber = a.PulpLumberAssay.spf_dimensional_lumber()
        source = ureg.Quantity(100, "USD/MBF")
        result = c.convert_commodity_quote(source, "USD/m^3", assay=lumber)
        self.assertAlmostEqual(result.magnitude, 100 / 2.359737216)

    def test_missing_density_and_invalid_inputs(self):
        assay = a.PulpLumberAssay.nbsk_softwood_pulp()
        invalid = {"product_type": ("invalid", None), "moisture_pct": (-1, 100, Decimal("NaN")), "dressed_recovery_pct": (0, 101), "bulk_density_kg_per_m3": (99, 2001), "unit_package_kg": (0, float("inf"))}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for bad in (-1, float("nan"), Decimal("Infinity")):
            with self.assertRaises(ValueError):
                assay.mbf_to_cubic_metres(bad)
            with self.assertRaises(ValueError):
                assay.cubic_metres_to_mbf(ureg.Quantity(bad, "m^3"))
            with self.assertRaises(ValueError):
                assay.mass_to_bdmt(ureg.Quantity(bad, "mt"))
            with self.assertRaises(ValueError):
                assay.packages_to_mass(bad)
        with self.assertRaises(ValueError):
            assay.admt_price_to_bdmt_price(Decimal("NaN"))
        missing = replace(assay, bulk_density_kg_per_m3=None)
        self.assertFalse(missing.as_context().funcs)
        with self.assertRaises(ValueError):
            missing.bulk_volume_to_mass(ureg.Quantity(1, "m^3"))
        with self.assertRaises(pint.DimensionalityError):
            a.convert_pulp_lumber(ureg.Quantity(1, "m^3"), "kg", assay=missing)
        with self.assertRaises(ValueError):
            a.convert_pulp_lumber(ureg.Quantity(1, "kg"), "kg", assay=assay, mode="invalid")


if __name__ == "__main__":
    unittest.main()
