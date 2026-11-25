"""Hydrogen storage, carrier yield, pricing and CI screening regressions."""

from dataclasses import replace
from decimal import Decimal
import unittest

import pint

from core.units.conversions import HydrogenAssay, create_hydrogen_context
from core.units.registry import ureg


class HydrogenAssayTests(unittest.TestCase):
    def test_all_presets_and_context_paths(self):
        names = ("gaseous_h2_ntp", "compressed_h2_350bar", "compressed_h2_700bar", "liquid_h2", "ammonia_carrier", "lohc_dbt", "green_h2_rfnbo", "blue_h2_atr")
        for name in names:
            assay = getattr(HydrogenAssay, name)()
            for basis in ("lhv", "hhv"):
                context = assay.as_context(basis=basis, context_name="test_hydrogen")
                self.assertEqual(context.name, "test_hydrogen")
                for magnitude in (2.5, Decimal("2.5")):
                    for source in ("m^3", "kg", "MWh", "Nm3_H2"):
                        for target in ("m^3", "kg", "MWh", "Nm3_H2"):
                            if source == target:
                                continue
                            with self.subTest(name=name, basis=basis, source=source, target=target, magnitude=magnitude):
                                quantity = ureg.Quantity(magnitude, source)
                                with self.assertRaises(pint.DimensionalityError):
                                    quantity.to(target)
                                result = quantity.to(target, context)
                                self.assertGreater(result.magnitude, 0)
                                self.assertAlmostEqual(float(result.to(source, context).magnitude), 2.5)
                                if isinstance(magnitude, Decimal):
                                    self.assertIsInstance(result.magnitude, Decimal)

    def test_reference_yields_and_storage_states(self):
        mass = ureg.Quantity(Decimal("1"), "metric_ton")
        ammonia = HydrogenAssay.ammonia_carrier()
        delivered = ammonia.carrier_to_h2_mass(mass)
        self.assertAlmostEqual(float(delivered.magnitude), 150.9175)
        self.assertAlmostEqual(float(ammonia.h2_to_carrier_mass(delivered).magnitude), 1)
        full_recovery = ammonia.carrier_to_h2_mass(mass, cracking_efficiency=1)
        self.assertAlmostEqual(float(full_recovery.magnitude), 177.55)
        lohc = HydrogenAssay.lohc_dbt()
        self.assertAlmostEqual(float(lohc.carrier_to_h2_mass(mass).magnitude), 55.8)
        for preset, expected_density in (("gaseous_h2_ntp", 0.08988), ("compressed_h2_350bar", 23.7), ("compressed_h2_700bar", 39.7), ("liquid_h2", 70.85)):
            assay = getattr(HydrogenAssay, preset)()
            self.assertAlmostEqual(float(ureg.Quantity(1, "m^3").to("kg", assay.as_context()).magnitude), expected_density)
        normal = HydrogenAssay.h2_mass_to_nm3(ureg.Quantity(Decimal("0.08988"), "kg"))
        self.assertEqual(normal.magnitude, 1)
        self.assertEqual(normal.units, ureg.Nm3_H2)
        self.assertAlmostEqual(float(HydrogenAssay.nm3_to_h2_mass(ureg.Quantity(1, "kNm3")).magnitude), 89.88)
        self.assertEqual(HydrogenAssay.nm3_to_h2_mass(ureg.Quantity(1, "m^3")).magnitude, 0.08988)

    def test_energy_electrolysis_and_price_parity(self):
        assay = HydrogenAssay.gaseous_h2_ntp()
        mass = ureg.Quantity(Decimal("1"), "kg")
        for basis, mj in (("lhv", 119.96), ("hhv", 141.88)):
            self.assertAlmostEqual(float(assay.mass_to_energy(mass, basis=basis).to("MJ").magnitude), mj)
            self.assertAlmostEqual(assay.kg_price_to_mwh_price(3, basis=basis), 3 * 3600 / mj, places=4)
            self.assertAlmostEqual(assay.mwh_price_to_kg_price(assay.kg_price_to_mwh_price(3, basis=basis), basis=basis), 3, places=4)
            self.assertAlmostEqual(assay.kg_price_to_mmbtu_price(3, basis=basis), 3 / float(ureg.Quantity(mj, "MJ").to("MMBtu_it").magnitude), places=4)
            self.assertAlmostEqual(assay.mmbtu_price_to_kg_price(assay.kg_price_to_mmbtu_price(3, basis=basis), basis=basis), 3, places=4)
        self.assertAlmostEqual(float(assay.electrolysis_power_required(mass).to("kWh").magnitude), 119.96 / 3.6 / 0.65)
        ammonia = HydrogenAssay.ammonia_carrier()
        price = ammonia.carrier_mt_price_to_h2_kg_price(500, cracking_cost_per_kg_h2=1.2)
        self.assertAlmostEqual(price, 500 / 150.9175 + 1.2, places=4)
        # Public price methods round to four decimals.
        self.assertAlmostEqual(ammonia.h2_kg_price_to_carrier_mt_price(price, cracking_cost_per_kg_h2=1.2), 500, delta=0.01)

    def test_regulatory_boundaries_and_explicit_year(self):
        assay = HydrogenAssay.gaseous_h2_ntp()
        self.assertIsNone(assay.us_45v_credit_usd_per_kg_h2)
        self.assertIsNone(assay.us_45v_credit_for_year(2026))
        self.assertIsNone(assay.is_eu_rfnbo_compliant)
        cases = ((0, 3, 0.656), (0.449999, 3, 0.656), (0.45, 1.002, 0.219), (1.499999, 1.002, 0.219), (1.5, 0.75, 0.164), (2.499999, 0.75, 0.164), (2.5, 0.6, 0.131), (4, 0.6, 0.131), (4.000001, 0, 0))
        for ci, nominal, current_base in cases:
            with self.subTest(ci=ci):
                stream = replace(assay, carbon_intensity_kg_co2e_per_kg_h2=ci)
                self.assertEqual(stream.us_45v_credit_usd_per_kg_h2, nominal)
                self.assertEqual(stream.us_45v_credit_for_year(2026), current_base)
                self.assertEqual(stream.us_45v_credit_for_year(2026, prevailing_wage_and_apprenticeship=True), round(current_base * 5, 3))
        with self.assertRaises(ValueError):
            assay.us_45v_credit_for_year(2027)
        threshold = 28.2 * 119.96 / 1000
        self.assertTrue(replace(assay, carbon_intensity_kg_co2e_per_kg_h2=threshold).is_eu_rfnbo_compliant)
        self.assertFalse(replace(assay, carbon_intensity_kg_co2e_per_kg_h2=threshold + 1e-6).is_eu_rfnbo_compliant)

    def test_blending_and_unknown_ci(self):
        green = HydrogenAssay.green_h2_rfnbo()
        blue = HydrogenAssay.blue_h2_atr()
        blend = green.blend_with(blue, self_h2_mass_fraction=0.75)
        self.assertAlmostEqual(blend.carbon_intensity_kg_co2e_per_kg_h2, 0.8)
        unknown = replace(blue, carbon_intensity_kg_co2e_per_kg_h2=None)
        self.assertIsNone(green.blend_with(unknown, self_h2_mass_fraction=0.5).carbon_intensity_kg_co2e_per_kg_h2)
        self.assertEqual(green.blend_with(unknown, self_h2_mass_fraction=1).carbon_intensity_kg_co2e_per_kg_h2, 0.4)
        with self.assertRaises(ValueError):
            HydrogenAssay.compressed_h2_350bar().blend_with(HydrogenAssay.compressed_h2_700bar(), self_h2_mass_fraction=0.5)
        # A blend just below a tier boundary must not round into the next tier.
        near = replace(green, carbon_intensity_kg_co2e_per_kg_h2=0.44999)
        self.assertEqual(near.blend_with(near, self_h2_mass_fraction=0.5).us_45v_credit_usd_per_kg_h2, 3)

    def test_invalid_specs_and_arguments(self):
        assay = HydrogenAssay.gaseous_h2_ntp()
        invalid = {"state": ("solid", ""), "density_kg_per_m3": (0, 2001, float("nan")), "h2_mass_fraction": (0, 1.1), "reference_temp_c": (-273.15, 201), "reference_pressure_bar": (0, float("inf")), "lhv_mj_per_kg_h2": (0, float("nan")), "hhv_mj_per_kg_h2": (100, float("inf")), "default_cracking_efficiency": (0, 1.1), "purity_pct": (0, 101), "carbon_intensity_kg_co2e_per_kg_h2": (-1, float("nan"))}
        for field, values in invalid.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                    replace(assay, **{field: value})
        for method in (assay.kg_price_to_mwh_price, assay.mwh_price_to_kg_price, assay.kg_price_to_mmbtu_price, assay.mmbtu_price_to_kg_price):
            with self.assertRaises(ValueError):
                method(1, basis="invalid")
            with self.assertRaises(ValueError):
                method(float("nan"))
        with self.assertRaises(ValueError):
            create_hydrogen_context(assay, basis="invalid")
        for efficiency in (0, 1.1, float("nan")):
            with self.assertRaises(ValueError):
                assay.carrier_to_h2_mass(ureg.Quantity(1, "kg"), cracking_efficiency=efficiency)
            with self.assertRaises(ValueError):
                assay.electrolysis_power_required(ureg.Quantity(1, "kg"), electrolyser_efficiency_lhv=efficiency)
        with self.assertRaises(TypeError):
            assay.volume_to_mass(1)
        with self.assertRaises(pint.DimensionalityError):
            assay.volume_to_mass(ureg.Quantity(1, "kg"))


if __name__ == "__main__":
    unittest.main()
