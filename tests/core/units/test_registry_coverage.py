"""Executable coverage checklist and independent published conversion factors."""

import unittest

import pint

import core.units as units


# A unit's presence does not certify a commodity specification or empirical model.
COVERAGE = {
    "mass": ("kg", "mt", "short_ton", "long_ton", "cwt", "grain"),
    "precious_metals": ("toz", "pennyweight", "tola", "tael_hk", "baht_gold", "don_korea", "carat"),
    "liquid_volume": ("gallon", "imperial_gallon", "oil_barrel", "liter", "kiloliter", "meter**3"),
    "dry_volume": ("dry_gallon", "bushel", "imperial_bushel", "peck", "hectoliter"),
    "gas_volume": ("scf", "Nm3", "Sm3", "Sm3_20C"),
    "energy": ("J", "Wh", "btu_it", "btu_th", "btu_us59", "cal_it", "calorie", "therm", "dth"),
    "energy_equivalents": ("boe", "toe", "tce"),
    "power": ("MW", "kVAR", "MVA", "MWth", "MWe"),
    "capacity": ("MW", "MWth", "MWe", "MVA", "MVAr", "Ah", "MWh"),
    "pressure": ("psia", "psig", "bara", "barg", "kpaa", "kpag", "inHg", "mmH2O"),
    "temperature": ("K", "degC", "degF", "degR", "delta_degC", "delta_degF"),
    "flow": ("kg/s", "m**3/hour", "mol/s", "MWh/day"),
    "viscosity": ("centipoise", "centistokes", "Pa*s", "SUS", "SFS"),
    "density": ("kg/m**3", "g/mL"),
    "concentration": ("ppm", "ppb", "parts_per_trillion", "mole_percent", "mass_percent", "volume_percent", "mg/L", "mol/L"),
    "chemical_strength": ("degree_Baume_heavy", "degree_Baume_light", "degree_Brix", "degree_Plato", "normality", "molarity"),
    "quality": ("percent", "basis_point"),
    "electrical_storage": ("Ah", "kAh", "C_rate", "Wh/kg", "Wh/liter"),
    "solar": ("W/m**2", "kWh/m**2"),
    "wind": ("m/s", "W/m**2", "percent"),
    "weather": ("HDD", "CDD", "HDH", "CDH", "millimeter", "m/s"),
    "agriculture": ("bushel_corn", "bushel_canola", "bushel_rice", "coffee_bag", "cotton_bale", "crate", "sack", "dozen"),
    "land": ("acre", "acre_us_survey", "hectare", "feddan_4200", "mu_china", "kg/hectare"),
    "forestry": ("board_foot", "cord", "meter**3", "admt", "bdmt"),
    "mining": ("wmt", "dmt", "dmtu", "contained_metal_tonne", "payable_metal_tonne"),
    "nuclear": ("kgU", "mtU", "SWU", "GWd/mtU"),
    "freight": ("TEU", "FEU", "tonne_mile", "short_ton_mile", "tonne_nautical_mile", "WS"),
    "shipping": ("dwt", "gross_tonnage", "net_tonnage", "displacement_tonne"),
    "environmental": ("EUA", "UKA", "VCU", "REC", "REGO", "RIN"),
    "financial": ("USD", "USD_cent", "basis_point", "percent", "price_tick"),
    "contracts": ("lot", "contract", "delivery_unit"),
    "water": ("acre_foot", "acre_foot_us_survey", "meter**3", "megaliter", "m**3/s"),
}


class RegistryCoverageTests(unittest.TestCase):
    def test_all_32_checklist_categories_have_explicit_representations(self):
        self.assertEqual(len(COVERAGE), 32)
        for category, expressions in COVERAGE.items():
            for expression in expressions:
                with self.subTest(category=category, unit=expression):
                    self.assertIsInstance(units.parse_unit(expression), pint.Unit)

    def test_nist_sp811_factors_at_published_precision(self):
        # Appendix B.8 factors are rounded except where the source marks exact.
        cases = (
            ("pound", "kg", .4535924), ("short_ton", "kg", 907.1847),
            ("long_ton", "kg", 1016.047), ("grain", "kg", 6.479891e-5),
            ("toz", "g", 31.10348), ("pennyweight", "g", 1.555174),
            ("carat", "g", .2), ("gallon", "L", 3.785412),
            ("imperial_gallon", "L", 4.54609), ("bushel", "L", 35.23907),
            ("peck", "L", 8.809768), ("hectare", "m**2", 10000),
            ("acre_us_survey", "m**2", 4046.873),
            ("acre_foot_us_survey", "m**3", 1233.489),
            ("Wh", "J", 3600), ("btu_it", "J", 1055.056),
            ("btu_th", "J", 1054.350), ("centipoise", "Pa*s", .001),
            ("centistokes", "m**2/s", 1e-6), ("atm", "Pa", 101325),
            ("bar", "Pa", 100000), ("mile", "m", 1609.344),
            ("nautical_mile", "m", 1852), ("Ah", "coulomb", 3600),
        )
        for source, target, expected in cases:
            with self.subTest(source=source):
                actual = float(units.Q_(1., source).to(target).magnitude)
                self.assertAlmostEqual(actual, expected, delta=max(1e-12, abs(expected) * 8e-7))

    def test_explicit_conventions_do_not_collapse(self):
        for source, target in (
            ("mole_percent", "mass_percent"), ("volume_percent", "mass_percent"),
            ("normality", "mol/L"), ("degree_Brix", "degree_Plato"),
            ("degree_Baume_heavy", "kg/m**3"), ("SUS", "centistokes"),
            ("SUS", "SFS"), ("gross_tonnage", "net_tonnage"),
            ("gross_tonnage", "kg"), ("gross_tonnage", "meter**3"),
            ("lot", "contract"), ("contract", "kg"), ("price_tick", "USD"),
            ("crate", "kg"),
        ):
            with self.subTest(source=source, target=target), self.assertRaises(pint.DimensionalityError):
                units.Q_(1, source).to(target)
        for token in ("ppt", "mu", "ton_mile"):
            with self.assertRaisesRegex(ValueError, "Ambiguous"):
                units.parse_quantity("1 " + token)

    def test_regional_and_volume_conventions_have_explicit_scaling(self):
        self.assertEqual(units.Q_(1, "don_korea").to("g").magnitude, 3.75)
        self.assertEqual(units.Q_(15, "mu_china").to("hectare").magnitude, 1)
        self.assertEqual(units.Q_(1, "feddan_4200").to("m**2").magnitude, 4200)
        self.assertAlmostEqual(units.Q_(1, "imperial_bushel").to("L").magnitude, 36.36872)
        self.assertAlmostEqual(units.Q_(1, "dry_gallon").to("bushel").magnitude, .125)
        self.assertNotEqual(units.Q_(1, "tonne_mile").to("tonne*meter").magnitude, units.Q_(1, "tonne_nautical_mile").to("tonne*meter").magnitude)
        self.assertAlmostEqual(units.Q_(1, "GWd/mtU").to("J/kgU").magnitude, 8.64e10)
