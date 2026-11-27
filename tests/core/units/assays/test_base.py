"""Independent mass balances and explicit multi-component settlements."""

from dataclasses import FrozenInstanceError
from decimal import Decimal as D
import unittest

import pint
import core.units as u


class SpecificationEngineTests(unittest.TestCase):
    def test_composition_normalization_is_explicit_and_balance_is_retained(self):
        values = {"Cu": D(".25"), "Zn": D(".15")}
        composition = u.ElementComposition(values)
        values["Cu"] = D(1)
        self.assertEqual(composition.unassayed_fraction, D(".6"))
        self.assertEqual(composition.fraction("Cu"), D(".25"))
        normalized = u.ElementComposition({"Cu": 25, "Zn": 15}, normalize=True)
        self.assertEqual(normalized.source_total, 40)
        self.assertEqual(normalized.fraction("Cu"), D(".625"))
        for invalid in ({"Cu": 25}, {"Cu": -1}, {"Cu": True}, {"Cu": float("nan")}):
            with self.subTest(invalid=invalid), self.assertRaises((ValueError, TypeError)):
                u.ElementComposition(invalid)
        with self.assertRaises(TypeError):
            composition.mass_fractions["Cu"] = D(1)
        grades = u.ElementComposition.from_grades({"Au": u.Q_(D(5), "g/metric_ton"), "Ag": u.Q_(D(500), "g/metric_ton")})
        self.assertEqual(grades.fraction("Au"), D(".000005"))
        self.assertEqual(grades.fraction("Ag"), D(".0005"))

    def test_multi_component_recovery_and_payability_keep_every_metal_independent(self):
        engine = u.ElementRecoveryEngine(u.ElementComposition({"Cu": D(".25"), "Zn": D(".15")}),
                                        {"Cu": D(".9"), "Zn": D(".8")}, {"Cu": D(".95"), "Zn": D(".9")}, {"Cu": D(".01")})
        dry = u.Q_(D(1000), "kg")
        self.assertEqual(engine.contained_mass(dry, component="Cu").magnitude, 250)
        self.assertEqual(engine.recovered_mass(dry, component="Cu").magnitude, 225)
        self.assertEqual(engine.payable_mass(dry, component="Cu").magnitude, D("204.25"))
        self.assertEqual(engine.payable_mass(dry, component="Zn").magnitude, 108)
        for stage in ("contained", "recovered", "payable"):
            output = getattr(engine, stage + "_mass")(dry, component="Cu")
            self.assertAlmostEqual(engine.feed_required(output, component="Cu", stage=stage).magnitude, 1000, places=20)
        missing = u.ElementRecoveryEngine(engine.composition)
        with self.assertRaisesRegex(ValueError, "recovery"):
            missing.recovered_mass(dry, component="Cu")
        with self.assertRaisesRegex(ValueError, "payability"):
            missing.payable_mass(dry, component="Cu")
        with self.assertRaises(ValueError):
            u.ElementRecoveryEngine(engine.composition, {"Au": 1})
        with self.assertRaisesRegex(ValueError, "basis"):
            engine.contained_mass(u.Q_(D(1), "wmt"), component="Cu")

    def test_settlement_dry_moisture_charge_and_impurity_accounting(self):
        assay = u.MultiMetalConcentrateAssay(name="Cu-Zn concentrate", moisture_fraction=D(".1"),
                    composition=u.ElementComposition({"Cu": D(".25"), "Zn": D(".15"), "As": D(".01")}),
                    recovery_fractions={"Cu": D(".9"), "Zn": D(".8")},
                    payable_fractions={"Cu": D(".95"), "Zn": D(".9")}, deduction_fractions={"Cu": D(".01")})
        result = assay.settle(u.Q_(D(1000), "kg"), {"Cu": u.Q_(D(10), "USD/kg"), "Zn": u.Q_(D(2), "USD/kg")},
                    currency="USD", treatment_charge=u.Q_(D(50), "USD/metric_ton"),
                    refining_charges={"Cu": u.Q_(D(".1"), "USD/kg")},
                    impurity_penalties=(u.ImpurityPenalty("As", D(".005"), u.Q_(D(20), "USD/metric_ton")),))
        self.assertEqual(result.gross_value.magnitude, D("2258.5"))
        self.assertEqual(result.deductions.magnitude, D("80.425"))
        self.assertEqual(result.net_value.magnitude, D("2178.075"))
        wet_result = assay.settle(u.Q_(D(1000), "kg"), {"Cu": u.Q_(D(10), "USD/kg"), "Zn": u.Q_(D(2), "USD/kg")}, currency="USD", mass_basis="received")
        self.assertEqual(wet_result.payable_kg["Cu"], D("183.825"))
        self.assertEqual(assay.wet_to_dry_mass(u.Q_(D(1000), "wmt")).to("kg").magnitude, D(900000))
        self.assertEqual(u.convert(u.Q_(D(1000), "wmt"), "dmt", assay=assay).magnitude, D(900))
        self.assertEqual(u.convert(u.Q_(D(1000), "USD/dmt"), "USD/wmt", assay=assay).magnitude, D(900))
        with self.assertRaisesRegex(ValueError, "explicit"):
            u.convert(u.Q_(D(1), "dmt"), "m^3", assay=assay)
        with self.assertRaises(pint.DimensionalityError):
            assay.settle(u.Q_(1, "kg"), {"Cu": u.Q_(1, "EUR/kg"), "Zn": u.Q_(1, "USD/kg")}, currency="USD")
        with self.assertRaisesRegex(ValueError, "Missing prices"):
            assay.settle(u.Q_(1, "kg"), {"Cu": u.Q_(1, "USD/kg")}, currency="USD")
        # Reading a result quantity cannot mutate the immutable settlement.
        result.net_value.ito("US_cent")
        self.assertEqual(result.net_value.units, u.ureg.USD)

    def test_stoichiometry_oxide_element_and_hydrate_parsing(self):
        carbonate = u.StoichiometricCompound("CaCO3")
        lime = carbonate.equivalent_mass(u.Q_(D(1000), "kg"), "CaO", element="Ca")
        co2 = carbonate.equivalent_mass(u.Q_(D(1000), "kg"), "CO2", element="C")
        self.assertAlmostEqual((lime + co2).magnitude, D(1000), places=20)
        hydrate = u.StoichiometricCompound("CaSO4(H2O)2")
        self.assertEqual(dict(hydrate.atoms), {"Ca": 1, "S": 1, "O": 6, "H": 4})
        oxide = u.StoichiometricCompound("Nd2O3")
        self.assertAlmostEqual(float(oxide.element_fraction("Nd")), 288.48 / 336.477, places=12)
        for formula in ("Xx2O3", "H0", "2H2O", "Ca(OH", "Ca()2", "CaSO4.2H2O", "NaCl)"):
            with self.subTest(formula=formula), self.assertRaises(ValueError):
                u.StoichiometricCompound(formula)

    def test_signed_quotes_and_same_dimension_mass_processing_stay_separate(self):
        assay = u.RefinedMetalAssay(name="cathode", metal="Cu", composition=u.ElementComposition({"Cu": D(".9999")}), density=u.Q_(D(8960), "kg/m^3"))
        for numeric in (1000., D(1000)):
            original = u.Q_(numeric, "kg")
            scaled = u.convert(original, "metric_ton", assay=assay)
            self.assertEqual(scaled.magnitude, 1)
            fine = u.convert_refined_metal(original, "kg", assay=assay, mode="fine_mass")
            self.assertAlmostEqual(float(fine.magnitude), 999.9)
            self.assertIsInstance(fine.magnitude, type(numeric))
        for price in (D(0), D(-100), D(100)):
            quote = u.Q_(price, "USD/metric_ton")
            volume_quote = u.convert(quote, "EUR/m^3", assay=assay, rates={"USD/EUR": D(".9")})
            self.assertAlmostEqual(volume_quote.magnitude, price * D("8.96") * D(".9"), places=20)
        with self.assertRaisesRegex(ValueError, "moisture"):
            assay.wet_to_dry_mass(u.Q_(1, "kg"))


if __name__ == "__main__":
    unittest.main()
