"""Tests for specification-dependent energy equivalents, ASTM D2161 viscosity, and chemical strength."""

from __future__ import annotations

from decimal import Decimal
import unittest

import pint

import core.units as units
from core.units.converters import chemical, energy, viscosity
from core.units.converters.viscosity import DomainError


class EnergyEquivalentsTests(unittest.TestCase):
    def test_gge_and_dge_specs_load_from_json_and_round_trip(self):
        self.assertIn("NIST_HB44", energy.GGE_SPECS)
        self.assertIn("DOE_AFDC", energy.GGE_SPECS)
        self.assertIn("EPA", energy.GGE_SPECS)
        self.assertIn("NIST_HB44", energy.DGE_SPECS)
        self.assertIn("DOE_AFDC", energy.DGE_SPECS)

        nist = energy.GGE_SPECS["NIST_HB44"]
        self.assertEqual(nist.energy_per_gallon_btu, 114118.8)
        self.assertEqual(nist.basis, "LHV")

        fuel_energy = units.Q_(Decimal("1141188"), "btu_it")
        gge = energy.energy_to_gge(fuel_energy, "NIST_HB44", basis="LHV")
        self.assertEqual(gge.value, Decimal("10"))
        self.assertEqual(gge.spec, "NIST_HB44")
        self.assertEqual(gge.basis, "LHV")
        self.assertEqual(gge.fuel_type, "GGE")
        self.assertIsInstance(gge.value, Decimal)

        restored = energy.gge_to_energy(gge, target_unit="btu_it")
        self.assertAlmostEqual(float(restored.magnitude), 1141188.0, places=9)
        self.assertIsInstance(restored.magnitude, Decimal)

        dge_energy = units.Q_(256976.0, "btu_it")
        dge = energy.energy_to_dge(dge_energy, "DOE_AFDC", basis="LHV")
        self.assertAlmostEqual(dge.value, 2.0, places=12)
        self.assertIsInstance(dge.value, float)
        restored_dge = energy.dge_to_energy(dge, target_unit="MJ")
        self.assertAlmostEqual(float(restored_dge.to("btu_it").magnitude), 256976.0, places=7)

    def test_gge_dge_basis_mismatch_refused_without_explicit_ratio(self):
        hhv_energy = units.Q_(Decimal("120000"), "btu_it")
        # NIST_HB44 is an LHV standard -> HHV input without ratio must fail
        with self.assertRaisesRegex(ValueError, "basis mismatch"):
            energy.energy_to_gge(hhv_energy, "NIST_HB44", basis="HHV")

        # Missing basis on generic energy unit must fail
        with self.assertRaisesRegex(ValueError, "basis"):
            energy.energy_to_gge(hhv_energy, "NIST_HB44")

        # Tagged MWh_lhv infers LHV automatically, while MWh_hhv conflicts with LHV spec
        tagged_lhv = energy.energy_to_gge(units.Q_(Decimal("1"), "MWh_lhv"), "EPA")
        self.assertGreater(tagged_lhv.value, Decimal("29"))
        with self.assertRaisesRegex(ValueError, "basis mismatch"):
            energy.energy_to_gge(units.Q_(Decimal("1"), "MWh_hhv"), "EPA")

        # Explicit physical LHV/HHV ratio in (0, 1] succeeds; inverted ratio > 1 fails
        with self.assertRaisesRegex(ValueError, r"\(0, 1\]"):
            energy.energy_to_gge(hhv_energy, "NIST_HB44", basis="HHV", basis_conversion_ratio=Decimal("1.10"))

        converted = energy.energy_to_gge(
            units.Q_(Decimal("114118.8") / Decimal("0.90"), "btu_it"),
            "NIST_HB44",
            basis="HHV",
            basis_conversion_ratio=Decimal("0.90"),
        )
        self.assertAlmostEqual(float(converted.value), 1.0, places=12)

        # Cross-type spec guard (DGE spec passed to energy_to_gge)
        with self.assertRaisesRegex(ValueError, "fuel_type"):
            energy.energy_to_gge(units.Q_(100000.0, "btu_it"), energy.DGE_SPECS["DOE_AFDC"], basis="LHV")

    def test_natural_gas_composition_and_volume_energy_conversions(self):
        comp = {"CH4": Decimal("0.94"), "C2H6": Decimal("0.04"), "N2": Decimal("0.02")}
        cv_iso15_hhv = energy.calculate_gas_calorific_value(
            comp, standard="ISO_6976", reference_conditions="iso15", basis="HHV"
        )
        cv_iso15_lhv = energy.calculate_gas_calorific_value(
            comp, standard="ISO_6976", reference_conditions="iso15", basis="LHV"
        )
        self.assertEqual(cv_iso15_hhv.units, units.parse_unit("MJ / Sm3"))
        self.assertGreater(cv_iso15_hhv.magnitude, cv_iso15_lhv.magnitude)
        self.assertTrue(37.0 < float(cv_iso15_hhv.magnitude) < 40.0)

        cv_gpa_1atm = energy.calculate_gas_calorific_value(
            {"CH4": 1.0}, standard="GPA_2172", reference_conditions="atm60f", basis="HHV"
        )
        cv_gpa_ferc = energy.calculate_gas_calorific_value(
            {"CH4": 1.0}, standard="GPA_2172", reference_conditions="ferc", basis="HHV"
        )
        self.assertGreater(cv_gpa_ferc.magnitude, cv_gpa_1atm.magnitude)
        self.assertTrue(1010.0 < cv_gpa_1atm.magnitude < 1015.0)

        # Volume -> energy & tuple unpacking + inverse round-trip
        vol = units.Q_(Decimal("1000"), "Sm3")
        res = energy.gas_volume_to_energy(vol, cv_iso15_hhv, basis="HHV", reference_conditions="iso15")
        energy_qty, basis_tag = res
        self.assertEqual(basis_tag, "HHV")
        self.assertEqual(res.reference_conditions, "iso15")
        self.assertIsInstance(energy_qty.magnitude, Decimal)

        restored_vol = energy.energy_to_gas_volume(
            energy_qty, cv_iso15_hhv, basis="HHV", reference_conditions="iso15", target_unit="Sm3"
        )
        self.assertAlmostEqual(float(restored_vol.magnitude), 1000.0, places=10)

        # Mismatched gas reference states (Nm3 at 0 °C vs MJ/Sm3 at 15 °C) must fail
        with self.assertRaises(pint.DimensionalityError):
            energy.gas_volume_to_energy(units.Q_(1000.0, "Nm3"), cv_iso15_hhv, basis="HHV")

        # Geometric m^3 without explicit reference_conditions must fail
        with self.assertRaisesRegex(ValueError, "reference_conditions"):
            energy.gas_volume_to_energy(units.Q_(100.0, "m^3"), units.Q_(38.0, "MJ / m^3"), basis="HHV")

        # Pure inert gas has zero calorific value -> must fail
        with self.assertRaisesRegex(ValueError, "zero calorific"):
            energy.calculate_gas_calorific_value({"N2": 1.0})


class ViscosityASTMTests(unittest.TestCase):
    def test_sus_low_and_high_regimes_and_temperature_correction(self):
        # Published ASTM D2161 Table 1 checkpoints at 100 °F
        self.assertAlmostEqual(viscosity.cst_to_sus(2.0, temp_f=100.0), 32.6026, places=2)
        self.assertAlmostEqual(viscosity.cst_to_sus(10.0, temp_f=100.0), 58.91, delta=0.15)

        # Exact asymptotic linear regime >= 70 cSt at 100 °F
        self.assertAlmostEqual(viscosity.cst_to_sus(100.0, temp_f=100.0), 463.24, places=10)
        self.assertAlmostEqual(viscosity.cst_to_sus(70.0, temp_f=100.0), 70.0 * 4.6324, places=10)

        # Temperature correction at 210 °F: factor = 1 + 0.000061 * (210 - 100) = 1.00671
        sus_100 = viscosity.cst_to_sus(15.0, temp_f=100.0)
        sus_210 = viscosity.cst_to_sus(15.0, temp_f=210.0)
        self.assertAlmostEqual(sus_210 / sus_100, 1.0 + 0.000061 * 110.0, places=12)

        # Exact round-trip inversion across both table and linear regimes and temperatures
        for nu in (1.81, 2.0, 4.5, 10.3, 32.7, 69.9, 70.0, 150.0, 500.0):
            for t_f in (100.0, 140.0, 210.0):
                with self.subTest(nu=nu, temp_f=t_f):
                    sus = viscosity.cst_to_sus(nu, temp_f=t_f)
                    back = viscosity.sus_to_cst(sus, temp_f=t_f)
                    self.assertAlmostEqual(back, nu, places=10)

    def test_sfs_table_and_asymptotic_regimes(self):
        # Published ASTM D2161 Table 3 checkpoints at 122 °F
        self.assertAlmostEqual(viscosity.cst_to_sfs(48.0), 25.11, delta=0.05)
        self.assertAlmostEqual(viscosity.cst_to_sfs(50.0), 26.03, delta=0.05)
        self.assertAlmostEqual(viscosity.cst_to_sfs(100.0), 48.63, delta=0.05)
        self.assertAlmostEqual(viscosity.cst_to_sfs(200.0), 94.77, delta=0.05)

        # Exact asymptotic linear regime >= 400 cSt at 122 °F
        self.assertAlmostEqual(viscosity.cst_to_sfs(400.0), 400.0 * 0.4717, places=10)
        self.assertAlmostEqual(viscosity.cst_to_sfs(1000.0), 471.7, places=10)

        for nu in (48.0, 55.5, 100.0, 235.0, 399.5, 400.0, 850.0):
            with self.subTest(nu=nu):
                sfs = viscosity.cst_to_sfs(nu, temp_f=units.Q_(50.0, "degC"))
                back = viscosity.sfs_to_cst(sfs, temp_f=122.0)
                self.assertAlmostEqual(back, nu, places=10)

    def test_domain_bounds_raise_domain_error(self):
        # SUS undefined below 1.81 cSt (~32.0–32.6 s)
        with self.assertRaises(DomainError):
            viscosity.cst_to_sus(1.80)
        with self.assertRaises(DomainError):
            viscosity.sus_to_cst(31.5)

        # SFS undefined below 48 cSt (~25.1 s at 122 °F) and at non-122 °F temperatures
        with self.assertRaises(DomainError):
            viscosity.cst_to_sfs(47.9)
        with self.assertRaises(DomainError):
            viscosity.sfs_to_cst(24.9)
        with self.assertRaises(DomainError):
            viscosity.cst_to_sfs(100.0, temp_f=100.0)

        # Temperature interval rejected; <= absolute zero raises DomainError
        with self.assertRaises(ValueError):
            viscosity.cst_to_sus(10.0, temp_f=units.Q_(10.0, "delta_degF"))
        with self.assertRaises(DomainError):
            viscosity.cst_to_sus(10.0, temp_f=-500.0)

    def test_two_stage_dynamic_kinematic_saybolt_pipeline(self):
        rho = units.Q_(Decimal("850"), "kg/m^3")
        mu = units.Q_(Decimal("85"), "centipoise")  # 85 cP / 0.85 g/cm^3 = 100 cSt

        kin = viscosity.dynamic_to_kinematic_viscosity(mu, rho)
        self.assertEqual(kin.units, units.ureg.centistokes)
        self.assertIsInstance(kin.magnitude, Decimal)
        self.assertAlmostEqual(float(kin.magnitude), 100.0, places=10)

        # Full two-stage cP -> SUS -> cP via convert_viscosity
        sus_q = viscosity.convert_viscosity(mu, "SUS", density=rho, temperature=units.Q_(100, "degF"))
        self.assertEqual(sus_q.units, units.ureg.SUS)
        self.assertIsInstance(sus_q.magnitude, Decimal)
        self.assertAlmostEqual(float(sus_q.magnitude), 463.24, places=9)

        mu_back = viscosity.convert_viscosity(
            sus_q, "centipoise", density=rho, temperature=units.Q_(100, "degF")
        )
        self.assertAlmostEqual(float(mu_back.magnitude), 85.0, places=9)
        self.assertIsInstance(mu_back.magnitude, Decimal)

        # Missing density or temperature on cross-stage conversion must fail
        with self.assertRaisesRegex(ValueError, "density"):
            viscosity.convert_viscosity(mu, "SUS", temperature=100.0)
        with self.assertRaisesRegex(ValueError, "temperature"):
            viscosity.convert_viscosity(kin, "SUS")


class ChemicalStrengthTests(unittest.TestCase):
    def test_normality_and_molarity_require_explicit_equivalence_factor(self):
        # Diprotic H2SO4 complete neutralisation (z = 2)
        m_qty = units.Q_(Decimal("1.5"), "mol/L")
        n_qty = chemical.molarity_to_normality(m_qty, 2)
        self.assertEqual(n_qty.units, units.ureg.normality)
        self.assertEqual(n_qty.magnitude, Decimal("3.0"))

        restored = chemical.normality_to_molarity(n_qty, 2, target_unit="mol/L")
        self.assertEqual(restored.magnitude, Decimal("1.5"))

        # Scalar mode
        self.assertAlmostEqual(chemical.molarity_to_normality(0.25, 5), 1.25)
        self.assertAlmostEqual(chemical.normality_to_molarity(1.25, 5), 0.25)

        # Reject <= 0, bool, or chemical formula strings
        for bad_z in (0, -2, True, "H2SO4"):
            with self.subTest(bad_z=bad_z), self.assertRaises((ValueError, TypeError)):
                chemical.molarity_to_normality(1.0, bad_z)  # type: ignore[arg-type]

    def test_heavy_and_light_baume_two_branch_hyperbolic_scales(self):
        # 66 °Bé commercial H2SO4 -> SG = 145 / (145 - 66) = 145 / 79 = 1.835443...
        sg_66 = chemical.baume_to_sg(Decimal("66"), scale="heavy")
        self.assertAlmostEqual(float(sg_66), 145.0 / 79.0, places=12)
        self.assertAlmostEqual(float(chemical.sg_to_baume(sg_66, scale="heavy")), 66.0, places=12)

        # Rational modulus 144.3
        sg_rat = chemical.baume_to_sg(50.0, scale="heavy", modulus="rational")
        self.assertAlmostEqual(sg_rat, 144.3 / (144.3 - 50.0), places=12)
        self.assertAlmostEqual(chemical.sg_to_baume(sg_rat, scale="heavy", modulus="rational"), 50.0, places=12)

        # Light Baumé (US 140 / SG - 130): SG = 0.80 -> 175 - 130 = 45 °Bé
        be_light = chemical.sg_to_baume(units.Q_(Decimal("0.80"), "dimensionless"), scale="light")
        self.assertEqual(be_light.units, units.ureg.degree_Baume_light)
        self.assertEqual(be_light.magnitude, Decimal("45"))
        sg_back = chemical.baume_to_sg(be_light, scale="light")
        self.assertEqual(sg_back.magnitude, Decimal("0.8"))

        # Strict domain errors at branch boundaries
        with self.assertRaises(ValueError):
            chemical.sg_to_baume(0.99, scale="heavy")
        with self.assertRaises(ValueError):
            chemical.sg_to_baume(1.01, scale="light")
        with self.assertRaises(ValueError):
            chemical.baume_to_sg(145.0, scale="heavy")
        with self.assertRaises(ValueError):
            chemical.baume_to_sg(9.5, scale="light")
        # Passing light Baumé Quantity into heavy branch must fail
        with self.assertRaises(pint.DimensionalityError):
            chemical.baume_to_sg(units.Q_(20.0, "degree_Baume_light"), scale="heavy")

    def test_plato_and_brix_20c_scales_and_divergence(self):
        # Plato Newton-Raphson round-trip across 0..30 °P
        for p in (0.0, 5.0, 12.0, 20.0, 30.0):
            with self.subTest(plato=p):
                sg = chemical.plato_to_sg(Decimal(str(p)))
                p_back = chemical.sg_to_plato(sg)
                self.assertIsInstance(sg, Decimal)
                self.assertAlmostEqual(float(p_back), p, places=11)

        # Brix ICUMSA 20/20 °C round-trip across full 0..85 °Bx range + NBS Table 109 checkpoints
        checkpoints = ((0.0, 1.00000), (10.0, 1.03998), (50.0, 1.23187), (70.0, 1.34956), (85.0, 1.44794))
        for bx, expected_sg in checkpoints:
            with self.subTest(brix=bx):
                sg_q = chemical.brix_to_sg(units.Q_(bx, "degree_Brix"))
                self.assertAlmostEqual(float(sg_q.magnitude), expected_sg, places=5)
                bx_back = chemical.sg_to_brix(sg_q)
                self.assertEqual(bx_back.units, units.ureg.degree_Brix)
                self.assertAlmostEqual(float(bx_back.magnitude), bx, places=10)

        # Interpolated non-integer Brix round-trips cleanly
        for bx in (12.34, 42.75, 67.8, 83.2):
            sg_val = chemical.brix_to_sg(bx)
            self.assertAlmostEqual(chemical.sg_to_brix(sg_val), bx, places=10)

        # Plato and Brix use distinct empirical curves and diverge slightly at 30%
        self.assertNotAlmostEqual(chemical.plato_to_sg(30.0), chemical.brix_to_sg(30.0), places=4)

        # Reference temperature guard (20 °C mandatory) and range bounds
        with self.assertRaisesRegex(ValueError, "20"):
            chemical.sg_to_brix(1.05, ref_temp_c=15.0)
        with self.assertRaisesRegex(ValueError, "20"):
            chemical.plato_to_sg(12.0, ref_temp_c=15.56)
        with self.assertRaises(ValueError):
            chemical.plato_to_sg(31.0)
        with self.assertRaises(ValueError):
            chemical.brix_to_sg(86.0)

        # Unified router test
        norm = chemical.convert_chemical_strength(units.Q_(Decimal("0.5"), "molarity"), "normality", equivalence_factor=3)
        self.assertEqual(norm, units.Q_(Decimal("1.5"), "normality"))
        brix_q = chemical.convert_chemical_strength(units.Q_(Decimal("1.08292"), "dimensionless"), "degree_Brix")
        self.assertAlmostEqual(float(brix_q.magnitude), 20.0, places=3)


if __name__ == "__main__":
    unittest.main()
