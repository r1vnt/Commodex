"""Numerical solver hardening tests for real-gas compressibility (SolverResult & failure taxonomy)."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
import math
import sys
from types import SimpleNamespace
import unittest
from unittest.mock import patch
import warnings

import core.units as units
from core.units import conversions


class ZFactorSolverHardeningTests(unittest.TestCase):
    """Verifies SolverResult structure, unrounded precision, failure taxonomy, and HEOS vs GERG-2008 identity."""

    def setUp(self) -> None:
        self.pressure = units.Q_(8.0, "MPa")
        self.temperature = units.Q_(310.0, "K")
        self.composition = {"CH4": 0.92, "C2H6": 0.05, "N2": 0.02, "CO2": 0.01}

    def test_solver_result_structure_immutability_and_unrounded_precision(self) -> None:
        res = units.solve_z_factor(
            self.pressure,
            self.temperature,
            relative_density=0.61,
            method="dak",
            tol=1e-12,
        )
        self.assertIsInstance(res, units.SolverResult)
        self.assertTrue(res.converged)
        self.assertGreaterEqual(res.iterations, 1)
        self.assertLessEqual(res.residual, 1e-11)
        self.assertGreaterEqual(res.residual, 0.0)
        self.assertEqual(res.method, "DAK")
        self.assertIsNone(res.backend_version)
        self.assertIsInstance(res.warnings, tuple)
        self.assertAlmostEqual(float(res), res.value, places=15)

        # Unrounded full IEEE-754 float precision (not truncated to 4-6 decimals)
        self.assertNotEqual(res.value, round(res.value, 6))

        # calculate_z_factor parity with return_result=False vs True
        scalar_z = units.calculate_z_factor(
            self.pressure,
            self.temperature,
            relative_density=0.61,
            method="dak",
            tol=1e-12,
        )
        detailed_z = units.calculate_z_factor(
            self.pressure,
            self.temperature,
            relative_density=0.61,
            method="dak",
            tol=1e-12,
            return_result=True,
        )
        self.assertEqual(scalar_z, res.value)
        self.assertEqual(detailed_z, res)

        # GasAssay.solve_z_factor and GasAssay.z_factor(return_result=True)
        assay = units.GasAssay.henry_hub()
        assay_res = assay.solve_z_factor(self.pressure, self.temperature, method="dak")
        self.assertIsInstance(assay_res, units.SolverResult)
        self.assertEqual(assay.z_factor(self.pressure, self.temperature, method="dak", return_result=True), assay_res)
        self.assertEqual(assay.z_factor(self.pressure, self.temperature, method="dak"), assay_res.value)

        with self.assertRaises(FrozenInstanceError):
            res.value = 1.0  # type: ignore[misc]

    def test_distinguishes_convergence_nonconvergence_domain_and_backend_failures(self) -> None:
        # 1. Failure to converge -> SolverConvergenceError (not Domain or Backend)
        with self.assertRaises(units.SolverConvergenceError) as conv_ctx:
            units.solve_z_factor(
                self.pressure,
                self.temperature,
                relative_density=0.62,
                method="dak",
                max_iter=1,
                tol=1e-16,
            )
        conv_err = conv_ctx.exception
        self.assertIsInstance(conv_err, ArithmeticError)
        self.assertNotIsInstance(conv_err, units.SolverDomainError)
        self.assertNotIsInstance(conv_err, units.SolverBackendError)
        self.assertIsNotNone(conv_err.result)
        assert conv_err.result is not None
        self.assertFalse(conv_err.result.converged)
        self.assertEqual(conv_err.result.iterations, 1)
        self.assertGreater(conv_err.result.residual, 1e-16)
        self.assertTrue(any("failed to converge" in w for w in conv_err.result.warnings))

        # Non-raising diagnostic mode returns SolverResult(converged=False, ...)
        nonconv = units.solve_z_factor(
            self.pressure,
            self.temperature,
            relative_density=0.62,
            method="dak",
            max_iter=1,
            tol=1e-16,
            raise_on_nonconvergence=False,
        )
        self.assertFalse(nonconv.converged)
        self.assertEqual(nonconv.iterations, 1)
        self.assertGreater(nonconv.residual, 1e-16)
        self.assertEqual(nonconv.method, "DAK")

        # 2. Out-of-domain inputs -> SolverDomainError (even when raise_on_nonconvergence=False)
        out_of_domain_cases = [
            # relative_density out of [0.50, 1.50]
            lambda: units.solve_z_factor(self.pressure, self.temperature, relative_density=0.35, method="dak", raise_on_nonconvergence=False),
            lambda: units.solve_z_factor(self.pressure, self.temperature, relative_density=1.80, method="dak"),
            # Negative absolute pressure & sub-zero Kelvin temperature
            lambda: units.solve_z_factor(units.Q_(-5.0, "psia"), self.temperature, method="dak"),
            lambda: units.solve_z_factor(self.pressure, units.Q_(-10.0, "K"), method="dak"),
            # Extreme pressure (P_pr > 30) & cryogenic sub-critical liquid state (T_pr < 0.70 or T_pr < 1 & P_pr >= 1)
            lambda: units.solve_z_factor(units.Q_(250.0, "MPa"), self.temperature, relative_density=0.60, method="dak"),
            lambda: units.solve_z_factor(units.Q_(5.0, "MPa"), units.Q_(120.0, "K"), relative_density=0.60, method="dak"),
            lambda: units.solve_z_factor(units.Q_(10.0, "MPa"), units.Q_(180.0, "K"), relative_density=0.60, method="dak"),
            # Non-hydrocarbon component in DAK
            lambda: units.solve_z_factor(self.pressure, self.temperature, method="dak", molar_composition={"H2": 1.0}),
            # Invalid composition sum or unknown component
            lambda: units.solve_z_factor(self.pressure, self.temperature, method="aga8", molar_composition={"CH4": 0.5}),
            lambda: units.solve_z_factor(self.pressure, self.temperature, method="heos", molar_composition={"XENON": 1.0}),
            # Sour-gas fractions out of [0, 1)
            lambda: units.solve_z_factor(self.pressure, self.temperature, method="dak", co2_mol_frac=0.7, h2s_mol_frac=0.4),
        ]
        for idx, case in enumerate(out_of_domain_cases):
            with self.subTest(domain_case=idx):
                with self.assertRaises(units.SolverDomainError) as dom_ctx:
                    case()
                self.assertIsInstance(dom_ctx.exception, ValueError)
                self.assertNotIsInstance(dom_ctx.exception, units.SolverConvergenceError)
                self.assertNotIsInstance(dom_ctx.exception, units.SolverBackendError)

        # 3. Backend failure -> SolverBackendError & SolverBackendUnavailableError
        with patch.dict(sys.modules, {"pyaga8": None}):
            with self.assertRaises(units.SolverBackendUnavailableError) as miss_ctx:
                units.solve_z_factor(self.pressure, self.temperature, method="aga8", molar_composition=self.composition)
            self.assertIsInstance(miss_ctx.exception, units.SolverBackendError)
            self.assertIsInstance(miss_ctx.exception, ImportError)
            self.assertNotIsInstance(miss_ctx.exception, units.SolverConvergenceError)
            self.assertNotIsInstance(miss_ctx.exception, units.SolverDomainError)

        def crashing_props(*_args):
            raise RuntimeError("Helmholtz flash failed to bracket vapor root")

        fake_cp = SimpleNamespace(__version__="8.1.0")
        fake_cp_mod = SimpleNamespace(PropsSI=crashing_props)
        with patch.dict(sys.modules, {"CoolProp": fake_cp, "CoolProp.CoolProp": fake_cp_mod}):
            with self.assertRaises(units.SolverBackendError) as back_ctx:
                units.solve_z_factor(self.pressure, self.temperature, method="heos", molar_composition=self.composition)
            self.assertNotIsInstance(back_ctx.exception, units.SolverConvergenceError)
            self.assertNotIsInstance(back_ctx.exception, units.SolverDomainError)

            # Non-physical / NaN backend return also raises SolverBackendError
            with patch.object(fake_cp_mod, "PropsSI", return_value=-0.25):
                with self.assertRaises(units.SolverBackendError):
                    units.solve_z_factor(self.pressure, self.temperature, method="gerg2008", molar_composition=self.composition)

        # pyaga8 0.1.18 success & runtime crash paths
        class FakeComposition:
            pass

        class FakeDetailOk:
            z = 0.87654321098

            def set_composition(self, _mix):
                pass

            def calc_density(self):
                pass

            def calc_properties(self):
                pass

        class FakeDetailCrash(FakeDetailOk):
            def calc_density(self):
                raise RuntimeError("AGA-8 bracket failure")

        fake_aga_ok = SimpleNamespace(__version__="0.1.18", Composition=FakeComposition, Detail=FakeDetailOk)
        with patch.dict(sys.modules, {"pyaga8": fake_aga_ok}):
            aga_res = units.solve_z_factor(
                self.pressure, self.temperature, method="aga8", molar_composition=self.composition
            )
            self.assertEqual(aga_res.method, "aga8")
            self.assertEqual(aga_res.backend_version, "0.1.18")
            self.assertEqual(aga_res.value, 0.87654321098)
            self.assertTrue(aga_res.converged)

        fake_aga_crash = SimpleNamespace(__version__="0.1.18", Composition=FakeComposition, Detail=FakeDetailCrash)
        with patch.dict(sys.modules, {"pyaga8": fake_aga_crash}):
            with self.assertRaises(units.SolverBackendError) as aga_err:
                units.solve_z_factor(
                    self.pressure, self.temperature, method="aga8", molar_composition=self.composition
                )
            self.assertNotIsInstance(aga_err.exception, units.SolverConvergenceError)
            self.assertNotIsInstance(aga_err.exception, units.SolverDomainError)

    def test_coolprop_heos_vs_gerg2008_explicit_identity_and_auto_fallback(self) -> None:
        recorded_fluids: list[str] = []

        def fake_props_si(output: str, p_key: str, p_pa: float, t_key: str, t_k: float, fluid: str) -> float:
            recorded_fluids.append(fluid)
            return 0.88123456789 if fluid.startswith("HEOS::") else 0.88234567891

        fake_cp_sub = SimpleNamespace(PropsSI=fake_props_si)
        fake_cp_pkg = SimpleNamespace(__version__="8.1.0", CoolProp=fake_cp_sub)

        with patch.dict(
            sys.modules,
            {"pyaga8": None, "CoolProp": fake_cp_pkg, "CoolProp.CoolProp": fake_cp_sub},
        ):
            # 1. Explicit HEOS:: reports method="HEOS" and CoolProp version 8.1.0
            heos_res = units.solve_z_factor(
                self.pressure, self.temperature, method="heos", molar_composition=self.composition
            )
            self.assertEqual(heos_res.method, "HEOS")
            self.assertEqual(heos_res.backend_version, "8.1.0")
            self.assertEqual(heos_res.value, 0.88123456789)
            self.assertTrue(recorded_fluids[-1].startswith("HEOS::"))

            # 2. Explicit GERG-2008 (both "gerg2008" and "gerg-2008") uses GERG2008:: and reports "GERG-2008"
            for alias in ("gerg2008", "GERG-2008"):
                gerg_res = units.solve_z_factor(
                    self.pressure, self.temperature, method=alias, molar_composition=self.composition
                )
                self.assertEqual(gerg_res.method, "GERG-2008")
                self.assertEqual(gerg_res.backend_version, "8.1.0")
                self.assertEqual(gerg_res.value, 0.88234567891)
                self.assertTrue(recorded_fluids[-1].startswith("GERG2008::"))

            # 3. auto fallback when pyaga8 is absent constructs HEOS:: and MUST report method="HEOS", NEVER GERG-2008
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                auto_res = units.solve_z_factor(
                    self.pressure, self.temperature, method="auto", molar_composition=self.composition
                )
            self.assertEqual(auto_res.method, "HEOS")
            self.assertNotEqual(auto_res.method, "GERG-2008")
            self.assertEqual(auto_res.backend_version, "8.1.0")
            self.assertEqual(auto_res.value, 0.88123456789)
            self.assertTrue(any("HEOS" in w and "not AGA-8 or strict GERG-2008" in w for w in auto_res.warnings))
            self.assertEqual(len(caught), 1)
            self.assertIn("HEOS", str(caught[0].message))

        # 4. auto fallback when BOTH pyaga8 and CoolProp are absent reports method="DAK" with warning
        with patch.dict(sys.modules, {"pyaga8": None, "CoolProp": None, "CoolProp.CoolProp": None}):
            with warnings.catch_warnings(record=True) as caught_dak:
                warnings.simplefilter("always")
                dak_auto = units.solve_z_factor(
                    self.pressure, self.temperature, relative_density=0.60, method="auto", molar_composition=self.composition
                )
            self.assertEqual(dak_auto.method, "DAK")
            self.assertIsNone(dak_auto.backend_version)
            self.assertTrue(dak_auto.converged)
            self.assertTrue(any("auto selected DAK" in w for w in dak_auto.warnings))
            self.assertEqual(len(caught_dak), 1)

    def test_ideal_and_zero_pressure_limits_and_dak_advisory_warnings(self) -> None:
        ideal_res = units.solve_z_factor(self.pressure, self.temperature, method="ideal")
        self.assertEqual(ideal_res.value, 1.0)
        self.assertTrue(ideal_res.converged)
        self.assertEqual(ideal_res.iterations, 0)
        self.assertEqual(ideal_res.residual, 0.0)
        self.assertEqual(ideal_res.method, "ideal")
        self.assertIsNone(ideal_res.backend_version)
        self.assertEqual(ideal_res.warnings, ())

        # Low-pressure standard state records low-P_pr advisory warning in SolverResult.warnings
        low_p_res = units.solve_z_factor(
            units.Q_(101.325, "kPa"),
            units.Q_(288.15, "K"),
            relative_density=0.60,
            method="dak",
        )
        self.assertTrue(low_p_res.converged)
        self.assertTrue(math.isclose(low_p_res.value, 0.998, abs_tol=0.003))
        self.assertTrue(any("below the primary DAK (1975) fitted floor of 0.20" in w for w in low_p_res.warnings))


if __name__ == "__main__":
    unittest.main()
