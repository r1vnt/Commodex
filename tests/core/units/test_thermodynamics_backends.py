"""Real optional EOS bindings, NIST checkpoints and dependency-free API checks.

Reference mixture, T/P and Z are from NIST's independent test vectors:
https://github.com/usnistgov/AGA8/blob/master/AGA8CODE/C/Detail_test_01.cpp
https://github.com/usnistgov/AGA8/blob/master/AGA8CODE/C/GERG2008_test_01.cpp
These tests never download dependencies. Real-backend tests run when installed.
"""

from dataclasses import replace
from decimal import Decimal
from importlib.util import find_spec
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import core.units as u


NIST_MIXTURE = dict(zip(
    ("CH4", "N2", "CO2", "C2H6", "C3H8", "iC4H10", "nC4H10", "iC5H12", "nC5H12",
     "C6H14", "C7H16", "C8H18", "C9H20", "C10H22", "H2", "O2", "CO", "H2O", "H2S", "He", "Ar"),
    (.77824, .02, .06, .08, .03, .0015, .003, .0005, .00165, .00215, .00088,
     .00024, .00015, .00009, .004, .005, .002, .0001, .0025, .007, .001),
))


class Aga8ApiTests(unittest.TestCase):
    def test_pyaga8_composition_and_state_contract_without_optional_dependencies(self):
        calls = []

        class Composition:
            pass

        class Detail:
            z = .9851778445897175

            def set_composition(self, composition):
                self.composition = composition
                if not isinstance(composition, Composition):
                    raise TypeError("Expected pyaga8.Composition")

            def calc_density(self):
                calls.append((self.pressure, self.temperature, vars(self.composition)))

            def calc_properties(self):
                calls.append("properties")

        with patch.dict("sys.modules", {"pyaga8": SimpleNamespace(Detail=Detail, Composition=Composition)}):
            result = u.calculate_z_factor(u.Q_(1, "MPa"), u.Q_(300, "K"), method="aga8", molar_composition={"CH4": .9, "N2": .1})
        self.assertEqual(result, Detail.z)
        self.assertEqual(calls[0][:2], (1000., 300.))
        self.assertEqual(calls[0][2]["methane"], .9)
        self.assertEqual(calls[0][2]["nitrogen"], .1)
        self.assertEqual(calls[0][2]["hydrogen_sulfide"], 0.)
        self.assertEqual(len(calls[0][2]), 21)
        self.assertEqual(calls[1], "properties")


class RealGasEosTests(unittest.TestCase):
    @unittest.skipUnless(find_spec("pyaga8"), "optional pyaga8 not installed")
    def test_detail_matches_nist_21_component_checkpoint(self):
        z = u.calculate_z_factor(u.Q_(50, "MPa"), u.Q_(400, "K"), method="aga8", molar_composition=NIST_MIXTURE)
        self.assertAlmostEqual(z, 1.173801364147326, delta=1e-8)

    @unittest.skipUnless(find_spec("CoolProp"), "optional CoolProp not installed")
    def test_gerg_matches_nist_21_component_checkpoint(self):
        z = u.calculate_z_factor(u.Q_(50, "MPa"), u.Q_(400, "K"), method="gerg2008", molar_composition=NIST_MIXTURE)
        self.assertAlmostEqual(z, 1.174690666383717, delta=1e-8)

    @unittest.skipUnless(find_spec("CoolProp"), "optional CoolProp not installed")
    def test_gas_liquid_dew_and_cargo_energy_use_actual_heos_and_gerg(self):
        composition = u.GasComposition({"CH4": Decimal(".9"), "C2H6": Decimal(".05"), "N2": Decimal(".05")})
        state = u.GasReferenceConditions.from_standard("iso15", compressibility_factor=None)
        liquid = u.LiquidState(u.Q_(110, "K"), u.Q_(150000, "Pa"))
        methane = u.GasComposition({"CH4": 1})
        for method in ("heos", "gerg2008"):
            with self.subTest(backend=method):
                backend = u.builtin_gas_backend(method, temperature_range=(u.Q_(90, "K"), u.Q_(400, "K")),
                                                pressure_range=(u.Q_(50000, "Pa"), u.Q_(5000000, "Pa")))
                properties = u.gas_properties(composition, state, backend=backend, air_compressibility=1)
                self.assertTrue(Decimal(".7") < properties.density.magnitude < Decimal(".8"))
                density = backend.liquid_density(methane, liquid)
                self.assertTrue(420 < density.magnitude < 430)
                warmer = backend.liquid_density(methane, replace(liquid, temperature=u.Q_(111, "K")))
                self.assertLess(warmer, density)
                rich = backend.liquid_density(u.GasComposition({"CH4": Decimal(".9"), "C2H6": Decimal(".1")}), liquid)
                self.assertNotEqual(rich, density)
                with self.assertRaisesRegex(ValueError, "not a single liquid"):
                    backend.liquid_density(composition, liquid)
                dew = u.gas_dew_point(composition, u.Q_(101325, "Pa"), kind="hydrocarbon", backend=backend)
                self.assertTrue(140 < dew.temperature.magnitude < 147)
                cargo = u.LNGCargo("test cargo", methane, liquid, density_backend=backend)
                volume = u.Q_(Decimal(100), "meter**3")
                vapor = cargo.gas_equivalent(volume, state, backend=backend, target_unit="Sm3")
                self.assertAlmostEqual(float(cargo.liquid_equivalent(vapor, state, backend=backend).magnitude), 100., places=9)
                gas_energy = u.composition_gas_energy(vapor, methane, state, backend=backend, basis="HHV", target_unit="MMBtu_it")
                self.assertAlmostEqual(float(gas_energy.magnitude / cargo.energy(volume, basis="HHV").magnitude), 1., places=11)
                with self.assertRaisesRegex(ValueError, "not a single gas"):
                    u.gas_density(methane, replace(state, temperature=liquid.temperature, absolute_pressure=liquid.absolute_pressure), backend=backend)


if __name__ == "__main__":
    unittest.main()
