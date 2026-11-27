"""Gas conversions agree across legacy, shared and master routes."""

from decimal import Decimal
import unittest

import core.units as u
from core.units import conversions as c


STATES = ("Nm3", "Sm3", "scf", "scf_tx", "scf_60F_14696psi_dry", "Sm3_20C", "m3_Groningen")


class GasRouteParityTests(unittest.TestCase):
    def assertRelativeClose(self, actual, expected):
        self.assertAlmostEqual(float(actual.magnitude) / float(expected.to(actual.units).magnitude), 1., places=11)

    def test_every_reference_amount_flow_and_quote_in_both_numeric_types(self):
        for density in (None, .62):
            gas = u.GasAssay(gcv_btu_per_scf=1030, ncv_btu_per_scf=930, relative_density=density)
            heat = u.heat_content_from_assay(gas)
            for basis in ("gcv", "ncv"):
                context = c.create_gas_heating_context(gas, basis=basis)
                for state in STATES:
                    for magnitude in (1000., Decimal(1000)):
                        with self.subTest(density=density, state=state, basis=basis, numeric=type(magnitude)):
                            quantity = u.Q_(magnitude, state)
                            with u.ureg.context(context):
                                expected = quantity.to("MWh")
                            for route in (
                                lambda q, target: u.convert(q, target, assay=gas, basis=basis),
                                lambda q, target: u.convert(q, target, assay=gas, heat_content=heat, basis=basis),
                                lambda q, target: u.convert_heat_content(q, target, assay=heat, basis=basis),
                            ):
                                energy = route(quantity, "MWh")
                                self.assertRelativeClose(energy, expected)
                                self.assertRelativeClose(route(energy, state), quantity)
                                self.assertIsInstance(energy.magnitude, type(magnitude))
                                flow = route(u.Q_(magnitude, state + "/day"), "MW")
                                self.assertRelativeClose(flow, u.Q_(expected.magnitude / 24, "MW"))
                                quote = u.Q_(magnitude / 20, "EUR/MWh")
                                price = route(quote, "EUR/" + state)
                                factor = expected.magnitude / magnitude
                                self.assertRelativeClose(price, u.Q_(quote.magnitude * factor, "EUR/" + state))
                                self.assertRelativeClose(route(price, "EUR/MWh"), quote)
                                self.assertIsInstance(price.magnitude, type(magnitude))

    def test_specific_energy_between_all_reference_pairs_preserves_supplied_value(self):
        for density in (None, .62):
            gas = u.GasAssay(gcv_btu_per_scf=1030, ncv_btu_per_scf=930, relative_density=density)
            heat = u.heat_content_from_assay(gas)
            for basis in ("gcv", "ncv"):
                context = c.create_gas_heating_context(gas, basis=basis)
                for source in STATES:
                    with u.ureg.context(context):
                        source_heat = u.Q_(Decimal(1), source).to("MJ").magnitude
                    self.assertRelativeClose(heat.energy_density("MJ/" + source, kind=basis), u.Q_(source_heat, "MJ/" + source))
                    for target in STATES:
                        with u.ureg.context(context):
                            target_volume = u.Q_(Decimal(1), target).to(source).magnitude
                        for magnitude in (Decimal(40), 40.):
                            with self.subTest(source=source, target=target, basis=basis, density=density, numeric=type(magnitude)):
                                original = u.Q_(magnitude, "MJ/" + source)
                                result = u.convert(original, "MJ/" + target, assay=gas, basis=basis)
                                self.assertRelativeClose(result, u.Q_(Decimal(str(magnitude)) * target_volume, "MJ/" + target))
                                restored = u.convert(result, original.units, assay=heat, basis=basis)
                                self.assertRelativeClose(restored, original)
                                self.assertIsInstance(result.magnitude, type(magnitude))

    def test_mass_and_volume_energy_density_routes_match(self):
        gas = u.GasAssay.from_composition({"CH4": Decimal(".9"), "N2": Decimal(".1")})
        heat = u.heat_content_from_assay(gas)
        for basis in ("gcv", "ncv"):
            mass_heat = heat.energy_density("MJ/kg", kind=basis)
            for state in STATES:
                with self.subTest(basis=basis, state=state):
                    volume_heat = heat.energy_density("MJ/" + state, kind=basis)
                    converted = u.convert(volume_heat, "MJ/kg", assay=heat, basis=basis)
                    self.assertRelativeClose(converted, mass_heat)

    def test_absent_reference_mapping_still_requires_explicit_density_state(self):
        heat = u.HeatContentAssay(ncv=u.Q_(40, "MJ/Nm3"))
        with self.assertRaisesRegex(ValueError, "metadata"):
            heat.energy_density("MJ/Sm3")


if __name__ == "__main__":
    unittest.main()
