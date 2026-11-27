"""Petroleum measurement identities, applicability and conservation contracts.

Engineering fixtures are synthetic measured-coefficient examples, not normative
API/ASTM checkpoints or validation of a standards implementation.
"""

from dataclasses import replace
from datetime import date
from decimal import Decimal
import math
import unittest

import pint

import core.units as u
from core.units import conversions


D = Decimal


class PetroleumMeasurementTests(unittest.TestCase):
    def setUp(self):
        self.reference = u.LiquidState(u.Q_(D(15), "degC"), u.Q_(D(101325), "Pa"))
        self.observed = u.LiquidState(u.Q_(D(35), "degC"), u.Q_(D(2000000), "Pa"))
        self.spec = u.LiquidProductSpec("example crude", "crude", self.reference, u.Q_(D(800), "kg/meter**3"), D(".01"))
        self.bounds = u.LiquidApplicability((u.Q_(0, "degC"), u.Q_(80, "degC")),
                                           (u.Q_(10000, "Pa"), u.Q_(5000000, "Pa")),
                                           (u.Q_(500, "kg/meter**3"), u.Q_(1200, "kg/meter**3")))
        self.method = u.constant_coefficient_liquid_method(product="crude", thermal_expansion=u.Q_(D(".001"), "1/kelvin"),
                                                          compressibility=u.Q_(D("1e-9"), "1/Pa"),
                                                          applicability=self.bounds, source="synthetic laboratory example")

    def test_ctl_cpl_ctpl_have_correct_physical_signs_and_reference_identity(self):
        factors = self.method.factors(self.spec, self.observed)
        self.assertAlmostEqual(float(factors.ctl), math.exp(-.02), places=14)
        self.assertAlmostEqual(float(factors.cpl), math.exp(.001898675), places=14)
        self.assertLess(factors.ctl, 1)
        self.assertGreater(factors.cpl, 1)
        self.assertEqual(factors.ctpl, factors.ctl * factors.cpl)
        baseline = self.method.factors(self.spec, self.reference)
        self.assertEqual((baseline.ctl, baseline.cpl, baseline.ctpl), (1, 1, 1))

    def test_volume_and_density_preserve_mass_and_round_trip_types(self):
        for magnitude in (D("123.456789012345"), 123.456789012345):
            with self.subTest(kind=type(magnitude)):
                volume = u.Q_(magnitude, "oil_barrel")
                standard = u.correct_liquid_volume(volume, spec=self.spec, observed=self.observed, method=self.method)
                restored = u.correct_liquid_volume(standard, spec=self.spec, observed=self.observed, method=self.method, to_reference=False)
                self.assertIsInstance(standard.magnitude, type(magnitude))
                self.assertAlmostEqual(float(restored.magnitude), float(magnitude), places=10)
                observed_density = u.liquid_density_at_state(self.spec, self.observed, method=self.method)
                observed_mass = volume.to("meter**3").magnitude * (observed_density.magnitude if isinstance(magnitude, D) else float(observed_density.magnitude))
                reference_mass = standard.to("meter**3").magnitude * self.spec.reference_density.magnitude if isinstance(magnitude, D) else standard.to("meter**3").magnitude * float(self.spec.reference_density.magnitude)
                self.assertAlmostEqual(float(observed_mass), float(reference_mass), places=8)

    def test_reference_density_solver_checks_residual_and_never_extrapolates(self):
        observed = u.liquid_density_at_state(self.spec, self.observed, method=self.method)
        density = u.infer_liquid_reference_density(observed.to("g/liter"), name=self.spec.name, product="crude",
                                                  reference=self.reference, observed=self.observed, method=self.method)
        self.assertAlmostEqual(float(density.magnitude), 800, places=8)
        self.assertIsInstance(density.magnitude, D)
        with self.assertRaisesRegex(ValueError, "convergence"):
            u.infer_liquid_reference_density(observed, name="test", product="crude", reference=self.reference,
                                            observed=self.observed, method=self.method, max_iterations=1)
        with self.assertRaisesRegex(ValueError, "bracketed"):
            u.infer_liquid_reference_density(u.Q_(2000, "kg/meter**3"), name="test", product="crude", reference=self.reference,
                                            observed=self.observed, method=self.method)

    def test_method_bounds_cover_reference_state_observed_state_and_density(self):
        for observed in (replace(self.observed, temperature=u.Q_(81, "degC")), replace(self.observed, absolute_pressure=u.Q_(6000000, "Pa"))):
            with self.subTest(observed=observed), self.assertRaisesRegex(ValueError, "applicability"):
                self.method.factors(self.spec, observed)
        for spec in (replace(self.spec, reference_density=u.Q_(1300, "kg/meter**3")),
                     replace(self.spec, reference=replace(self.reference, temperature=u.Q_(-1, "degC")))):
            with self.subTest(spec=spec), self.assertRaisesRegex(ValueError, "applicability"):
                self.method.factors(spec, self.observed)
        for temperature in (0, 80):
            self.method.factors(self.spec, replace(self.observed, temperature=u.Q_(temperature, "degC")))

    def test_product_specific_algorithms_and_validation_metadata_are_explicit(self):
        with self.assertRaisesRegex(ValueError, "does not cover"):
            self.method.factors(replace(self.spec, product="gasoline"), self.observed)
        for product in ("lpg", "ngl", "denatured_ethanol", "ethanol_gasoline"):
            with self.subTest(product=product), self.assertRaisesRegex(ValueError, "dedicated"):
                u.constant_coefficient_liquid_method(product=product, thermal_expansion=u.Q_(.001, "1/K"),
                                                     compressibility=u.Q_(1e-9, "1/Pa"), applicability=self.bounds, source="test")
        with self.assertRaisesRegex(ValueError, "validation reference"):
            replace(self.method, standard="Example standard, edition 1")
        wrapped = replace(self.method, standard="Example standard, edition 1", validation_reference="synthetic test fixture only")
        self.assertEqual(wrapped.factors(self.spec, self.observed), self.method.factors(self.spec, self.observed))
        with self.assertRaisesRegex(ValueError, "non-unity"):
            replace(self.method, calculator=lambda spec, observed: u.LiquidCorrectionFactors(D(".9"), D(1))).factors(self.spec, self.reference)
        with self.assertRaises(TypeError):
            replace(self.method, calculator=lambda spec, observed: (.9, 1)).factors(self.spec, self.observed)

    def test_lpg_ngl_require_their_own_method_and_saturation_states(self):
        reference = u.LiquidState(u.Q_(15, "degC"), u.Q_(1000000, "Pa"), u.Q_(700000, "Pa"))
        observed = u.LiquidState(u.Q_(35, "degC"), u.Q_(2000000, "Pa"), u.Q_(1500000, "Pa"))
        for product in ("lpg", "ngl"):
            spec = replace(self.spec, product=product, reference=reference)
            method = replace(self.method, method_id="example dedicated volatile-liquid backend", products=frozenset({product}))
            self.assertGreater(method.factors(spec, observed).ctpl, 0)
            with self.assertRaisesRegex(ValueError, "saturation pressure"):
                method.factors(spec, replace(observed, saturation_pressure=None))
            with self.assertRaisesRegex(ValueError, "saturation pressure"):
                replace(spec, reference=self.reference)
        with self.assertRaisesRegex(ValueError, "below"):
            u.LiquidState(u.Q_(35, "degC"), u.Q_(1000000, "Pa"), u.Q_(1500000, "Pa"))

    def test_alcohol_methods_require_known_composition_and_composition_bounds(self):
        with self.assertRaisesRegex(ValueError, "ethanol volume fraction"):
            replace(self.spec, product="ethanol_gasoline")
        spec = replace(self.spec, product="ethanol_gasoline", ethanol_volume_fraction=D(".1"))
        method = replace(self.method, method_id="example E10 laboratory backend", products=frozenset({"ethanol_gasoline"}))
        with self.assertRaisesRegex(ValueError, "composition applicability"):
            method.factors(spec, self.observed)
        method = replace(method, applicability=replace(self.bounds, ethanol_fraction_range=(D(".05"), D(".15"))))
        self.assertGreater(method.factors(spec, self.observed).ctpl, 0)
        with self.assertRaisesRegex(ValueError, "Ethanol fraction"):
            method.factors(replace(spec, ethanol_volume_fraction=D(".2")), self.observed)
        with self.assertRaisesRegex(ValueError, "dedicated ethanol"):
            replace(self.spec, ethanol_volume_fraction=D(".1"))

    def test_inventory_separates_free_water_and_suspended_sediment_water(self):
        result = u.calculate_petroleum_inventory(u.Q_(D(100), "meter**3"), free_water_volume=u.Q_(D(10000), "liter"),
                                                bsw_fraction=D(".02"), bsw_basis="standard_volume", spec=self.spec,
                                                observed=self.observed, method=self.method)
        self.assertEqual(result.tov.magnitude, 100)
        self.assertAlmostEqual(float(result.free_water.magnitude), 10, places=12)
        self.assertAlmostEqual(float(result.gov.magnitude), 90, places=12)
        self.assertEqual(result.gsv.magnitude, result.gov.magnitude * result.factors.ctpl)
        self.assertEqual(result.nsv.magnitude, result.gsv.magnitude * D(".98"))
        self.assertEqual(result.nsv + result.standard_sediment_and_water, result.gsv)
        self.assertEqual(result.net_physical_mass.magnitude, result.nsv.magnitude * 800)
        self.assertEqual((result.method_id, result.method_version), (self.method.method_id, self.method.version))

    def test_unknown_analytical_basis_or_invalid_inventory_never_defaults_to_zero(self):
        arguments = dict(free_water_volume=u.Q_(0, "meter**3"), bsw_fraction=D(0), bsw_basis="standard_volume",
                         spec=self.spec, observed=self.observed, method=self.method)
        for changes in (dict(free_water_volume=u.Q_(101, "meter**3")), dict(bsw_fraction=D(1)),
                        dict(bsw_fraction=D("NaN")), dict(bsw_fraction=True), dict(bsw_basis="observed_volume"), dict(bsw_basis="mass")):
            with self.subTest(changes=changes), self.assertRaises((TypeError, ValueError)):
                u.calculate_petroleum_inventory(u.Q_(100, "meter**3"), **{**arguments, **changes})
        for missing in ("free_water_volume", "bsw_fraction", "bsw_basis"):
            arguments_without = dict(arguments)
            del arguments_without[missing]
            with self.subTest(missing=missing), self.assertRaises(TypeError):
                u.calculate_petroleum_inventory(u.Q_(100, "meter**3"), **arguments_without)

    def test_air_buoyancy_round_trip_uses_air_and_reference_weight_densities(self):
        conditions = dict(liquid_density=u.Q_(800, "kg/meter**3"), air_density=u.Q_(D("1.2"), "kg/meter**3"),
                          weight_density=u.Q_(8000, "kg/meter**3"))
        factor = (1 - D("1.2") / 800) / (1 - D("1.2") / 8000)
        for magnitude in (D(1000), 1000.):
            mass = u.Q_(magnitude, "kg")
            apparent = u.physical_to_apparent_mass(mass, **conditions)
            physical = u.apparent_to_physical_mass(apparent, **conditions)
            self.assertIsInstance(apparent.magnitude, type(magnitude))
            self.assertAlmostEqual(float(apparent.magnitude), float(D(1000) * factor), places=10)
            self.assertAlmostEqual(float(physical.magnitude), 1000, places=10)
        self.assertEqual(u.physical_to_apparent_mass(u.Q_(D(1000), "kg"), **{**conditions, "air_density": u.Q_(0, "kg/meter**3")}).magnitude, 1000)
        with self.assertRaises(ValueError):
            u.physical_to_apparent_mass(u.Q_(1, "kg"), **{**conditions, "air_density": u.Q_(900, "kg/meter**3")})

    def test_liquid_states_reject_gauge_datum_temperature_intervals_and_nonfinite(self):
        for changes in (dict(temperature=u.Q_(20, "delta_degC")), dict(temperature=u.Q_(0, "kelvin")),
                        dict(temperature=u.Q_(float("nan"), "kelvin")), dict(absolute_pressure=u.Q_(10, "psig")),
                        dict(absolute_pressure=u.Q_(0, "Pa")), dict(absolute_pressure=u.Q_(float("inf"), "Pa"))):
            with self.subTest(changes=changes), self.assertRaises((ValueError, TypeError)):
                replace(self.reference, **changes)
        with self.assertRaises(pint.DimensionalityError):
            u.correct_liquid_volume(u.Q_(1, "kg"), spec=self.spec, observed=self.observed, method=self.method)

    def test_frozen_inputs_cannot_be_changed_through_in_place_quantity_conversion(self):
        original = self.method.factors(self.spec, self.observed)
        self.observed.temperature.ito("degF")
        self.observed.absolute_pressure.ito("bar")
        self.spec.reference_density.ito("kg/liter")
        self.bounds.pressure_range[1].ito("bar")
        self.assertEqual(self.method.factors(self.spec, self.observed), original)
        self.assertEqual(self.observed.temperature.units, u.ureg.kelvin)
        self.assertEqual(self.spec.reference_density.magnitude, 800)

    def test_master_router_round_trips_amounts_flows_and_prices_without_double_correction(self):
        parameters = dict(assay=self.spec, liquid_method=self.method, liquid_observed_state=self.observed)
        factor = self.method.factors(self.spec, self.observed).ctpl
        amount = u.convert(u.Q_(D(100), "meter**3"), "kg", **parameters)
        self.assertEqual(amount.magnitude, 100 * 800 * factor)
        self.assertAlmostEqual(float(u.convert(amount, "meter**3", **parameters).magnitude), 100, places=10)
        flow = u.convert(u.Q_(D(100), "meter**3/hour"), "kg/second", **parameters)
        self.assertAlmostEqual(float(flow.magnitude), float(D(100) * 800 * factor / 3600), places=10)
        price = u.convert(u.Q_(D(100), "USD/meter**3"), "EUR/kg", rates={"EUR": D(".9")}, **parameters)
        self.assertAlmostEqual(float(price.magnitude), float(D(90) / (800 * factor)), places=12)
        restored = u.convert(price, "USD/meter**3", rates={"EUR": D(".9")}, **parameters)
        self.assertAlmostEqual(float(restored.magnitude), 100, places=10)
        same_basis = u.convert(u.Q_(D(100), "liter"), "meter**3", **parameters)
        self.assertAlmostEqual(float(same_basis.magnitude), .1, places=12)

    def test_master_router_rejects_incompatible_or_half_supplied_correction_options(self):
        for parameters in (dict(liquid_method=self.method), dict(liquid_observed_state=self.observed),
                           dict(assay=self.spec, liquid_method=self.method, liquid_observed_state=self.observed, observed_temp=20),
                           dict(assay=u.NGLAssay.propane(), liquid_method=self.method, liquid_observed_state=self.observed)):
            with self.subTest(parameters=parameters), self.assertRaises((ValueError, TypeError)):
                u.convert(u.Q_(D(100), "meter**3"), "kg", **parameters)
        with self.assertRaisesRegex(ValueError, "volume endpoint"):
            u.convert(u.Q_(D(100), "kg"), "gram", assay=self.spec, liquid_method=self.method, liquid_observed_state=self.observed)

    def test_crude_bridge_separates_legacy_temperature_estimate_from_full_factor_path(self):
        crude = u.CrudeAssay.brent()
        spec = crude.as_liquid_spec()
        self.assertEqual(spec.reference_density.magnitude, D(str(crude.density_kg_per_m3)))
        factors = crude.correction_factors(self.observed, method=self.method)
        corrected = crude.correct_volume_to_60f(u.Q_(D(10), "oil_barrel"), u.Q_(35, "degC"),
                                               bsw_percentage=D(2), method=self.method, observed_pressure=u.Q_(2000000, "Pa"))
        self.assertEqual(corrected.magnitude, D(10) * factors.ctpl * D(".98"))
        result = u.convert(u.Q_(D(10), "meter**3"), "kg", assay=crude, liquid_method=self.method, liquid_observed_state=self.observed)
        self.assertAlmostEqual(float(result.magnitude), float(10 * factors.ctpl * spec.reference_density.magnitude), places=7)
        self.assertEqual(crude.volume_correction_factor(60), 1)
        self.assertAlmostEqual(crude.density_in_air_kg_per_m3, (crude.density_kg_per_m3 - 1.2) / (1 - 1.2 / 8000))
        with self.assertRaises(ValueError):
            crude.correct_volume_to_60f(u.Q_(1, "oil_barrel"), 60, observed_pressure=u.Q_(101325, "Pa"))
        with self.assertRaises(ValueError):
            crude.volume_correction_factor(-500)

    def test_distillate_requires_a_specific_correction_family(self):
        assay = u.CrudeAssay(density_kg_per_m3=820, product_group="distillate")
        with self.assertRaisesRegex(ValueError, "exact correction class"):
            assay.as_liquid_spec()
        self.assertEqual(assay.as_liquid_spec(product="jet").product, u.LiquidProduct.JET)
        method = replace(self.method, products=frozenset({"jet"}))
        self.assertGreater(assay.correction_factors(self.observed, method=method).ctpl, 0)
        with self.assertRaises(ValueError):
            assay.blend_with(u.CrudeAssay.brent(), self_volume_fraction=.5)


class TankCalibrationTests(unittest.TestCase):
    def setUp(self):
        points = tuple(u.TankCalibrationPoint(u.Q_(level, "meter"), u.Q_(volume, "meter**3"), u.Q_(D(".1"), "meter**3"))
                       for level, volume in ((0, 0), (1, 10), (2, 30)))
        self.table = u.TankCalibrationTable("tank-1", "bottom-datum", "certificate-1", u.Q_(20, "degC"), points,
                                            D(1), valid_from="2026-01-01", valid_to="2027-01-01")
        self.parameters = dict(datum_id="bottom-datum", level_standard_uncertainty=u.Q_(D(".01"), "meter"), measurement_date="2026-11-25")

    def test_interpolation_and_uncertainty_propagation_preserve_certificate_metadata(self):
        reading = self.table.volume_at(u.Q_(D(".5"), "meter"), **self.parameters)
        self.assertEqual(reading.volume.magnitude, 5)
        self.assertAlmostEqual(float(reading.standard_uncertainty.magnitude), math.sqrt(.1**2 + .1**2), places=12)
        self.assertEqual((reading.tank_id, reading.datum_id, reading.certificate_id), ("tank-1", "bottom-datum", "certificate-1"))
        independent = replace(self.table, endpoint_correlation=D(0)).volume_at(u.Q_(D(".5"), "meter"), **self.parameters)
        self.assertLess(independent.standard_uncertainty, reading.standard_uncertainty)

    def test_knot_slope_is_conservative_and_calibration_bounds_never_extrapolate(self):
        reading = self.table.volume_at(u.Q_(1, "meter"), **self.parameters)
        self.assertAlmostEqual(float(reading.standard_uncertainty.magnitude), math.sqrt(.1**2 + .2**2), places=12)
        for level in (-.1, 2.1):
            with self.subTest(level=level), self.assertRaises(ValueError):
                self.table.volume_at(u.Q_(level, "meter"), **self.parameters)
        with self.assertRaises(ValueError):
            self.table.volume_at(u.Q_(1, "meter"), **{**self.parameters, "datum_id": "ullage"})
        for when in (None, "2025-12-31", "2027-01-01"):
            with self.subTest(when=when), self.assertRaises(ValueError):
                self.table.volume_at(u.Q_(1, "meter"), **{**self.parameters, "measurement_date": when})

    def test_invalid_calibration_values_and_correlations_are_rejected(self):
        with self.assertRaises(ValueError):
            replace(self.table, endpoint_correlation=D("1.01"))
        with self.assertRaises(ValueError):
            replace(self.table, points=tuple(reversed(self.table.points)))
        with self.assertRaises(ValueError):
            u.TankCalibrationPoint(u.Q_(1, "meter"), u.Q_(10, "meter**3"), u.Q_(-1, "meter**3"))


class PetroleumBlendingTests(unittest.TestCase):
    def setUp(self):
        self.reference = u.LiquidState(u.Q_(60, "degF"), u.Q_(101325, "Pa"))
        self.light = u.LiquidProductSpec("light", "crude", self.reference, u.Q_(D(800), "kg/meter**3"), D(".005"))
        self.heavy = u.LiquidProductSpec("heavy", "crude", self.reference, u.Q_(D(1000), "kg/meter**3"), D(".03"))
        self.ideal = u.BlendVolumeModel("explicit additive volumes", "engineering assumption", u.ideal_blend_volume)
        self.components = (u.LiquidBlendComponent(self.light, u.Q_(D(1), "meter**3")),
                           u.LiquidBlendComponent(self.heavy, u.Q_(D(1000), "kg")))

    def test_mass_volume_sulfur_and_api_blending_conserve_their_actual_bases(self):
        result = u.blend_petroleum(self.components, volume_model=self.ideal)
        self.assertEqual(result.total_mass.magnitude, 1800)
        self.assertEqual(result.blend_volume.magnitude, 2)
        self.assertEqual(result.reference_density.magnitude, 900)
        self.assertEqual(result.sulfur_mass_fraction * result.total_mass.magnitude, D(34))
        api_expected = D("141.5") / (D(900) / D("999.016")) - D("131.5")
        self.assertEqual(result.api_gravity, api_expected)
        arithmetic_api = (u.CrudeAssay(density_kg_per_m3=800).api_gravity + u.CrudeAssay(density_kg_per_m3=1000).api_gravity) / 2
        self.assertNotAlmostEqual(float(result.api_gravity), arithmetic_api, places=5)

    def test_non_additive_volume_model_changes_density_without_changing_mass_or_sulfur(self):
        model = u.BlendVolumeModel("measured 1 percent contraction", "synthetic experiment",
                                   lambda components: u.ideal_blend_volume(components) * D(".99"))
        result = u.blend_petroleum(self.components, volume_model=model)
        self.assertEqual(result.blend_volume.magnitude, D("1.98"))
        self.assertEqual(result.total_mass.magnitude, 1800)
        self.assertGreater(result.reference_density.magnitude, 900)
        self.assertAlmostEqual(float(result.sulfur_mass_fraction * 1800), 34, places=12)
        with self.assertRaises(ValueError):
            u.blend_petroleum(self.components, volume_model=replace(model, calculator=lambda components: u.Q_(-1, "meter**3")))

    def test_nonlinear_viscosity_requires_a_common_temperature_and_explicit_method(self):
        components = (replace(self.components[0], kinematic_viscosity=u.Q_(D(10), "centistokes"), viscosity_temperature=u.Q_(50, "degC")),
                      replace(self.components[1], kinematic_viscosity=u.Q_(D(100), "centistokes"), viscosity_temperature=u.Q_(50, "degC")))
        with self.assertRaises(ValueError):
            u.blend_petroleum(components, volume_model=self.ideal)
        result = u.blend_petroleum(components, volume_model=self.ideal, viscosity_method="refutas_mass")
        reference_vbn = (800 * conversions._viscosity_blending_number(10) + 1000 * conversions._viscosity_blending_number(100)) / 1800
        self.assertAlmostEqual(float(result.kinematic_viscosity.magnitude), conversions._viscosity_from_vbn(reference_vbn), places=10)
        self.assertNotAlmostEqual(float(result.kinematic_viscosity.magnitude), (800 * 10 + 1000 * 100) / 1800)
        with self.assertRaises(ValueError):
            u.blend_petroleum((components[0], replace(components[1], viscosity_temperature=u.Q_(40, "degC"))), volume_model=self.ideal, viscosity_method="refutas_mass")
        with self.assertRaises(ValueError):
            u.blend_petroleum((replace(components[0], kinematic_viscosity=u.Q_(D(".2"), "centistokes")), components[1]), volume_model=self.ideal, viscosity_method="refutas_mass")

    def test_blending_never_guesses_missing_quality_reference_state_or_alcohol_contraction(self):
        unknown = replace(self.components[1], spec=replace(self.heavy, sulfur_mass_fraction=None))
        self.assertIsNone(u.blend_petroleum((self.components[0], unknown), volume_model=self.ideal).sulfur_mass_fraction)
        other_reference = replace(self.reference, temperature=u.Q_(15, "degC"))
        with self.assertRaises(ValueError):
            u.blend_petroleum((self.components[0], replace(self.components[1], spec=replace(self.heavy, reference=other_reference))), volume_model=self.ideal)
        ethanol = replace(self.light, product="denatured_ethanol", ethanol_volume_fraction=D(".95"))
        with self.assertRaisesRegex(ValueError, "non-ideal"):
            u.blend_petroleum((replace(self.components[0], spec=ethanol), self.components[1]), volume_model=self.ideal)
        nonstandard = u.blend_petroleum((replace(self.components[0], spec=replace(self.light, reference=other_reference)),), volume_model=self.ideal)
        with self.assertRaisesRegex(ValueError, "60°F"):
            _ = nonstandard.api_gravity
        # Bubble point is composition-specific, not part of the common T/P basis.
        first = replace(self.light, product="lpg", reference=replace(self.reference, absolute_pressure=u.Q_(2000000, "Pa"), saturation_pressure=u.Q_(700000, "Pa")))
        second = replace(self.heavy, product="ngl", reference=replace(first.reference, saturation_pressure=u.Q_(1200000, "Pa")))
        result = u.blend_petroleum((u.LiquidBlendComponent(first, u.Q_(1, "meter**3")), u.LiquidBlendComponent(second, u.Q_(1, "meter**3"))), volume_model=self.ideal)
        self.assertIsNone(result.reference.saturation_pressure)
        self.assertIsNotNone(result.components[0].spec.reference.saturation_pressure)


if __name__ == "__main__":
    unittest.main()
