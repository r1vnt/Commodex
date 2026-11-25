"""Arbitrary gas states, water accounting and independent combustion metadata."""

from dataclasses import FrozenInstanceError, replace
from decimal import Decimal
import unittest
from unittest.mock import patch

import pint

import core.units as units
from core.units import conversions as c


class GasReferenceConditionTests(unittest.TestCase):
    def state(self, temperature=300, pressure=200, z="0.8", water="0"):
        return c.GasReferenceConditions(units.Q_(Decimal(temperature), "K"), units.Q_(Decimal(pressure), "kPa"), Decimal(z), Decimal(water))

    def test_arbitrary_nonideal_states_and_numeric_round_trips(self):
        source, target = self.state(), self.state(360, 100, ".96")
        for amount in (100., Decimal(100)):
            with patch.object(c, "calculate_z_factor", side_effect=AssertionError("No gas model may be guessed")):
                original = units.Q_(amount, "meter**3")
                result = c.convert_gas_reference_volume(original, "liter", source=source, target=target)
                restored = c.convert_gas_reference_volume(result, original.units, source=target, target=source)
            self.assertAlmostEqual(float(result.magnitude), 288000)
            self.assertAlmostEqual(float(restored.magnitude), 100)
            self.assertIsInstance(result.magnitude, type(amount))
            self.assertIsInstance(restored.magnitude, type(amount))

    def test_wet_dry_conservation_uses_nonwater_moles(self):
        wet, dry = self.state(z="1", water=".02"), self.state(z="1")
        volume = units.Q_(Decimal(100), "meter**3")
        result = c.convert_gas_reference_volume(volume, "meter**3", source=wet, target=dry)
        self.assertEqual(result.magnitude, Decimal(98))
        self.assertEqual(c.convert_gas_reference_volume(result, volume.units, source=dry, target=wet).magnitude, 100)
        target_wet = self.state(z="1", water=".01")
        result = c.convert_gas_reference_volume(volume, volume.units, source=wet, target=target_wet)
        self.assertEqual(result.magnitude, Decimal(100) * Decimal(".98") / Decimal(".99"))
        with self.assertRaisesRegex(ValueError, "wet/dry"):
            c.convert_gas_reference_volume(volume, volume.units, source=wet, target=dry, conserve="total_gas")

    def test_unknown_water_and_compressibility_are_not_ideal_dry_defaults(self):
        known = self.state()
        unknown_water = replace(known, water_mole_fraction=None)
        unknown_z = replace(known, compressibility_factor=None)
        volume = units.Q_(1., "meter**3")
        with self.assertRaisesRegex(ValueError, "water_mole_fraction"):
            c.convert_gas_reference_volume(volume, volume.units, source=unknown_water, target=known)
        with self.assertRaisesRegex(ValueError, "compressibility_factor"):
            c.convert_gas_reference_volume(volume, volume.units, source=unknown_z, target=known)
        fixed_moles = c.convert_gas_reference_volume(volume, volume.units, source=unknown_water, target=unknown_water, conserve="total_gas")
        self.assertEqual(fixed_moles.magnitude, 1)

    def test_absolute_states_and_field_validation(self):
        valid = self.state()
        for field, values in {
            "temperature": (units.Q_(0, "K"), units.Q_(-1, "K"), units.Q_(10, "delta_degC"), units.Q_(float("nan"), "K")),
            "absolute_pressure": (units.Q_(1, "psig"), units.Q_(1, "kpag"), units.Q_(0, "Pa"), units.Q_(float("inf"), "Pa")),
            "compressibility_factor": (Decimal(0), Decimal(-1), Decimal("NaN"), True),
            "water_mole_fraction": (Decimal(1), Decimal(-1), Decimal("Infinity"), True),
            "standard_id": ("", 123),
        }.items():
            for value in values:
                with self.subTest(field=field, value=value), self.assertRaises((ValueError, TypeError, pint.PintError)):
                    replace(valid, **{field: value})
        offset_state = replace(valid, temperature=units.Q_(0, "degC"))
        self.assertEqual(offset_state.temperature.magnitude, Decimal("273.15"))
        with self.assertRaises(FrozenInstanceError):
            valid.standard_id = "modified"
        original_temperature = units.Q_(Decimal(300), "K")
        state = replace(valid, temperature=original_temperature)
        original_temperature.ito("degC")
        self.assertEqual(state.temperature.magnitude, 300)

    def test_combustion_temperature_is_independent_and_not_a_metering_factor(self):
        metered = self.state(273, 100, "1")
        combustion_15 = replace(metered, combustion_reference_temperature=units.Q_(15, "degC"))
        combustion_25 = replace(metered, combustion_reference_temperature=units.Q_(25, "degC"))
        volume = units.Q_(Decimal(10), "meter**3")
        result = c.convert_gas_reference_volume(volume, volume.units, source=combustion_15, target=combustion_25)
        self.assertEqual(result.magnitude, 10)
        self.assertNotEqual(combustion_15.temperature, combustion_15.combustion_reference_temperature)

    def test_named_standard_bridges_and_contradictions(self):
        arbitrary = self.state(320, 250, ".95")
        for key, (tagged, geometry, _, _) in c._GAS_REFERENCE_PRESETS.items():
            with self.subTest(standard=key):
                named = c.GasReferenceConditions.from_standard(key, compressibility_factor=Decimal(".99"))
                original = units.Q_(Decimal(100), geometry)
                tagged_result = c.convert_gas_reference_volume(original, tagged, source=arbitrary, target=named)
                restored = c.convert_gas_reference_volume(tagged_result, geometry, source=named, target=arbitrary)
                self.assertAlmostEqual(float(restored.magnitude), 100, places=10)
                with self.assertRaisesRegex(ValueError, "contradicts"):
                    c.convert_gas_reference_volume(original, tagged, source=arbitrary, target=replace(named, temperature=units.Q_(999, "K")))
        with self.assertRaisesRegex(ValueError, "Groningen"):
            c.GasReferenceConditions.from_standard("groningen")
        with self.assertRaisesRegex(ValueError, "calorific"):
            c.convert_gas_reference_volume(units.Q_(1, "m3_Groningen"), "meter**3", source=arbitrary, target=arbitrary)

    def test_invalid_volume_and_mutated_quantity_members_are_rejected(self):
        source = self.state()
        for quantity in (units.Q_(-1, "meter**3"), units.Q_(float("nan"), "meter**3"), units.Q_(1, "kg")):
            with self.assertRaises(ValueError):
                c.convert_gas_reference_volume(quantity, "meter**3", source=source, target=source)
        source.temperature._magnitude = Decimal(-1)
        with self.assertRaises(ValueError):
            c.convert_gas_reference_volume(units.Q_(1, "meter**3"), "meter**3", source=source, target=source)

    def test_public_exports_and_display_format_do_not_change_state_parsing(self):
        self.assertIs(units.GasReferenceConditions, c.GasReferenceConditions)
        self.assertIs(units.convert_gas_reference_volume, c.convert_gas_reference_volume)
        previous = units.ureg.formatter.default_format
        try:
            units.ureg.formatter.default_format = "L"
            self.assertEqual(self.state().temperature.magnitude, 300)
        finally:
            units.ureg.formatter.default_format = previous

    def test_master_router_requires_complete_unambiguous_reference_request(self):
        source, target = self.state(), self.state(360, 100, ".96")
        original = units.Q_(Decimal(100), "meter**3")
        routed = units.convert(original, "meter**3", source_reference=source, target_reference=target)
        self.assertEqual(routed.magnitude, 288)
        with self.assertRaisesRegex(ValueError, "Both"):
            units.convert(original, "meter**3", source_reference=source)
        with self.assertRaisesRegex(ValueError, "cannot be combined"):
            units.convert(original, "meter**3", source_reference=source, target_reference=target, assay=c.GasAssay(gcv_btu_per_scf=1000))
        with self.assertRaisesRegex(ValueError, "requires"):
            units.convert(original, "meter**3", gas_conservation="total_gas")
