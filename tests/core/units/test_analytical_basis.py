"""Tests for the single universal moisture conversion engine and assay delegation."""

from __future__ import annotations

from decimal import Decimal as D
import unittest
from unittest.mock import patch

import pint

import core.units as u
from core.units import agricultural_specs as ag
from core.units import conversions as conv
from core.units import moisture as m
from core.units.registry import ForeignRegistryError, Q_


class UniversalMoistureFormulaTests(unittest.TestCase):
    def test_wet_to_dry_and_two_moisture_level_identities(self) -> None:
        # 1. m_dry = m_wet * (1 - w)
        w = D("0.1234567890123456789")
        wet_q = Q_(D("1000"), "wmt")
        dry_q = m.wet_to_dry_mass(wet_q, w)
        self.assertEqual(dry_q.units, u.ureg.dmt)
        self.assertEqual(dry_q.magnitude, D("1000") * (D("1") - w))

        restored_wet = m.dry_to_wet_mass(dry_q, w)
        self.assertEqual(restored_wet.units, u.ureg.wmt)
        self.assertAlmostEqual(restored_wet.magnitude, D("1000"), places=22)

        # 2. m_2 = m_1 * (1 - w_1) / (1 - w_2) & transitive dry-matter conservation
        w1, w2, w3 = D("0.22"), D("0.14"), D("0.08")
        m1 = Q_(D("500"), "metric_ton")
        m2 = m.convert_moisture_mass(m1, source_moisture=w1, target_moisture=w2)
        m3 = m.convert_moisture_mass(m2, source_moisture=w2, target_moisture=w3)
        m1_to_3 = m.convert_moisture_mass(m1, source_moisture=w1, target_moisture=w3)

        # Dry matter m_i * (1 - w_i) is strictly invariant across all moisture levels
        dm1 = m1.magnitude * (D("1") - w1)
        dm2 = m2.magnitude * (D("1") - w2)
        dm3 = m3.magnitude * (D("1") - w3)
        self.assertAlmostEqual(dm1, dm2, places=22)
        self.assertAlmostEqual(dm2, dm3, places=22)
        self.assertAlmostEqual(m3.magnitude, m1_to_3.magnitude, places=22)

        # Intensive property conservation: v_1 * m_1 == v_2 * m_2
        v1 = D("240.50")
        v2 = m.rebase_moisture_intensive(v1, source_moisture=w1, target_moisture=w2)
        self.assertAlmostEqual(v1 * m1.magnitude, v2 * m2.magnitude, places=20)

    def test_validation_and_semantic_tag_guards(self) -> None:
        for bad in (-D("0.01"), D("1.0"), D("1.5"), float("nan"), float("inf"), True):
            with self.subTest(bad=bad), self.assertRaises((ValueError, TypeError)):
                m.validate_moisture_fraction(bad)  # type: ignore[arg-type]

        # Reject already-dry semantic unit as wet_to_dry source
        with self.assertRaises(pint.DimensionalityError):
            m.wet_to_dry_mass(Q_(D("100"), "dmt"), moisture_pct=D("10"))
        # Reject already-wet semantic unit as dry_to_wet source
        with self.assertRaises(pint.DimensionalityError):
            m.dry_to_wet_mass(Q_(D("100"), "wmt"), moisture_pct=D("10"))

        foreign = pint.UnitRegistry().Quantity(100, "kg")
        with self.assertRaises(ForeignRegistryError):
            m.wet_to_dry_mass(foreign, moisture_pct=10)


class UniversalAssayDelegationTests(unittest.TestCase):
    def test_every_agricultural_biomass_coal_and_ore_assay_uses_universal_moisture(self) -> None:
        """Prove via spy wrapping that every domain assay delegates to core.units.moisture."""
        coal = conv.CoalAssay.newcastle_6000()
        biomass = conv.BiomassAssay.forestry_wood_chips()
        petcoke = conv.PetcokeAssay.usgc_fuel_grade_65s()
        ore = conv.OreAssay.iron_ore_62fe()
        battery = conv.BatteryMaterialAssay.spodumene_sc6()
        bulk_ore = conv.BulkOreAssay("iron", D("0.08"), {"Fe": D("0.62")})
        grain = ag.GrainAssay.corn_56lb(actual_moisture_pct=20.0)
        pulp = ag.PulpLumberAssay.nbsk_softwood_pulp()
        rubber = ag.RubberAssay.field_latex_60drc()
        honey = ag.HoneyApicultureAssay.bulk_polyfloral_export()
        wheat_q = u.WheatQualitySpec(moisture_pct=D("15.0"), protein_pct=D("12.5"))
        multi_conc = u.MultiMetalConcentrateAssay(
            name="Cu",
            moisture_fraction=D("0.09"),
            composition=u.ElementComposition({"Cu": D("0.28")}),
        )

        cases = [
            ("CoalAssay.wmt_to_dmt", lambda: coal.wmt_to_dmt(100.0), 90.0),
            ("CoalAssay.dmt_to_wmt", lambda: coal.dmt_to_wmt(90.0), 100.0),
            ("CoalAssay.mass_at_moisture", lambda: coal.mass_at_moisture(Q_(D("100"), "mt"), D("25")).magnitude, D("120")),
            ("BiomassAssay.wmt_to_odt", lambda: biomass.wmt_to_odt(Q_(D("100"), "mt")).magnitude, D("65")),
            ("BiomassAssay.odt_to_wmt", lambda: biomass.odt_to_wmt(Q_(D("65"), "mt")).magnitude, D("100")),
            ("BiomassAssay.mass_at_moisture", lambda: biomass.mass_at_moisture(Q_(D("100"), "mt"), D("18.75")).magnitude, D("80")),
            ("PetcokeAssay.wmt_to_dmt", lambda: petcoke.wmt_to_dmt(Q_(D("100"), "wmt")).magnitude, D("92")),
            ("PetcokeAssay.dmt_to_wmt", lambda: petcoke.dmt_to_wmt(Q_(D("92"), "dmt")).magnitude, D("100")),
            ("OreAssay.wmt_to_dmt", lambda: ore.wmt_to_dmt(Q_(D("100"), "wmt")).magnitude, D("92")),
            ("OreAssay.dmt_to_wmt", lambda: ore.dmt_to_wmt(Q_(D("92"), "dmt")).magnitude, D("100")),
            ("OreAssay.mass_at_moisture", lambda: ore.mass_at_moisture(Q_(D("100"), "mt"), D("20")).magnitude, D("115")),
            ("BatteryMaterialAssay.wmt_to_dmt", lambda: battery.wmt_to_dmt(Q_(D("100"), "wmt")).magnitude, D("92")),
            ("BulkOreAssay.wmt_to_dmt", lambda: bulk_ore.wmt_to_dmt(D("100")), D("92")),
            ("convert_wet_to_dry_mass", lambda: conv.convert_wet_to_dry_mass(Q_(D("100"), "wmt"), D("8")).magnitude, D("92")),
            ("convert_dry_to_wet_mass", lambda: conv.convert_dry_to_wet_mass(Q_(D("92"), "dmt"), D("8")).magnitude, D("100")),
            ("GrainAssay.apply_moisture_deduction", lambda: grain.apply_moisture_deduction(Q_(D("84.5"), "mt")).magnitude, D("80")),
            ("GrainAssay.mass_at_moisture", lambda: grain.mass_at_moisture(Q_(D("80"), "mt"), D("20"), source_moisture_pct=D("15.5")).magnitude, D("84.5")),
            ("PulpLumberAssay.mass_to_bdmt", lambda: pulp.mass_to_bdmt(Q_(D("100"), "mt")).magnitude, D("90")),
            ("PulpLumberAssay.mass_to_admt", lambda: pulp.mass_to_admt(Q_(D("90"), "mt"), source_moisture_pct=D("0")).magnitude, D("100")),
            ("RubberAssay.wet_to_dry_mass", lambda: rubber.wet_to_dry_mass(Q_(D("100"), "mt")).magnitude, D("60")),
            ("HoneyApicultureAssay.normalise", lambda: honey.normalise_to_base_moisture(Q_(D("82"), "mt")).magnitude, D("82.5")),
            ("WheatQualitySpec.normalized_mass", lambda: wheat_q.normalized_mass(Q_(D("86.5"), "mt"), base_moisture_pct=D("13.5")).magnitude, D("85")),
            ("MultiMetalConcentrateAssay.wet_to_dry", lambda: multi_conc.wet_to_dry_mass(Q_(D("100"), "kg")).magnitude, D("91")),
        ]

        orig_convert = m.convert_moisture_mass
        for label, fn, expected in cases:
            with self.subTest(case=label):
                with patch("core.units.moisture.convert_moisture_mass", wraps=orig_convert) as spy:
                    # Also patch direct imported references in consumer modules to point to the spy
                    with (
                        patch.object(conv, "_universal_convert_moisture_mass", spy),
                        patch.object(ag, "_universal_convert_moisture_mass", spy),
                        patch.object(u.agricultural_quality, "_universal_convert_moisture_mass", spy),
                    ):
                        actual = fn()
                        self.assertGreaterEqual(spy.call_count, 1, f"{label} did not invoke convert_moisture_mass!")
                        self.assertAlmostEqual(float(actual), float(expected), places=9)


if __name__ == "__main__":
    unittest.main()
