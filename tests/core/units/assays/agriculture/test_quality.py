"""Unit tests for physical agricultural quality specifications and normalisation."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from decimal import Decimal as D
import unittest

import pint

import core.units as u
from core.units import agricultural_specs as ag
from core.units.agricultural_quality import (
    BarleyQualitySpec,
    CocoaQualitySpec,
    CoffeeQualitySpec,
    CornQualitySpec,
    CottonQualitySpec,
    MilkQualitySpec,
    QualityThreshold,
    RiceQualitySpec,
    RubberQualitySpec,
    SoybeanQualitySpec,
    SugarQualitySpec,
    TimberQualitySpec,
    VegetableOilQualitySpec,
    WheatQualitySpec,
    WoolQualitySpec,
    normalize_mass_to_moisture,
    quality_spec_from_assay,
    rebase_constituent_pct,
)
from core.units.materials import MaterialSpec, material_spec_from_assay
from core.units.registry import ForeignRegistryError, Q_


class GrainAndOilseedQualityTests(unittest.TestCase):
    def test_wheat_protein_basis_pirten_falling_number_and_test_weight(self) -> None:
        # 12.5% protein on 12.0% US FGIS moisture basis -> dry basis = 12.5 * 100 / 88 = 14.204545...%
        hrw = WheatQualitySpec.from_winchester_test_weight(
            moisture_pct=D("13.5"),
            protein_pct=D("12.5"),
            protein_moisture_basis_pct=D("12.0"),
            falling_number_s=D("350"),
            wet_gluten_pct=D("28.0"),
            dry_gluten_pct=D("10.0"),
            gluten_index=D("90.0"),
            test_weight_lb_per_bu=D("60.0"),
        )
        self.assertAlmostEqual(hrw.protein_at_moisture_basis(D("0")), D("12.5") * D("100") / D("88"), places=18)
        self.assertAlmostEqual(hrw.test_weight_lb_per_winchester_bu, D("60.0"), places=16)
        self.assertGreater(hrw.test_weight_kg_per_hl, D("77.0"))

        # Blend sprout-damaged wheat (FN=150s => LN=60) with sound wheat (FN=350s => LN=20)
        # Equal dry mass => LN_blend = 40 => FN_blend = 6000/40 + 50 = 200s (strictly < 250s arithmetic mean!)
        sprouted = WheatQualitySpec(
            moisture_pct=D("13.5"),
            protein_pct=D("11.5"),
            protein_moisture_basis_pct=D("12.0"),
            falling_number_s=D("150"),
            test_weight=hrw.test_weight,
        )
        blended = hrw.blend_with(sprouted, self_mass_fraction=D("0.5"))
        self.assertAlmostEqual(blended.falling_number_s, D("200"), places=16)

        # Benchmark comparison aligns protein moisture basis automatically (CWRS 13.5% MB vs US 12.0% MB)
        cwrs_bench = WheatQualitySpec(
            moisture_pct=D("13.5"),
            protein_pct=hrw.protein_at_moisture_basis(D("13.5")),
            protein_moisture_basis_pct=D("13.5"),
        )
        deltas = hrw.compare_to_benchmark(cwrs_bench)
        self.assertAlmostEqual(deltas["protein_pct"], D("0"), places=18)

    def test_corn_bcfm_aflatoxin_and_clean_mass_normalisation(self) -> None:
        yc2 = CornQualitySpec(
            moisture_pct=D("16.5"),
            broken_kernels_pct=D("2.2"),
            foreign_material_pct=D("1.3"),
            damaged_kernels_pct=D("4.0"),
            heat_damaged_pct=D("0.1"),
            aflatoxin_ppb=D("15.0"),
            test_weight_lb_per_bu=D("56.0"),
        )
        self.assertEqual(yc2.bcfm_pct, D("3.5"))
        self.assertTrue(yc2.meets_aflatoxin_limit(D("20")))
        self.assertFalse(yc2.meets_aflatoxin_limit(D("10")))

        # Deduct FM above 1.0% (0.3% dockage => 997 mt) then shrink from 16.5% to 15.5% moisture
        norm = yc2.normalized_clean_mass(
            Q_(D("1000"), "mt"),
            base_moisture_pct=D("15.5"),
            deductible_fm_threshold_pct=D("1.0"),
        )
        expected = D("1000") * D("0.997") * (D("83.5") / D("84.5"))
        self.assertAlmostEqual(norm.magnitude, expected, places=18)

    def test_soybeans_profat_moisture_rebasing_and_blending(self) -> None:
        brazil_soy = SoybeanQualitySpec(
            moisture_pct=D("13.5"),
            protein_pct=D("35.5"),
            oil_content_pct=D("19.5"),
            basis_moisture_pct=D("13.0"),
            damaged_beans_pct=D("2.0"),
            foreign_material_pct=D("1.0"),
        )
        self.assertEqual(brazil_soy.profat_pct(), D("55.0"))
        self.assertAlmostEqual(
            brazil_soy.profat_pct(target_moisture_pct=D("0")),
            D("55.0") * D("100") / D("87"),
            places=18,
        )
        oil_mt = brazil_soy.contained_oil_mass(Q_(D("1000"), "mt"))
        self.assertAlmostEqual(oil_mt.magnitude, D("1000") * (D("19.5") * D("86.5") / D("87.0")) / D("100"), places=16)

    def test_rice_milling_yield_and_broken_basis_consistency(self) -> None:
        # US rough rice 56/70 => 56% head rice of paddy, 70% total milled => 20% broken of milled
        rough = RiceQualitySpec(
            milling_recovery_pct=D("70.0"),
            head_rice_yield_pct=D("56.0"),
            head_rice_basis="of_paddy",
        )
        self.assertAlmostEqual(rough.broken_rice_pct, D("20.0"), places=16)
        self.assertAlmostEqual(rough.head_rice_pct_of_milled, D("80.0"), places=16)

        parts = rough.milling_breakdown(Q_(D("1000"), "mt"))
        self.assertEqual(parts["milled_total"].magnitude, D("700.0"))
        self.assertEqual(parts["head_rice"].magnitude, D("560.0"))
        self.assertEqual(parts["broken_rice"].magnitude, D("140.0"))
        self.assertEqual(parts["hull_and_bran_loss"].magnitude, D("300.0"))
        self.assertEqual(rough.required_paddy_for_head_rice(Q_(D("560"), "mt")).magnitude, D("1000"))

    def test_barley_malting_audit_and_extract_solids(self) -> None:
        malting = BarleyQualitySpec(
            moisture_pct=D("13.5"),
            protein_pct=D("10.8"),
            germination_energy_pct=D("97.0"),
            protein_dry_basis=True,
            plump_kernels_pct=D("92.0"),
            fine_extract_db_pct=D("81.0"),
        )
        audit = malting.meets_malting_thresholds(min_plump_kernels_pct=D("90.0"))
        self.assertTrue(audit.passed)
        self.assertEqual(audit.failed_attributes, ())

        feed_barley = BarleyQualitySpec(
            moisture_pct=D("15.0"),
            protein_pct=D("13.2"),
            germination_energy_pct=D("88.0"),
        )
        failed = feed_barley.meets_malting_thresholds()
        self.assertFalse(failed.passed)
        self.assertEqual(set(failed.failed_attributes), {"germination_energy_pct", "protein_db_pct", "moisture_pct"})
        self.assertEqual(failed.deviations["germination_energy_pct"], D("-7.0"))
        self.assertEqual(failed.deviations["protein_db_pct"], D("1.7"))

        # 100 mt @ 13.5% moisture = 86.5 mt dry * 81% extract = 70.065 mt extract solids
        self.assertEqual(malting.extractable_malt_solids_mass(Q_(D("100"), "mt")).magnitude, D("70.0650"))


class SoftsFibreDairyOilAndTimberQualityTests(unittest.TestCase):
    def test_coffee_cocoa_and_sugar_quality_normalisation(self) -> None:
        coffee = CoffeeQualitySpec(
            screen_size=D("17"),
            moisture_pct=D("11.5"),
            primary_defects_count=D("0"),
            secondary_defects_count=D("8"),
            processing_loss_pct=D("18.0"),
            water_activity=D("0.58"),
        )
        self.assertAlmostEqual(coffee.screen_aperture_mm, D("17") * D("25.4") / D("64"), places=18)
        self.assertEqual(coffee.full_defect_equivalents(), D("1.6"))
        self.assertEqual(coffee.green_to_roasted_mass(Q_(D("100"), "kg")).magnitude, D("82.0"))
        self.assertEqual(coffee.roasted_to_green_mass(Q_(D("82"), "kg")).magnitude, D("100"))

        cocoa = CocoaQualitySpec(
            bean_count_per_100g=D("92"),
            fermentation_pct=D("82"),
            fat_content_pct=D("54"),
            moisture_pct=D("7.5"),
            shell_pct=D("12.5"),
            defective_beans_pct=D("2.0"),
        )
        self.assertEqual(cocoa.nib_yield_pct, D("80.0"))
        self.assertEqual(cocoa.bean_to_nib_mass(Q_(D("1000"), "kg")).magnitude, D("800.0"))
        self.assertEqual(cocoa.extractable_butter_mass(Q_(D("1000"), "kg"), press_fat_recovery_pct=D("90")).magnitude, D("388.80"))

        vhp_sugar = SugarQualitySpec(
            polarisation_deg_z=D("99.3"),
            icumsa_colour_iu=D("600"),
            brix_pct=D("99.6"),
            moisture_pct=D("0.10"),
            ash_pct=D("0.12"),
        )
        self.assertAlmostEqual(vhp_sugar.apparent_purity_pct, D("99.3") / D("99.6") * D("100"), places=16)
        # 1,000 mt of 99.3° Pol VHP sugar = 1000 * 99.3 / 96.0 = 1,034.375 mt 96° Pol equivalent
        self.assertEqual(vhp_sugar.normalize_to_reference_pol(Q_(D("1000"), "mt")).magnitude, D("1034.375"))

    def test_cotton_and_rubber_quality_normalisation(self) -> None:
        cotton = CottonQualitySpec.from_staple_32nds(
            staple_32nds_inch=D("36"),  # 1-1/8 inch = 28.575 mm
            micronaire=D("4.0"),
            strength_gpt=D("30.5"),
            colour_grade="31",
            leaf_grade=3,
            moisture_regain_pct=D("7.0"),
        )
        self.assertEqual(cotton.staple_length_mm, D("28.575"))
        self.assertEqual(cotton.staple_inches, D("1.125"))
        # Bale weighed dry at 7.0% regain conditioned to standard 8.5% regain: 100 * 108.5 / 107.0
        cond = cotton.conditioned_invoice_mass(Q_(D("107"), "kg"), standard_regain_pct=D("8.5"))
        self.assertEqual(cond.magnitude, D("108.5"))

        # Harmonic micronaire blending: 50/50 blend of Mic 3.0 and Mic 6.0 -> 4.0 (not 4.5!)
        c_fine = CottonQualitySpec(staple_length_mm=D("28"), micronaire=D("3.0"), strength_gpt=D("29"))
        c_coarse = CottonQualitySpec(staple_length_mm=D("28"), micronaire=D("6.0"), strength_gpt=D("29"))
        self.assertEqual(c_fine.blend_with(c_coarse, self_mass_fraction=D("0.5")).micronaire, D("4"))

        # ISO 2000 TSR screening
        rubber = RubberQualitySpec(
            dirt_pct=D("0.12"),
            ash_pct=D("0.70"),
            volatile_matter_pct=D("0.60"),
            mooney_viscosity_ml1_4_100c=D("75"),
            nitrogen_pct=D("0.45"),
        )
        self.assertTrue(rubber.evaluate_iso2000_grade("TSR20").passed)
        tsr5_screen = rubber.evaluate_iso2000_grade("TSR5")
        self.assertFalse(tsr5_screen.passed)
        self.assertEqual(set(tsr5_screen.failed_attributes), {"dirt_pct", "ash_pct"})

    def test_milk_ecm_pearson_square_and_vegetable_oil_stoichiometry(self) -> None:
        whole = MilkQualitySpec(
            butterfat_pct=D("4.0"),
            protein_pct=D("3.3"),
            lactose_pct=D("4.8"),
        )
        skim = MilkQualitySpec(
            butterfat_pct=D("0.1"),
            protein_pct=D("3.5"),
            lactose_pct=D("5.0"),
            density_kg_per_m3=D("1035.0"),
        )
        self.assertEqual(whole.milk_solids_pct, D("7.3"))
        self.assertEqual(whole.solids_not_fat_pct, D("8.82"))
        # Gaines 4% FCM for 4.0% fat milk equals exact milk mass!
        self.assertEqual(whole.fat_corrected_milk_4pct(Q_(D("1000"), "kg")).magnitude, D("1000.00"))
        # Pearson Square standardisation to 2.05% fat (exact midpoint between 4.0% and 0.1%)
        w_whole, w_skim = whole.pearson_square_fat_weights(skim, target_butterfat_pct=D("2.05"))
        self.assertEqual(w_whole, D("0.5"))
        self.assertEqual(w_skim, D("0.5"))
        standardised = whole.blend_with(skim, self_mass_fraction=w_whole)
        self.assertEqual(standardised.butterfat_pct, D("2.05"))

        # Crude Palm Oil (palmitic FFA basis) vs Soybean Oil (oleic FFA basis)
        cpo = VegetableOilQualitySpec(
            free_fatty_acid_pct=D("4.5"),
            peroxide_value_meq_per_kg=D("5.0"),
            iodine_value_g_per_100g=D("52.0"),
            moisture_and_impurities_pct=D("0.20"),
            ffa_basis="palmitic",
        )
        # Acid value = 0.045 * 56105.6 / 256.42
        self.assertAlmostEqual(cpo.acid_value_mg_koh_per_g, D("0.045") * D("56105.6") / D("256.42"), places=16)
        self.assertAlmostEqual(cpo.ffa_on_basis("palmitic"), D("4.5"), places=16)

    def test_wool_iwto_blending_and_timber_ps20_nominal_vs_actual(self) -> None:
        fine_merino = WoolQualitySpec(
            micron_diameter_um=D("18.0"),
            vegetable_matter_pct=D("1.2"),
            clean_yield_pct=D("65.0"),
            staple_length_mm=D("85.0"),
        )
        broad_merino = WoolQualitySpec(
            micron_diameter_um=D("22.0"),
            vegetable_matter_pct=D("2.0"),
            clean_yield_pct=D("70.0"),
            staple_length_mm=D("95.0"),
        )
        blended_wool = fine_merino.blend_with(broad_merino, self_mass_fraction=D("0.5"))
        self.assertEqual(blended_wool.clean_yield_pct, D("67.50"))
        self.assertEqual(blended_wool.vegetable_matter_pct, D("1.60"))
        self.assertGreater(blended_wool.micron_diameter_um, D("19.5"))
        self.assertLess(blended_wool.micron_diameter_um, D("20.1"))

        # DOC PS 20 North American 2x4 S4S Dry: nominal 2"x4" (8 in^2) vs actual 1.5"x3.5" (5.25 in^2) -> 0.65625
        spf_2x4 = TimberQualitySpec.from_north_american_nominal(2, 4, species="SPF", moisture_content_pct=D("19.0"))
        self.assertEqual(spf_2x4.actual_to_nominal_area_ratio, D("0.65625"))
        # 1 piece of 2x4 x 12 ft = (2 * 4 * 12) / 12 = 8.0 nominal board feet
        self.assertEqual(spf_2x4.piece_nominal_board_feet(Q_(D("12"), "foot")), D("8"))
        # Round-trip MBF <-> actual m^3
        actual_m3 = spf_2x4.nominal_mbf_to_actual_volume(D("10"))
        self.assertAlmostEqual(actual_m3.magnitude, D("10") * D("2.359737216") * D("0.65625"), places=18)
        self.assertAlmostEqual(spf_2x4.actual_volume_to_nominal_mbf(actual_m3), D("10"), places=18)
        # Oven-dry 19% MC <-> wet-basis MC (19 / 119 * 100)
        self.assertAlmostEqual(spf_2x4.wet_basis_moisture_pct, D("1900") / D("119"), places=18)

    def test_material_spec_bridges_assay_adapters_and_validation(self) -> None:
        all_specs = (
            WheatQualitySpec(moisture_pct=D("13"), protein_pct=D("12")),
            CornQualitySpec(moisture_pct=D("15"), broken_kernels_pct=D("2"), foreign_material_pct=D("1")),
            SoybeanQualitySpec(moisture_pct=D("13"), protein_pct=D("35"), oil_content_pct=D("19")),
            RiceQualitySpec(milling_recovery_pct=D("70"), broken_rice_pct=D("5"), head_rice_basis="of_milled"),
            BarleyQualitySpec(moisture_pct=D("13.5"), protein_pct=D("11"), germination_energy_pct=D("96")),
            CoffeeQualitySpec(screen_size=D("16"), moisture_pct=D("11.5")),
            CocoaQualitySpec(bean_count_per_100g=D("95"), fermentation_pct=D("80"), fat_content_pct=D("52")),
            SugarQualitySpec(polarisation_deg_z=D("96"), icumsa_colour_iu=D("800")),
            CottonQualitySpec.from_staple_32nds(staple_32nds_inch=D("36"), micronaire=D("4.2"), strength_gpt=D("29")),
            RubberQualitySpec(dirt_pct=D("0.1"), ash_pct=D("0.6"), volatile_matter_pct=D("0.5"), mooney_viscosity_ml1_4_100c=D("70")),
            MilkQualitySpec(butterfat_pct=D("4.0"), protein_pct=D("3.3"), lactose_pct=D("4.8")),
            VegetableOilQualitySpec(free_fatty_acid_pct=D("0.5"), peroxide_value_meq_per_kg=D("3"), iodine_value_g_per_100g=D("130")),
            WoolQualitySpec(micron_diameter_um=D("19"), vegetable_matter_pct=D("1.2"), clean_yield_pct=D("66")),
            TimberQualitySpec.from_north_american_nominal(2, 6),
        )
        for qspec in all_specs:
            with self.subTest(spec=type(qspec).__name__):
                mspec = material_spec_from_assay(qspec)
                self.assertIsInstance(mspec, MaterialSpec)
                self_deltas = qspec.compare_to_benchmark(qspec)
                for delta_val in self_deltas.values():
                    self.assertAlmostEqual(delta_val, D("0"), places=18)

        # Existing agricultural *Assay -> *QualitySpec adapters
        for assay in (
            ag.GrainAssay.wheat_60lb(),
            ag.GrainAssay.corn_56lb(),
            ag.GrainAssay.soybeans_60lb(),
            ag.GrainAssay.barley_48lb(),
            ag.RiceAssay.thai_white_5pct_broken(),
            ag.SoftsAssay.raw_sugar_no11(),
            ag.SoftsAssay.green_coffee_60kg_bag(),
            ag.SoftsAssay.cocoa_beans_62_5kg_bag(),
            ag.SoftsAssay.us_cotton_480lb_bale(),
            ag.LivestockDairyAssay.whole_milk_standard(),
            ag.OliveOilAssay.extra_virgin_evoo(),
            ag.RubberAssay.tsr20_standard(),
            ag.WoolAssay.australian_merino_19micron(),
            ag.PulpLumberAssay.spf_dimensional_lumber(),
        ):
            with self.subTest(assay=assay.name):
                adapted = quality_spec_from_assay(assay)
                self.assertIsInstance(adapted.to_material_spec(), MaterialSpec)

        # Immutability & foreign-registry rejection
        w = WheatQualitySpec(moisture_pct=D("13"), protein_pct=D("12"))
        with self.assertRaises(FrozenInstanceError):
            w.protein_pct = D("14")  # type: ignore[misc]
        foreign_ureg = pint.UnitRegistry()
        with self.assertRaises(ForeignRegistryError):
            WheatQualitySpec(
                moisture_pct=D("13"),
                protein_pct=D("12"),
                test_weight=foreign_ureg.Quantity(78, "kg / hectoliter"),
            )
        with self.assertRaises(ForeignRegistryError):
            w.normalized_mass(foreign_ureg.Quantity(100, "metric_ton"))


if __name__ == "__main__":
    unittest.main()
