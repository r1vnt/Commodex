"""Unit tests for composable MaterialComposition, MaterialProperties, and ProcessingSpecification."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from decimal import Decimal as D
import unittest

import pint

import core.units as u
from core.units import agricultural_specs as ag
from core.units import conversions as conv
from core.units.materials import (
    CompositionBasis,
    MaterialComposition,
    MaterialProperties,
    MaterialSpec,
    ProcessingSpecification,
    create_material_context,
    material_spec_from_assay,
)
from core.units.registry import ForeignRegistryError, Q_


class MaterialCompositionTests(unittest.TestCase):
    def test_analytical_basis_round_trips_ar_ad_db_daf(self) -> None:
        # Thermal coal / concentrate style proximate + elemental composition on AS_RECEIVED basis:
        # moisture = 10%, ash = 9% AR (=> 10% DB), fixed_carbon = 54% AR (=> 60% DB => 66.666...% DAF), S = 0.9% AR (=> 1% DB)
        comp_ar = MaterialComposition(
            components={"fixed_carbon": D("0.54"), "volatile_matter": D("0.261")},
            basis=CompositionBasis.AS_RECEIVED,
            moisture_fraction=D("0.10"),
            ash_fraction=D("0.09"),
            impurities={"S": D("0.009")},
        )

        comp_db = comp_ar.to_basis("dry")
        self.assertEqual(comp_db.basis, CompositionBasis.DRY)
        self.assertAlmostEqual(comp_db.components["fixed_carbon"], D("0.60"), places=18)
        self.assertAlmostEqual(comp_db.ash_fraction, D("0.10"), places=18)
        self.assertAlmostEqual(comp_db.impurities["S"], D("0.01"), places=18)

        comp_daf = comp_db.to_basis("daf")
        self.assertEqual(comp_daf.basis, CompositionBasis.DRY_ASH_FREE)
        self.assertAlmostEqual(comp_daf.components["fixed_carbon"], D("2") / D("3"), places=18)

        comp_ad = comp_daf.to_basis(CompositionBasis.AIR_DRIED, target_moisture_fraction=D("0.05"))
        self.assertEqual(comp_ad.basis, CompositionBasis.AIR_DRIED)
        self.assertAlmostEqual(comp_ad.components["fixed_carbon"], D("0.57"), places=18)
        self.assertAlmostEqual(comp_ad.ash_fraction, D("0.095"), places=18)

        restored_ar = comp_ad.to_basis(CompositionBasis.AS_RECEIVED, target_moisture_fraction=D("0.10"))
        self.assertAlmostEqual(restored_ar.components["fixed_carbon"], D("0.54"), places=18)
        self.assertAlmostEqual(restored_ar.ash_fraction, D("0.09"), places=18)
        self.assertAlmostEqual(restored_ar.impurities["S"], D("0.009"), places=18)

    def test_wet_dry_and_contained_material_calculations(self) -> None:
        cu_conc = MaterialComposition(
            components={"Cu": D("0.28"), "Au": D("0.000010")},
            basis=CompositionBasis.DRY,
            moisture_fraction=D("0.08"),
            impurities={"As": D("0.003")},
        )

        # 1,000 wmt @ 8% moisture -> 920 dmt
        dmt = cu_conc.wet_to_dry_mass(Q_(D("1000"), "wmt"))
        self.assertEqual(dmt.units, u.parse_unit("dmt"))
        self.assertEqual(dmt.magnitude, D("920.00"))
        self.assertEqual(cu_conc.dry_to_wet_mass(dmt).magnitude, D("1000"))

        # Contained Cu from wet mass (1000 wmt -> 920 dmt * 0.28 = 257.6 mt Cu)
        cu_from_wet = cu_conc.contained_mass(Q_(D("1000"), "wmt"), "Cu")
        self.assertEqual(cu_from_wet.magnitude, D("257.6000"))
        cu_from_dry = cu_conc.contained_mass(dmt, "Cu")
        self.assertEqual(cu_from_dry.magnitude, D("257.6000"))

        # Inverse: required wet and dry mass for 257.6 mt contained Cu
        req_dry = cu_conc.required_mass_for_contained(cu_from_wet, "Cu", output_basis="dry")
        self.assertEqual(req_dry.units, u.parse_unit("dmt"))
        self.assertEqual(req_dry.magnitude, D("920"))

        req_wet = cu_conc.required_mass_for_contained(cu_from_wet, "Cu", output_basis="ar")
        self.assertEqual(req_wet.units, u.parse_unit("wmt"))
        self.assertEqual(req_wet.magnitude, D("1000"))

        # Master router integration with MaterialComposition
        self.assertEqual(u.convert(Q_(D("1000"), "wmt"), "dmt", assay=cu_conc).magnitude, D("920.00"))
        self.assertEqual(
            u.convert(Q_(D("920"), "dmt"), "contained_metal_tonne", assay=cu_conc, component="Cu").magnitude,
            D("257.60"),
        )

    def test_basis_aware_blending_conserves_wet_dry_water_and_metal_mass(self) -> None:
        lot_a = MaterialComposition(
            components={"Cu": D("0.30")},
            basis=CompositionBasis.DRY,
            moisture_fraction=D("0.10"),
            impurities={"As": D("0.004")},
        )
        lot_b = MaterialComposition(
            components={"Cu": D("0.20")},
            basis=CompositionBasis.DRY,
            moisture_fraction=D("0.05"),
            impurities={"As": D("0.001")},
        )

        # Blend 600 wmt of Lot A with 400 wmt of Lot B
        w_a, w_b = D("600"), D("400")
        blended = lot_a.blend_with(lot_b, self_weight=w_a, other_weight=w_b, weight_basis="as_received")

        total_wmt = Q_(w_a + w_b, "wmt")
        expected_dmt = lot_a.wet_to_dry_mass(Q_(w_a, "wmt")).magnitude + lot_b.wet_to_dry_mass(Q_(w_b, "wmt")).magnitude
        self.assertAlmostEqual(blended.wet_to_dry_mass(total_wmt).magnitude, expected_dmt, places=18)

        expected_cu = (
            lot_a.contained_mass(Q_(w_a, "wmt"), "Cu").magnitude
            + lot_b.contained_mass(Q_(w_b, "wmt"), "Cu").magnitude
        )
        self.assertAlmostEqual(blended.contained_mass(total_wmt, "Cu").magnitude, expected_cu, places=18)


class MaterialPropertiesTests(unittest.TestCase):
    def test_density_porosity_purity_moisture_and_molar_conversions(self) -> None:
        props = MaterialProperties(
            density=Q_(D("2500"), "kg / meter ** 3"),
            porosity=D("0.40"),
            molar_mass=Q_(D("100"), "g / mole"),
            viscosity=Q_(D("10"), "centipoise"),
            moisture=D("0.10"),
            purity=D("0.80"),
        )
        # Effective bulk density = 2500 * (1 - 0.40) = 1500 kg/m^3
        self.assertEqual(props.effective_bulk_density.to("kg/m^3").magnitude, D("1500.00"))
        self.assertEqual(props.implied_porosity, D("0.40"))

        # True vs bulk volume <-> mass
        self.assertEqual(props.volume_to_mass(Q_(D("2"), "m^3"), target_unit="mt").magnitude, D("5.000"))
        self.assertEqual(props.volume_to_mass(Q_(D("2"), "m^3"), target_unit="mt", use_bulk=True).magnitude, D("3.000"))
        self.assertEqual(props.mass_to_volume(Q_(D("3"), "mt"), target_unit="m^3", use_bulk=True).magnitude, D("2"))

        # Purity gross <-> pure
        self.assertEqual(props.gross_to_pure(Q_(D("100"), "kg")).magnitude, D("80.00"))
        self.assertEqual(props.pure_to_gross(Q_(D("80"), "kg")).magnitude, D("100"))

        # Moisture wet <-> dry
        self.assertEqual(props.wet_to_dry(Q_(D("100"), "wmt")).magnitude, D("90.0"))
        self.assertEqual(props.dry_to_wet(Q_(D("90"), "dmt")).magnitude, D("100"))

        # Molar mass gross vs active_only
        self.assertEqual(props.mass_to_moles(Q_(D("10"), "kg")).magnitude, D("100"))
        self.assertEqual(props.mass_to_moles(Q_(D("10"), "kg"), active_only=True).magnitude, D("80.0"))
        self.assertEqual(props.moles_to_mass(Q_(D("80"), "mole"), active_only=True).magnitude, D("10"))

        # Dynamic <-> kinematic viscosity via true density (10 cP / 2.5 g/cm^3 = 4 cSt)
        self.assertAlmostEqual(props.kinematic_viscosity("cSt").magnitude, D("4"), places=12)
        self.assertAlmostEqual(props.dynamic_viscosity("cP").magnitude, D("10"), places=12)

        # Pint Context & Master Router
        ctx = create_material_context(props, active_only=True)
        with u.ureg.context(ctx):
            self.assertEqual(Q_(D("10"), "kg").to("mole").magnitude, D("80.0"))
        self.assertEqual(u.convert(Q_(D("10"), "kg"), "mole", assay=props, active_only=True).magnitude, D("80.0"))
        self.assertEqual(u.convert(Q_(D("100"), "wmt"), "dmt", assay=props).magnitude, D("90.0"))

    def test_price_normalisation_and_harmonic_refutas_blending(self) -> None:
        caustic = MaterialProperties(
            density=Q_(D("1525"), "kg / meter ** 3"),
            molar_mass=Q_(D("40"), "g / mole"),
            purity=D("0.50"),
            moisture=D("0.50"),
        )
        # $300/mt 50% solution -> $600/mt 100% pure NaOH
        pure_px = caustic.normalize_price(
            Q_(D("300"), "USD / metric_tonne"),
            from_purity_basis="gross",
            to_purity_basis="pure",
        )
        self.assertEqual(pure_px.magnitude, D("600"))

        # Tagged wet -> dry quote + FX conversion ($90/wmt @ 10% moisture -> $100/dmt -> 90 EUR/dmt at 0.9 EUR/USD)
        ore_props = MaterialProperties(moisture=D("0.10"))
        eur_dmt_px = ore_props.normalize_price(
            Q_(D("90"), "USD / wmt"),
            "EUR / dmt",
            rates={"EUR": D("0.9")},
        )
        self.assertEqual(eur_dmt_px.units, u.parse_unit("EUR / dmt"))
        self.assertEqual(eur_dmt_px.magnitude, D("90.0"))

        # Harmonic density & molar-mass + Refutas VBN viscosity blending
        p1 = MaterialProperties(
            density=Q_(D("800"), "kg / meter ** 3"),
            molar_mass=Q_(D("100"), "g / mole"),
            viscosity=Q_(D("5"), "cSt"),
        )
        p2 = MaterialProperties(
            density=Q_(D("1000"), "kg / meter ** 3"),
            molar_mass=Q_(D("200"), "g / mole"),
            viscosity=Q_(D("500"), "cSt"),
        )
        blend = p1.blend_with(p2, self_mass_fraction=D("0.5"))
        # 1 / (0.5/800 + 0.5/1000) = 888.888... kg/m^3 (conserves total ideal volume)
        v_sum = p1.mass_to_volume(Q_(D("500"), "kg")).magnitude + p2.mass_to_volume(Q_(D("500"), "kg")).magnitude
        self.assertAlmostEqual(blend.mass_to_volume(Q_(D("1000"), "kg")).magnitude, v_sum, places=18)
        # 1 / (0.5/100 + 0.5/200) = 133.333... g/mol (conserves total moles)
        mol_sum = p1.mass_to_moles(Q_(D("500"), "kg")).magnitude + p2.mass_to_moles(Q_(D("500"), "kg")).magnitude
        self.assertAlmostEqual(blend.mass_to_moles(Q_(D("1000"), "kg")).magnitude, mol_sum, places=18)
        # Refutas viscosity is strictly below arithmetic mean (252.5 cSt) and above 5 cSt
        visc_cst = blend.kinematic_viscosity("cSt").magnitude
        self.assertGreater(visc_cst, D("15"))
        self.assertLess(visc_cst, D("60"))


class ProcessingAndMaterialSpecTests(unittest.TestCase):
    def test_processing_yields_gpw_margins_and_blending(self) -> None:
        crush = ProcessingSpecification(
            input_material="soybeans",
            output_yields={"meal": D("0.79"), "oil": D("0.19")},
            recovery_rates={"meal": D("1.0"), "oil": D("0.98")},
            process_losses={"hull_and_moisture_loss": D("0.02")},
        )
        self.assertEqual(crush.total_output_yield, D("0.98"))
        self.assertEqual(crush.total_process_loss, D("0.02"))
        self.assertEqual(crush.unaccounted_mass_fraction, D("0"))

        # Forward & inverse mass balance
        oil_out = crush.input_to_output(Q_(D("1000"), "mt"), "oil")
        self.assertEqual(oil_out.magnitude, D("186.2000"))
        req_beans = crush.output_to_input(oil_out, "oil")
        self.assertEqual(req_beans.magnitude, D("1000"))

        outputs = crush.all_outputs(Q_(D("100"), "mt"))
        self.assertEqual(outputs["meal"].magnitude, D("79.00"))
        self.assertEqual(outputs["oil"].magnitude, D("18.6200"))

        # Gross Processing Worth & Crush Margin
        prices = {
            "meal": Q_(D("400"), "USD / metric_tonne"),
            "oil": Q_(D("1000"), "USD / metric_tonne"),
        }
        # GPW = 0.79 * 400 + 0.1862 * 1000 = 316 + 186.2 = 502.2 USD/mt
        gpw = crush.gross_processing_worth(prices)
        self.assertEqual(gpw.magnitude, D("502.2000"))
        margin = crush.processing_margin(
            Q_(D("440"), "USD / metric_tonne"),
            prices,
            processing_cost=Q_(D("20"), "USD / metric_tonne"),
        )
        self.assertEqual(margin.magnitude, D("42.2000"))

        # Processing blend conserves recovered outputs
        crush_b = ProcessingSpecification(
            input_material="soybeans",
            output_yields={"meal": D("0.78"), "oil": D("0.20")},
            recovery_rates={"meal": D("1.0"), "oil": D("0.95")},
            process_losses={"hull_and_moisture_loss": D("0.02")},
        )
        blended_proc = crush.blend_with(crush_b, self_input_fraction=D("0.6"))
        expected_oil = (
            crush.input_to_output(Q_(D("600"), "mt"), "oil").magnitude
            + crush_b.input_to_output(Q_(D("400"), "mt"), "oil").magnitude
        )
        self.assertAlmostEqual(
            blended_proc.input_to_output(Q_(D("1000"), "mt"), "oil").magnitude,
            expected_oil,
            places=18,
        )

    def test_unified_material_spec_routing_and_specialised_assay_adapters(self) -> None:
        spec = MaterialSpec(
            name="Copper Concentrate",
            properties=MaterialProperties(
                bulk_density=Q_(D("1800"), "kg / meter ** 3"),
                moisture=D("0.08"),
            ),
            composition=MaterialComposition(
                components={"Cu": D("0.30")},
                basis=CompositionBasis.DRY,
            ),
            processing=ProcessingSpecification(
                input_material="Copper Concentrate",
                output_yields={"Cu": D("0.30")},
                recovery_rates={"Cu": D("0.965")},
            ),
        )
        # Moisture auto-synchronised between MaterialProperties and MaterialComposition
        self.assertEqual(spec.composition.moisture_fraction, D("0.08"))
        # Router: volume <-> mass, wmt <-> dmt <-> contained <-> payable, quotes, and thermal override
        self.assertEqual(u.convert(Q_(D("10"), "m^3"), "mt", assay=spec).magnitude, D("18"))
        self.assertEqual(u.convert(Q_(D("1000"), "wmt"), "dmt", assay=spec).magnitude, D("920.00"))
        self.assertEqual(u.convert(Q_(D("1000"), "wmt"), "contained_metal_tonne", assay=spec).magnitude, D("276.0000"))
        self.assertEqual(u.convert(Q_(D("1000"), "wmt"), "payable_metal_tonne", assay=spec).magnitude, D("266.3400000"))
        self.assertEqual(
            u.convert(Q_(D("180"), "USD / m^3"), "USD / mt", assay=spec).magnitude,
            D("100"),
        )
        self.assertEqual(
            u.convert(Q_(D("1"), "mt"), "GJ", assay=spec, ncv=Q_(D("20"), "MJ/kg")).magnitude,
            D("20"),
        )

        # Non-destructive adapters across ALL existing specialised assay presets + composite assays
        from tests.core.units.test_agricultural_specs import agricultural_presets
        from tests.core.units.test_registry_compatibility import assay_presets

        all_assays = [
            *(pair for pair in assay_presets()),
            *(pair for pair in agricultural_presets()),
            ("BulkOreAssay", conv.BulkOreAssay("fe_fines", D("0.08"), {"Fe": D("0.62")})),
            (
                "RefineryCrudeAssay",
                conv.RefineryCrudeAssay(
                    conv.CrudeAssay.brent(),
                    (conv.DistillationCut("gasoline", 30, 180, 35, 38, 0.74), conv.DistillationCut("resid", 180, 650, 65, 62, 0.92)),
                ),
            ),
            (
                "MultiMetalConcentrateAssay",
                u.MultiMetalConcentrateAssay(
                    name="Cu-Zn",
                    moisture_fraction=D("0.09"),
                    composition=u.ElementComposition({"Cu": D("0.25"), "Zn": D("0.10")}),
                    recovery_fractions={"Cu": D("0.95")},
                ),
            ),
        ]
        self.assertGreater(len(all_assays), 100)
        for label, assay in all_assays:
            with self.subTest(assay=label):
                adapted = material_spec_from_assay(assay)
                self.assertIsInstance(adapted, MaterialSpec)
                self.assertIs(MaterialSpec.from_assay(adapted), adapted)

        # Verify numerical parity between ChemicalAssay.caustic_soda_50pct() and its MaterialSpec
        caustic_assay = conv.ChemicalAssay.caustic_soda_50pct()
        caustic_spec = MaterialSpec.from_assay(caustic_assay)
        self.assertAlmostEqual(
            float(u.convert(Q_(D("1000"), "kg"), "mole", assay=caustic_spec, active_only=True).magnitude),
            float(u.convert(Q_(D("1000"), "kg"), "mole", assay=caustic_assay, active_only=True).magnitude),
            places=9,
        )

    def test_strict_input_validation_immutability_and_foreign_registry_rejection(self) -> None:
        # Immutability
        comp = MaterialComposition(components={"Fe": D("0.62")})
        with self.assertRaises(FrozenInstanceError):
            comp.basis = CompositionBasis.AS_RECEIVED  # type: ignore[misc]
        with self.assertRaises(TypeError):
            comp.components["Fe"] = D("0.50")  # type: ignore[index]

        # Closure > 1, overlapping component/impurity, boolean/NaN scalars, moisture == 1
        with self.assertRaises(ValueError):
            MaterialComposition(components={"Fe": D("0.70"), "SiO2": D("0.35")}, basis="dry")
        with self.assertRaises(ValueError):
            MaterialComposition(components={"S": D("0.02")}, impurities={"S": D("0.01")})
        with self.assertRaises(TypeError):
            MaterialComposition(components={"Fe": True})  # type: ignore[dict-item]
        with self.assertRaises(ValueError):
            MaterialComposition(components={"Fe": float("nan")})
        with self.assertRaises(ValueError):
            MaterialComposition(components={"Fe": D("0.5")}, moisture_fraction=D("1.0"))

        # Bulk density > true density & inconsistent porosity
        with self.assertRaises(ValueError):
            MaterialProperties(
                density=Q_(D("1500"), "kg/m^3"),
                bulk_density=Q_(D("1800"), "kg/m^3"),
            )
        with self.assertRaises(ValueError):
            MaterialProperties(
                density=Q_(D("2000"), "kg/m^3"),
                bulk_density=Q_(D("1500"), "kg/m^3"),
                porosity=D("0.50"),
            )

        # Processing mass balance violation (yields + losses > 1)
        with self.assertRaises(ValueError):
            ProcessingSpecification(
                input_material="crude",
                output_yields={"gasoline": D("0.60"), "distillate": D("0.35")},
                process_losses={"flare": D("0.10")},
            )

        # Conflicting moisture between MaterialProperties and MaterialComposition
        with self.assertRaises(ValueError):
            MaterialSpec(
                name="bad",
                properties=MaterialProperties(moisture=D("0.10")),
                composition=MaterialComposition(components={"Cu": D("0.25")}, moisture_fraction=D("0.08")),
            )

        # Foreign-registry pint.Quantity rejection
        foreign_ureg = pint.UnitRegistry()
        with self.assertRaises(ForeignRegistryError):
            MaterialProperties(density=foreign_ureg.Quantity(1000, "kg / meter ** 3"))
        valid_props = MaterialProperties(density=Q_(D("1000"), "kg / meter ** 3"))
        with self.assertRaises(ForeignRegistryError):
            valid_props.volume_to_mass(foreign_ureg.Quantity(1, "meter ** 3"))


if __name__ == "__main__":
    unittest.main()
