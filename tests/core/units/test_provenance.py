"""Tests for ConversionResult audit metadata, uncertainty propagation with covariance, and IEEE 754 precision."""

from __future__ import annotations

from dataclasses import fields
from decimal import Decimal
import unittest

import core.units as cu
from core.units import (
    BARREL_TO_M3,
    BiomassAssay,
    ConversionClassification,
    ConversionResult,
    CrudeAssay,
    EthanolFeedstockAssay,
    GasAssay,
    GasReferenceConditions,
    GrainAssay,
    HeatContentAssay,
    IEEE754_FLOAT64_DECIMAL_DIGITS,
    IEEE754_FLOAT64_MANTISSA_BITS,
    IEEE754_FLOAT64_MIN_SIGNIFICANT_DIGITS,
    LivestockDairyAssay,
    MissingReferenceConditionsError,
    OreAssay,
    PulpLumberAssay,
    Q_,
    ReferenceConditions,
    SoftsAssay,
    StandardsApplicabilityError,
    WATER_DENSITY_60F_KG_M3,
    convert,
    convert_actual_to_standard_gas_volume,
    convert_chemical_strength,
    convert_currency,
    convert_fx,
    convert_gas_energy,
    convert_gas_reference_volume,
    convert_heat_content,
    convert_hydrocarbon,
    convert_livestock_dairy,
    convert_pressure,
    convert_viscosity,
    format_ieee754,
    pad_ieee754_decimal,
    propagate_callable_uncertainty,
    propagate_uncertainty,
)


class ConversionAuditTests(unittest.TestCase):
    """Verify auditable ConversionResult reconstruction and dataclass contracts."""

    def test_conversion_result_exact_field_layout_and_brent_audit_trail(self) -> None:
        expected_fields = (
            "quantity",
            "original_quantity",
            "commodity_id",
            "conversion_path",
            "assay_id",
            "assay_version",
            "reference_conditions",
            "method",
            "standard",
            "standard_version",
            "registry_version",
            "registry_hash",
            "estimated",
            "uncertainty",
            "warnings",
        )
        self.assertEqual(tuple(f.name for f in fields(ConversionResult)), expected_fields)

        brent = CrudeAssay.brent()
        result = convert(
            Q_(1000, "bbl"),
            "metric_ton",
            assay=brent,
            return_metadata=True,
            strict=True,
        )
        self.assertIsInstance(result, ConversionResult)
        self.assertEqual(result.original_quantity.units, Q_(1, "bbl").units)
        self.assertEqual(result.original_quantity.magnitude, 1000)
        self.assertEqual(result.quantity.units, cu.parse_unit("metric_ton"))

        expected_mt = 1000.0 * BARREL_TO_M3 * (141.5 / (38.3 + 131.5)) * WATER_DENSITY_60F_KG_M3 / 1000.0
        self.assertAlmostEqual(float(result.quantity.magnitude), expected_mt, places=13)

        # Reconstruct why the conversion produced its result
        self.assertEqual(result.commodity_id, "crude_oil")
        self.assertEqual(result.assay_id, "CrudeAssay.brent")
        self.assertIsNotNone(result.assay_version)
        assert result.assay_version is not None
        self.assertTrue(result.assay_version.startswith("2026-11-27+sha256:"))
        self.assertEqual(result.classification, ConversionClassification.ASSAY_DEPENDENT.value)
        self.assertEqual(result.standard, "API MPMS Ch. 11.1 / ASTM D1250")
        self.assertEqual(result.standard_version, "60F_101.325kPa")
        self.assertEqual(result.registry_version, cu.REGISTRY_VERSION)
        self.assertEqual(result.registry_hash, cu.get_registry_definition_hash(cu.ureg))
        self.assertGreaterEqual(len(result.conversion_path), 2)
        self.assertEqual(result.conversion_path[0], "bbl")
        self.assertEqual(result.conversion_path[-1], "metric_ton")

        ref = result.reference_conditions
        self.assertIsInstance(ref, ReferenceConditions)
        assert ref is not None
        self.assertAlmostEqual(float(ref.temperature.to("degF").magnitude), 60.0, places=12)
        self.assertAlmostEqual(float(ref.pressure.to("Pa").magnitude), 101325.0, places=9)
        self.assertAlmostEqual(float(ref.density.to("kg/m^3").magnitude), brent.density_kg_per_m3, places=12)
        self.assertAlmostEqual(float(ref.parameters["api_gravity"]), 38.3, places=12)
        self.assertTrue(result.estimated)  # Reference benchmark value_kind == 'indicative_benchmark'
        self.assertIsNone(result.uncertainty)

        # Mutating returned quantity copies never corrupts the frozen ConversionResult snapshot
        q_view = result.quantity
        q_view.ito("kilogram")
        self.assertEqual(result.quantity.units, cu.parse_unit("metric_ton"))
        d_view = ref.density
        d_view.ito("g/cm^3")
        self.assertEqual(result.reference_conditions.density.units, cu.parse_unit("kg/m^3"))

        # Changing assay density deterministically changes assay_version hash
        custom_crude = CrudeAssay(name="Cargo #19", api_gravity=39.1)
        custom_res = convert(Q_(1000, "bbl"), "metric_ton", assay=custom_crude, return_metadata=True, strict=True)
        self.assertNotEqual(custom_res.assay_version, result.assay_version)
        self.assertFalse(custom_res.estimated)

        audit_dict = result.to_dict()
        self.assertEqual(audit_dict["assay_id"], "CrudeAssay.brent")
        self.assertEqual(audit_dict["registry_version"], cu.REGISTRY_VERSION)

    def test_all_nine_uncertainty_classifications_match_domain_taxonomy(self) -> None:
        # 1. kg -> tonne: Exact
        r_mass = convert(Q_(Decimal("2500"), "kg"), "metric_ton", return_metadata=True, strict=True)
        self.assertEqual(r_mass.classification, ConversionClassification.EXACT.value)
        self.assertFalse(r_mass.estimated)
        self.assertEqual(r_mass.uncertainty, Decimal("0"))
        self.assertEqual(format(r_mass.uncertainty, "f"), "0.0000000000000000")

        # 2. MWh -> GJ: Exact
        r_energy = convert(Q_(Decimal("10"), "MWh"), "GJ", return_metadata=True, strict=True)
        self.assertEqual(r_energy.classification, ConversionClassification.EXACT.value)
        self.assertFalse(r_energy.estimated)
        self.assertEqual(r_energy.uncertainty, Decimal("0"))

        # 3. Petroleum barrel -> litre: Exact
        r_vol = convert(Q_(Decimal("1"), "bbl"), "liter", return_metadata=True, strict=True)
        self.assertEqual(r_vol.classification, ConversionClassification.EXACT.value)
        self.assertFalse(r_vol.estimated)
        self.assertAlmostEqual(float(r_vol.quantity.magnitude), 158.987294928, places=12)

        # 4. USD -> GBP: Market-dependent
        r_fx = convert(
            Q_(Decimal("100"), "USD"),
            "GBP",
            rates={"GBP": Decimal("0.7850")},
            return_metadata=True,
            strict=True,
        )
        self.assertEqual(r_fx.classification, ConversionClassification.MARKET_DEPENDENT.value)
        self.assertFalse(r_fx.estimated)
        self.assertIsNone(r_fx.uncertainty)
        self.assertEqual(r_fx.reference_conditions.parameters["fx_rate"], Decimal("0.7850"))

        # 5. Barrel of crude -> tonne: Assay-dependent
        r_crude = convert(
            Q_(Decimal("1000"), "bbl"),
            "metric_ton",
            assay=CrudeAssay.brent(),
            return_metadata=True,
            strict=True,
        )
        self.assertEqual(r_crude.classification, ConversionClassification.ASSAY_DEPENDENT.value)

        # 6. LNG tonne -> MMBtu: Composition-dependent
        lng_assay = GasAssay(name="Sabine Pass LNG", gcv_btu_per_scf=1030.0, relative_density=0.58)
        r_lng = convert(
            Q_(Decimal("100"), "metric_ton"),
            "MMBtu",
            assay=lng_assay,
            basis="hhv",
            return_metadata=True,
            strict=True,
        )
        self.assertEqual(r_lng.classification, ConversionClassification.COMPOSITION_DEPENDENT.value)
        self.assertEqual(r_lng.commodity_id, "lng")

        # 7. Live cattle -> carcass mass: Estimated or measured yield
        cattle = LivestockDairyAssay.live_cattle_standard()
        r_cattle = convert_livestock_dairy(
            Q_(Decimal("50"), "metric_ton"),
            "metric_ton",
            assay=cattle,
            mode="live_to_carcass",
            return_metadata=True,
            strict=True,
        )
        self.assertEqual(r_cattle.classification, ConversionClassification.YIELD_DEPENDENT.value)
        self.assertEqual(r_cattle.quantity.magnitude, Decimal("31"))

        # 8. Gas actual volume -> standard volume: State and model-dependent
        r_gas = convert_actual_to_standard_gas_volume(
            Q_(Decimal("1000"), "m^3"),
            Q_(Decimal("5000"), "kPa"),
            Q_(Decimal("290"), "K"),
            target_standard="iso15",
            method="dak",
            return_metadata=True,
            strict=True,
        )
        self.assertEqual(r_gas.classification, ConversionClassification.STATE_AND_MODEL_DEPENDENT.value)
        self.assertTrue(r_gas.estimated)

        # 9. Biomass tonnes -> MWh: Moisture and calorific-dependent
        pellets = BiomassAssay.enplus_a1_pellets()
        r_bio = convert(
            Q_(Decimal("100"), "metric_ton"),
            "MWh",
            assay=pellets,
            return_metadata=True,
            strict=True,
        )
        self.assertEqual(r_bio.classification, ConversionClassification.MOISTURE_AND_CALORIFIC_DEPENDENT.value)
        self.assertGreater(float(r_bio.quantity.magnitude), 400.0)

    def test_strict_mode_rejects_approximate_fallbacks_and_accepts_explicit_specs(self) -> None:
        # Legacy temperature-only VCF rejected in strict mode
        brent = CrudeAssay.brent()
        with self.assertRaises(StandardsApplicabilityError):
            convert(Q_(1000, "bbl"), "metric_ton", assay=brent, observed_temp=85.0, strict=True)

        # Unqualified neutral pressure unit in datum shift rejected in strict mode
        with self.assertRaises(MissingReferenceConditionsError):
            convert_pressure(Q_(100.0, "psig"), "psi", strict=True)
        # Explicit psia succeeds in strict mode
        ok_p = convert_pressure(Q_(0.0, "psig"), "psia", strict=True, return_metadata=True)
        self.assertAlmostEqual(float(ok_p.quantity.magnitude), cu.STANDARD_ATMOSPHERE_PSI, places=12)

        # Auto Z-factor EOS selection rejected in strict mode without explicit Z or named method
        with self.assertRaises(MissingReferenceConditionsError):
            convert_actual_to_standard_gas_volume(
                Q_(100.0, "m^3"), Q_(2000.0, "kPa"), Q_(300.0, "K"), strict=True
            )
        # Explicit Z factors succeed in strict mode with estimated=False
        exact_gas = convert_actual_to_standard_gas_volume(
            Q_(100.0, "m^3"),
            Q_(2000.0, "kPa"),
            Q_(300.0, "K"),
            z_actual=0.965,
            z_standard=0.998,
            strict=True,
            return_metadata=True,
        )
        self.assertFalse(exact_gas.estimated)

        # Latent-water NCV approximation rejected in strict mode
        est_hc = HeatContentAssay(
            gcv=Q_(Decimal("25"), "MJ/kg"),
            hydrogen_pct=5.0,
            moisture_pct=10.0,
            method="latent_water",
        )
        with self.assertRaises(StandardsApplicabilityError):
            convert_heat_content(Q_(10.0, "kg"), "MJ", assay=est_hc, basis="lhv", strict=True)


class UncertaintyPropagationTests(unittest.TestCase):
    """Verify independent and covariance-aware first-order Taylor error propagation."""

    def test_propagate_uncertainty_independent_and_correlated_exact_formulas(self) -> None:
        # Product f(V, rho) = V * rho at V = 100, rho = 0.85 -> y = 85
        # J_V = 0.85, J_rho = 100; sigma_V = 2.0, sigma_rho = 0.01
        # Independent: var = (0.85*2)^2 + (100*0.01)^2 = 1.7^2 + 1.0^2 = 2.89 + 1.00 = 3.89
        indep = propagate_uncertainty(
            sensitivities={"V": Decimal("0.85"), "rho": Decimal("100")},
            standard_uncertainties={"V": Decimal("2.0"), "rho": Decimal("0.01")},
            value=Decimal("85"),
        )
        self.assertAlmostEqual(float(indep.variance), 3.89, places=14)
        self.assertAlmostEqual(float(indep.standard_uncertainty), 3.89 ** 0.5, places=14)
        self.assertEqual(indep.covariance_contribution, Decimal("0"))

        # Positive correlation r = +0.5 -> cov_term = 2 * 1.7 * 1.0 * 0.5 = +1.70 -> var = 5.59
        pos_corr = propagate_uncertainty(
            sensitivities={"V": Decimal("0.85"), "rho": Decimal("100")},
            standard_uncertainties={"V": Decimal("2.0"), "rho": Decimal("0.01")},
            value=Decimal("85"),
            correlations={("V", "rho"): Decimal("0.5")},
        )
        self.assertAlmostEqual(float(pos_corr.variance), 5.59, places=14)
        self.assertAlmostEqual(float(pos_corr.covariance_contribution), 1.70, places=14)
        self.assertGreater(pos_corr.standard_uncertainty, indep.standard_uncertainty)

        # Equivalent explicit covariance Cov(V, rho) = 0.5 * 2.0 * 0.01 = 0.01
        pos_cov = propagate_uncertainty(
            sensitivities={"V": Decimal("0.85"), "rho": Decimal("100")},
            standard_uncertainties={"V": Decimal("2.0"), "rho": Decimal("0.01")},
            value=Decimal("85"),
            covariance={("V", "rho"): Decimal("0.01")},
        )
        self.assertEqual(pos_cov.standard_uncertainty, pos_corr.standard_uncertainty)

        # Perfect negative correlation on equal-impact variables cancels variance to zero
        zero_cancel = propagate_uncertainty(
            sensitivities={"x1": Decimal("1"), "x2": Decimal("1")},
            standard_uncertainties={"x1": Decimal("3"), "x2": Decimal("3")},
            value=Decimal("10"),
            correlations={("x1", "x2"): Decimal("-1.0")},
        )
        self.assertEqual(zero_cancel.standard_uncertainty, Decimal("0"))
        self.assertEqual(format(zero_cancel.standard_uncertainty, "f"), "0.0000000000000000")

        # Invalid correlation / non-PSD / Cauchy-Schwarz violations rejected
        with self.assertRaises(ValueError):
            propagate_uncertainty({"x": 1, "y": 1}, {"x": 1, "y": 1}, correlations={("x", "y"): Decimal("1.2")})
        with self.assertRaises(ValueError):
            propagate_uncertainty({"x": 1, "y": 1}, {"x": 1, "y": 2}, covariance={("x", "y"): Decimal("2.5")})
        # 3-variable correlation matrix with all pairwise r = -0.8 is not PSD: sum_var = 3 - 6*0.8 = -1.8 < 0
        with self.assertRaises(ValueError):
            propagate_uncertainty(
                {"a": 1, "b": 1, "c": 1},
                {"a": 1, "b": 1, "c": 1},
                correlations={("a", "b"): -0.8, ("b", "c"): -0.8, ("a", "c"): -0.8},
            )

    def test_end_to_end_correlated_uncertainty_in_convert(self) -> None:
        brent = CrudeAssay.brent()
        # m = V_m3 * rho / 1000
        res_uncorr = convert(
            Q_(Decimal("1000"), "bbl"),
            "metric_ton",
            assay=brent,
            return_metadata=True,
            uncertainties={"quantity": Decimal("5"), "density": Decimal("2.5")},
        )
        res_neg_corr = convert(
            Q_(Decimal("1000"), "bbl"),
            "metric_ton",
            assay=brent,
            return_metadata=True,
            uncertainties={"quantity": Decimal("5"), "density": Decimal("2.5")},
            correlations={("quantity", "density"): Decimal("-0.6")},
        )
        res_pos_corr = convert(
            Q_(Decimal("1000"), "bbl"),
            "metric_ton",
            assay=brent,
            return_metadata=True,
            uncertainties={"quantity": Decimal("5"), "density": Decimal("2.5")},
            correlations={("quantity", "density"): Decimal("0.6")},
        )
        assert res_uncorr.uncertainty is not None
        assert res_neg_corr.uncertainty is not None
        assert res_pos_corr.uncertainty is not None

        # Analytical check: m = V * mt_per_bbl = V * BARREL_TO_M3 * rho / 1000
        j_v = BARREL_TO_M3 * brent.density_kg_per_m3 / 1000.0
        j_rho = 1000.0 * BARREL_TO_M3 / 1000.0
        expected_uncorr = ((j_v * 5.0) ** 2 + (j_rho * 2.5) ** 2) ** 0.5
        expected_neg = ((j_v * 5.0) ** 2 + (j_rho * 2.5) ** 2 + 2.0 * (-0.6) * (j_v * 5.0) * (j_rho * 2.5)) ** 0.5
        expected_pos = ((j_v * 5.0) ** 2 + (j_rho * 2.5) ** 2 + 2.0 * (0.6) * (j_v * 5.0) * (j_rho * 2.5)) ** 0.5

        self.assertAlmostEqual(float(res_uncorr.uncertainty), expected_uncorr, places=9)
        self.assertAlmostEqual(float(res_neg_corr.uncertainty), expected_neg, places=9)
        self.assertAlmostEqual(float(res_pos_corr.uncertainty), expected_pos, places=9)
        self.assertLess(res_neg_corr.uncertainty, res_uncorr.uncertainty)
        self.assertGreater(res_pos_corr.uncertainty, res_uncorr.uncertainty)

        # FX uncertainty propagation (both amount and FX rate)
        fx_res = convert_currency(
            Q_(Decimal("1000"), "USD"),
            "GBP",
            rates={"GBP": Decimal("0.80")},
            return_metadata=True,
            uncertainties={"quantity": Decimal("10"), "fx_rate": Decimal("0.005")},
        )
        # y = Q * r = 800; dy/dQ = 0.80, dy/dr = 1000 -> sigma = sqrt((8)^2 + (5)^2) = sqrt(89)
        self.assertAlmostEqual(float(fx_res.uncertainty), 89.0 ** 0.5, places=11)

        # Callable uncertainty propagation on nonlinear API gravity conversion
        api_prop = propagate_callable_uncertainty(
            lambda p: convert(Q_(p["bbl"], "bbl"), "metric_ton", api_gravity=float(p["api"])),
            inputs={"bbl": Decimal("1000"), "api": Decimal("38.3")},
            standard_uncertainties={"bbl": Decimal("2"), "api": Decimal("0.2")},
        )
        self.assertGreater(api_prop.standard_uncertainty, Decimal("0"))
        self.assertLess(api_prop.sensitivities["api"], Decimal("0"))  # Higher API -> lighter crude -> fewer tonnes

        # Secondary converter metadata support (viscosity & chemical strength)
        visc_res = convert_viscosity(
            Q_(Decimal("12.5"), "centipoise"),
            "centistokes",
            density=Q_(Decimal("850"), "kg/m^3"),
            return_metadata=True,
        )
        self.assertEqual(visc_res.classification, ConversionClassification.ASSAY_DEPENDENT.value)

        chem_res = convert_chemical_strength(
            Q_(Decimal("2.0"), "molarity"),
            "normality",
            equivalence_factor=2,
            return_metadata=True,
        )
        self.assertEqual(chem_res.classification, ConversionClassification.COMPOSITION_DEPENDENT.value)
        self.assertEqual(chem_res.quantity.magnitude, Decimal("4"))


class IEEE754PrecisionAndTrailingZerosTests(unittest.TestCase):
    """Verify IEEE 754 binary64 15-17 sig-digit accuracy and full trailing-zero padding."""

    def test_ieee754_constants_and_trailing_zero_padding(self) -> None:
        self.assertEqual(IEEE754_FLOAT64_MANTISSA_BITS, 53)
        self.assertEqual(IEEE754_FLOAT64_MIN_SIGNIFICANT_DIGITS, 15)
        self.assertEqual(IEEE754_FLOAT64_DECIMAL_DIGITS, 17)

        # 100.00 -> 17 significant digits with full trailing zeros: 100.00000000000000
        d100 = pad_ieee754_decimal(Decimal("100.00"))
        self.assertEqual(d100, Decimal("100"))
        self.assertEqual(str(d100), "100.00000000000000")
        self.assertEqual(format_ieee754(100.0), "100.00000000000000")
        self.assertEqual(format_ieee754(100), "100.00000000000000")
        self.assertEqual(format_ieee754(0), "0.0000000000000000")
        self.assertEqual(format(pad_ieee754_decimal(0), "f"), "0.0000000000000000")

        # High-precision Decimal (>17 sig digits) is never truncated
        high_prec = Decimal("10.123456789123456789")
        self.assertEqual(str(pad_ieee754_decimal(high_prec)), "10.123456789123456789")

        # Decimal conversions retain full trailing zeros to >=17 significant digits
        q_out = convert(Q_(Decimal("100.00"), "metric_ton"), "metric_ton")
        self.assertEqual(str(q_out.magnitude), "100.00000000000000")

        q_scaled = convert(Q_(Decimal("100000"), "kilogram"), "metric_ton")
        self.assertEqual(str(q_scaled.magnitude), "100.00000000000000")

        q_fx = convert_fx(Q_(Decimal("100"), "USD"), "EUR", {"USD/EUR": Decimal("1.00")})
        self.assertEqual(str(q_fx.magnitude), "100.00000000000000")

    def test_agricultural_and_physical_helpers_retain_full_ieee754_precision(self) -> None:
        # Previously truncated to 4 decimals via round(..., 4) / .quantize(Decimal('0.0001'))
        sugar_px = SoftsAssay.cents_per_lb_to_usd_per_mt(18.5)
        expected_px = (18.5 / 100.0) * (1000.0 / 0.45359237)
        self.assertAlmostEqual(sugar_px, expected_px, places=12)

        lumber = PulpLumberAssay.nbsk_softwood_pulp()
        mbf_dec = lumber.m3_price_to_mbf_price(Decimal("100"))
        self.assertEqual(str(mbf_dec), "235.97372160000000")

        corn = GrainAssay.corn_56lb(actual_moisture_pct=20.0)
        adj_px = corn.adjust_price_for_moisture(Decimal("200"))
        self.assertAlmostEqual(adj_px, Decimal("200") * Decimal("80") / Decimal("84.5"), places=16)

        etoh = EthanolFeedstockAssay.us_corn_ethanol()
        be_gal = etoh.breakeven_ethanol_usd_per_gallon(4.5, 150.0)
        self.assertAlmostEqual(be_gal, 3.225 / 2.85, places=15)

        # Ore moisture conversion retains full trailing zeros on exact Decimal result
        ore = OreAssay.iron_ore_62fe()
        dry_q = ore.wmt_to_dmt(Q_(Decimal("100"), "wmt"))
        self.assertEqual(dry_q.magnitude, Decimal("92"))
        self.assertEqual(str(dry_q.magnitude), "92.000000000000000")


if __name__ == "__main__":
    unittest.main()
