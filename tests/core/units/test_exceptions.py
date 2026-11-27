from __future__ import annotations

from decimal import Decimal
import subprocess
import sys
import unittest

import core.units as cu
from core.units import (
    AmbiguousUnitError,
    BatteryMaterialAssay,
    BulkOreAssay,
    CarbonDioxideAssay,
    CommodityUnitError,
    CompositionBasis,
    ConversionConvergenceError,
    ForeignRegistryError,
    GrainAssay,
    HeatContentAssay,
    IncompatibleAnalyticalBasisError,
    IncompatibleCommodityError,
    InvalidPhysicalStateError,
    MaterialComposition,
    MaterialProperties,
    MaterialSpecification,
    MissingAssayError,
    MissingFXRateError,
    MissingReferenceConditionsError,
    OilseedCrushAssay,
    OreAssay,
    PoultryEggAssay,
    Q_,
    RegistryCompatibilityError,
    RegistryFrozenError,
    RegistryVersionMismatchError,
    SolverConvergenceError,
    SolverDomainError,
    SolverError,
    StandardsApplicabilityError,
    UnsupportedConversionError,
    assert_registry_compatible,
    convert,
    convert_bulk_ore,
    convert_chemical_strength,
    convert_commodity_quote,
    convert_currency,
    convert_grain,
    convert_moisture_mass,
    convert_specification,
    convert_viscosity,
    cst_to_sfs,
    cst_to_sus,
    energy_to_gge,
    gas_volume_to_energy,
    gge_to_energy,
    parse_quote,
    parse_unit,
    quantity_from_dict,
    quantity_to_dict,
    resolve_currency,
    sg_to_brix,
    sg_to_plato,
    triangulate_fx_rate,
    ureg,
)
from core.units.converters.petroleum import (
    LiquidApplicability,
    LiquidProduct,
    LiquidProductSpec,
    LiquidState,
    constant_coefficient_liquid_method,
    infer_liquid_reference_density,
)
from core.units.converters.viscosity import DomainError as ViscosityDomainError


class TestDomainExceptionHierarchy(unittest.TestCase):
    """Verify MRO, package re-exports, and end-to-end raising across core.units."""

    def test_all_requested_exceptions_exported_and_inherit_commodity_unit_error(self) -> None:
        expected_classes = (
            AmbiguousUnitError,
            MissingAssayError,
            MissingReferenceConditionsError,
            IncompatibleCommodityError,
            IncompatibleAnalyticalBasisError,
            UnsupportedConversionError,
            InvalidPhysicalStateError,
            ConversionConvergenceError,
            MissingFXRateError,
            RegistryVersionMismatchError,
            RegistryFrozenError,
            StandardsApplicabilityError,
        )
        self.assertTrue(issubclass(CommodityUnitError, Exception))
        self.assertIs(getattr(cu, "CommodityUnitError"), CommodityUnitError)
        for exc_cls in expected_classes:
            with self.subTest(exc=exc_cls.__name__):
                self.assertIs(getattr(cu, exc_cls.__name__), exc_cls)
                self.assertTrue(issubclass(exc_cls, CommodityUnitError))

        # Existing specialized exceptions also integrate cleanly into CommodityUnitError
        self.assertTrue(issubclass(ForeignRegistryError, CommodityUnitError))
        self.assertTrue(issubclass(RegistryCompatibilityError, RegistryVersionMismatchError))
        self.assertTrue(issubclass(SolverError, CommodityUnitError))
        self.assertTrue(issubclass(SolverConvergenceError, ConversionConvergenceError))
        self.assertTrue(issubclass(SolverDomainError, InvalidPhysicalStateError))
        self.assertTrue(issubclass(ViscosityDomainError, StandardsApplicabilityError))

        # Backward-compatible Python builtin inheritance
        self.assertTrue(issubclass(RegistryFrozenError, RuntimeError))
        self.assertTrue(issubclass(ConversionConvergenceError, ArithmeticError))
        self.assertTrue(issubclass(ConversionConvergenceError, ValueError))
        for val_cls in (
            AmbiguousUnitError,
            MissingAssayError,
            MissingReferenceConditionsError,
            IncompatibleCommodityError,
            IncompatibleAnalyticalBasisError,
            UnsupportedConversionError,
            InvalidPhysicalStateError,
            MissingFXRateError,
            RegistryVersionMismatchError,
            StandardsApplicabilityError,
        ):
            self.assertTrue(issubclass(val_cls, ValueError))

    def test_ambiguous_unit_error(self) -> None:
        # Registry strictly blocked ambiguous tokens
        for bad_unit in ("ppt", "mu", "ton_mile"):
            with self.subTest(unit=bad_unit):
                with self.assertRaises(AmbiguousUnitError) as ctx:
                    parse_unit(bad_unit)
                self.assertIsInstance(ctx.exception, CommodityUnitError)
                self.assertIsInstance(ctx.exception, ValueError)

        # Quote parser ambiguous token without source conventions
        for quote_str in ("12.50 c/lb", "100 USD/ton", "50 USD/gal", "2500 USD/oz"):
            with self.subTest(quote=quote_str):
                with self.assertRaises(AmbiguousUnitError):
                    parse_quote(quote_str)

        # Currency resolver ambiguous symbol without currency_hint
        with self.assertRaises(AmbiguousUnitError):
            resolve_currency("$")

    def test_missing_assay_error(self) -> None:
        # Cross-dimensional quote without assay
        with self.assertRaises(MissingAssayError):
            convert_commodity_quote(Q_(75.0, "USD / oil_barrel"), "USD / metric_ton")

        # BulkOreAssay wet/dry without moisture
        with self.assertRaises(MissingAssayError):
            BulkOreAssay(grades={"Fe": Decimal("0.62")}).wmt_to_dmt(Decimal("100"))

        # OreAssay volume without bulk density
        ore_no_rho = OreAssay(name="Test Ore", primary_element="Fe", grade_pct=62.0, moisture_pct=8.0)
        with self.assertRaises(MissingAssayError):
            ore_no_rho.wmt_to_volume(Q_(100.0, "wmt"))

        # MaterialProperties density conversion without density
        with self.assertRaises(MissingAssayError):
            MaterialProperties().volume_to_mass(Q_(1.0, "m^3"))

        # HeatContentAssay without GCV or NCV
        with self.assertRaises(MissingAssayError):
            HeatContentAssay()

        # Chemical strength molarity <-> normality without equivalence_factor
        with self.assertRaises(MissingAssayError):
            convert_chemical_strength(Q_(1.0, "molarity"), "normality")

        # Dynamic <-> kinematic viscosity without density
        with self.assertRaises(MissingAssayError):
            convert_viscosity(Q_(10.0, "centipoise"), "centistokes")

    def test_missing_reference_conditions_error(self) -> None:
        # Geometric gas volume to energy without reference_conditions
        with self.assertRaises(MissingReferenceConditionsError):
            gas_volume_to_energy(Q_(100.0, "m^3"), Q_(38.0, "MJ / m^3"), basis="HHV")

        # Partial source_reference without target_reference in master router
        ref = cu.GasReferenceConditions.from_standard("iso15")
        with self.assertRaises(MissingReferenceConditionsError):
            convert(Q_(100.0, "m^3"), "Sm3", source_reference=ref)

        # Historical currency without explicit as_of date
        with self.assertRaises(MissingReferenceConditionsError):
            resolve_currency("HRK")

        # SUS viscosity conversion without temperature
        with self.assertRaises(MissingReferenceConditionsError):
            convert_viscosity(Q_(50.0, "centistokes"), "SUS")

        # LPG reference state without saturation pressure
        state = LiquidState(Q_(15.0, "degC"), Q_(101.325, "kPa"))
        with self.assertRaises(MissingReferenceConditionsError):
            LiquidProductSpec("Propane", LiquidProduct.LPG, state, Q_(507.0, "kg/m^3"))

    def test_incompatible_commodity_error(self) -> None:
        # Blending iron ore with copper concentrate
        fe_ore = OreAssay.iron_ore_62fe()
        cu_conc = OreAssay.copper_concentrate()
        with self.assertRaises(IncompatibleCommodityError):
            fe_ore.blend_with(cu_conc, self_wet_mass_fraction=0.5)

        # Non-lithium battery material LCE conversion
        niso4 = BatteryMaterialAssay.nickel_sulphate_hexahydrate()
        with self.assertRaises(IncompatibleCommodityError):
            niso4.dmt_to_lce(Q_(10.0, "dmt"))

        # Passing DGE GallonEquivalent into gge_to_energy
        dge = cu.energy_to_dge(Q_(1.0, "MMBtu_it"), "NIST_HB44", basis="LHV")
        with self.assertRaises(IncompatibleCommodityError):
            gge_to_energy(dge)

        # Egg method called on broiler poultry assay
        broiler = PoultryEggAssay.broiler_chicken_standard()
        with self.assertRaises(IncompatibleCommodityError):
            broiler.eggs_to_shell_mass(12)

        # CO2 purity mismatch across assays
        food_co2 = CarbonDioxideAssay(
            name="Food CO2",
            state="gas",
            purity_fraction=Decimal("0.999"),
            density=Q_(1.98, "kg/m^3"),
            reference_temperature=Q_(273.15, "K"),
            reference_pressure=Q_(101.325, "kPa"),
        )
        ind_co2 = CarbonDioxideAssay(
            name="Industrial CO2",
            state="gas",
            purity_fraction=Decimal("0.990"),
            density=Q_(1.98, "kg/m^3"),
            reference_temperature=Q_(273.15, "K"),
            reference_pressure=Q_(101.325, "kPa"),
        )
        with self.assertRaises(IncompatibleCommodityError):
            food_co2.equivalent_volume(Q_(100.0, "m^3"), target=ind_co2)

    def test_incompatible_analytical_basis_error(self) -> None:
        # Energy basis mismatch (HHV -> LHV GGE without basis_conversion_ratio)
        with self.assertRaises(IncompatibleAnalyticalBasisError):
            energy_to_gge(Q_(1.0, "MMBtu_it"), "NIST_HB44", basis="HHV")

        # Dry-basis HeatContentAssay declaring nonzero moisture
        with self.assertRaises(IncompatibleAnalyticalBasisError):
            HeatContentAssay(ncv=Q_(25.0, "MJ/kg"), basis="DB", moisture_pct=10.0)

        # Conflicting semantic target unit tag in moisture conversion (wet_to_dry targeting wmt)
        with self.assertRaises(IncompatibleAnalyticalBasisError):
            cu.wet_to_dry_mass(Q_(100.0, "kg"), moisture_pct=10.0, target_unit="wmt")

        # MaterialComposition non-gravimetric projection (MOLAR -> DRY)
        comp = MaterialComposition(components={"CH4": Decimal("0.95")}, basis=CompositionBasis.MOLAR)
        with self.assertRaises(IncompatibleAnalyticalBasisError):
            comp.to_basis(CompositionBasis.DRY)

    def test_unsupported_conversion_error(self) -> None:
        # Unsupported processing mode in agricultural converter
        corn = GrainAssay.corn_56lb()
        with self.assertRaises(UnsupportedConversionError):
            convert_grain(Q_(100.0, "metric_ton"), "metric_ton", assay=corn, mode="nonexistent_mode")

        # Unsupported mode in convert_specification
        spec = MaterialSpecification(name="Silica Sand", density=Q_(1600, "kg/m^3"))
        with self.assertRaises(UnsupportedConversionError):
            convert_specification(Q_(1.0, "m^3"), "kg", assay=spec, mode="unknown_mode")

        # Heat-content overrides supplied to a non-energy conversion
        with self.assertRaises(UnsupportedConversionError):
            convert(Q_(100.0, "kg"), "metric_ton", gcv=Q_(25.0, "MJ/kg"))

        # Mismatched mass exponents in convert_bulk_ore
        with self.assertRaises(UnsupportedConversionError):
            convert_bulk_ore(Q_(100.0, "wmt"), "1 / dmt", moisture_fraction=Decimal("0.08"))

    def test_invalid_physical_state_error(self) -> None:
        # Negative physical mass in master router
        with self.assertRaises(InvalidPhysicalStateError):
            convert(Q_(-5.0, "metric_ton"), "kilogram")

        # Out-of-range moisture fraction >= 1.0
        with self.assertRaises(InvalidPhysicalStateError):
            convert_moisture_mass(100.0, source_moisture=1.0, target_moisture=0.1)

        # NCV > GCV thermodynamic violation
        with self.assertRaises(InvalidPhysicalStateError):
            HeatContentAssay(gcv=Q_(20.0, "MJ/kg"), ncv=Q_(22.0, "MJ/kg"))

        # Liquid pressure below saturation pressure
        with self.assertRaises(InvalidPhysicalStateError):
            LiquidState(
                temperature=Q_(20.0, "degC"),
                absolute_pressure=Q_(100.0, "kPa"),
                saturation_pressure=Q_(250.0, "kPa"),
            )

        # Real-gas solver out-of-domain physical state (negative pressure)
        with self.assertRaises(InvalidPhysicalStateError) as ctx:
            cu.solve_z_factor(Q_(-50.0, "kPa"), Q_(300.0, "K"), 0.60)
        self.assertIsInstance(ctx.exception, SolverDomainError)

    def test_conversion_convergence_error(self) -> None:
        ref = LiquidState(Q_(15.0, "degC"), Q_(101.325, "kPa"))
        obs = LiquidState(Q_(35.0, "degC"), Q_(101.325, "kPa"))
        app = LiquidApplicability(
            temperature_range=(Q_(0.0, "degC"), Q_(60.0, "degC")),
            pressure_range=(Q_(50.0, "kPa"), Q_(500.0, "kPa")),
            density_range=(Q_(800.0, "kg/m^3"), Q_(900.0, "kg/m^3")),
        )
        method = constant_coefficient_liquid_method(
            product=LiquidProduct.CRUDE,
            thermal_expansion=Q_("0.0008 / K"),
            compressibility=Q_("1e-9 / Pa"),
            applicability=app,
            source="test",
        )
        # Unbracketed measured density (far outside [800, 900] kg/m^3)
        with self.assertRaises(ConversionConvergenceError) as ctx:
            infer_liquid_reference_density(
                Q_(650.0, "kg/m^3"),
                name="Light Condensate",
                product=LiquidProduct.CRUDE,
                reference=ref,
                observed=obs,
                method=method,
            )
        self.assertIsInstance(ctx.exception, CommodityUnitError)
        self.assertIsInstance(ctx.exception, ArithmeticError)
        self.assertIsInstance(ctx.exception, ValueError)

    def test_missing_fx_rate_error(self) -> None:
        # convert_currency missing target rate
        with self.assertRaises(MissingFXRateError):
            convert_currency(Q_(100.0, "USD"), "EUR", rates={"GBP": Decimal("0.79")})

        # triangulate_fx_rate disconnected currency graph
        with self.assertRaises(MissingFXRateError):
            triangulate_fx_rate({"USD/EUR": Decimal("0.92"), "JPY/KRW": Decimal("9.0")}, "USD", "JPY")

        # convert_commodity_quote cross-currency without rates
        with self.assertRaises(MissingFXRateError):
            convert_commodity_quote(Q_(80.0, "USD / oil_barrel"), "EUR / oil_barrel")

        # MaterialProperties.normalize_price cross-currency without fx_rate
        with self.assertRaises(MissingFXRateError):
            MaterialProperties().normalize_price(Q_(100.0, "USD / metric_ton"), "EUR / metric_ton")

    def test_registry_version_mismatch_error(self) -> None:
        with self.assertRaises(RegistryVersionMismatchError) as ctx:
            assert_registry_compatible(expected_registry_version="999.0.0")
        self.assertIsInstance(ctx.exception, RegistryCompatibilityError)
        self.assertIsInstance(ctx.exception, CommodityUnitError)

        payload = quantity_to_dict(Q_(42.0, "oil_barrel"), include_metadata=True)
        payload["definition_hash"] = "0" * 64
        with self.assertRaises(RegistryVersionMismatchError):
            quantity_from_dict(payload)

    def test_registry_frozen_error_in_isolated_subprocess(self) -> None:
        code = """
import core.units as cu
cu.freeze_registry()
assert cu.is_registry_frozen()
try:
    cu.register_currency("XTS9", decimals=2)
except cu.RegistryFrozenError as exc:
    assert isinstance(exc, cu.CommodityUnitError)
    assert isinstance(exc, RuntimeError)
else:
    raise AssertionError("Expected RegistryFrozenError on register_currency")

try:
    cu.ureg.define("forbidden_custom_unit = 123 * meter")
except cu.RegistryFrozenError as exc:
    assert isinstance(exc, cu.CommodityUnitError)
else:
    raise AssertionError("Expected RegistryFrozenError on ureg.define")
"""
        proc = subprocess.run(
            [sys.executable, "-c", code],
            capture_output=True,
            text=True,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, msg=f"STDOUT:\n{proc.stdout}\nSTDERR:\n{proc.stderr}")
        self.assertFalse(cu.is_registry_frozen())

    def test_standards_applicability_error(self) -> None:
        # ASTM D2161 below SUS minimum kinematic viscosity (1.81 cSt)
        with self.assertRaises(StandardsApplicabilityError) as ctx:
            cst_to_sus(1.2)
        self.assertIsInstance(ctx.exception, ViscosityDomainError)

        # ASTM D2161 SFS outside 122 °F reference temperature
        with self.assertRaises(StandardsApplicabilityError):
            cst_to_sfs(100.0, temp_f=100.0)

        # ASBC Plato above 30 °P ceiling
        with self.assertRaises(StandardsApplicabilityError):
            sg_to_plato(1.18)

        # ICUMSA Brix non-20 °C hydrometer reference
        with self.assertRaises(StandardsApplicabilityError):
            sg_to_brix(1.05, ref_temp_c=15.0)

        # Vegetable oil VCF outside [-20, 100] °C applicability
        cpo = OilseedCrushAssay.crude_palm_oil_cpo()
        with self.assertRaises(StandardsApplicabilityError):
            cpo.volume_correction_factor(145.0)


if __name__ == "__main__":
    unittest.main()
