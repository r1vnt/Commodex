"""All 36 additional sector specifications and their public conversion APIs."""

from dataclasses import replace
from decimal import Decimal as D
import inspect
import re
import unittest

import pint
import core.units as u
from core.units import material_specs as metals, chemical_specs as chemicals, energy_environment_specs as energy


def gas_inputs(component="N2"):
    backend = u.builtin_gas_backend("ideal", temperature_range=(u.Q_(270, "K"), u.Q_(400, "K")),
                                    pressure_range=(u.Q_(50000, "Pa"), u.Q_(200000, "Pa")))
    return dict(gas_composition=u.GasComposition({component: 1}),
                reference=u.GasReferenceConditions.from_standard("iso15", compressibility_factor=D(1)), backend=backend)


def examples():
    density = u.Q_(D(2000), "kg/m^3")
    composition = u.ElementComposition({"Cu": D(".25"), "Zn": D(".15"), "Fe": D(".5")})
    shared = dict(name="declared test material", density=density, composition=composition, moisture_fraction=D(".1"))
    yield u.RefinedMetalAssay(**shared, metal="Cu")
    yield u.MetalScrapAssay(**shared, scrap_grade="shredded")
    yield u.SteelProductAssay(**shared, product="HRC", thickness=u.Q_(2, "mm"), width=u.Q_(1, "m"))
    yield u.DirectReducedIronAssay(name="HBI", density=density, total_iron_fraction=D(".92"), metallisation_fraction=D(".94"))
    yield u.IronPelletAssay(**shared, product="pellet")
    yield u.MultiMetalConcentrateAssay(**shared)
    yield u.RareEarthAssay(name="NdPr oxide", density=density, oxide_fractions={"Nd2O3": D(".6"), "Pr6O11": D(".2")})
    yield u.CriticalMineralAssay(**shared, commodity="tungsten", reported_compound="WO3")
    yield u.PreciousMetalConcentrateAssay(**shared)
    yield u.RefractoryMaterialAssay(**shared, loss_on_ignition_fraction=D(".05"))
    yield u.IndustrialMineralAssay(**shared, mineral="barite")
    yield u.CementClinkerAssay(**shared, clinker_fraction=D(".8"))
    yield u.AggregateAssay(**shared)
    yield u.GemstoneAssay(**shared)
    yield u.BatteryRecyclingAssay(**shared)
    yield u.PolymerAssay(name="PET", density=u.Q_(D(1380), "kg/m^3"), polymer="PET")
    yield u.RecycledPolymerAssay(name="rPET", density=u.Q_(D(1300), "kg/m^3"), polymer="PET", contamination_fraction=D(".1"), wash_recovery_fraction=D(".9"), melt_recovery_fraction=D(".95"))
    yield u.IndustrialGasAssay(name="nitrogen", **gas_inputs())
    yield u.RefrigerantAssay(name="declared blend", density=u.Q_(D(1100), "kg/m^3"), blend_mass_fractions={"R32": D(".5"), "R125": D(".5")}, component_molar_masses={"R32": u.Q_(52, "g/mol"), "R125": u.Q_(120, "g/mol")})
    solution = dict(name="declared liquid", solute_formula="NaOH", active_mass_fraction=D(".5"), density=u.Q_(D(1530), "kg/m^3"), reference_temperature=u.Q_(20, "degC"))
    yield u.ChemicalSolutionAssay(**solution)
    yield u.SolventAssay(**{**solution, "solute_formula": "C3H8O", "active_mass_fraction": D(".99"), "density": u.Q_(D(790), "kg/m^3")})
    yield u.SulphurProductAssay(name="molten S", density=u.Q_(D(1800), "kg/m^3"), state="molten", sulphur_mass_fraction=D(".999"))
    yield u.ChlorAlkaliAssay(name="salt", density=density, salt_purity_fraction=D(".99"), salt_conversion_fraction=D(".9"))
    yield u.AmmoniaAssay(name="aqueous NH3", density=u.Q_(D(910), "kg/m^3"), active_mass_fraction=D(".25"), state="aqueous")
    yield u.GraphiteAssay(name="flake", density=density, fixed_carbon_fraction=D(".95"), product_form="flake")
    yield u.ExplorationBrineAssay(name="Li brine", density=u.Q_(D(1200), "kg/m^3"), concentrations={"Li": u.Q_(D(500), "mg/liter")}, extraction_recoveries={"Li": D(".8")})
    yield u.BiogasAssay(name="biogas", **gas_inputs("CH4"))
    yield u.BiomethaneAssay(name="RNG", **gas_inputs("CH4"))
    heat = u.HeatContentAssay(gcv=u.Q_(D(22), "MJ/kg"), ncv=u.Q_(D(20), "MJ/kg"))
    yield u.SyntheticFuelAssay(name="e-diesel", fuel="e-diesel", density=u.Q_(D(800), "kg/m^3"), heat_content=heat)
    yield u.WasteDerivedFuelAssay(name="RDF", fuel="RDF", density=u.Q_(D(350), "kg/m^3"), moisture_fraction=D(".2"), inert_dry_fraction=D(".15"), heat_content=heat)
    yield u.CarbonDioxideAssay(name="liquid CO2", state="liquid", purity_fraction=D(".999"), density=u.Q_(D(1000), "kg/m^3"), reference_temperature=u.Q_(250, "K"), reference_pressure=u.Q_(2, "MPa"))
    yield u.MethanolFuelAssay(name="marine methanol", active_mass_fraction=D(".99"), density=u.Q_(D(790), "kg/m^3"), heat_content=heat)
    h2_heat = u.HeatContentAssay(gcv=u.Q_(D(142), "MJ/kg"), ncv=u.Q_(D(120), "MJ/kg"))
    yield u.HydrogenCarrierAssay(name="LOHC", carrier="LOHC", density=u.Q_(D(850), "kg/m^3"), h2_mass_fraction=D(".06"), cracking_recovery_fraction=D(".9"), hydrogen_heat_content=h2_heat)
    yield u.EnergyStorageSpecification(name="battery", charge_efficiency=D(".95"), discharge_efficiency=D(".9"), standing_loss_fraction_per_hour=D(0))
    yield u.EnvironmentalUnitSpecification(scheme="example-scheme", vintage="2026", certificate_type="energy", issuer="example issuer", quantity_per_certificate=u.Q_(D(1), "MWh"))
    yield u.WaterAssay(name="saline water", density=u.Q_(D(1025), "kg/m^3"), salinity_mass_fraction=D(".035"), dissolved_solids_fraction=D(".04"))


class SectorSpecificationTests(unittest.TestCase):
    def test_all_36_models_named_apis_exports_and_round_trips(self):
        models = list(examples())
        declared = {value for module in (metals, chemicals, energy) for value in vars(module).values()
                    if inspect.isclass(value) and value.__name__ in module.__all__ and value.__name__.endswith(("Assay", "Specification"))}
        self.assertEqual(len(declared), 36)
        self.assertEqual({type(model) for model in models}, declared)
        for module in (metals, chemicals, energy):
            self.assertEqual(len(module.__all__), len(set(module.__all__)))
            for name in module.__all__:
                self.assertIs(getattr(u, name), getattr(module, name))
        for model in models:
            stem = re.sub(r"(?<!^)(?=[A-Z])", "_", type(model).__name__.removesuffix("Assay").removesuffix("Specification")).lower()
            converter = getattr(u, "convert_" + stem)
            factory = getattr(u, "create_" + stem + "_context")
            self.assertIsInstance(factory(model), pint.Context)
            for numeric in (1.25, D("1.25")):
                if isinstance(model, u.EnergyStorageSpecification):
                    original, target = u.Q_(numeric, "MWh"), "MMBtu"
                elif isinstance(model, u.EnvironmentalUnitSpecification):
                    original, target = u.Q_(numeric, "dimensionless"), "MWh"
                else:
                    original, target = u.Q_(numeric, "m^3"), "kg"
                with self.subTest(model=type(model).__name__, numeric=type(numeric)):
                    output = converter(original, target, assay=model)
                    restored = u.convert(output, original.units, assay=model)
                    self.assertAlmostEqual(float(restored.magnitude), float(numeric), places=10)
                    self.assertIsInstance(output.magnitude, type(numeric))
                    self.assertIsInstance(restored.magnitude, type(numeric))

    def test_refined_metals_and_all_polymer_families(self):
        for symbol in ("Cu", "Zn", "Ni", "Pb", "Sn", "Al"):
            assay = u.RefinedMetalAssay(name=symbol, metal=symbol, composition=u.ElementComposition({symbol: D(".999")}))
            fine = assay.fine_mass(u.Q_(D(1000), "kg"))
            self.assertEqual(fine.magnitude, 999)
            self.assertEqual(assay.gross_mass_from_fine(fine).magnitude, 1000)
        for family in ("PE", "HDPE", "LDPE", "LLDPE", "PP", "PVC", "PET", "PS", "ABS"):
            assay = u.PolymerAssay(name=family, polymer=family)
            self.assertEqual(assay.pellet_mass(u.Q_(D(100), "kg"), process_yield_fraction=D(".98")).magnitude, 98)

    def test_scrap_dri_beneficiation_and_black_mass_balances(self):
        composition = u.ElementComposition({"Cu": D(".8"), "Al": D(".1")})
        scrap = u.MetalScrapAssay(name="scrap", scrap_grade="copper", composition=composition, recovery_fractions={"Cu": D(".9"), "Al": D(".8")})
        quantity = u.Q_(D(1000), "kg")
        self.assertEqual(scrap.unrecovered_mass(quantity).magnitude, 200)
        dri = u.DirectReducedIronAssay(name="HBI", total_iron_fraction=D(".92"), metallisation_fraction=D(".94"))
        self.assertEqual(dri.metallic_iron_mass(quantity).magnitude, D("864.8"))
        self.assertEqual((dri.metallic_iron_mass(quantity) + dri.oxide_bound_iron_mass(quantity)).magnitude, 920)
        pellet = u.IronPelletAssay(name="fines", product="fines", composition=u.ElementComposition({"Fe": D(".5")}))
        product = pellet.beneficiated_product_mass(quantity, target_iron_fraction=D(".65"), iron_recovery_fraction=D(".9"))
        self.assertAlmostEqual(product.magnitude * D(".65"), 450, places=20)
        battery = u.BatteryRecyclingAssay(name="black mass", composition=u.ElementComposition({"Li": D(".05"), "Ni": D(".2"), "Co": D(".1")}), recovery_fractions={"Li": D(".8"), "Ni": D(".9"), "Co": D(".95")})
        black = battery.black_mass_from_pack(quantity, black_mass_yield_fraction=D(".4"))
        self.assertEqual(black.magnitude, 400)
        self.assertEqual(battery.recovered_mass(black, component="Li").magnitude, 16)
        self.assertEqual(battery.recovered_mass(black, component="Ni").magnitude, 72)

    def test_steel_geometry_coating_and_gemstone_carat_recovery(self):
        flat = u.SteelProductAssay(name="HRC", product="HRC", density=u.Q_(D(7850), "kg/m^3"), width=u.Q_(D(1), "m"), thickness=u.Q_(D(2), "mm"), coating_mass_per_area=u.Q_(D(275), "g/m^2"))
        self.assertEqual(flat.piece_mass(u.Q_(D(10), "m")).magnitude, 157)
        self.assertEqual(flat.coating_mass(u.Q_(D(20), "m^2")).magnitude, D("5.5"))
        round_product = replace(flat, product="rebar", width=None, thickness=None, diameter=u.Q_(D(10), "mm"))
        self.assertAlmostEqual(float(round_product.piece_mass(u.Q_(D(1), "m")).magnitude), .6165375583, places=9)
        gems = u.GemstoneAssay(name="rough")
        polished = gems.polished_mass(u.Q_(D(100), "carat"), recovery_fraction=D(".4"))
        self.assertEqual(polished.to("carat").magnitude, 40)
        self.assertEqual(gems.rough_mass(polished, recovery_fraction=D(".4")).to("carat").magnitude, 100)

    def test_rare_earth_critical_mineral_and_graphite_conversions(self):
        assay = u.RareEarthAssay(name="TREO", oxide_fractions={"Nd2O3": D(".4"), "Pr6O11": D(".1"), "Dy2O3": D(".05"), "Tb4O7": D(".01")})
        sample = u.Q_(D(1000), "kg")
        self.assertEqual(assay.treo_mass(sample).magnitude, 560)
        self.assertGreater(assay.ndpr_element_mass(sample).magnitude, 400)
        self.assertLess(assay.ndpr_element_mass(sample).magnitude, 500)
        tungsten = u.CriticalMineralAssay(name="WO3", commodity="tungsten", reported_compound="WO3")
        self.assertAlmostEqual(float(tungsten.reported_to_element(sample, element="W").magnitude), 1000 * 183.84 / 231.837, places=9)
        graphite = u.GraphiteAssay(name="flake", fixed_carbon_fraction=D(".95"), product_form="flake")
        result = graphite.purified_mass(sample, target_fraction=D(".9995"), recovery_fraction=D(".9"))
        self.assertAlmostEqual(result.magnitude * D(".9995"), 855, places=20)

    def test_calcination_cement_refractory_and_aggregate_inverses(self):
        cement = u.CementClinkerAssay(name="cement", clinker_fraction=D(".75"))
        sample = u.Q_(D(1000), "kg")
        solids = cement.calcined_solid_mass(sample, limestone_purity_fraction=D(".95"), conversion_fraction=D(".9"))
        gas = cement.calcination_co2_mass(sample, limestone_purity_fraction=D(".95"), conversion_fraction=D(".9"))
        self.assertEqual((solids + gas).magnitude, 1000)
        self.assertAlmostEqual(cement.clinker_from_cement(cement.cement_from_clinker(sample)).magnitude, 1000, places=20)
        refractory = u.RefractoryMaterialAssay(name="magnesia", loss_on_ignition_fraction=D(".1"))
        self.assertAlmostEqual(refractory.uncalcined_mass(refractory.calcined_mass(sample)).magnitude, 1000, places=20)
        aggregate = u.AggregateAssay(name="sand")
        volume = u.Q_(D(100), "m^3")
        compacted = aggregate.compacted_volume(volume, volume_reduction_fraction=D(".15"))
        self.assertEqual(compacted.magnitude, 85)
        self.assertEqual(aggregate.loose_volume(compacted, volume_reduction_fraction=D(".15")).magnitude, 100)

    def test_solution_active_molarity_dilution_and_legacy_routes(self):
        solution = u.ChemicalSolutionAssay(name="caustic", solute_formula="NaOH", active_mass_fraction=D(".5"), density=u.Q_(D(1530), "kg/m^3"), reference_temperature=u.Q_(20, "degC"))
        original = u.Q_(D(1000), "kg")
        self.assertEqual(solution.active_mass(original).magnitude, 500)
        self.assertEqual(solution.diluent_required(original, target_fraction=D(".25")).magnitude, 1000)
        self.assertEqual(solution.diluted_mass(original, target_fraction=D(".25")).magnitude, 2000)
        self.assertAlmostEqual(float(solution.molarity.magnitude), 765 / 39.997, places=10)
        moles = solution.active_to_moles(solution.active_mass(original))
        self.assertAlmostEqual(solution.moles_to_active(moles).magnitude, 500, places=20)
        legacy = solution.to_chemical_assay()
        for numeric in (1.25, D("1.25")):
            volume = u.Q_(numeric, "m^3")
            self.assertAlmostEqual(float(u.convert(volume, "kg", assay=solution).magnitude), float(u.convert(volume, "kg", assay=legacy, active_only=False).magnitude), places=9)
            self.assertAlmostEqual(float(solution.active_mass(solution.volume_to_mass(volume)).magnitude), float(legacy.volume_to_active_mass(volume).to("kg").magnitude), places=9)
            self.assertAlmostEqual(float(u.convert(volume, "mol", assay=solution).magnitude), float(u.convert(volume, "mol", assay=legacy, active_only=True).magnitude), places=7)
        measured_legacy = replace(legacy, molecular_weight_g_per_mol=40.)
        adapted = u.ChemicalSolutionAssay.from_chemical_assay(measured_legacy)
        self.assertEqual(adapted.molar_mass.magnitude, 40)
        self.assertEqual(u.convert(u.Q_(D(1000), "kg"), "kmol", assay=adapted).magnitude, D("12.5"))
        solvent = u.SolventAssay(name="IPA", solute_formula="C3H8O", active_mass_fraction=D(".99"), density=u.Q_(D(790), "kg/m^3"), reference_temperature=u.Q_(20, "degC"))
        petro = solvent.to_petrochem_assay()
        self.assertEqual(u.convert(u.Q_(D(1), "m^3"), "kg", assay=solvent).magnitude, u.convert(u.Q_(D(1), "m^3"), "kg", assay=petro).magnitude)
        self.assertEqual(u.SolventAssay.from_petrochem_assay(petro).active_mass_fraction, D(".99"))

    def test_refrigerants_brines_and_chlor_alkali_stoichiometry(self):
        blend = u.RefrigerantAssay(name="HFC blend", blend_mass_fractions={"R32": D(".5"), "R125": D(".5")}, component_molar_masses={"R32": u.Q_(D(52), "g/mol"), "R125": u.Q_(D(120), "g/mol")})
        self.assertAlmostEqual(blend.molar_mass.magnitude, D(1) / (D(".5") / 52 + D(".5") / 120), places=20)
        self.assertAlmostEqual(sum(blend.mole_fractions.values()), 1, places=25)
        brine = u.ExplorationBrineAssay(name="Li", concentrations={"Li": u.Q_(D(500), "mg/liter")}, extraction_recoveries={"Li": D(".8")})
        volume = u.Q_(D(1000), "m^3")
        self.assertAlmostEqual(brine.dissolved_mass(volume, component="Li").magnitude, 500, places=10)
        output = brine.recoverable_mass(volume, component="Li")
        self.assertAlmostEqual(output.magnitude, 400, places=10)
        self.assertAlmostEqual(brine.brine_required(output, component="Li").magnitude, 1000, places=20)
        chlor = u.ChlorAlkaliAssay(name="cell", salt_purity_fraction=D(1), salt_conversion_fraction=D(1))
        salt = u.Q_(D(1000), "kg")
        products = sum((chlor.co_product_mass(salt, product=p, recovery_fraction=D(1)).magnitude for p in ("Cl2", "NaOH", "H2")), D(0))
        self.assertAlmostEqual(products, D(1000) + chlor.water_consumed(salt).magnitude, places=20)
        sulphur = u.SulphurProductAssay(name="S", state="solid", sulphur_mass_fraction=D(1))
        acid = sulphur.acid_equivalent(u.Q_(D(32.06), "kg"), acid_mass_fraction=D(1))
        self.assertAlmostEqual(acid.magnitude, D("98.072"), places=10)

    def test_gas_upgrading_methane_and_total_mass_are_conserved(self):
        inputs = gas_inputs()
        wet = u.GasComposition({"CH4": D(".588"), "CO2": D(".392"), "H2O": D(".02")})
        feed = u.BiogasAssay(name="wet biogas", **{**inputs, "gas_composition": wet, "reference": replace(inputs["reference"], water_mole_fraction=D(".02"))})
        result = feed.upgrade(u.Q_(D(1000), "m^3"), methane_recovery_fraction=D(".95"), target_methane_fraction=D(".98"), target_reference=inputs["reference"], target_backend=inputs["backend"], target_unit="Sm3")
        self.assertEqual(result.product_composition.water_mole_fraction, 0)
        self.assertEqual(result.product_composition.mole_fractions["CH4"], D(".98"))
        self.assertAlmostEqual((result.product_mass + result.offgas_mass + result.removed_water_mass).magnitude, result.feed_mass.magnitude, places=20)
        input_ch4_moles = result.feed_mass.magnitude / wet.molar_mass.to("kg/mol").magnitude * D(".588")
        output_ch4_moles = result.product_mass.magnitude / result.product_composition.molar_mass.to("kg/mol").magnitude * D(".98")
        self.assertAlmostEqual(output_ch4_moles / input_ch4_moles, D(".95"), places=25)
        rng = u.BiomethaneAssay(name="RNG", **gas_inputs("CH4"))
        self.assertTrue(all(rng.injection_quality(minimum_methane=D(".98"), maximum_co2=D(".02"), maximum_h2s=D(".00001"), maximum_water=D(".0001")).values()))
        nitrogen = u.IndustrialGasAssay(name="N2", **inputs)
        tagged_mass = u.convert(u.Q_(D(1000), "Sm3"), "kg", assay=nitrogen)
        self.assertGreater(tagged_mass.magnitude, 1000)
        explicit = u.convert_industrial_gas(tagged_mass, "Sm3", assay=nitrogen, mode="mass_to_gas_volume")
        self.assertAlmostEqual(explicit.magnitude, D(1000), places=20)
        moles = u.convert(u.Q_(D(1000), "Sm3"), "kmol", assay=nitrogen)
        self.assertAlmostEqual(float(u.convert(moles, "Sm3", assay=nitrogen).magnitude), 1000, places=9)
        with self.assertRaises(ValueError):
            u.convert(u.Q_(D(1), "Nm3"), "kg", assay=nitrogen)

    def test_fuels_and_hydrogen_carriers_reuse_heat_content(self):
        heat = u.HeatContentAssay(gcv=u.Q_(D(22), "MJ/kg"), ncv=u.Q_(D(20), "MJ/kg"))
        fuel = u.SyntheticFuelAssay(name="e-diesel", fuel="e-diesel", density=u.Q_(D(800), "kg/m^3"), heat_content=heat)
        self.assertEqual(u.convert(u.Q_(D(100), "kg"), "MJ", assay=fuel).magnitude, 2000)
        self.assertEqual(u.convert(u.Q_(D(1), "m^3"), "MJ", assay=fuel, basis="hhv").magnitude, 17600)
        self.assertEqual(u.heat_content_from_assay(fuel).density.magnitude, 800)
        explicit = u.convert(u.Q_(D(1), "kg"), "MJ", assay=u.AggregateAssay(name="declared sample"), ncv=u.Q_(D(18), "MJ/kg"))
        self.assertEqual(explicit.magnitude, 18)
        price = u.convert(u.Q_(D(100), "USD/metric_ton"), "EUR/GJ", assay=fuel, rates={"USD/EUR": D(".9")})
        self.assertEqual(price.magnitude, D("4.5"))
        h2 = u.HeatContentAssay(gcv=u.Q_(D(142), "MJ/kg"), ncv=u.Q_(D(120), "MJ/kg"))
        carrier = u.HydrogenCarrierAssay(name="LOHC", carrier="LOHC", h2_mass_fraction=D(".06"), cracking_recovery_fraction=D(".9"), hydrogen_heat_content=h2)
        quantity = u.Q_(D(1000), "kg")
        self.assertEqual(carrier.delivered_hydrogen(quantity).magnitude, 54)
        self.assertEqual(carrier.carrier_required(carrier.delivered_hydrogen(quantity)).magnitude, 1000)
        self.assertEqual(u.convert(quantity, "MJ", assay=carrier).magnitude, 6480)
        rdf = u.WasteDerivedFuelAssay(name="RDF", fuel="RDF", moisture_fraction=D(".2"), inert_dry_fraction=D(".15"))
        self.assertEqual(rdf.combustible_dry_mass(quantity).magnitude, 680)

    def test_co2_state_equivalence_and_water_mass_balance(self):
        liquid = u.CarbonDioxideAssay(name="CO2", state="liquid", purity_fraction=D(".999"), density=u.Q_(D(1000), "kg/m^3"), reference_temperature=u.Q_(250, "K"), reference_pressure=u.Q_(2, "MPa"))
        gas = replace(liquid, state="gas", density=u.Q_(D(2), "kg/m^3"), reference_temperature=u.Q_(288, "K"), reference_pressure=u.Q_(101325, "Pa"))
        self.assertEqual(liquid.equivalent_volume(u.Q_(D(1), "m^3"), target=gas).magnitude, 500)
        water = u.WaterAssay(name="water", density=u.Q_(D(1025), "kg/m^3"), salinity_mass_fraction=D(".035"), dissolved_solids_fraction=D(".04"))
        quantity = u.Q_(D(1000), "kg")
        self.assertEqual(water.salt_mass(quantity).magnitude, 35)
        self.assertEqual((water.dissolved_solids_mass(quantity) + water.solvent_water_mass(quantity)).magnitude, 1000)
        self.assertEqual(water.dilution_water_required(quantity, target_salinity=D(".0175")).magnitude, 1000)

    def test_storage_losses_and_certificate_identity(self):
        storage = u.EnergyStorageSpecification(name="battery", charge_efficiency=D(".95"), discharge_efficiency=D(".9"), standing_loss_fraction_per_hour=D(".01"), capacity=u.Q_(D(100), "MWh"))
        quantity = u.Q_(D(100), "MWh")
        output = storage.output_energy(quantity, duration=u.Q_(D(2), "hour"))
        self.assertEqual(output.magnitude, D("83.79855"))
        self.assertAlmostEqual(storage.input_required(output, duration=u.Q_(D(2), "hour")).magnitude, 100, places=20)
        with self.assertRaises(ValueError):
            storage.output_energy(quantity)
        with self.assertRaises(ValueError):
            storage.stored_energy(u.Q_(D(1000), "MWh"))
        cert = u.EnvironmentalUnitSpecification(scheme="example", vintage="2026", issuer="issuer", certificate_type="energy", quantity_per_certificate=u.Q_(D(1), "MWh"))
        small = replace(cert, quantity_per_certificate=u.Q_(D(1), "kWh"))
        self.assertEqual(cert.transfer_to(u.Q_(D(1), "dimensionless"), small).magnitude, 1000)
        allowance = replace(cert, scheme="example allowance", certificate_type="allowance", quantity_per_certificate=u.Q_(D(1), "EUA"))
        self.assertEqual(u.convert(u.Q_(D(10), "dimensionless"), "EUA", assay=allowance).magnitude, 10)
        with self.assertRaises(pint.DimensionalityError):
            u.convert(u.Q_(D(10), "dimensionless"), "metric_ton", assay=allowance)
        with self.assertRaises(ValueError):
            cert.transfer_to(u.Q_(D(1), "dimensionless"), replace(small, vintage="2025"))

    def test_missing_process_data_and_invalid_fractions_fail_explicitly(self):
        with self.assertRaises(ValueError):
            u.RefinedMetalAssay(name="unknown", metal="Cu")
        with self.assertRaises(ValueError):
            u.PolymerAssay(name="unknown", polymer="XYZ")
        with self.assertRaises(ValueError):
            u.WaterAssay(name="bad", salinity_mass_fraction=D(".1"), dissolved_solids_fraction=D(".05"))
        with self.assertRaises(ValueError):
            u.EnergyStorageSpecification(name="bad", charge_efficiency=95, discharge_efficiency=D(".9"))
        with self.assertRaises(ValueError):
            u.ExplorationBrineAssay(name="bad", concentrations={"Li": u.Q_(D(2000), "kg/m^3")}, density=u.Q_(D(1000), "kg/m^3"))
        with self.assertRaises(ValueError):
            u.RecycledPolymerAssay(name="bad", polymer="PET", contamination_fraction=-1, wash_recovery_fraction=1, melt_recovery_fraction=1)
        with self.assertRaisesRegex(ValueError, "calorific"):
            u.convert(u.Q_(D(1), "kg"), "MWh", assay=u.AggregateAssay(name="gravel"))

    def test_quantity_metadata_is_immutable_and_dataclass_replacement_remains_valid(self):
        original = u.Q_(D(500), "mg/liter")
        brine = u.ExplorationBrineAssay(name="brine", concentrations={"Li": original})
        original.ito("kg/m^3")
        brine.concentrations["Li"].ito("g/liter")
        clone = replace(brine, name="same sample")
        self.assertAlmostEqual(clone.concentrations["Li"].magnitude, D(".5"), places=12)
        self.assertEqual(clone.concentrations["Li"].units, u.ureg.kg / u.ureg.meter ** 3)
        refrigerant = u.RefrigerantAssay(name="blend", blend_mass_fractions={"R32": D(1)}, component_molar_masses={"R32": u.Q_(D(52), "g/mol")})
        refrigerant.component_molar_masses["R32"].ito("kg/mol")
        self.assertEqual(replace(refrigerant, name="clone").molar_mass.magnitude, D(52))
        denomination = u.Q_(D(1), "MWh")
        certificate = u.EnvironmentalUnitSpecification(scheme="scheme", vintage="2026", issuer="issuer", certificate_type="energy", quantity_per_certificate=denomination)
        denomination.ito("kWh")
        self.assertEqual(certificate.quantity_per_certificate.magnitude, 1)
        self.assertEqual(certificate.quantity_per_certificate.units, u.ureg.MWh)


if __name__ == "__main__":
    unittest.main()
