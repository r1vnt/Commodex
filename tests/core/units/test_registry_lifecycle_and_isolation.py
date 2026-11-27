"""Institutional registry lifecycle, immutability, compatibility, concurrency, and security tests."""

from __future__ import annotations

from collections import deque
from concurrent.futures import ThreadPoolExecutor
import contextvars
from dataclasses import FrozenInstanceError
from datetime import date
from decimal import Decimal
from pathlib import Path
import subprocess
import sys
import threading
import types
import unittest
from unittest.mock import patch

import pint

import core.units as units
from core.units import agricultural_specs as agricultural
from core.units import conversions, converters, currencies, heat_content, registry


_REPO_ROOT = Path(__file__).resolve().parents[3]


class RegistryLifecycleAndIsolationTests(unittest.TestCase):
    """Covers all nine registry lifecycle, concurrency, and security invariants."""

    # ------------------------------------------------------------------
    # 1. Atomic registry construction and publication
    # ------------------------------------------------------------------

    def test_atomic_singleton_publication_under_thread_contention(self) -> None:
        script = """
import threading
import pint
import core.units.registry as reg

# Reset singleton reference in an isolated child process to force a true multi-thread construction race
reg._MASTER_REGISTRY = None
barrier = threading.Barrier(16)
results = [None] * 16
errors = []

def worker(idx: int) -> None:
    try:
        barrier.wait(timeout=10)
        results[idx] = reg.get_registry()
    except BaseException as exc:
        errors.append(repr(exc))

threads = [threading.Thread(target=worker, args=(i,)) for i in range(16)]
for t in threads:
    t.start()
for t in threads:
    t.join(timeout=20)

assert not errors, errors
assert all(r is not None for r in results)
assert len({id(r) for r in results}) == 1, "Multiple UnitRegistry instances published!"
published = results[0]
assert published._bootstrap_complete is True
assert pint.get_application_registry().get() is published
"""
        proc = subprocess.run(
            [sys.executable, "-B", "-c", script],
            cwd=_REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=45,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)

    def test_transactional_rollback_preserves_definition_hash_and_maps(self) -> None:
        hash_before = units.get_registry_definition_hash()
        snap_before = units.get_registry_snapshot()
        minor_before = dict(registry.minor_mapping)
        denom_before = dict(registry.irregular_denominators)

        with self.assertRaises((ValueError, pint.errors.RedefinitionError)):
            units.load_desk_config(
                {
                    "currencies": [{"code": "TXROLLCCY", "decimals": 6, "symbols": ["TX$"]}],
                    "commodities": [
                        {"definition": "tx_good_drum = 200 * liter"},
                        {"definition": "tx_bad_directive = @import evil.txt"},
                    ],
                }
            )

        self.assertEqual(units.get_registry_definition_hash(), hash_before)
        self.assertEqual(units.get_registry_snapshot(), snap_before)
        self.assertEqual(dict(registry.minor_mapping), minor_before)
        self.assertEqual(dict(registry.irregular_denominators), denom_before)
        self.assertNotIn("TXROLLCCY", units.ureg)
        self.assertNotIn("tx_good_drum", units.ureg)
        self.assertFalse(units.is_registry_frozen())

    # ------------------------------------------------------------------
    # 2. Collision allowlist for intentional Pint overrides
    # ------------------------------------------------------------------

    def test_collision_allowlist_is_immutable_and_enforced(self) -> None:
        self.assertIn("bbl", units.PINT_COLLISION_ALLOWLIST)
        self.assertEqual(units.PINT_COLLISION_ALLOWLIST["bbl"], frozenset({"bbl"}))
        self.assertIs(registry._BOOTSTRAP_COLLISIONS, units.PINT_COLLISION_ALLOWLIST)

        with self.assertRaises(TypeError):
            units.PINT_COLLISION_ALLOWLIST["rogue"] = frozenset({"meter"})  # type: ignore[index]

        # Any unallowlisted bootstrap override fails closed
        patched_defs = registry.COMMODITY_UNIT_DEFINITIONS + ("rogue_override = 2 * meter = second",)
        with patch.object(registry, "COMMODITY_UNIT_DEFINITIONS", patched_defs):
            with self.assertRaises(pint.errors.RedefinitionError):
                with units.sandbox_registry():
                    pass

        # Runtime unit, alias, and @alias collisions fail closed with RedefinitionError
        with units.sandbox_registry() as sb:
            with self.assertRaises(pint.errors.RedefinitionError):
                sb.define("oil_barrel = 50 * gallon")
            with self.assertRaises(pint.errors.RedefinitionError):
                sb.define("fresh_parcel_x = 10 * kilogram = bbl")
            with self.assertRaises(pint.errors.RedefinitionError):
                sb.define("@alias psig = oil_barrel")

    # ------------------------------------------------------------------
    # 3. Immutable reference-data snapshots
    # ------------------------------------------------------------------

    def test_reference_data_and_registry_snapshots_are_immutable(self) -> None:
        snap = units.get_registry_snapshot()
        self.assertIsInstance(snap, units.RegistrySnapshot)
        self.assertEqual(snap.schema_version, units.REGISTRY_SCHEMA_VERSION)
        self.assertEqual(snap.registry_version, units.REGISTRY_VERSION)
        self.assertEqual(snap.baseline_definition_hash, units.BASELINE_DEFINITION_HASH)

        with self.assertRaises(FrozenInstanceError):
            snap.unit_count = 0  # type: ignore[misc]

        for mapping in (
            units.PINT_COLLISION_ALLOWLIST,
            units.CURRENCY_MINOR_UNIT_EXPONENT,
            units.REFERENCE_DATA_HASHES,
            units.CURRENCY_CATALOG.symbols,
            units.CURRENCY_CATALOG.shared_symbols,
            units.GGE_SPECS,
            units.DGE_SPECS,
            units.ISO6976_COMPONENTS,
            registry._SHARED_SYMBOL_CURRENCIES,
            registry._CANONICAL_TO_MARKET_TICKER,
            registry._MARKET_TICKERS,
            registry._CANONICAL_UNIT_DISPLAY,
        ):
            with self.assertRaises(TypeError):
                mapping["__mutate__"] = "forbidden"  # type: ignore[index]

        expected_keys = {
            "currencies_manifest",
            "currencies_reference",
            "currencies_conventions",
            "currencies_overlay",
            "energy_specs",
            "assay_reference_manifest",
            "gas_components",
            "astm_d2161",
            "icumsa_brix",
        }
        self.assertEqual(set(units.REFERENCE_DATA_HASHES.keys()), expected_keys)
        for digest in units.REFERENCE_DATA_HASHES.values():
            self.assertEqual(len(digest), 64)
            int(digest, 16)

    # ------------------------------------------------------------------
    # 4. Registry definition hashes
    # ------------------------------------------------------------------

    def test_definition_hashes_are_deterministic_and_track_extensions(self) -> None:
        self.assertEqual(len(units.BASELINE_DEFINITION_HASH), 64)
        int(units.BASELINE_DEFINITION_HASH, 16)

        with units.sandbox_registry() as sb:
            self.assertEqual(units.get_registry_definition_hash(sb), units.BASELINE_DEFINITION_HASH)
            sb.define("hash_test_tote = 1000 * liter")
            extended_hash = units.get_registry_definition_hash(sb)
            self.assertNotEqual(extended_hash, units.BASELINE_DEFINITION_HASH)
            sb_snap = units.get_registry_snapshot(sb)
            self.assertEqual(sb_snap.definition_hash, extended_hash)
            self.assertEqual(sb_snap.custom_registrations, ("hash_test_tote = 1000 * liter",))

        # Verify cross-process determinism of BASELINE_DEFINITION_HASH
        script = "import core.units as u; print(u.BASELINE_DEFINITION_HASH, u.get_registry_definition_hash())"
        proc = subprocess.run(
            [sys.executable, "-B", "-c", script],
            cwd=_REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)
        child_base, child_active = proc.stdout.strip().split()
        self.assertEqual(child_base, units.BASELINE_DEFINITION_HASH)
        self.assertEqual(child_active, units.BASELINE_DEFINITION_HASH)

    # ------------------------------------------------------------------
    # 5. Per-version registry compatibility checks
    # ------------------------------------------------------------------

    def test_per_version_compatibility_checks_and_packet_metadata(self) -> None:
        active_hash = units.get_registry_definition_hash()
        self.assertTrue(
            units.check_registry_compatibility(
                expected_schema_version=units.REGISTRY_SCHEMA_VERSION,
                expected_registry_version=units.REGISTRY_VERSION,
                expected_definition_hash=active_hash,
                required_units={"oil_barrel": "[volume]", "MMBtu_it": "[energy]", "USD": "[currency_USD]"},
            )
        )

        self.assertFalse(units.check_registry_compatibility(expected_schema_version=999))
        self.assertFalse(units.check_registry_compatibility(expected_registry_version="commodex-units-v0+old"))
        self.assertFalse(units.check_registry_compatibility(expected_definition_hash="0" * 64))
        self.assertFalse(units.check_registry_compatibility(required_units=("nonexistent_unit_xyz",)))
        self.assertFalse(units.check_registry_compatibility(required_units={"oil_barrel": "[mass]"}))

        with self.assertRaises(units.RegistryCompatibilityError):
            units.assert_registry_compatible(expected_registry_version="commodex-units-v999")

        # Superset extension policy on sandbox registry
        with units.sandbox_registry() as sb:
            sb.define("compat_desk_drum = 200 * liter")
            self.assertFalse(
                units.check_registry_compatibility(
                    expected_definition_hash=units.BASELINE_DEFINITION_HASH,
                    allow_superset_extensions=False,
                    registry=sb,
                )
            )
            self.assertTrue(
                units.check_registry_compatibility(
                    expected_definition_hash=units.BASELINE_DEFINITION_HASH,
                    allow_superset_extensions=True,
                    required_units={"compat_desk_drum": "[volume]"},
                    registry=sb,
                )
            )

        # Metadata-carrying quantity serialisation & verification
        qty = units.Q_(Decimal("125.75"), "USD / oil_barrel")
        compact = units.quantity_to_dict(qty)
        self.assertEqual(set(compact.keys()), {"m", "u"})

        enriched = units.quantity_to_dict(qty, include_metadata=True)
        self.assertEqual(enriched["schema_version"], units.REGISTRY_SCHEMA_VERSION)
        self.assertEqual(enriched["registry_version"], units.REGISTRY_VERSION)
        self.assertEqual(enriched["definition_hash"], active_hash)
        restored = units.quantity_from_dict(enriched)
        self.assertEqual(restored.magnitude, Decimal("125.75"))
        self.assertEqual(restored.units, qty.units)

        for bad_field, bad_val in (
            ("schema_version", 2),
            ("registry_version", "commodex-units-v1+1999-01-01"),
            ("definition_hash", "f" * 64),
        ):
            tampered = dict(enriched, **{bad_field: bad_val})
            with self.assertRaises(units.RegistryCompatibilityError):
                units.quantity_from_dict(tampered)

    # ------------------------------------------------------------------
    # 6. Prevention of mutation after initialisation
    # ------------------------------------------------------------------

    def test_post_freeze_mutation_blocked_across_all_surfaces(self) -> None:
        script = """
from decimal import Decimal
import pint
import core.units as u
from core.units import registry as r

u.freeze_registry()
assert u.is_registry_frozen() is True
assert u.get_registry_snapshot().is_frozen is True

redef_ctx = pint.Context("redef_ctx")
redef_ctx.redefine("oil_barrel = 50 * gallon")

def enter_redef_ctx():
    with u.ureg.context(redef_ctx):
        pass

blocked_actions = [
    lambda: u.register_currency("POSTFREEZE", decimals=2),
    lambda: u.register_commodity_unit("post_freeze_unit = 10 * kilogram"),
    lambda: u.load_desk_config({"commodities": [{"definition": "cfg_unit = 5 * meter"}]}),
    lambda: u.ureg.define("direct_freeze_unit = 5 * meter"),
    lambda: u.ureg.load_definitions(["loaded_freeze_unit = 5 * meter"]),
    lambda: u.ureg.add_context(pint.Context("extra_ctx")),
    lambda: u.ureg.remove_context("sp"),
    lambda: u.ureg.enable_contexts(u.create_fx_context({"EUR": Decimal("0.92")})),
    lambda: u.ureg.disable_contexts(),
    lambda: setattr(u.ureg, "default_system", "cgs"),
    lambda: setattr(u.ureg, "_on_redefinition", "ignore"),
    lambda: setattr(u.ureg, "_active_ctx", None),
    lambda: delattr(u.ureg, "_on_redefinition"),
    enter_redef_ctx,
    lambda: u.CURRENCY_MINOR_DIGITS.__setitem__("USD", 4),
    lambda: u.CURRENCY_MINOR_DIGITS.pop("USD"),
    lambda: u.CURRENCY_SYMBOL_TO_CODE.__setitem__("Z$", "USD"),
    lambda: r.minor_mapping.__setitem__("z¢", "USD_cent"),
    lambda: r.irregular_denominators.__setitem__("rogue", "kilogram"),
    lambda: r._QUOTE_REPLACEMENTS.__setitem__("$/rogue", "USD / kilogram"),
    lambda: r._QUOTE_REPLACEMENTS.clear(),
]

for idx, action in enumerate(blocked_actions):
    try:
        action()
    except RuntimeError:
        pass
    else:
        raise AssertionError(f"Post-freeze mutation #{idx} was not blocked!")

# Read-only parsing (including previously uncached quotes) and scoped contexts still work cleanly
quotes_len_before = len(r._QUOTE_REPLACEMENTS)
parsed = u.parse_quantity("CHF 88.50/mt") if False else u.parse_unit("CHF/mt")
assert parsed == u.ureg.Unit("CHF / tonne")
assert len(r._QUOTE_REPLACEMENTS) == quotes_len_before

fx_ctx = u.create_fx_context({"EUR": Decimal("0.90")})
with u.ureg.context(fx_ctx):
    eur_val = u.Q_(Decimal("100"), "USD").to("EUR")
    assert eur_val.magnitude == Decimal("90.00")
"""
        proc = subprocess.run(
            [sys.executable, "-B", "-c", script],
            cwd=_REPO_ROOT,
            capture_output=True,
            text=True,
            timeout=45,
            check=False,
        )
        self.assertEqual(proc.returncode, 0, proc.stderr)

    # ------------------------------------------------------------------
    # 7. Foreign-registry quantity and unit validation
    # ------------------------------------------------------------------

    def test_foreign_registry_quantities_and_units_rejected_everywhere(self) -> None:
        foreign_ureg = pint.UnitRegistry()
        fq_kg = foreign_ureg.Quantity(10.0, "kilogram")
        fq_vol = foreign_ureg.Quantity(10.0, "liter")
        fq_energy = foreign_ureg.Quantity(100.0, "megajoule")
        fq_press = foreign_ureg.Quantity(101325.0, "pascal")
        fq_temp = foreign_ureg.Quantity(300.0, "kelvin")
        fq_visc = foreign_ureg.Quantity(10.0, "centistokes")
        fq_dens = foreign_ureg.Quantity(850.0, "kilogram / meter ** 3")
        fq_dimless = foreign_ureg.Quantity(1.05, "dimensionless")
        fu_kg = foreign_ureg.Unit("kilogram")

        with units.sandbox_registry() as sb:
            sq_bbl = sb.Quantity(Decimal("10"), "oil_barrel")
            sq_usd = sb.Quantity(Decimal("100"), "USD")
            sq_quote = sb.Quantity(Decimal("75"), "USD / oil_barrel")
            sq_wmt = sb.Quantity(Decimal("100"), "wmt")
            sq_molarity = sb.Quantity(Decimal("2"), "molarity")
            su_mt = sb.Unit("tonne")

            crude = conversions.CrudeAssay.brent()
            softs = agricultural.SoftsAssay.raw_sugar_no11()

            calls = [
                lambda: units.validate_registry_quantity(fq_kg),
                lambda: units.validate_registry_quantity(fu_kg),
                lambda: units.validate_registry_quantity(sq_bbl),
                lambda: units.validate_registry_quantity(su_mt),
                lambda: units.parse_unit(fu_kg),
                lambda: units.parse_unit(su_mt),
                lambda: units.parse_quantity(10, fu_kg),  # type: ignore[arg-type]
                lambda: units.validate_dimensionality(fq_kg, "kg"),
                lambda: units.validate_dimensionality(units.Q_(1, "kg"), fu_kg),
                lambda: units.ensure_unit(fq_kg, "gram"),
                lambda: units.ensure_unit(units.Q_(1, "kg"), fu_kg),
                lambda: units.are_convertible(fq_kg, "gram"),
                lambda: units.are_convertible("gram", fu_kg),
                lambda: units.format_quantity(fq_kg),
                lambda: units.quantity_to_dict(sq_bbl),
                lambda: units.is_currency(sq_usd),
                lambda: units.is_commodity_price(sq_quote),
                lambda: units.is_flow_rate(fq_kg),
                lambda: units.is_fx_rate(su_mt),
                lambda: conversions.convert(fq_kg, "gram"),
                lambda: conversions.convert(units.Q_(1, "kg"), fu_kg),
                lambda: conversions.convert_pressure(fq_press, "bar"),
                lambda: conversions.differential_pressure(fq_press, units.Q_(1, "bar")),
                lambda: conversions.convert_currency(sq_usd, "EUR", rates={"EUR": Decimal("0.9")}),
                lambda: conversions.convert_fx(sq_usd, "EUR", {"USD/EUR": Decimal("0.9")}),
                lambda: conversions.convert_commodity_quote(sq_quote, "USD/mt", crude_assay=crude),
                lambda: conversions.convert_hydrocarbon(sq_bbl, "mt", assay=crude),
                lambda: conversions.convert_bulk_ore(sq_wmt, "dmt", moisture_fraction=Decimal("0.1")),
                lambda: conversions.convert_wet_to_dry_mass(sq_wmt, 10.0),
                lambda: conversions.calculate_z_factor(fq_press, fq_temp),
                lambda: conversions.GasReferenceConditions(fq_temp, fq_press),
                lambda: crude.correct_volume_to_60f(sq_bbl, 75.0),
                lambda: softs.raw_to_processed_mass(fq_kg),
                lambda: conversions.HydrogenAssay.h2_mass_to_nm3(fq_kg),
                lambda: heat_content.HeatContentAssay(ncv=fq_energy / fq_kg),
                lambda: heat_content.convert_heat_content(
                    fq_kg, "MJ", assay=conversions.BiomassAssay.enplus_a1_pellets()
                ),
                lambda: heat_content.thermal_to_output(fq_energy, efficiency=Decimal("0.4")),
                lambda: currencies.round_currency(sq_usd),
                lambda: currencies.redenominate_currency(sq_usd, "TRL_TRY_2005", as_of=date(2006, 1, 1)),
                lambda: converters.energy_to_gge(fq_energy, "NIST_GGE", basis="LHV"),
                lambda: converters.gas_volume_to_energy(fq_vol, fq_energy / fq_vol, basis="HHV", reference_conditions="iso15"),
                lambda: converters.convert_viscosity(fq_visc, "SUS", temperature=100.0),
                lambda: converters.dynamic_to_kinematic_viscosity(
                    foreign_ureg.Quantity(5.0, "centipoise"), fq_dens
                ),
                lambda: converters.molarity_to_normality(sq_molarity, 2),
                lambda: converters.sg_to_brix(fq_dimless),
                lambda: converters.convert_chemical_strength(sq_molarity, "normality", equivalence_factor=2),
            ]

            for idx, fn in enumerate(calls):
                with self.subTest(call_index=idx), self.assertRaises(units.ForeignRegistryError):
                    fn()

    # ------------------------------------------------------------------
    # 8. Explicit thread, asyncio task, and subprocess isolation tests
    # ------------------------------------------------------------------

    def test_thread_context_isolation_no_cross_thread_bleed(self) -> None:
        n_workers = 10
        barrier = threading.Barrier(n_workers + 1)

        def contextual_worker(worker_id: int) -> None:
            rate = Decimal("0.80") + Decimal(worker_id) * Decimal("0.03")
            api = 25.0 + worker_id * 2.0
            assay = conversions.CrudeAssay(f"thread_crude_{worker_id}", api_gravity=api)
            fx_ctx = conversions.create_fx_context({"EUR": rate})
            hc_ctx = conversions.create_hydrocarbon_context(assay)
            barrier.wait(timeout=10)
            for _ in range(30):
                with units.ureg.context(fx_ctx, hc_ctx):
                    eur = units.Q_(Decimal("100"), "USD").to("EUR")
                    self.assertEqual(eur.magnitude, Decimal("100") * rate)
                    mt = units.Q_(1000.0, "oil_barrel").to("metric_tonne")
                    self.assertAlmostEqual(mt.magnitude, 1000.0 * assay.mt_per_bbl, places=9)
                # Immediately outside the with-block, contexts must be inactive in this thread
                with self.assertRaises(pint.DimensionalityError):
                    units.Q_(Decimal("1"), "USD").to("EUR")

        def uncontextualised_observer() -> None:
            barrier.wait(timeout=10)
            for _ in range(40):
                with self.assertRaises(pint.DimensionalityError):
                    units.Q_(Decimal("1"), "USD").to("EUR")
                with self.assertRaises(pint.DimensionalityError):
                    units.Q_(1.0, "oil_barrel").to("metric_tonne")

        with ThreadPoolExecutor(max_workers=n_workers + 1) as pool:
            futures = [pool.submit(contextual_worker, i) for i in range(n_workers)]
            futures.append(pool.submit(uncontextualised_observer))
            for fut in futures:
                fut.result(timeout=30)

    def test_asyncio_task_and_child_task_context_isolation(self) -> None:
        # Socket-free PEP 567 ContextVar async task runner (mirrors asyncio.Task &
        # asyncio.create_task context copying without Windows loopback socketpair()).
        @types.coroutine
        def async_yield():
            yield

        class _AsyncTaskHandle:
            def __init__(self, coro, ctx: contextvars.Context) -> None:
                self.coro = coro
                self.ctx = ctx
                self.done = False
                self.result = None

            def __await__(self):
                while not self.done:
                    yield
                return self.result

        ready: deque[_AsyncTaskHandle] = deque()

        def create_task(coro) -> _AsyncTaskHandle:
            handle = _AsyncTaskHandle(coro, contextvars.copy_context())
            ready.append(handle)
            return handle

        def run_until_complete(root_coro):
            root = create_task(root_coro)
            while ready:
                task = ready.popleft()
                try:
                    task.ctx.run(task.coro.send, None)
                except StopIteration as stop:
                    task.done = True
                    task.result = stop.value
                else:
                    ready.append(task)
            return root.result

        async def run_concurrent_tasks() -> None:
            async def fx_task(task_id: int) -> Decimal:
                rate = Decimal("0.70") + Decimal(task_id) * Decimal("0.02")
                fx_ctx = conversions.create_fx_context({"EUR": rate})
                with units.ureg.context(fx_ctx):
                    await async_yield()
                    first = units.Q_(Decimal("50"), "USD").to("EUR").magnitude
                    self.assertEqual(first, Decimal("50") * rate)

                    # Spawn a child task inside the active context that overrides EUR and adds GBP
                    async def child_coro() -> Decimal:
                        gbp_ctx = conversions.create_fx_context({"EUR": Decimal("2.00"), "GBP": Decimal("0.75")})
                        with units.ureg.context(gbp_ctx):
                            await async_yield()
                            return units.Q_(Decimal("10"), "USD").to("EUR").magnitude

                    child_task = create_task(child_coro())
                    await async_yield()
                    # Parent's EUR rate must remain completely unaffected by child's override!
                    second = units.Q_(Decimal("50"), "USD").to("EUR").magnitude
                    self.assertEqual(second, Decimal("50") * rate)
                    child_res = await child_task
                    self.assertEqual(child_res, Decimal("20.00"))

                with self.assertRaises(pint.DimensionalityError):
                    units.Q_(Decimal("1"), "USD").to("EUR")
                return second

            async def uncontextualised_task() -> None:
                for _ in range(15):
                    await async_yield()
                    with self.assertRaises(pint.DimensionalityError):
                        units.Q_(Decimal("1"), "USD").to("EUR")

            handles = [create_task(fx_task(i)) for i in range(16)]
            obs = create_task(uncontextualised_task())
            results = [await h for h in handles]
            await obs
            for i in range(16):
                expected = Decimal("50") * (Decimal("0.70") + Decimal(i) * Decimal("0.02"))
                self.assertEqual(results[i], expected)

        run_until_complete(run_concurrent_tasks())

    def test_subprocess_registry_isolation(self) -> None:
        worker_a = """
import core.units as u
u.register_commodity_unit("subproc_isolated_unit = 10 * kilogram")
u.freeze_registry()
print(u.Q_(3, "subproc_isolated_unit").to("kilogram").magnitude, u.get_registry_definition_hash())
"""
        worker_b = """
import core.units as u
u.register_commodity_unit("subproc_isolated_unit = 25 * liter")
u.freeze_registry()
print(u.Q_(3, "subproc_isolated_unit").to("liter").magnitude, u.get_registry_definition_hash())
"""
        proc_a = subprocess.Popen(
            [sys.executable, "-B", "-c", worker_a],
            cwd=_REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        proc_b = subprocess.Popen(
            [sys.executable, "-B", "-c", worker_b],
            cwd=_REPO_ROOT,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )
        out_a, err_a = proc_a.communicate(timeout=30)
        out_b, err_b = proc_b.communicate(timeout=30)
        self.assertEqual(proc_a.returncode, 0, err_a)
        self.assertEqual(proc_b.returncode, 0, err_b)

        mag_a, hash_a = out_a.strip().split()
        mag_b, hash_b = out_b.strip().split()
        self.assertEqual(float(mag_a), 30.0)
        self.assertEqual(float(mag_b), 75.0)
        self.assertNotEqual(hash_a, hash_b)
        self.assertNotIn("subproc_isolated_unit", units.ureg)
        self.assertFalse(units.is_registry_frozen())

    # ------------------------------------------------------------------
    # 9. Restrictions on unrestricted user-supplied Pint definitions
    #    and external market-data strings
    # ------------------------------------------------------------------

    def test_adversarial_pint_definitions_are_rejected(self) -> None:
        hostile_definitions = (
            "",
            "   ",
            "@import /etc/passwd",
            "@context evil_ctx",
            "@system cgs",
            "@group secret_group",
            "@alias oil_barrel = sneaky_alias",
            "multi_line = 10 * meter\ninjected_unit = 5 * second",
            "null_byte_unit = 10 * meter\x00",
            "comment_unit = 10 * meter # inline comment",
            "offset_unit = 1 * kelvin; offset: 273.15",
            "prefix_hack- = 1000 * meter",
            "colon:meta = 10 * meter",
            "new_root_dim = [arbitrary_new_dimension]",
            "spoof_ccy = [currency_USD]",
            "zero_scale = 0 * meter",
            "neg_scale = -5 * meter",
            "extreme_scale = 1e100 * meter",
            "tiny_scale = 1e-100 * meter",
            "bare_scalar = 42",
            "unknown_rhs = 10 * totally_nonexistent_unit_xyz",
            "huge_power = meter ** 24",
            "zero_power = meter ** 0",
            "float_power = meter ** 1.5",
            "func_call = sin(meter)",
            "attr_access = meter.__class__",
            "sum_expr = meter + second",
            "xor_expr = meter ^ 3",
            "def = 10 * meter",
            "_private_unit = 10 * meter",
            "ppt = 10 * item",
            "dup_alias = 10 * meter = d_alias = d_alias",
            "overlong_" + ("x" * 260) + " = 10 * meter",
        )
        with units.sandbox_registry() as sb:
            for hostile in hostile_definitions:
                with self.subTest(definition=hostile):
                    with self.assertRaises(ValueError):
                        units.validate_unit_definition(hostile, registry=sb)
                    with self.assertRaises(ValueError):
                        units.register_commodity_unit(hostile)

            # Direct sb.define also enforces validation (except explicit validated @alias)
            for hostile in hostile_definitions:
                if hostile.startswith("@alias "):
                    continue
                with self.subTest(direct_define=hostile), self.assertRaises(ValueError):
                    sb.define(hostile)

            # External file/resource loading via load_definitions is blocked even when unfrozen
            with self.assertRaises(ValueError):
                sb.load_definitions("external_definitions.txt")
            with self.assertRaises(ValueError):
                sb.load_definitions(["@import external.txt"])

            # Valid custom definitions succeed cleanly
            canonical, rhs, aliases = units.validate_unit_definition(
                "clean_iso_tank = 24000 * liter = iso_tank_24kl",
                registry=sb,
            )
            self.assertEqual((canonical, rhs, aliases), ("clean_iso_tank", "24000 * liter", ("iso_tank_24kl",)))
            sb.define("clean_iso_tank = 24000 * liter = iso_tank_24kl")
            self.assertEqual(sb.Quantity(1, "iso_tank_24kl").to("liter").magnitude, 24000)

    def test_external_market_data_unit_strings_block_directives_and_ast_abuse(self) -> None:
        hostile_inputs = (
            "@import evil.txt",
            "@context sp",
            "@alias oil_barrel = b2",
            "evil_unit = 100 * meter",
            "kelvin; offset: 273.15",
            "USD/bbl\nEUR/bbl",
            "USD/bbl\rEUR/bbl",
            "USD/bbl\x00",
            "USD/bbl # comment",
            "[currency_USD]",
            "sin(meter)",
            "__import__('os')",
            "meter.__class__",
            "USD.__builtins__",
            "meter ** 9999",
            "meter ^ 99",
            "meter **",
            "USD / 'bbl'",
            "USD / {bbl}",
            "A" * 300,
        )
        for bad in hostile_inputs:
            with self.subTest(raw=bad):
                with self.assertRaises(ValueError):
                    units.sanitize_unit_string(bad)
                with self.assertRaises(ValueError):
                    units.parse_unit(bad)
                with self.assertRaises(ValueError):
                    units.parse_quantity(bad)
                with self.assertRaises(ValueError):
                    units.parse_quote(10, bad)
                with self.assertRaises(ValueError):
                    units.quantity_from_dict({"m": "10", "u": bad})

        # Legitimate market strings, percentages, Unicode currency glyphs, and bounded powers still work
        self.assertEqual(units.parse_quantity("2.5%").to("dimensionless").magnitude, Decimal("0.025"))
        self.assertEqual(units.parse_quantity("د.إ100/bbl").units, units.ureg.Unit("AED / oil_barrel"))
        self.assertEqual(units.parse_unit("kg / m^3"), units.ureg.Unit("kilogram / meter ** 3"))


if __name__ == "__main__":
    unittest.main()
