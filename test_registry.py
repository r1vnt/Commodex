"""Regression coverage for registry bootstrap and market quote parsing."""

from concurrent.futures import ThreadPoolExecutor
from decimal import Decimal
from importlib.util import find_spec
import json
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
import unittest
from unittest.mock import patch

import pint

from core.units import registry as units
from core.units.config import load_desk_config


class RegistryTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        units.verify_registry_integrity()

    def test_annual_capacity(self):
        annual_volume = units.Q_(1, "bcma") * units.Q_(365, "day")
        self.assertAlmostEqual(annual_volume.to("BSm3_15C_101325kPa_dry").magnitude, 1)

    def test_bootstrap_restores_strict_redefinitions(self):
        with self.assertRaises(pint.errors.RedefinitionError):
            units.ureg.define("trading_year_365 = 366 * day")
        with self.assertRaises(pint.errors.RedefinitionError):
            units.register_commodity_unit("regression_collision = meter = oil_barrel")

    def test_runtime_currency_and_aliases(self):
        units.register_currency("BTC", decimals=8, aliases=("XBT", "bitcoin"))
        for code in ("BTC", "btc", "XBT", "BiTcOiN"):
            self.assertEqual(units.sanitize_unit_string(f"{code}/MWh"), "BTC / megawatt_hour")
            self.assertEqual(units.parse_quantity(1, f"{code}/MWh").units, units.ureg.Unit("BTC / megawatt_hour"))

    def test_liberian_symbol(self):
        self.assertEqual(units.parse_quantity("L$100").units, units.ureg.LRD)

    def test_default_currency_overrides_cached_quote(self):
        self.assertEqual(units.sanitize_unit_string("$/bbl"), "USD / oil_barrel")
        self.assertEqual(units.sanitize_unit_string("$/bbl", "CAD"), "CAD / oil_barrel")
        self.assertEqual(units.sanitize_unit_string("$/bbl"), "USD / oil_barrel")
        with self.assertRaises(ValueError):
            units.sanitize_unit_string("$/bbl", "EUR")

    def test_leading_decimal_quotes(self):
        for raw, magnitude in (("$.50/bbl", ".50"), ("-$.85/gal", "-.85"), ("$.5e2/bbl", "50")):
            self.assertEqual(units.parse_quantity(raw).magnitude, Decimal(magnitude))

    def test_pressure_carbon_intensity_bullion_and_grid_units(self):
        self.assertEqual(units.parse_quantity("14.73 psia").to("psi").magnitude, Decimal("14.73"))
        self.assertEqual(units.parse_quantity("600 psig").to("psi").magnitude, 600)
        self.assertEqual(units.parse_quantity("1 bara").to("bar").magnitude, 1)
        self.assertEqual(units.parse_quantity("1 barg").to("bar").magnitude, 1)

        carbon_intensity = units.Q_(35.5, "CI")
        self.assertEqual(carbon_intensity.dimensionality, (units.ureg.gram / units.ureg.megajoule).dimensionality)
        self.assertEqual(units.Q_(1, "kgCO2e_per_MWh").to("kgCO2e / megawatt_hour").magnitude, 1)

        self.assertAlmostEqual(units.Q_(1, "kilobar_gold").to("toz").magnitude, 32.1507465, places=6)
        self.assertEqual(units.Q_(1, "ten_tola_bar").to("tola").magnitude, 10)
        self.assertEqual(units.Q_(1, "good_delivery_gold").to("toz").magnitude, 400)
        self.assertEqual(units.Q_(1, "good_delivery_silver").to("toz").magnitude, 1000)

        self.assertEqual(units.Q_(1, "MW_per_minute").to("kW_per_minute").magnitude, 1000)
        self.assertEqual(units.Q_(1, "C_rate").to("1 / hour").magnitude, 1)

    def test_latex_quantity_format(self):
        quantity = units.Q_(75.5, "USD / oil_barrel")
        self.assertEqual(
            units.format_quantity(quantity, style="latex", precision=2),
            r"75.50\,\frac{\mathrm{USD}}{\mathrm{bbl}}",
        )

    def test_standalone_unit_parser(self):
        self.assertEqual(units.parse_unit("$/MMBtu"), units.ureg.Unit("USD / MMBtu_it"))
        self.assertEqual(units.parse_unit("c/bu"), units.ureg.Unit("USD_cent / bushel"))
        parsed = units.parse_unit(units.ureg.Unit("oil_barrel"))
        self.assertEqual(parsed, units.ureg.Unit("oil_barrel"))
        with self.assertRaises(ValueError):
            units.parse_unit("")

    def test_fx_rate_classification(self):
        for expression in ("EUR / USD", "GBP / USD", "USD / JPY"):
            self.assertTrue(units.is_fx_rate(units.parse_unit(expression)))
            self.assertTrue(units.is_fx_rate(units.Q_(1, expression)))

        self.assertFalse(units.is_fx_rate(units.Q_(1, "USD")))
        self.assertFalse(units.is_fx_rate(units.Q_(1, "USD / oil_barrel")))
        self.assertFalse(units.is_fx_rate("EUR / USD"))

    def test_dimensional_convertibility(self):
        self.assertTrue(units.are_convertible("bbl", "gallon"))
        self.assertTrue(units.are_convertible("cwt", "mt"))
        self.assertTrue(units.are_convertible(units.Q_(2, "bbl"), units.ureg.Unit("gallon")))
        self.assertFalse(units.are_convertible("bbl", "MWh"))
        self.assertFalse(units.are_convertible("USD", "EUR"))
        self.assertFalse(units.are_convertible("not_a_unit", "gallon"))

    def test_cache_miss_is_resolved_once(self):
        cache = units.QuoteReplacementsMap()
        with patch.object(units, "_resolve_quote_expression", return_value="USD / meter") as resolve:
            with ThreadPoolExecutor(max_workers=8) as pool:
                values = list(pool.map(lambda _: cache["test_quote"], range(32)))
            self.assertEqual(values, ["USD / meter"] * 32)
            resolve.assert_called_once_with("test_quote")

    def test_registry_lock_is_reentrant(self):
        lock = units._REGISTRY_LOCK
        self.assertTrue(lock.acquire())
        try:
            self.assertTrue(lock.acquire(blocking=False))
            lock.release()
        finally:
            lock.release()

    def test_load_desk_config_from_mapping_and_json(self):
        load_desk_config(
            {
                "currencies": [
                    {"code": "TSTCFG", "decimals": 4, "aliases": ["TSTC"]},
                ],
                "commodities": [
                    {"definition": "config_test_parcel = 3 * kilogram = ctp"},
                ],
            },
            freeze_after_load=False,
        )
        self.assertEqual(units.Q_(1, "config_test_parcel").to("kilogram").magnitude, 3)
        self.assertEqual(units.sanitize_unit_string("TSTC/MWh"), "TSTCFG / megawatt_hour")
        self.assertFalse(units.is_registry_frozen())

        config = {
            "commodities": [
                {"definition": "json_config_test_unit = 5 * meter"},
            ],
        }
        with TemporaryDirectory() as directory:
            config_path = Path(directory) / "desk.json"
            config_path.write_text(json.dumps(config), encoding="utf-8")
            load_desk_config(config_path, freeze_after_load=False)
        self.assertEqual(units.Q_(1, "json_config_test_unit").to("meter").magnitude, 5)

    @unittest.skipUnless(find_spec("yaml"), "PyYAML is optional")
    def test_load_desk_config_from_yaml(self):
        with TemporaryDirectory() as directory:
            config_path = Path(directory) / "desk.yaml"
            config_path.write_text(
                "commodities:\n  - definition: yaml_config_test_unit = 7 * meter\n",
                encoding="utf-8",
            )
            load_desk_config(config_path, freeze_after_load=False)
        self.assertEqual(units.Q_(1, "yaml_config_test_unit").to("meter").magnitude, 7)

    def test_load_desk_config_freezes_after_bootstrap(self):
        code = (
            "from core.units.config import load_desk_config; "
            "from core.units.registry import is_registry_frozen, register_commodity_unit; "
            "load_desk_config({'commodities': [{'definition': 'boot_config_unit = 2 * meter'}]}); "
            "assert is_registry_frozen(); "
            "\ntry: register_commodity_unit('late_config_unit = 3 * meter')"
            "\nexcept RuntimeError: pass"
            "\nelse: raise AssertionError('registry accepted registration after freeze')"
        )
        result = subprocess.run(
            [sys.executable, "-c", code],
            cwd=Path(__file__).parent,
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(result.returncode, 0, result.stderr)


if __name__ == "__main__":
    unittest.main()
