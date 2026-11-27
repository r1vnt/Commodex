"""Commodex config implementation and public interfaces."""
from __future__ import annotations

from pathlib import Path
from typing import Any
import json

def load_desk_config(source: str | Path | dict[str, Any], *, freeze_after_load: bool = True,) -> None:
    """Register desk currencies and commodities from a mapping or JSON/YAML file.

    YAML support is optional and requires PyYAML. When ``freeze_after_load`` is
    true, the registry is sealed after all configured entries are registered.
    """
    from core.units.registry import _registry_transaction
    from core.units.registry import freeze_registry
    from core.units.registry import register_commodity_unit
    from core.units.registry import register_currency
    if isinstance(source, (str, Path)):
        file_path = Path(source)
        if not file_path.exists():
            raise FileNotFoundError(f"Desk configuration file not found: {file_path}")

        raw_text = file_path.read_text(encoding="utf-8")
        suffix = file_path.suffix.lower()
        if suffix == ".json":
            payload = json.loads(raw_text)
        elif suffix in {".yaml", ".yml"}:
            try:
                import yaml
            except ImportError as error:
                raise ImportError("Loading YAML desk configuration requires the optional PyYAML package.") from error
            payload = yaml.safe_load(raw_text)
        else:
            raise ValueError(f"Unsupported config format: {file_path.suffix} (use YAML or JSON)")
    elif isinstance(source, dict):
        payload = source
    else:
        raise TypeError(f"Expected file path or dict, got {type(source).__name__}")

    if not isinstance(payload, dict):
        raise ValueError("Desk configuration must contain a mapping at its root.")

    currencies = payload.get("currencies", [])
    commodities = payload.get("commodities", [])
    if not isinstance(currencies, list):
        raise ValueError("'currencies' must be a list.")
    if not isinstance(commodities, list):
        raise ValueError("'commodities' must be a list.")
    for index, currency in enumerate(currencies):
        if not isinstance(currency, dict) or not isinstance(currency.get("code"), str):
            raise ValueError(f"Currency at index {index} must be a mapping with a string 'code'.")
    for index, commodity in enumerate(commodities):
        if not isinstance(commodity, dict) or not isinstance(commodity.get("definition"), str):
            raise ValueError(f"Commodity at index {index} must be a mapping with a string 'definition'.")

    with _registry_transaction():
        for currency in currencies:
            register_currency( code=currency["code"], decimals=currency.get("decimals", 2), symbols=tuple(currency.get("symbols", ())), aliases=tuple(currency.get("aliases", ())), sub_units=currency.get("sub_units"), quote_aliases=currency.get("quote_aliases"),)

        for commodity in commodities:
            register_commodity_unit(definition=commodity["definition"], quote_aliases=commodity.get("quote_aliases"),)

        if freeze_after_load:
            freeze_registry()


__all__ = ('load_desk_config',)
_PUBLIC_ALIASES = {'register_currency': ('registry', 'register_currency'), 'freeze_registry': ('registry', 'freeze_registry'), '_registry_transaction': ('registry', '_registry_transaction'), 'register_commodity_unit': ('registry', 'register_commodity_unit')}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
