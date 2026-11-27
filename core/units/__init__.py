"""Commodex unit, conversion and assay public APIs."""
from importlib import import_module as _import_module
import sys as _sys
import types as _types
import json as _json
from pathlib import Path as _Path

class _CompatibilityModule(_types.ModuleType):
    def __getattribute__(self, name):
        bindings = super().__getattribute__("__dict__").get("_bindings", {})
        if name in bindings:
            module, symbol = bindings[name]
            return getattr(_import_module("core.units." + module), symbol)
        return super().__getattribute__(name)
    def __getattr__(self, name):
        try:
            module, symbol = self._bindings[name]
        except KeyError:
            raise AttributeError(name) from None
        return getattr(_import_module("core.units." + module), symbol)
    def __setattr__(self, name, value):
        bindings = self.__dict__.get("_bindings", {})
        if name in bindings:
            module, symbol = bindings[name]
            setattr(_import_module("core.units." + module), symbol, value)
        else:
            super().__setattr__(name, value)
    def __delattr__(self, name):
        if name in self._bindings:
            module, symbol = self._bindings[name]
            delattr(_import_module("core.units." + module), symbol)
        else:
            super().__delattr__(name)

_LAYOUT = _json.loads((_Path(__file__).with_name("reference_data") / "aliases/python_imports.json").read_text(encoding="utf-8"))
_COMPATIBILITY = _LAYOUT["compatibility"]
_ORIGINAL_EXPORTS = _LAYOUT["original_exports"]
for _old, _bindings in _COMPATIBILITY.items():
    _module = _CompatibilityModule(__name__ + "." + _old)
    _module._bindings = _bindings
    _module.__all__ = tuple(_ORIGINAL_EXPORTS.get(_old, ()))
    if _old == "converters":
        _module.__path__ = []
    _sys.modules[_module.__name__] = _module
    if "." not in _old:
        globals()[_old] = _module
_MODULES = _LAYOUT["modules"]
for _name in _MODULES:
    _import_module(__name__ + "." + _name)
_EXPORTS = _LAYOUT["exports"]
for _name, (_owner, _symbol) in _EXPORTS.items():
    globals()[_name] = getattr(_import_module(__name__ + "." + _owner), _symbol)
__all__ = tuple(_EXPORTS)
for _old, _bindings in _COMPATIBILITY.items():
    _module = _sys.modules[__name__ + "." + _old]
    for _name, (_owner, _symbol) in _bindings.items():
        _module.__dict__[_name] = getattr(_import_module(__name__ + "." + _owner), _symbol)
_AGGREGATES = _LAYOUT["aggregates"]
for _old, _bindings in _AGGREGATES.items():
    _module = _import_module(__name__ + "." + _old)
    for _name, (_owner, _symbol) in _bindings.items():
        if _name not in vars(_module):
            setattr(_module, _name, getattr(_import_module(__name__ + "." + _owner), _symbol))
