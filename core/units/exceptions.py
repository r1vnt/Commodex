"""Commodex exceptions implementation and public interfaces."""
from __future__ import annotations
from typing import Final

class CommodityUnitError(Exception):
    """Base exception for all domain-specific unit, assay, quote, FX, and registry errors."""

class AmbiguousUnitError(CommodityUnitError, ValueError):
    """Raised when an ambiguous unit, market ticker, quote token, or currency glyph lacks context."""

class MissingAssayError(CommodityUnitError, ValueError):
    """Raised when a conversion requires an assay, density, moisture, grade, or calorific specification."""

class MissingReferenceConditionsError(CommodityUnitError, ValueError):
    """Raised when a gas, liquid, calorific, or dated convention lacks explicit reference conditions."""

class IncompatibleCommodityError(CommodityUnitError, ValueError):
    """Raised when blending, transferring, or converting across incompatible commodities or grades."""

class IncompatibleAnalyticalBasisError(CommodityUnitError, ValueError):
    """Raised when analytical bases (e.g. AR/AD/DB/DAF, LHV/HHV, wet/dry, gross/pure) conflict."""

class UnsupportedConversionError(CommodityUnitError, ValueError):
    """Raised when a requested conversion mode, route, or option combination is unsupported."""

class InvalidPhysicalStateError(CommodityUnitError, ValueError):
    """Raised when physical quantities or thermodynamic states violate physical bounds."""

class ConversionConvergenceError(CommodityUnitError, ArithmeticError, ValueError):
    """Raised when an iterative numerical EOS, root-finding, or property inversion fails to converge."""

class MissingFXRateError(CommodityUnitError, ValueError):
    """Raised when cross-currency conversion lacks required FX fixings or a triangulation path."""

class RegistryVersionMismatchError(CommodityUnitError, ValueError):
    """Raised when registry schema version, registry version, definition hash, or unit contracts mismatch."""

class RegistryFrozenError(CommodityUnitError, RuntimeError):
    """Raised when attempting to mutate or register definitions on a frozen UnitRegistry."""

class StandardsApplicabilityError(CommodityUnitError, ValueError):
    """Raised when an empirical correlation or measurement standard is applied outside its valid scope."""

__all__ = ('AmbiguousUnitError', 'CommodityUnitError', 'ConversionConvergenceError', 'IncompatibleAnalyticalBasisError', 'IncompatibleCommodityError', 'InvalidPhysicalStateError', 'MissingAssayError', 'MissingFXRateError', 'MissingReferenceConditionsError', 'RegistryFrozenError', 'RegistryVersionMismatchError', 'StandardsApplicabilityError', 'UnsupportedConversionError')
_PUBLIC_ALIASES = {}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)
