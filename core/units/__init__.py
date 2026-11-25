"""Public unit registry, conversion, assay and desk-configuration APIs.

CleanFuelAssay is the lightweight multi-carrier thermodynamic benchmark.
HydrogenAssay is the detailed H2/storage/carrier specification with recovery,
electrolyser demand and 45V/RFNBO screening. Choose the latter when those
fields matter; matching carrier names do not imply matching assay semantics.
Crude/products require explicit calorific data for heat-content conversions.

Pressure datum shifts use convert()/convert_pressure(), not raw Pint .to().
"""

from .currencies import *
from .currencies import __all__ as _currency_exports
from .registry import *
from .registry import __all__ as _registry_exports
from .conversions import *
from .conversions import __all__ as _conversion_exports
from .agricultural_specs import *
from .agricultural_specs import __all__ as _agricultural_exports
from .heat_content import *
from .heat_content import __all__ as _heat_exports
from .config import load_desk_config
from .config import __all__ as _config_exports

__all__ = tuple(dict.fromkeys(
    _currency_exports + _registry_exports + _conversion_exports + _agricultural_exports + _heat_exports + _config_exports
))
