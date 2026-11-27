"""Commodex blending implementation and public interfaces."""
from __future__ import annotations



__all__ = ('BlendVolumeModel', 'LiquidBlendComponent', 'blend_heat_content', 'blend_petroleum', 'ideal_blend_volume')
_PUBLIC_ALIASES = {'BlendVolumeModel': ('assays.energy.petroleum', 'BlendVolumeModel'), 'blend_heat_content': ('processing', 'blend_heat_content'), 'LiquidBlendComponent': ('assays.energy.petroleum', 'LiquidBlendComponent'), 'blend_petroleum': ('assays.energy.petroleum', 'blend_petroleum'), 'ideal_blend_volume': ('assays.energy.petroleum', 'ideal_blend_volume')}

def __getattr__(name):
    from importlib import import_module
    try:
        module, symbol = _PUBLIC_ALIASES[name]
    except KeyError:
        raise AttributeError(name) from None
    return getattr(import_module("core.units." + module), symbol)

# Static declarations for lazy public interfaces.
from typing import TYPE_CHECKING, Any, Callable
if TYPE_CHECKING:
    from core.units.assays.energy.petroleum import BlendVolumeModel
    from core.units.assays.energy.petroleum import LiquidBlendComponent
    from core.units.processing import blend_heat_content
    from core.units.assays.energy.petroleum import blend_petroleum
    from core.units.assays.energy.petroleum import ideal_blend_volume
