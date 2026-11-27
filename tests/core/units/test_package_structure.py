"""The requested physical layout must contain implementations, not leaf facades."""
import importlib
import pickle
from pathlib import Path

from tests.core.units._test_layout import check


def test_exact_library_layout_and_documentation():
    check(Path(__file__).resolve().parents[3])


def test_assays_are_defined_in_their_canonical_modules():
    import core.units as units
    targets = {
        "CrudeAssay": "assays.energy.petroleum", "RefineryCrudeAssay": "assays.energy.refining",
        "GasAssay": "assays.energy.natural_gas", "CoalAssay": "assays.energy.coal",
        "GrainAssay": "assays.agriculture.grains", "WheatQualitySpec": "assays.agriculture.grains",
        "PulseAssay": "assays.agriculture.pulses", "OreAssay": "assays.metals.ores",
        "RefinedMetalAssay": "assays.metals.refined_metals", "ChemicalAssay": "assays.industrial.chemicals",
        "IndustrialGasAssay": "assays.industrial.industrial_gases",
        "CarbonDioxideAssay": "assays.environmental.carbon_dioxide", "WaterAssay": "assays.environmental.water",
    }
    for name, module in targets.items():
        cls = getattr(units, name)
        assert cls.__module__ == "core.units." + module
        assert getattr(importlib.import_module(cls.__module__), name) is cls
    assert "class CrudeAssay" not in (Path(units.__file__).parent / "conversions.py").read_text(encoding="utf-8")
    for name in units.__all__:
        value=getattr(units,name)
        if isinstance(value,type) and name.endswith('Assay') and name != 'HeatContentAssay':
            assert value.__module__.startswith('core.units.assays.'),(name,value.__module__)


def test_legacy_imports_and_pickle_class_references_keep_working():
    import core.units as units
    from core.units.conversions import CrudeAssay
    from core.units.material_specs import RefinedMetalAssay
    assert CrudeAssay is units.CrudeAssay
    assert RefinedMetalAssay is units.RefinedMetalAssay
    assert pickle.loads(b"ccore.units.conversions\nCrudeAssay\n.") is CrudeAssay
    assert pickle.loads(pickle.dumps(CrudeAssay.brent())) == CrudeAssay.brent()
