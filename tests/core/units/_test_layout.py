"""Verify the physical package layout against the requested Commodex tree."""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
ROOT_MODULES = "__init__ registry config exceptions types dimensions context validation precision conversions agricultural_specs heat_content quote_parser currencies reference_conditions analytical_basis composition blending processing thermodynamics logistics_units provenance standards uncertainty serialisation batch".split()
SECTORS = {
    "energy": "petroleum refining natural_gas lng ngl coal petcoke bitumen biomass biofuels synthetic_fuels hydrogen nuclear waste_fuels".split(),
    "metals": "ores concentrates refined_metals precious_metals ferroalloys alloys steel scrap battery_materials rare_earths critical_minerals".split(),
    "agriculture": "grains rice pulses oilseeds vegetable_oils softs livestock dairy poultry seafood forestry fibres rubber feedstocks fruit specialty".split(),
    "industrial": "petrochemicals chemicals fertilisers polymers industrial_gases construction_materials industrial_minerals".split(),
    "environmental": "carbon_dioxide water certificate_units".split(),
}
REFERENCE_CATEGORIES = {"units", "currencies", "aliases", "standards", "assays"}


def expected_files():
    files = {name + ".py" for name in ROOT_MODULES} | {"assays/__init__.py", "assays/base.py"}
    for sector, names in SECTORS.items():
        files |= {f"assays/{sector}/{name}.py" for name in ["__init__", *names]}
    return files


def check(root=ROOT):
    package = root / "core/units"
    actual = {path.relative_to(package).as_posix() for path in package.rglob("*.py")}
    assert actual == expected_files(), {"missing": expected_files() - actual, "unexpected": actual - expected_files()}
    assert {path.name for path in (package / "reference_data").iterdir() if path.is_dir() and path.name != '__pycache__'} == REFERENCE_CATEGORIES
    assert {path.name for path in package.iterdir() if path.is_dir() and path.name != "__pycache__"} == {"assays", "reference_data"}
    assert not list(package.rglob('*.md'))
    assert not any((root / directory).exists() for directory in ('docs', 'examples', 'requirements', 'tools'))
    assert (root/'tests/core/units').is_dir()
    assert not any(path.name.startswith(('test_', '_test', '_family_contracts')) for path in package.rglob('*') if path.is_file())
    print(f"Package layout verified: {len(actual)} runtime modules; all tests outside core/ in top-level tests/.")


if __name__ == "__main__":
    check()
