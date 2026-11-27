"""Pinned reference integrity, constructor compatibility and installed resources."""
from datetime import date
from dataclasses import fields
import importlib
import json
from pathlib import Path
import pkgutil
import shutil
import subprocess

import pytest

import core.units as units
from core.units.standards import (REFERENCE_DIRECTORY, REFERENCE_VERSION,
    reference_record, reference_records, _reference_document)


DIRECTORY = REFERENCE_DIRECTORY
CATEGORIES = tuple(path.stem for path in (DIRECTORY / "assays" / REFERENCE_VERSION).glob("*.json"))


def test_git_checkout_preserves_checksum_protected_reference_bytes(repository_root):
    if not shutil.which('git') or not (repository_root/'.git').exists():
        pytest.skip('Git checkout attributes apply only to repository checkouts')
    paths=sorted(path.relative_to(repository_root).as_posix()
                 for path in DIRECTORY.rglob('*') if path.suffix in {'.json','.xml'})
    result=subprocess.run(['git','-c','safe.directory='+repository_root.as_posix(),
                           'check-attr','-z','--stdin','text'],input=('\0'.join(paths)+'\0').encode(),
                          cwd=repository_root,capture_output=True,timeout=15)
    assert result.returncode == 0,result.stderr.decode(errors='replace')
    fields=result.stdout.decode().rstrip('\0').split('\0')
    attributes={fields[index]:fields[index+2] for index in range(0,len(fields),3)}
    assert set(attributes) == set(paths)
    assert all(value == 'unset' for value in attributes.values()),attributes


@pytest.mark.parametrize("category", CATEGORIES)
def test_every_preset_loads_and_preserves_constructor_arguments(category):
    records = reference_records(category)
    assert records
    for record in records:
        classname, method = record.id.split(".")
        cls = getattr(units, classname)
        actual = getattr(cls, method)()
        args, kwargs = record.constructor_arguments()
        expected = cls(*args, **kwargs)
        # Public fields of both objects must retain the migrated benchmark.
        assert all(getattr(actual, field.name) == getattr(expected, field.name) for field in fields(actual))
        assert record.source["citation"]
        assert record.value_kind == "indicative_benchmark"
        assert record.uncertainty["status"] == "unknown"
        with pytest.raises(TypeError):
            record.source["citation"] = "changed"
        with pytest.raises(ValueError, match="validity is unknown"):
            reference_record(category, record.id, as_of=date(2026, 11, 27))


def test_benchmark_data_is_defensive_and_brent_values_are_compatible():
    record = reference_record("crude_oil", "CrudeAssay.brent")
    args, kwargs = record.constructor_arguments()
    kwargs["api_gravity"] = 0
    assert record.constructor_arguments()[1]["api_gravity"] == 38.3
    assert units.CrudeAssay.brent().api_gravity == pytest.approx(38.3)
    assert units.CrudeAssay.brent().sulfur_percentage == .4


def test_parameterised_presets_keep_explicit_caller_overrides():
    assert units.GrainAssay.corn_56lb(actual_moisture_pct=20).actual_moisture_pct == 20
    assert units.UraniumAssay.leu_4_5pct(tails_u235_pct=.3).tails_u235_pct == .3
    olive = units.OliveOilAssay.lampante_crude(free_fatty_acid_pct=4, refining_yield_pct=93)
    assert olive.free_fatty_acid_pct == 4
    assert olive.refining_yield_pct == 93
    assert "4%" in olive.name


def test_reference_manifest_and_all_documents_are_checked():
    from hashlib import sha256
    manifest = json.loads((DIRECTORY / "standards" / REFERENCE_VERSION / "manifest.json").read_text(encoding="utf-8"))
    for relative, digest in manifest["files"].items():
        section, name = relative.split("/")
        assert sha256((DIRECTORY / section / REFERENCE_VERSION / name).read_bytes()).hexdigest() == digest
        assert _reference_document(section, Path(name).stem)["version"] == REFERENCE_VERSION
    with pytest.raises(ValueError):
        reference_records("../crude_oil")
    with pytest.raises(ValueError):
        reference_records("crude_oil", version="../2026-11-27")


def test_corruption_is_rejected(tmp_path, monkeypatch):
    import core.units.standards as reference
    version = "test-corruption"
    directory = tmp_path
    (directory / "assays" / version).mkdir(parents=True)
    (directory / "standards" / version).mkdir(parents=True)
    (directory / "assays" / version / "crude_oil.json").write_text('{}', encoding="utf-8")
    (directory / "standards" / version / "manifest.json").write_text(json.dumps({"schema_version": 1,
        "version": version, "files": {"assays/crude_oil.json": "0" * 64}}), encoding="utf-8")
    monkeypatch.setattr(reference, "REFERENCE_DIRECTORY", tmp_path)
    with pytest.raises(ValueError, match="integrity"):
        reference_records("crude_oil", version=version)


def test_sector_interfaces_preserve_class_identity():
    import core.units.assays as assays
    for module in pkgutil.walk_packages(assays.__path__, assays.__name__ + "."):
        interface = importlib.import_module(module.name)
        for name in getattr(interface, "__all__", ()):
            assert getattr(interface, name) is getattr(units, name)


def test_recommended_standards_never_assert_compliance():
    from core.units.standards import recommended_standards
    standards = recommended_standards()
    assert len(standards) == 20
    assert all(record["compliance_verified"] is False for record in standards)
    assert all(record["source"]["edition"] is None for record in standards)


def test_modular_interfaces_batch_and_local_precision():
    names = ("types", "dimensions", "context", "validation", "precision", "quote_parser",
        "reference_conditions", "analytical_basis", "composition", "blending", "processing",
        "thermodynamics", "logistics_units", "provenance", "serialisation", "batch", "uncertainty")
    for name in names:
        importlib.import_module("core.units." + name)
    from decimal import getcontext
    from core.units.precision import decimal_precision
    from core.units.batch import convert_many
    original = getcontext().prec
    with decimal_precision(40):
        assert getcontext().prec == 40
    assert getcontext().prec == original
    assert [q.magnitude for q in convert_many((units.Q_(1, "kg"), units.Q_(2, "kg")), "g")] == [1000, 2000]
