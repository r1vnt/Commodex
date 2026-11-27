"""Static, dependency and release-evidence checks."""
import ast
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import pytest


def test_all_sources_compile(repository_root):
    for root in (repository_root/'core',repository_root/'tests'):
        for path in root.rglob('*.py'): compile(path.read_text(encoding='utf-8'),str(path),'exec')
    from tests.core.units._test_layout import check
    check(repository_root)


def test_critical_lint_rules_for_imports_names_syntax_and_control_flow(repository_root):
    pytest.importorskip('ruff',reason='Install commodex[dev] for the static lint gate')
    result=subprocess.run([sys.executable,'-B','-m','ruff','check','core/units','tests','--select','E9,F63,F7,F82','--output-format','concise'],cwd=repository_root,capture_output=True,text=True,timeout=60)
    assert result.returncode == 0,result.stdout+result.stderr


def test_numeric_and_storage_boundary_type_checks(repository_root):
    pytest.importorskip('mypy',reason='Install commodex[dev] for gradual boundary type checks')
    result=subprocess.run([sys.executable,'-B','-m','mypy','--follow-imports=skip','--ignore-missing-imports','--cache-dir',str(repository_root/'artifacts/mypy-cache'),'core/units/precision.py','core/units/batch.py','core/units/serialisation.py'],cwd=repository_root,capture_output=True,text=True,timeout=60)
    assert result.returncode == 0,result.stdout+result.stderr


def test_high_confidence_dead_and_unreachable_code_detection(repository_root):
    pytest.importorskip('vulture',reason='Install the dev vulture tool for high-confidence dead-code analysis')
    from tests.core.units._test_layout import expected_files
    result=subprocess.run([sys.executable,'-B','-m','vulture',*('core/units/'+path for path in expected_files()),'--min-confidence','100'],cwd=repository_root,capture_output=True,text=True,timeout=60)
    assert result.returncode == 0,result.stdout+result.stderr


def test_cyclomatic_complexity_analysis_records_real_function_metrics(repository_root):
    pytest.importorskip('radon',reason='Install the dev radon tool for cyclomatic-complexity reporting')
    from tests.core.units._test_layout import expected_files
    result=subprocess.run([sys.executable,'-B','-m','radon','cc',*('core/units/'+path for path in expected_files()),'-j'],cwd=repository_root,capture_output=True,text=True,timeout=60)
    assert result.returncode == 0,result.stderr
    values=json.loads(result.stdout)
    functions=[item for entries in values.values() for item in entries]
    assert functions and all(item['complexity'] >= 1 for item in functions)
    (repository_root/'artifacts/complexity.json').write_text(json.dumps(values,indent=2)+'\n',encoding='utf-8')


def test_public_api_family_and_semantic_path_coverage_are_declared(repository_root):
    import core.units as units
    assert len(units.__all__) == len(set(units.__all__))
    assert len(units.__all__) > 500
    assert all(hasattr(units,name) for name in units.__all__)
    from tests.core.units._test_layout import expected_files
    base=repository_root/'tests/core/units'
    missing=[]
    for path in expected_files():
        if path.endswith('__init__.py'): continue
        source=Path(path)
        target=base/source.parent/('test_'+source.name)
        if not target.is_file(): missing.append(target.relative_to(repository_root).as_posix())
    assert not missing,missing


@pytest.mark.network_audit
def test_pinned_dependency_vulnerability_audit(repository_root,tmp_path):
    pytest.importorskip('pip_audit',reason='Install commodex[dev] for dependency advisory auditing')
    import tomllib
    packages=tomllib.loads((repository_root/'pylock.toml').read_text(encoding='utf-8'))['packages']
    requirements=tmp_path/'locked-requirements.txt'
    requirements.write_text('\n'.join(item['name']+'=='+item['version'] for item in packages if 'version' in item)+'\n',encoding='utf-8')
    result=subprocess.run([sys.executable,'-B','-m','pip_audit','--no-deps','--disable-pip','-r',str(requirements),'--format','json'],cwd=repository_root,capture_output=True,text=True,timeout=180)
    (repository_root/'artifacts/dependency-audit.json').write_text(result.stdout,encoding='utf-8')
    assert result.returncode == 0,result.stdout+result.stderr
