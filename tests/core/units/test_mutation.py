"""Actual source mutants run in isolated package copies, never the live library."""
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import xml.etree.ElementTree as ET
import pytest

pytestmark=pytest.mark.mutation
MUTANTS=[
    ('arithmetic_operator','processing.py','current -= loss','current += loss','test_loss_arithmetic'),
    ('conversion_factor','conversions.py','BARREL_TO_M3: Final[float] = 0.158987294928','BARREL_TO_M3: Final[float] = 0.17','test_barrel_factor'),
    ('reciprocal_transformation','processing.py','factor = eta if source_role == "thermal" else 1 / eta','factor = eta if source_role == "thermal" else eta','test_inverse_generation'),
    ('comparison_operator','processing.py','if self.output_basis not in {"gross", "net"}:','if self.output_basis in {"gross", "net"}:','test_valid_output_basis'),
    ('boundary_condition','heat_content.py','if result > 1:','if result >= 1:','test_efficiency_one_is_valid'),
    ('conditional_logic','processing.py','if source_role != target_role:','if source_role == target_role:','test_generation_conditional'),
    ('boolean_logic','processing.py','if (self.efficiency is None) == (self.heat_rate is None):','if (self.efficiency is None) != (self.heat_rate is None):','test_exactly_one_process_parameter'),
    ('exception_removal','processing.py','raise ValueError("Declare gross or net generation output basis.")','pass','test_invalid_output_basis'),
    ('validation_removal','processing.py','if value <= 0:','if value < 0:','test_zero_cop_rejected'),
    ('rounding_policy','currencies.py','increment, rounding = Decimal(1).scaleb(-record.minor_digits), ROUND_HALF_EVEN','increment, rounding = Decimal(1).scaleb(-record.minor_digits), "ROUND_HALF_UP"','test_currency_half_even_rounding'),
    ('solver_convergence','thermodynamics.py','if abs(f_val) <= conv_threshold:','if abs(f_val) >= conv_threshold:','test_solver_residual'),
    ('currency_conversion','conversions.py','factor = cross_rate_dec * source_scale / target_scale','factor = source_scale / cross_rate_dec / target_scale','test_fx_direction'),
    ('conservation_equation','processing.py','rejected = energy - gross - heat','rejected = energy + gross + heat','test_chp_conservation'),
    ('registry_protection','registry.py','_REGISTRY_FROZEN = True','_REGISTRY_FROZEN = False','test_frozen_registry_is_reported'),
]


@pytest.fixture(scope='session')
def mutation_results(tmp_path_factory,repository_root):
    temporary=tmp_path_factory.mktemp('mutation-campaign')
    results={}
    for name,file,needle,replacement,oracle in MUTANTS:
        target=temporary/name
        shutil.copytree(repository_root/'core',target/'core',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        path=target/'core/units'/file
        text=path.read_text(encoding='utf-8')
        assert needle in text,'Stale mutation target: '+name
        path.write_text(text.replace(needle,replacement,1),encoding='utf-8')
        junit=target/'result.xml'
        environment=dict(os.environ)
        environment['PYTHONPATH']=os.pathsep.join([str(target),str(repository_root),environment.get('PYTHONPATH','')])
        code='import core.units; from pathlib import Path; assert Path(core.units.__file__).is_relative_to(Path('+repr(str(target))+')); import pytest; raise SystemExit(pytest.main('+repr([str(repository_root/'tests/core/units/test_mutation_oracles.py')+'::'+oracle,'-q','--tb=short','--junitxml='+str(junit),'--basetemp='+str(target/'pytest-temp'),'-o','cache_dir='+str(target/'cache')])+'))'
        completed=subprocess.run([sys.executable,'-B','-c',code],cwd=target,env=environment,capture_output=True,text=True,timeout=75)
        assert junit.exists(),'Invalid mutant/import failure '+name+': '+completed.stdout+completed.stderr
        tree=ET.parse(junit)
        errors=sum(int(element.get('errors','0')) for element in tree.getroot().findall('testsuite'))
        failures=sum(int(element.get('failures','0')) for element in tree.getroot().findall('testsuite'))
        assert errors == 0,'Mutation harness error for '+name+': '+completed.stdout+completed.stderr
        results[name]={'killed':completed.returncode == 1 and failures > 0,'returncode':completed.returncode,'failures':failures}
        print('Mutation '+name+': '+('killed' if results[name]['killed'] else 'SURVIVED'),flush=True)
    report={'scope':'seeded source mutants only','mutants':results,'score':sum(value['killed'] for value in results.values())/len(results)}
    output=repository_root/'artifacts/mutation-results.json'
    output.parent.mkdir(exist_ok=True)
    output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    return report


@pytest.mark.parametrize('name',[item[0] for item in MUTANTS])
def test_seeded_source_mutant_is_killed(name,mutation_results):
    assert mutation_results['mutants'][name]['killed'],name


def test_seeded_mutation_score_release_gate(mutation_results):
    assert mutation_results['score'] >= .90
    assert len(mutation_results['mutants']) == 14
