"""Distribution metadata, optional absence and opt-in installation acceptance."""
import importlib.metadata
import json
import os
import shutil
from pathlib import Path
import subprocess
import sys
import tarfile
import tomllib
import zipfile
import pytest
import core.units as u
from tests.core.units._test_layout import expected_files


def test_metadata_and_current_supported_python_pint_versions(repository_root):
    from packaging.specifiers import SpecifierSet
    metadata=tomllib.loads((repository_root/'pyproject.toml').read_text(encoding='utf-8'))['project']
    assert metadata['name'] == 'commodex'
    assert sys.version.split()[0] in SpecifierSet(metadata['requires-python'])
    import pint
    assert pint.__version__ in SpecifierSet(metadata['dependencies'][0].removeprefix('pint'))


def test_import_without_optional_backends_and_no_network_side_effects():
    program="""
import importlib.abc, sys, socket
blocked={'numpy','pandas','pyarrow','pint_pandas','yaml','pyaga8','CoolProp','uncertainties'}
class Block(importlib.abc.MetaPathFinder):
    def find_spec(self,fullname,path=None,target=None):
        if fullname.split('.')[0] in blocked: raise ModuleNotFoundError(fullname)
sys.meta_path.insert(0,Block())
def forbidden(*args,**kwargs): raise AssertionError('import attempted networking')
socket.create_connection=forbidden
import core.units as u
assert u.convert(u.Q_(1,'kg'),'g').magnitude == 1000
assert u.calculate_z_factor(u.Q_(100,'kPa'),u.Q_(300,'K'),method='ideal') == 1
"""
    result=subprocess.run([sys.executable,'-B','-c',program],capture_output=True,text=True,timeout=25)
    assert result.returncode == 0,result.stderr


@pytest.mark.parametrize('order',['ascending','descending','permuted'])
def test_all_module_import_orders_are_cycle_safe(order,repository_root):
    names=sorted('core.units.'+name[:-3].replace('/','.') for name in expected_files() if not name.endswith('__init__.py'))
    if order == 'descending': names.reverse()
    elif order == 'permuted':
        import random
        random.Random(444).shuffle(names)
    code='import importlib\n'+''.join('importlib.import_module('+repr(name)+')\n' for name in names)+'import core.units as u\nassert u.get_registry() is u.ureg\n'
    result=subprocess.run([sys.executable,'-B','-c',code],capture_output=True,text=True,timeout=30)
    assert result.returncode == 0,result.stderr


@pytest.mark.release_install
def test_clean_wheel_and_source_distribution_installation(tmp_path,repository_root):
    pytest.importorskip('build',reason='Install the dev build tool for source-distribution acceptance')
    output=tmp_path/'dist'
    clean_source=tmp_path/'source'
    clean_source.mkdir()
    for directory in ('core','tests'):
        shutil.copytree(repository_root/directory,clean_source/directory,ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
    for filename in ('pyproject.toml','pylock.toml','MANIFEST.in','requirements.txt'):
        shutil.copyfile(repository_root/filename,clean_source/filename)
    build=subprocess.run([sys.executable,'-B','-m','build','--no-isolation','--outdir',str(output)],cwd=clean_source,capture_output=True,text=True,timeout=180)
    assert build.returncode == 0,build.stdout+build.stderr
    wheel=next(output.glob('*.whl'))
    source=next(output.glob('*.tar.gz'))
    with zipfile.ZipFile(wheel) as archive:
        files={name.removeprefix('core/units/') for name in archive.namelist() if name.startswith('core/units/') and name.endswith('.py')}
        assert files == expected_files()
        assert not any(name.startswith('tests/') or '/test_' in name or '/_test' in name for name in archive.namelist())
        assert 'core/units/reference_data/aliases/python_imports.json' in archive.namelist()
    with tarfile.open(source) as archive:
        assert any(name.endswith('/tests/core/units/_test_checklist.json') for name in archive.getnames())
    for artifact,label in ((wheel,'wheel'),(source,'sdist')):
        target=tmp_path/label
        installed=subprocess.run([sys.executable,'-B','-m','pip','install','--no-deps','--no-build-isolation','--target',str(target),str(artifact)],capture_output=True,text=True,timeout=120)
        assert installed.returncode == 0,installed.stdout+installed.stderr
        code = 'import sys; sys.path.insert(0,'+repr(str(target))+'); from pathlib import Path; import core.units as u; assert Path(u.__file__).is_relative_to(Path('+repr(str(target))+')); u.verify_registry_integrity(); assert u.CrudeAssay.brent().sulfur_percentage == .4; from core.units.standards import reference_records; assert len(reference_records("natural_gas")) == 9'
        checked=subprocess.run([sys.executable,'-B','-c',code],cwd=tmp_path,capture_output=True,text=True,timeout=30)
        assert checked.returncode == 0,checked.stdout+checked.stderr


@pytest.mark.release_install
def test_editable_installation_resolves_back_to_source(tmp_path,repository_root):
    target=tmp_path/'editable'
    result=subprocess.run([sys.executable,'-B','-m','pip','install','--no-deps','--no-build-isolation','--target',str(target),'-e',str(repository_root)],capture_output=True,text=True,timeout=120)
    assert result.returncode == 0,result.stdout+result.stderr
    code='import site; site.addsitedir('+repr(str(target))+'); import core.units as u; from pathlib import Path; assert Path(u.__file__).resolve() == Path('+repr(str(repository_root/'core/units/__init__.py'))+').resolve()'
    result=subprocess.run([sys.executable,'-B','-c',code],cwd=tmp_path,capture_output=True,text=True,timeout=30)
    assert result.returncode == 0,result.stderr
