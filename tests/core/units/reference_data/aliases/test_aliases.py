"""Legacy imports are constrained to the canonical implementation namespace."""
import importlib
import json
from pathlib import Path
import core.units as u


def test_import_compatibility_bindings_resolve_to_real_canonical_symbols():
    root=Path(u.__file__).parent
    data=json.loads((root/'reference_data/aliases/python_imports.json').read_text(encoding='utf-8'))
    assert data['schema_version'] == 1
    for bindings in [data['exports'],*data['compatibility'].values()]:
        for module,name in bindings.values():
            assert '..' not in module and '/' not in module and '\\' not in module
            path=root/(module.replace('.','/')+'.py')
            assert path.is_file()
            assert hasattr(importlib.import_module('core.units.'+module),name)
