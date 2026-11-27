"""Every frozen unit definition is loadable and has unambiguous provenance fields."""
from pathlib import Path
import json
import core.units as u
from core.units.standards import REFERENCE_DIRECTORY,REFERENCE_VERSION


def test_all_versioned_unit_definitions_and_aliases_are_registered():
    path=REFERENCE_DIRECTORY/'units'/REFERENCE_VERSION/'registry_units.json'
    records=json.loads(path.read_text(encoding='utf-8'))['records']
    assert len(records) == 333
    for record in records:
        assert record['source']['citation'] and record['units']
        assert 'valid_from' in record and 'valid_to' in record and 'uncertainty' in record
        assert u.ureg.Unit(record['id'])
        pieces=[value.strip() for value in record['definition'].split('=')]
        for alias in pieces[2:]:
            assert u.ureg.Unit(alias) == u.ureg.Unit(record['id'])
