"""Configuration schemas, format parity, atomicity and cache invalidation."""
import json
from pathlib import Path
import pytest
import pint
import core.units as u
from core.units.registry import _registry_transaction


@pytest.mark.parametrize('extension',['mapping','json','yaml'])
def test_configuration_formats_agree_and_rollback_outer_transaction(extension,tmp_path):
    payload={'commodities':[{'definition':'checklist_parcel = 3 * kilogram'}]}
    if extension == 'mapping': source=payload
    else:
        source=tmp_path/('desk.'+extension)
        if extension == 'json': source.write_text(json.dumps(payload),encoding='utf-8')
        else:
            yaml=pytest.importorskip('yaml')
            source.write_text(yaml.safe_dump(payload),encoding='utf-8')
    before=u.get_registry_definition_hash()
    class Rollback(Exception): pass
    with pytest.raises(Rollback):
        with _registry_transaction():
            u.load_desk_config(source,freeze_after_load=False)
            assert u.Q_(2,'checklist_parcel').to('kg').magnitude == 6
            assert u.get_registry_definition_hash() != before
            raise Rollback()
    assert u.get_registry_definition_hash() == before
    with pytest.raises(pint.UndefinedUnitError): u.parse_unit('checklist_parcel')


@pytest.mark.parametrize('payload',[[],{'currencies':{}},{'commodities':'kg'},{'currencies':[{}]},{'commodities':[{'definition':42}]}])
def test_configuration_schema_rejects_invalid_shapes_without_partial_registration(payload):
    before=u.get_registry_definition_hash()
    with pytest.raises((TypeError,ValueError)):
        u.load_desk_config(payload,freeze_after_load=False)
    assert u.get_registry_definition_hash() == before


def test_duplicate_and_conflicting_definitions_are_atomic():
    before=u.get_registry_definition_hash()
    with pytest.raises((ValueError,pint.PintError)):
        u.load_desk_config({'commodities':[{'definition':'checklist_duplicate = 2 * kg'},{'definition':'checklist_duplicate = 3 * kg'}]},freeze_after_load=False)
    assert u.get_registry_definition_hash() == before
