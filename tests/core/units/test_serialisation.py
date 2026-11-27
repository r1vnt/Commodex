"""Scalar schema, trusted pickle and spawned-process compatibility contracts."""
from decimal import Decimal as D
import json
import multiprocessing
import pickle
import pytest
import core.units as u


def _restore_in_child(payload, connection):
    value = u.quantity_from_dict(payload)
    connection.send(u.quantity_to_dict(value, include_metadata=True))
    connection.close()


@pytest.mark.parametrize("unit", ["kg","USD/MWh_th","wmt/day","EUA","heat_rate","kpag"])
def test_quantity_unit_decimal_json_and_trusted_pickle_roundtrips(unit):
    original = u.Q_(D("1.234567890123456789012345"), unit)
    payload = u.quantity_to_dict(original, include_metadata=True)
    restored = u.quantity_from_dict(json.loads(json.dumps(payload)))
    assert restored == original
    assert isinstance(restored.magnitude, D)
    assert pickle.loads(pickle.dumps(original)) == original


def test_assay_and_reference_condition_trusted_pickle_roundtrip():
    for original in (u.CrudeAssay.brent(),u.GasReferenceConditions.from_standard("iso15")):
        assert pickle.loads(pickle.dumps(original)) == original


def test_cross_process_spawn_quantity_schema_and_registry_hash():
    context = multiprocessing.get_context("spawn")
    parent, child = context.Pipe(duplex=False)
    payload = u.quantity_to_dict(u.Q_(D("123.45"),"MWh_e"),include_metadata=True)
    process = context.Process(target=_restore_in_child,args=(payload,child))
    process.start()
    try:
        assert parent.poll(20), "spawned deserialisation timed out"
        assert parent.recv() == payload
    finally:
        process.join(20)
        if process.is_alive(): process.terminate(); process.join(5)
        parent.close(); child.close()
    assert process.exitcode == 0


@pytest.mark.parametrize("change", [{"schema_version":999},{"definition_hash":"0"*64},{"registry_version":"unknown"}])
def test_unknown_or_incompatible_schema_and_registry_metadata_are_rejected(change):
    payload = u.quantity_to_dict(u.Q_(D(1),"kg"),include_metadata=True)
    payload.update(change)
    with pytest.raises((ValueError, u.RegistryCompatibilityError)):
        u.quantity_from_dict(payload)
