"""Currency snapshots expose historical and current provenance separately."""
from core.units.currencies import CURRENCY_CATALOG,DATA_DIRECTORY


def test_currency_reference_resources_remain_available_after_packaging_moves():
    for name in ('manifest.json','reference.json','overlay.json','conventions.json','list-one.xml','list-three.xml'):
        assert (DATA_DIRECTORY/name).is_file()
    assert any(record.status == 'historical' for record in CURRENCY_CATALOG.records)
    assert any(record.status == 'current' for record in CURRENCY_CATALOG.records)
