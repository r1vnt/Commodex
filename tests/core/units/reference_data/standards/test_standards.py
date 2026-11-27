"""Standards inventory coverage and explicit limits on compliance claims."""
from core.units.standards import recommended_standards


def test_all_recommended_standards_are_provenance_records_not_compliance_certificates():
    values=recommended_standards()
    assert len(values) == 20
    assert all(record['status'] == 'recommended_reference' for record in values)
    assert all(record['compliance_verified'] is False for record in values)
    assert all(record['source']['citation'] for record in values)
