"""Benchmark identities, version selection and reference argument ownership."""
import pytest
import core.units as u


@pytest.mark.parametrize('category',['crude_oil','refined_products','natural_gas','coal','alternative_fuels','metals','chemicals','agriculture'])
def test_benchmarks_do_not_claim_measured_uncertainty_or_unknown_historical_validity(category):
    for record in u.reference_records(category):
        assert record.value_kind == 'indicative_benchmark'
        assert record.uncertainty['status'] == 'unknown'
        assert record.source['citation']
        with pytest.raises(ValueError,match='validity'):
            u.reference_record(category,record.id,as_of='2020-01-01')
