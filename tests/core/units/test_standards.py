"""Selected public independent checkpoints and a strict official-fixture runner."""
from decimal import Decimal as D
import hashlib
import json
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import pytest
import core.units as u

FIXTURES=Path(__file__).parent/'reference_data/standards/public_factors.json'
PUBLIC=json.loads(FIXTURES.read_text(encoding='utf-8'))
EXTERNAL_STANDARDS=['ISO 80000','ISO 13443','ISO 6976','AGA Report No. 8','GERG-2008','API MPMS 11.1','ASTM D1250','API MPMS 11.2.4','GPA 8217','ISO 6578','GIIGNL','ISO 1928','ISO 18125','ISO 18134','ISO 712','USDA','LBMA','LPPM','ISO 668','JCGM 100','CME','ICE','LME']


@pytest.mark.parametrize('case',PUBLIC['factors'],ids=lambda case:case['source']+'-'+case['unit'])
def test_public_nist_and_eia_conversion_checkpoints(case):
    actual=u.Q_(D(1),case['unit']).to(case['target']).magnitude
    assert abs(D(str(actual))-D(case['value'])) <= abs(D(case['value']))*D(case['rtol'])
    assert PUBLIC['sources'][case['source']].startswith('https://')


@pytest.mark.parametrize('case',PUBLIC['defining_constants'],ids=lambda case:case['unit'])
def test_bipm_si_defining_constant_checkpoints(case):
    assert float(u.Q_(1,case['unit']).to(case['target']).magnitude) == pytest.approx(float(case['value']),rel=2e-15)


def test_iso4217_snapshot_matches_original_six_xml_independently():
    from core.units.currencies import DATA_DIRECTORY,CURRENCY_CATALOG
    raw=ET.fromstring((DATA_DIRECTORY/'list-one.xml').read_bytes())
    by_code={record.code:record for record in CURRENCY_CATALOG.records if record.status == 'current'}
    count=0
    for entry in raw.iter('CcyNtry'):
        code=entry.findtext('Ccy')
        if not code: continue
        record=by_code[code]
        assert record.numeric_code == entry.findtext('CcyNbr').zfill(3)
        digits=entry.findtext('CcyMnrUnts')
        assert record.minor_digits == (int(digits) if digits and digits.isdigit() else None)
        count+=1
    assert count > 150


@pytest.mark.standards_reference
@pytest.mark.parametrize('standard',EXTERNAL_STANDARDS)
def test_supplied_official_reference_cases_are_versioned_and_executed(standard,request):
    root=request.config.getoption('--standards-fixtures')
    if not root:
        pytest.skip(standard+': applicable edition and independently reviewed official fixture not supplied; no compliance claim')
    path=Path(root)/(re.sub('[^a-z0-9]+','_',standard.lower()).strip('_')+'.json')
    if not path.exists(): pytest.skip(standard+': official fixture file absent')
    fixture=json.loads(path.read_text(encoding='utf-8'))
    assert fixture['standard'] == standard and fixture['edition'] and fixture['source'] and fixture['cases']
    assert fixture['source'].startswith('https://')
    for case in fixture['cases']:
        api=case['api']
        if api == 'quantity':
            result=u.Q_(D(case['value']),case['unit']).to(case['target'])
        elif api == 'convert':
            result=u.convert(u.Q_(D(case['value']),case['unit']),case['target'],**case.get('kwargs',{}))
        elif api == 'gas_z':
            result=u.Q_(u.calculate_z_factor(u.Q_(D(case['pressure_kpa']),'kPa'),u.Q_(D(case['temperature_k']),'K'),**case['kwargs']),'dimensionless')
        else: pytest.fail('Unsupported official-fixture adapter: '+api)
        actual=D(str(result.to(case.get('expected_unit',str(result.units))).magnitude))
        expected=D(case['expected'])
        assert abs(actual-expected) <= D(case.get('atol','0'))+abs(expected)*D(case.get('rtol','1e-10'))
