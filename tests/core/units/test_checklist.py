"""No requested checklist item or mirrored source module may silently disappear."""
from collections import Counter
import json
from pathlib import Path
import subprocess
import sys


def test_complete_34_section_700_item_checklist_has_existing_suite_targets(repository_root):
    path=repository_root/'tests/core/units/_test_checklist.json'
    data=json.loads(path.read_text(encoding='utf-8'))
    assert [section['id'] for section in data['sections']] == list(range(1,35))
    rows=[item for section in data['sections'] for item in section['items']]
    assert len(rows) == 700
    assert len({row['id'] for row in rows}) == 700
    assert set(row['status'] for row in rows) <= {'suite','fixture','gate','unsupported'}
    for row in rows:
        assert row['targets'] and row['reason']
        assert all((repository_root/target).exists() for target in row['targets']),row


def test_no_test_modules_remain_at_repository_root(repository_root):
    assert not list(repository_root.glob('test_*.py'))
    assert (repository_root/'tests/core/units/assays/energy/test_petroleum.py').is_file()
    assert (repository_root/'tests/core/units/reference_data/standards/test_standards.py').is_file()
    assert not list((repository_root/'core').rglob('test_*.py'))
