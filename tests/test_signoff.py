"""The human handoff must fail validation when a displayed result goes stale."""
import json
from pathlib import Path

import pytest

from scripts.check_data_card_numbers import signoff_table_checks


ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("old,new", [
    (None, None),
    ("**3,965 compounds**", "**3,966 compounds**"),
    ("0.970 [0.857, 1.088]", "0.970 [0.858, 1.088]"),
    ("| Shared exact scaffolds | 268 | 0 |", ""),
    ("| Shared exact scaffolds | 268 | 0 |",
     "| Shared exact scaffolds | 268 | 0 |\n| Shared exact scaffolds | 268 | 0 |"),
])
def test_signoff_tables_detect_stale_missing_and_duplicate_results(old, new):
    text = (ROOT / "GATE2_SIGNOFF.md").read_text()
    if old is not None:
        assert old in text
        text = text.replace(old, new)
    logs = [json.loads((ROOT / path).read_text()) for path in (
        "results/curation/CHEMBL284/curation_log.json",
        "results/splits/CHEMBL284/audit.json",
        "results/leakage_probe/CHEMBL284/probe.json")]
    checks = signoff_table_checks(text, *logs)
    assert len(checks) == 15
    assert all(ok for _, ok in checks) == (old is None)
