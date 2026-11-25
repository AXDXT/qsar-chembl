"""The human handoff must fail validation when a displayed result goes stale."""
import json
from pathlib import Path

import pytest

from scripts.check_data_card_numbers import gate3_feature_checks, signoff_table_checks


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


@pytest.mark.parametrize("old,new", [
    (None, None),
    ("| Rows featurised | 3,965 |", "| Rows featurised | 3,964 |"),
    ("| Featurisation failures | 0 |", "| Featurisation failures | 1 |"),
    ("| Occupied bits at 2,048 | 2,028 |", "| Occupied bits at 2,048 | 2,048 |"),
    ("| Morgan bits (radius 2, nonchiral) | 2,048 |",
     "| Morgan bits (radius 2, nonchiral) | 1,024 |"),
    ("| Largest absolute r among descriptors, random training rows | 0.980 (HeavyAtomCount / MolWt) |",
     "| Largest absolute r among descriptors, random training rows | 0.981 (HeavyAtomCount / MolWt) |"),
    ("| Compounds losing ≥1 identifier at 4,096 | 905 |", ""),
    ("| Descriptor columns | 12 |", "| Descriptor columns | 12 |\n| Descriptor columns | 12 |"),
])
def test_gate3_feature_figures_detect_stale_missing_and_duplicate_rows(old, new):
    text = (ROOT / "DATA_CARD.md").read_text()
    if old is not None:
        assert old in text
        text = text.replace(old, new)
    manifest = json.loads((ROOT / "data/features/CHEMBL284/manifest.json").read_text())
    checks = gate3_feature_checks(text, manifest)
    assert len(checks) == 16
    assert all(ok for _, ok in checks) == (old is None)
