"""Gate 4 search summaries must remain traceable to completed training-only runs."""
import hashlib
import json
from pathlib import Path
import re

import pytest

from scripts.check_data_card_numbers import (gate4_combined_checks, gate4_control_checks,
                                            load_gate4_combined, load_gate4_controls)


ROOT = Path(__file__).resolve().parents[1]
DIRECTORY = ROOT / "results/gate4_training/CHEMBL284"
CARD = (ROOT / "DATA_CARD.md").read_text()
LOADERS = {"descriptors": load_gate4_controls, "combined": load_gate4_combined}
CHECKS = {"descriptors": gate4_control_checks, "combined": gate4_combined_checks}
ROWS = [(block, line) for block, prefix, label in (
            ("descriptors", "Control", "Tuned descriptor control"),
            ("combined", "Combined", "Tuned descriptors+Morgan"))
        for line in CARD.splitlines()
        if (line.startswith(f"| {prefix} ")
            and not line.startswith(f"| {prefix} protocol / result "))
        or line.startswith((f"| scaffold / {label} (OOF) |",
                            f"| random / {label} (OOF) |"))]


@pytest.fixture(scope="module")
def runs():
    return {block: loader(DIRECTORY) for block, loader in LOADERS.items()}


@pytest.mark.parametrize("block, count", [("descriptors", 25), ("combined", 27)])
def test_completed_searches_match_every_displayed_result(runs, block, count):
    checks = CHECKS[block](CARD, runs[block])
    if block == "combined" and len(runs[block]) < 2:
        assert all(ok for _, ok in checks)
        assert all(b != block for b, _ in ROWS)
        return
    assert len(checks) == sum(b == block for b, _ in ROWS) == count
    assert all(ok for _, ok in checks)


@pytest.mark.parametrize("block, row", ROWS)
def test_every_search_cell_and_numeric_component_is_guarded(runs, block, row):
    cells = [cell.strip() for cell in row.strip("|").split("|")]
    # Mutate each value cell, including nonnumeric choices such as sqrt and False.
    for column in range(1, len(cells)):
        wrong = cells.copy()
        wrong[column] = "wrong"
        changed = CARD.replace(row, "| " + " | ".join(wrong) + " |")
        assert not all(ok for _, ok in CHECKS[block](changed, runs[block]))
        # Also move each number by the smallest displayed precision. This catches
        # an interval endpoint or one fold choice drifting inside a compound cell.
        for match in re.finditer(r"-?\d[\d,]*(?:\.\d+)?", cells[column]):
            token = match.group()
            precision = len(token.split(".")[1]) if "." in token else 0
            nudged = f"{float(token.replace(',', '')) + 10 ** -precision:.{precision}f}"
            wrong = cells.copy()
            wrong[column] = cells[column][:match.start()] + nudged + cells[column][match.end():]
            changed = CARD.replace(row, "| " + " | ".join(wrong) + " |")
            assert not all(ok for _, ok in CHECKS[block](changed, runs[block]))


@pytest.mark.parametrize("block, row", ROWS)
@pytest.mark.parametrize("duplicate", [False, True])
def test_each_missing_or_duplicate_search_row_fails(runs, block, row, duplicate):
    changed = CARD.replace(row, row + "\n" + row if duplicate else "")
    assert not all(ok for _, ok in CHECKS[block](changed, runs[block]))


@pytest.mark.parametrize("partial", [{}, {"scaffold": None}, {"random": None}])
def test_partial_search_rejects_published_scores(partial):
    assert all(ok for _, ok in gate4_combined_checks("Combined searches pending.", partial))
    for row in ("| Combined forest trees | 200 | 200 |",
                "| scaffold / Tuned descriptors+Morgan (OOF) | 0.99 |"):
        assert not all(ok for _, ok in gate4_combined_checks(row, partial))


@pytest.fixture(scope="module")
def completed_sources(tmp_path_factory):
    sources = {block: DIRECTORY / "scaffold" / block for block in LOADERS}
    if not (sources["combined"] / "report.json").exists():
        # Exercise the same combined report integrity checks before the real search finishes.
        from src import nested_cv
        from test_gate4_combined import _run

        root = tmp_path_factory.mktemp("finished-combined")
        with pytest.MonkeyPatch.context() as patch:
            patch.setattr(nested_cv, "RF_SPEC", {**nested_cv.RF_SPEC, "n_estimators": 2})
            _run(root)
        sources["combined"] = root / "combined"
    return sources


@pytest.mark.parametrize("corrupt", ["manifest", "report", "identity", "name"])
@pytest.mark.parametrize("block", LOADERS)
def test_search_loader_refuses_unchecksummed_changes(tmp_path, corrupt, block, completed_sources):
    path = tmp_path / "scaffold" / block
    path.mkdir(parents=True)
    for name in ("manifest.json", "report.json"):
        (path / name).write_bytes((completed_sources[block] / name).read_bytes())
    if corrupt == "manifest":
        name = "manifest.json"
        value = json.loads((path / name).read_text())
        value["seed"] += 1
    else:
        name = "report.json"
        value = json.loads((path / name).read_text())
        if corrupt == "report":
            value["result"]["oof_metrics"]["rmse"] += 0.001
        else:
            value[corrupt] = "wrong"
    (path / name).write_text(json.dumps(value))
    with pytest.raises(ValueError, match="checksum/identity"):
        LOADERS[block](tmp_path)


@pytest.mark.parametrize("corrupt", ["partition", "heldout_evaluations"])
@pytest.mark.parametrize("block", LOADERS)
def test_search_loader_refuses_rechecksummed_heldout_results(tmp_path, corrupt, block, completed_sources):
    def digest(value):
        return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                         allow_nan=False).encode()).hexdigest()

    path = tmp_path / "scaffold" / block
    path.mkdir(parents=True)
    manifest = json.loads((completed_sources[block] / "manifest.json").read_text())
    record = json.loads((completed_sources[block] / "report.json").read_text())
    if corrupt == "partition":
        manifest["rows"][0]["partition"] = "test"
    else:
        record["result"]["heldout_evaluations"] = 1
    record["identity"] = record["result"]["identity_sha256"] = digest(manifest)
    record["sha256"] = digest(record["result"])
    (path / "manifest.json").write_text(json.dumps(manifest))
    (path / "report.json").write_text(json.dumps(record))
    with pytest.raises(ValueError, match="training-only"):
        LOADERS[block](tmp_path)
