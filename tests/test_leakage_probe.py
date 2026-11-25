"""Offline leakage-probe checks: train-only prediction, uncertainty and input provenance."""
import csv
import json
import math

import numpy as np
import pytest
from rdkit import DataStructs

from src import leakage_probe as probe, retrieve, splits
from test_splits import _curated_fixture


def _bits(*indices):
    fp = DataStructs.ExplicitBitVect(8)
    for i in indices:
        fp.SetBit(i)
    return fp


def _toy():
    rows = [{"inchikey": key, "scaffold": core, "pIC50_median": label, "n_records": "1"}
            for key, core, label in [("A", "core1", "4"), ("B", "core2", "8"),
                                     ("C", "core3", "99"), ("D", "core4", "-99")]]
    return rows, [_bits(0, 1), _bits(0, 2), _bits(0), _bits(0, 1)]


def test_predictions_use_training_labels_only_and_average_all_maximum_ties():
    rows, fps = _toy()
    assignment = ["train", "train", "test", "test"]
    result = probe.predict_partition(rows, fps, assignment, scaffold_disjoint=True)
    assert [p["predicted_pIC50"] for p in result] == [6., 4.]
    assert [p["maximum_similarity_ties"] for p in result] == [2, 1]
    assert result[0]["nearest_tanimoto"] == .5
    rows[2]["pIC50_median"], rows[3]["pIC50_median"] = "0", "0"
    changed = probe.predict_partition(rows, fps, assignment)
    assert [p["predicted_pIC50"] for p in changed] == [6., 4.]
    # Reordering training rows must not select a different label among tied neighbours.
    order = [1, 0, 2, 3]
    reversed_train = probe.predict_partition([rows[i] for i in order], [fps[i] for i in order], assignment)
    assert [p["predicted_pIC50"] for p in reversed_train] == [6., 4.]


def test_prediction_refuses_identity_overlap_or_scaffold_overlap():
    rows, fps = _toy()
    rows[2]["inchikey"] = "A"
    with pytest.raises(AssertionError, match="InChIKey leakage"):
        probe.predict_partition(rows, fps, ["train", "train", "test", "test"])
    rows[2]["inchikey"] = "C"
    rows[2]["scaffold"] = "core1"
    with pytest.raises(AssertionError, match="Scaffold leakage"):
        probe.predict_partition(rows, fps, ["train", "train", "test", "test"], True)


def test_metrics_and_bootstrap_are_scaffold_grouped_not_individual_row_resampling():
    predictions = [{"scaffold": core, "observed_pIC50": y, "predicted_pIC50": y + e}
                   for core, y, e in [("A", 0., 0.), ("A", 2., 0.), ("B", 4., 3.)]]
    metrics = probe.point_metrics(predictions)
    assert metrics == pytest.approx({"rmse": np.sqrt(3), "mae": 1., "r2": -.125})
    intervals = probe.clustered_intervals(predictions, repeats=100)
    assert intervals == probe.clustered_intervals(predictions, repeats=100)
    assert intervals["scaffold_groups"] == 2
    # Whole-group draws can select both A, both B, or A+B. They preserve the two A rows.
    assert intervals["metrics"]["rmse"]["ci95"] == [0., 3.]
    assert intervals["metrics"]["r2"]["valid_draws"] < 100  # B+B has constant y


def test_constant_targets_have_undefined_r2_not_fabricated_perfect_score():
    predictions = [{"scaffold": "", "observed_pIC50": 1., "predicted_pIC50": 2.}] * 3
    assert probe.point_metrics(predictions) == {"rmse": 1., "mae": 1., "r2": None}
    interval = probe.clustered_intervals(predictions, repeats=10)
    assert interval["scaffold_groups"] == 1
    assert interval["metrics"]["r2"] == {"ci95": None, "valid_draws": 0}


@pytest.mark.parametrize("predictions", [[], [{"observed_pIC50": np.nan, "predicted_pIC50": 1.}]])
def test_invalid_metric_input_is_refused(predictions):
    with pytest.raises(ValueError, match="nonempty finite"):
        probe.point_metrics(predictions)


def _inputs(tmp_path):
    raw, curated = _curated_fixture(tmp_path)
    split_dir = tmp_path / "splits"
    splits.split_target(raw, curated, split_dir)
    return raw, curated, split_dir


def test_end_to_end_is_byte_reproducible_checksummed_and_preserves_inputs(tmp_path):
    raw, curated, split_dir = _inputs(tmp_path)
    before = {p: p.read_bytes() for directory in (raw, curated, split_dir)
              for p in directory.iterdir() if not p.name.startswith("._")}
    out = tmp_path / "probe"
    report = probe.run_probe(raw, curated, split_dir, out)
    first = {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert report == probe.run_probe(raw, curated, split_dir, out)
    assert first == {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert all(p.read_bytes() == contents for p, contents in before.items())
    for filename, meta in report["files"].items():
        assert retrieve._sha256(out / filename) == meta["sha256"]
    for name in ("random", "scaffold"):
        with (out / f"{name}_predictions.csv").open() as fh:
            predictions = list(csv.DictReader(fh))
        assert len(predictions) == report["splits"][name]["test_compounds"]
        assert {r["inchikey"] for r in predictions} == {
            r["inchikey"] for r in probe._read_csv(split_dir / f"{name}_assignments.csv")
            if r["partition"] == "test"}


@pytest.mark.parametrize("corrupt", ["raw", "curated", "split", "provenance", "assignment_labels"])
def test_corrupt_inputs_or_tampered_test_labels_fail_before_output(tmp_path, corrupt):
    raw, curated, split_dir = _inputs(tmp_path)
    if corrupt == "raw":
        (raw / retrieve.RAW_TARGET).write_text("{}")
    elif corrupt == "curated":
        with (curated / "compounds.csv").open("a") as fh:
            fh.write("corrupted\n")
    elif corrupt == "split":
        with (split_dir / "random_assignments.csv").open("a") as fh:
            fh.write("corrupted\n")
    else:
        audit = json.loads((split_dir / "audit.json").read_text())
        if corrupt == "provenance":
            audit["raw_manifest_sha256"] = "wrong"
        else:
            # Even recomputing the CSV checksum must not make changed held-out labels valid.
            path = split_dir / "random_assignments.csv"
            assignments = probe._read_csv(path)
            assignments[0]["pIC50_median"] = "100"
            splits._write_csv(path, assignments, list(assignments[0]))
            audit["files"][path.name]["sha256"] = retrieve._sha256(path)
        (split_dir / "audit.json").write_text(json.dumps(audit))
    with pytest.raises(ValueError, match="checksum|provenance|disagree"):
        probe.run_probe(raw, curated, split_dir, tmp_path / "probe")
    assert not (tmp_path / "probe").exists()


def test_outputs_cannot_overwrite_any_source_tree(tmp_path):
    raw, curated, split_dir = _inputs(tmp_path)
    for dest in (raw, raw / "derived", curated, split_dir, split_dir / "derived", tmp_path):
        with pytest.raises(ValueError, match="separate"):
            probe.run_probe(raw, curated, split_dir, dest)


def test_saved_dpp4_predictions_match_heldout_labels_and_reported_errors():
    """Check the actual coursework artifacts independently of the metric helper."""
    out = probe.RESULTS / "leakage_probe" / "CHEMBL284"
    report = json.loads((out / "probe.json").read_text())
    split_dir = probe.RESULTS / "splits" / "CHEMBL284"
    assert report["split_audit_sha256"] == retrieve._sha256(split_dir / "audit.json")
    for filename, metadata in report["files"].items():
        assert retrieve._sha256(out / filename) == metadata["sha256"]
    for name in ("random", "scaffold"):
        predictions = probe._read_csv(out / f"{name}_predictions.csv")
        assignment = probe._read_csv(split_dir / f"{name}_assignments.csv")
        train_keys = {r["inchikey"] for r in assignment if r["partition"] == "train"}
        test_labels = {r["inchikey"]: float(r["pIC50_median"]) for r in assignment
                       if r["partition"] == "test"}
        assert {r["inchikey"] for r in predictions} == set(test_labels)
        assert len(predictions) == len(test_labels)
        for row in predictions:
            assert float(row["observed_pIC50"]) == test_labels[row["inchikey"]]
            assert row["nearest_train_inchikey"] in train_keys
        errors = [float(r["predicted_pIC50"]) - test_labels[r["inchikey"]] for r in predictions]
        metrics = report["splits"][name]["metrics"]
        assert metrics["rmse"] == pytest.approx(math.sqrt(math.fsum(e * e for e in errors) / len(errors)))
        assert metrics["mae"] == pytest.approx(math.fsum(abs(e) for e in errors) / len(errors))
