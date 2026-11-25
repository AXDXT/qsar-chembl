"""Saved-data provenance and train-only Gate 4 pipeline boundary guards."""
import json

import numpy as np
import pytest

from src import gate4_training as training
from src import nested_cv, retrieve, splits
from test_baselines import _cached_inputs


def test_adapter_exposes_only_saved_training_rows_and_aligned_descriptors(tmp_path, monkeypatch):
    inputs = _cached_inputs(tmp_path)

    def forbidden(*args, **kwargs):
        pytest.fail("Training adapter must not regenerate splits or evaluate held-out rows")

    monkeypatch.setattr(splits, "make_partitions", forbidden)
    monkeypatch.setattr(training, "run_nested_cv", forbidden)
    prepared = training.load_training_inputs(*inputs)
    keys, descriptors, _ = training.load_features(inputs[-1])
    for split, data in prepared["splits"].items():
        saved = training._read_csv(inputs[2] / f"{split}_assignments.csv")
        expected = [r for r in saved if r["partition"] == "train"]
        assert data["rows"] == expected
        assert {r["inchikey"] for r in data["rows"]}.isdisjoint(
            r["inchikey"] for r in saved if r["partition"] == "test")
        np.testing.assert_array_equal(data["descriptors"], descriptors[
            [keys.index(r["inchikey"]) for r in expected]])
        assert data["source_checksums"]["assignments"] == retrieve._sha256(
            inputs[2] / f"{split}_assignments.csv")


@pytest.mark.parametrize("which", [0, 1, 2, 3])
def test_tampered_upstream_payload_fails_before_any_search(tmp_path, monkeypatch, which):
    inputs = _cached_inputs(tmp_path)
    names = [retrieve.RAW_TARGET, "compounds.csv", "scaffold_assignments.csv", "descriptors.csv"]
    path = inputs[which] / names[which]
    with path.open("a") as stream:
        stream.write("tampered\n")
    monkeypatch.setattr(training, "run_nested_cv", lambda *a, **k: pytest.fail("Search started"))
    cache = tmp_path / "cache"
    with pytest.raises(ValueError, match="checksum|provenance"):
        training.run_descriptor_control(*inputs, cache)
    assert not cache.exists()


def _rewrite_assignment(inputs, split, mutate):
    path = inputs[2] / f"{split}_assignments.csv"
    rows = training._read_csv(path)
    mutate(rows)
    splits._write_csv(path, rows, list(rows[0]))
    audit_path = inputs[2] / "audit.json"
    audit = json.loads(audit_path.read_text())
    audit["files"][path.name]["sha256"] = retrieve._sha256(path)
    audit_path.write_text(json.dumps(audit))
    feature_path = inputs[3] / retrieve.MANIFEST
    feature = json.loads(feature_path.read_text())
    feature["source"]["split_audit_sha256"] = retrieve._sha256(audit_path)
    feature_path.write_text(json.dumps(feature))


@pytest.mark.parametrize("corrupt", ["labels", "order", "scaffold", "partition", "counts"])
def test_alignment_guards_survive_rechecksummed_assignments(tmp_path, corrupt):
    inputs = _cached_inputs(tmp_path)

    def mutate(rows):
        if corrupt == "labels":
            rows[0]["pIC50_median"] = "10000"
        elif corrupt == "order":
            rows.reverse()
        elif corrupt == "scaffold":
            rows[0]["scaffold"] = "different-saved-core"
        elif corrupt == "partition":
            rows[0]["partition"] = "validation"
        else:
            next(r for r in rows if r["partition"] == "train")["partition"] = "test"

    _rewrite_assignment(inputs, "scaffold", mutate)
    with pytest.raises((ValueError, AssertionError)):
        training.load_training_inputs(*inputs)


def test_source_changes_during_loading_are_refused(tmp_path, monkeypatch):
    inputs = _cached_inputs(tmp_path)
    actual = training.load_features

    def changed(directory):
        result = actual(directory)
        with (inputs[1] / retrieve.MANIFEST).open("a") as stream:
            stream.write(" ")
        return result

    monkeypatch.setattr(training, "load_features", changed)
    with pytest.raises(ValueError, match="changed during"):
        training.load_training_inputs(*inputs)


def test_timestamp_only_rebuild_keeps_search_identity_but_changed_data_does_not(tmp_path):
    inputs = _cached_inputs(tmp_path)
    original = training.load_training_inputs(*inputs)
    path = inputs[1] / retrieve.MANIFEST
    manifest = json.loads(path.read_text())
    manifest["curated_utc"] = "a different run timestamp"
    path.write_text(json.dumps(manifest))
    rebuilt = training.load_training_inputs(*inputs)
    for split in original["splits"]:
        assert original["splits"][split]["source_checksums"] == rebuilt["splits"][split]["source_checksums"]
    # Change a held-out label consistently: it still cannot enter training,
    # but the frozen input identity must change so a cache cannot hide the edit.
    saved = training._read_csv(inputs[2] / "scaffold_assignments.csv")
    key = next(r["inchikey"] for r in saved if r["partition"] == "test")
    compounds_path = inputs[1] / "compounds.csv"
    compounds = training._read_csv(compounds_path)
    next(r for r in compounds if r["inchikey"] == key)["pIC50_median"] = "1000"
    splits._write_csv(compounds_path, compounds, list(compounds[0]))
    manifest["files"]["compounds.csv"]["sha256"] = retrieve._sha256(compounds_path)
    path.write_text(json.dumps(manifest))
    # Updating just curation is insufficient: the existing split/feature provenance fails.
    with pytest.raises(ValueError, match="provenance"):
        training.load_training_inputs(*inputs)


def test_real_nested_control_runs_both_saved_training_pools_and_resumes(tmp_path, monkeypatch):
    inputs = _cached_inputs(tmp_path)
    monkeypatch.setattr(nested_cv, "RF_SPEC", {**nested_cv.RF_SPEC, "n_estimators": 2})
    actual = nested_cv.run_nested_cv
    calls = []

    def small_search(rows, values, cache, **kwargs):
        assert all(r["partition"] == "train" for r in rows)
        assert kwargs["feature_spec"]["columns"] == list(training.DESCRIPTOR_NAMES)
        calls.append((cache, kwargs))
        return actual(rows, values, cache, **kwargs, folds=2, bootstrap_repeats=8,
                      grid={"max_features": [1.0], "min_samples_leaf": [1]})

    monkeypatch.setattr(training, "run_nested_cv", small_search)
    cache = tmp_path / "cache"
    before = {p: p.read_bytes() for directory in inputs for p in directory.iterdir()
              if not p.name.startswith("._")}
    reports = training.run_descriptor_control(*inputs, cache, n_jobs=2)
    assert list(reports) == ["scaffold", "random"]
    assert all(report["heldout_evaluations"] == 0 for report in reports.values())
    assert all(path == cache / kwargs["feature_spec"]["saved_split"] / "descriptors"
               for path, kwargs in calls)
    monkeypatch.setattr(nested_cv, "_fit", lambda *a, **k: pytest.fail("Cache retrained"))
    assert reports == training.run_descriptor_control(*inputs, cache, n_jobs=1)
    assert all(path.read_bytes() == contents for path, contents in before.items())


def test_cache_cannot_overlap_inputs_and_invalid_split_requests_fail(tmp_path):
    inputs = _cached_inputs(tmp_path)
    for cache in (*inputs, inputs[-1] / "search", tmp_path):
        with pytest.raises(ValueError, match="separate"):
            training.run_descriptor_control(*inputs, cache)
    for selected in ((), ("invalid",), ("random", "random")):
        with pytest.raises(ValueError, match="unique saved split"):
            training.run_descriptor_control(*inputs, tmp_path / "cache", splits=selected)


def test_cli_forwards_saved_split_worker_count_and_explicit_cache(tmp_path, monkeypatch):
    seen = {}

    def spy(*args, **kwargs):
        seen["args"], seen["kwargs"] = args, kwargs
        return {"scaffold": {"identity_sha256": "fixture"}}

    monkeypatch.setattr(training, "run_descriptor_control", spy)
    cache = tmp_path / "clean-build"
    assert training.main(["chembl284", "--split", "scaffold", "--n-jobs", "4",
                          "--cache-dir", str(cache)]) == 0
    assert seen["args"][-1] == cache
    assert seen["args"][0] == training.DATA_RAW / "CHEMBL284"
    assert seen["kwargs"] == {"splits": ("scaffold",), "n_jobs": 4}
