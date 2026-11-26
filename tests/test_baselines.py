"""Offline baseline guards: train-only fitting, saved splits and held-out scoring."""
import json

import numpy as np
import pytest

from src import baselines, features, retrieve, splits
from test_features import _inputs, _read_csv


def _toy():
    rows = [{"inchikey": str(i), "scaffold": f"core{i}", "pIC50_median": str(y)}
            for i, y in enumerate([4., 8., 100., -100.])]
    x = np.arange(4 * len(features.DESCRIPTOR_NAMES), dtype=float).reshape(4, -1)
    return rows, x, ["train", "train", "test", "test"]


def test_fixed_specs_are_logged_in_full_and_mean_is_a_real_training_mean():
    models = baselines.make_models()
    assert models["mean"].get_params()["strategy"] == "mean"
    params = models["descriptor_rf"].get_params()
    assert all(params[name] == value for name, value in baselines.RF_SPEC.items())
    assert params["random_state"] == baselines.SEED
    assert params["n_jobs"] == 1 and params["n_estimators"] == 200
    assert params["max_features"] == 1.0 and params["min_samples_leaf"] == 1


def test_held_out_labels_cannot_affect_either_fitted_model_or_its_predictions():
    rows, x, assignment = _toy()
    original = baselines.predict_partition(rows, x, assignment, True)
    assert {p["predicted_pIC50"] for p in original["mean"]} == {6.}
    rows[2]["pIC50_median"], rows[3]["pIC50_median"] = "9999", "9999"
    changed = baselines.predict_partition(rows, x, assignment, True)
    for name in original:
        assert [p["predicted_pIC50"] for p in changed[name]] == [
            p["predicted_pIC50"] for p in original[name]]
        assert [p["partition"] for p in changed[name]] == assignment


def test_no_transform_is_fitted_on_test_and_estimators_receive_only_raw_train_descriptors(monkeypatch):
    rows, x, assignment = _toy()
    actual_models = baselines.make_models()
    seen = {}

    class Spy:
        def __init__(self, name, model):
            self.name, self.model = name, model

        def fit(self, values, labels):
            seen[self.name] = (values.copy(), labels.copy())
            self.model.fit(values, labels)

        def predict(self, values):
            # An unnoticed scaling/imputation/filtering pass would change this array.
            np.testing.assert_array_equal(values, x)
            return self.model.predict(values)

    monkeypatch.setattr(baselines, "make_models", lambda: {
        name: Spy(name, model) for name, model in actual_models.items()})
    baselines.predict_partition(rows, x, assignment)
    for values, labels in seen.values():
        np.testing.assert_array_equal(values, x[:2])
        np.testing.assert_array_equal(labels, [4., 8.])
    assert set(seen) == set(actual_models)


def test_changing_test_features_cannot_change_training_predictions():
    rows, x, assignment = _toy()
    original = baselines.predict_partition(rows, x, assignment)
    x[2:] = -1e9
    changed = baselines.predict_partition(rows, x, assignment)
    for model in original:
        assert [p["predicted_pIC50"] for p in original[model][:2]] == [
            p["predicted_pIC50"] for p in changed[model][:2]]


@pytest.mark.parametrize("corrupt", ["shape", "nonfinite", "labels", "identity", "scaffold", "partition"])
def test_invalid_training_or_evaluation_inputs_fail(corrupt):
    rows, x, assignment = _toy()
    if corrupt == "shape":
        x = x[:, :-1]
    elif corrupt == "nonfinite":
        x[2, 0] = np.inf
    elif corrupt == "labels":
        rows[2]["pIC50_median"] = "nan"
    elif corrupt == "identity":
        rows[2]["inchikey"] = rows[0]["inchikey"]
    elif corrupt == "scaffold":
        rows[2]["scaffold"] = rows[0]["scaffold"]
    else:
        assignment = ["train"] * 4
    with pytest.raises((ValueError, AssertionError)):
        baselines.predict_partition(rows, x, assignment, True)


def _cached_inputs(tmp_path):
    raw, curated, split_dir = _inputs(tmp_path)
    feature_dir = tmp_path / "features"
    features.build_features(raw, curated, split_dir, feature_dir)
    return raw, curated, split_dir, feature_dir


def test_pipeline_reads_saved_splits_without_regenerating_and_is_byte_reproducible(tmp_path, monkeypatch):
    inputs = _cached_inputs(tmp_path)

    def forbidden(*args, **kwargs):
        pytest.fail("Baseline code must never regenerate a split")

    monkeypatch.setattr(splits, "make_partitions", forbidden)
    before = {p: p.read_bytes() for directory in inputs for p in directory.iterdir()
              if not p.name.startswith("._")}
    out = tmp_path / "baselines"
    report = baselines.run_baselines(*inputs, out)
    first = {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert report == baselines.run_baselines(*inputs, out)
    assert first == {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert all(p.read_bytes() == contents for p, contents in before.items())
    assert report["features"] == list(features.DESCRIPTOR_NAMES)
    for filename, meta in report["files"].items():
        assert retrieve._sha256(out / filename) == meta["sha256"]
    _assert_scores_match_partitions(report, out, inputs[2])


def _independent_metrics(rows):
    y = np.array([float(r["observed_pIC50"]) for r in rows])
    p = np.array([float(r["predicted_pIC50"]) for r in rows])
    sse = np.square(p - y).sum()
    denominator = np.square(y - y.mean()).sum()
    return {"r2": 1 - sse / denominator if denominator > 0 else None,
            "rmse": np.sqrt(sse / len(y)), "mae": np.abs(p - y).mean()}


def _assert_scores_match_partitions(report, out, split_dir):
    for split, result in report["splits"].items():
        saved = _read_csv(split_dir / f"{split}_assignments.csv")
        for model, scores in result["models"].items():
            values = _read_csv(out / f"{split}_{model}_predictions.csv")
            assert [r["inchikey"] for r in values] == [r["inchikey"] for r in saved]
            assert [r["partition"] for r in values] == [r["partition"] for r in saved]
            for side in ("train", "test"):
                selected = [r for r in values if r["partition"] == side]
                assert len(selected) == result[f"{side}_compounds"]
                assert {r["inchikey"] for r in selected} == {
                    r["inchikey"] for r in saved if r["partition"] == side}
                assert scores[side] == pytest.approx(_independent_metrics(selected))
            if model == "mean":
                mean = np.mean([float(r["pIC50_median"]) for r in saved if r["partition"] == "train"])
                assert {float(r["predicted_pIC50"]) for r in values} == {mean}


@pytest.mark.parametrize("corrupt", ["raw", "curated", "cache", "split", "regenerated", "labels", "order"])
def test_inconsistent_inputs_fail_before_results_are_written(tmp_path, corrupt):
    raw, curated, split_dir, feature_dir = _cached_inputs(tmp_path)
    if corrupt in ("raw", "curated", "cache", "split"):
        path = {"raw": raw / retrieve.RAW_TARGET, "curated": curated / "compounds.csv",
                "cache": feature_dir / "descriptors.csv",
                "split": split_dir / "random_assignments.csv"}[corrupt]
        with path.open("a") as fh:
            fh.write("tampered\n")
    else:
        # Even an updated assignment checksum cannot silently replace the cached split.
        path = split_dir / "random_assignments.csv"
        saved = _read_csv(path)
        if corrupt == "regenerated":
            a = next(r for r in saved if r["partition"] == "train")
            b = next(r for r in saved if r["partition"] == "test")
            a["partition"], b["partition"] = "test", "train"
        elif corrupt == "labels":
            saved[0]["pIC50_median"] = "100"
        else:
            saved.reverse()
        splits._write_csv(path, saved, list(saved[0]))
        audit = json.loads((split_dir / "audit.json").read_text())
        audit["files"][path.name]["sha256"] = retrieve._sha256(path)
        (split_dir / "audit.json").write_text(json.dumps(audit))
        if corrupt != "regenerated":
            # Pass provenance checks deliberately to exercise the content alignment guard.
            manifest = json.loads((feature_dir / retrieve.MANIFEST).read_text())
            manifest["source"]["split_audit_sha256"] = retrieve._sha256(split_dir / "audit.json")
            (feature_dir / retrieve.MANIFEST).write_text(json.dumps(manifest))
    out = tmp_path / "baselines"
    with pytest.raises(ValueError, match="checksum|provenance|row order|labels"):
        baselines.run_baselines(raw, curated, split_dir, feature_dir, out)
    assert not out.exists()


def test_output_cannot_overwrite_any_input_tree(tmp_path):
    inputs = _cached_inputs(tmp_path)
    for dest in (*inputs, inputs[-1] / "derived", tmp_path):
        with pytest.raises(ValueError, match="separate"):
            baselines.run_baselines(*inputs, dest)


def test_logged_dpp4_metrics_score_the_saved_test_keys_instead_of_training_rows():
    out = baselines.RESULTS / "baselines" / "CHEMBL284"
    report = json.loads((out / "baselines.json").read_text())
    split_dir = baselines.RESULTS / "splits" / "CHEMBL284"
    _assert_scores_match_partitions(report, out, split_dir)
    for split, result in report["splits"].items():
        assert report["source"]["assignment_sha256"][split] == retrieve._sha256(
            split_dir / f"{split}_assignments.csv")
        for model, scores in result["models"].items():
            test = [r for r in _read_csv(out / f"{split}_{model}_predictions.csv")
                    if r["partition"] == "test"]
            for r in test:
                r["observed_pIC50"] = float(r["observed_pIC50"])
                r["predicted_pIC50"] = float(r["predicted_pIC50"])
            assert scores["test_bootstrap"] == baselines.clustered_intervals(test)
            assert scores["train"]["rmse"] != scores["test"]["rmse"]
