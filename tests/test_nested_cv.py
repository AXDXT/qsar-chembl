"""Training-only nested-CV isolation and resumable-cache regression guards."""
import copy
import json

import numpy as np
import pytest
from sklearn.model_selection import KFold

from src import nested_cv as cv


@pytest.fixture(autouse=True)
def small_forest(monkeypatch):
    # Exercise real sklearn fitting with fewer trees only in these toy tests.
    monkeypatch.setattr(cv, "RF_SPEC", {**cv.RF_SPEC, "n_estimators": 3})


def _toy():
    rows = [{"inchikey": f"key{i}", "scaffold": f"core{i // 2}" if i > 1 else "",
             "pIC50_median": str(4 + i / 5), "partition": "train"} for i in range(20)]
    x = np.column_stack([np.arange(20), np.arange(20) % 3]).astype(float)
    return rows, x


def _run(tmp_path, rows=None, x=None, **kwargs):
    default_rows, default_x = _toy()
    settings = {"feature_spec": {"block": "toy", "columns": ["identity", "remainder"]},
                "source_checksums": {"assignments": "fixture-sha256"},
                "grid": {"max_features": [.5, 1.], "min_samples_leaf": [1]},
                "bootstrap_repeats": 20}
    settings.update(kwargs)
    return cv.run_nested_cv(default_rows if rows is None else rows,
                            default_x if x is None else x, tmp_path, **settings)


def test_five_by_five_grouping_every_fit_is_isolated_and_outer_predictions_are_oof(tmp_path, monkeypatch):
    rows, x = _toy()
    actual = cv._fit
    seen, active = [], {}
    forest = cv.RandomForestRegressor

    class AuditedForest(forest):
        def fit(self, values, labels):
            np.testing.assert_array_equal(values, x[active["fold"]["train"]])
            np.testing.assert_array_equal(labels, [float(rows[i]["pIC50_median"])
                                                  for i in active["fold"]["train"]])
            return super().fit(values, labels)

        def predict(self, values):
            np.testing.assert_array_equal(values, x[active["fold"]["validation"]])
            return super().predict(values)

    monkeypatch.setattr(cv, "RandomForestRegressor", AuditedForest)

    def spy(rows, x, y, fold, parameters):
        train, validation = fold["train"], fold["validation"]
        assert set(train).isdisjoint(validation)
        assert {rows[i]["scaffold"] for i in train}.isdisjoint(
            rows[i]["scaffold"] for i in validation)
        seen.append(copy.deepcopy(fold))
        active["fold"] = fold
        return actual(rows, x, y, fold, parameters)

    monkeypatch.setattr(cv, "_fit", spy)
    report = _run(tmp_path / "cache", rows, x)
    manifest = json.loads((tmp_path / "cache" / "manifest.json").read_text())
    assert len(report["outer_folds"]) == 5
    assert len(seen) == 5 * (2 * 5 + 1) + 2 * 5
    for outer_number, outer in enumerate(manifest["outer_plan"]):
        assert len(manifest["inner_plan"][outer_number]) == 5
        for inner in manifest["inner_plan"][outer_number]:
            assert set(inner["train"] + inner["validation"]) == set(outer["train"])
            assert set(inner["train"] + inner["validation"]).isdisjoint(outer["validation"])
        # First ten calls tune inside this outer training pool, then one scores its validation.
        assert seen[outer_number * 11: outer_number * 11 + 10] == (
            manifest["inner_plan"][outer_number] * 2)
        assert seen[outer_number * 11 + 10] == outer
        result = report["outer_folds"][outer_number]
        candidates = result["selection"]["candidates"]
        best = min(range(len(candidates)), key=lambda i: candidates[i]["mean_validation_rmse"])
        assert result["selection"]["selected_parameters"] == candidates[best]["parameters"]
    predictions = report["oof_predictions"]
    assert len(predictions) == len(rows)
    assert {p["inchikey"] for p in predictions} == {r["inchikey"] for r in rows}
    for p in predictions:
        i = next(i for i, row in enumerate(rows) if row["inchikey"] == p["inchikey"])
        assert i in manifest["outer_plan"][p["outer_fold"]]["validation"]
        assert i not in manifest["outer_plan"][p["outer_fold"]]["train"]
    assert report["oof_metrics"] == cv.point_metrics(predictions)
    assert report["oof_bootstrap"] == cv.clustered_intervals(predictions, repeats=20)
    assert report["heldout_evaluations"] == 0
    assert report["timing"]["n_jobs"] == report["timing"]["search_workers"] == 1
    assert sum(r["outer_folds"] for r in report["selection_stability"]["counts"]) == 5
    assert report["timing"]["cumulative_fit_predict_seconds"] > 0


def test_changing_outer_validation_labels_cannot_change_its_inner_selection(tmp_path):
    rows, x = _toy()
    original = _run(tmp_path / "original", rows, x)
    plan = json.loads((tmp_path / "original" / "manifest.json").read_text())["outer_plan"]
    changed = copy.deepcopy(rows)
    for i in plan[0]["validation"]:
        changed[i]["pIC50_median"] = "10000"
    altered = _run(tmp_path / "altered", changed, x)
    before, after = original["outer_folds"][0], altered["outer_folds"][0]
    assert before["selection"]["selected_parameters"] == after["selection"]["selected_parameters"]
    assert [r["fold_metrics"] for r in before["selection"]["candidates"]] == [
        r["fold_metrics"] for r in after["selection"]["candidates"]]
    assert [p["predicted_pIC50"] for p in original["oof_predictions"] if p["outer_fold"] == 0] == [
        p["predicted_pIC50"] for p in altered["oof_predictions"] if p["outer_fold"] == 0]
    assert before["metrics"] != after["metrics"]


def test_interrupted_search_resumes_each_completed_fit_and_complete_rerun_never_fits(tmp_path, monkeypatch):
    actual, calls = cv._fit, []

    def interrupted(*args):
        calls.append(1)
        if len(calls) == 4:
            raise KeyboardInterrupt("simulated turn cutoff")
        return actual(*args)

    monkeypatch.setattr(cv, "_fit", interrupted)
    with pytest.raises(KeyboardInterrupt):
        _run(tmp_path / "cache")
    assert len(list((tmp_path / "cache").glob("outer*-fold*.json"))) == 3
    assert not (tmp_path / "cache" / "report.json").exists()
    resumed = _run(tmp_path / "cache")
    assert len(calls) == 66  # 65 completed fits plus the interrupted one, no repeats

    def forbidden(*args):
        pytest.fail("A completed cached run must never refit")

    monkeypatch.setattr(cv, "_fit", forbidden)
    assert resumed == _run(tmp_path / "cache")
    assert not list((tmp_path / "cache").glob(".pending-*"))


@pytest.mark.parametrize("changed", ["features", "labels", "groups", "keys", "row_order", "source",
                                      "feature_spec", "grid", "trees", "code", "bootstrap", "version"])
def test_cache_never_silently_reuses_changed_inputs_or_configuration(tmp_path, monkeypatch, changed):
    _run(tmp_path / "cache")
    rows, x = _toy()
    settings = {}
    if changed == "features":
        x[0, 0] += 1
    elif changed == "labels":
        rows[0]["pIC50_median"] = "99"
    elif changed == "groups":
        rows[0]["scaffold"] = "new-scaffold"
    elif changed == "keys":
        rows[0]["inchikey"] = "new-key"
    elif changed == "row_order":
        rows, x = rows[::-1], x[::-1]
    elif changed == "source":
        settings["source_checksums"] = {"assignments": "changed-sha256"}
    elif changed == "feature_spec":
        settings["feature_spec"] = {"block": "different-feature-definition"}
    elif changed == "grid":
        settings["grid"] = {"max_features": [.5], "min_samples_leaf": [1]}
    elif changed == "trees":
        monkeypatch.setattr(cv, "RF_SPEC", {**cv.RF_SPEC, "n_estimators": 4})
    elif changed == "code":
        monkeypatch.setattr(cv, "_sha256", lambda path: "changed-code-sha256")
    elif changed == "bootstrap":
        settings["bootstrap_repeats"] = 21
    else:
        monkeypatch.setattr(cv.sklearn, "__version__", "changed-version")
    with pytest.raises(ValueError, match="cache inputs/settings/code changed"):
        _run(tmp_path / "cache", rows, x, **settings)


def test_tampered_cache_result_is_rejected(tmp_path):
    _run(tmp_path / "cache")
    path = tmp_path / "cache" / "report.json"
    record = json.loads(path.read_text())
    record["result"]["oof_metrics"]["rmse"] = 0
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="checksum/identity mismatch"):
        _run(tmp_path / "cache")


def test_random_inner_folds_are_refused_before_fitting_or_creating_cache(tmp_path, monkeypatch):
    class Ungrouped:
        def __init__(self, n_splits):
            self.splitter = KFold(n_splits=n_splits, shuffle=True, random_state=cv.SEED)

        def split(self, indices, groups):
            return self.splitter.split(indices)

    actual = cv.GroupKFold
    calls = []

    def only_inner_is_wrong(n_splits):
        calls.append(1)
        return actual(n_splits=n_splits) if len(calls) == 1 else Ungrouped(n_splits)

    monkeypatch.setattr(cv, "GroupKFold", only_inner_is_wrong)
    with pytest.raises(AssertionError, match="Scaffold leakage"):
        _run(tmp_path / "cache")
    assert not (tmp_path / "cache").exists()


@pytest.mark.parametrize("invalid", ["heldout", "missing_partition", "duplicate", "nonfinite",
                                      "shape", "too_few_groups", "grid_axis"])
def test_invalid_inputs_abort_before_creating_cache(tmp_path, invalid):
    rows, x = _toy()
    settings = {}
    if invalid == "heldout":
        rows[0]["partition"] = "test"
    elif invalid == "missing_partition":
        del rows[0]["partition"]
    elif invalid == "duplicate":
        rows[0]["inchikey"] = rows[1]["inchikey"]
    elif invalid == "nonfinite":
        x[0, 0] = np.nan
    elif invalid == "shape":
        x = x[:-1]
    elif invalid == "too_few_groups":
        for row in rows:
            row["scaffold"] = "one-group"
    else:
        settings["grid"] = {"n_estimators": [1]}
    with pytest.raises(ValueError):
        _run(tmp_path / "cache", rows, x, **settings)
    assert not (tmp_path / "cache").exists()


def test_default_grid_fixes_baseline_tree_count_and_only_searches_two_axes():
    assert len(list(cv.ParameterGrid(cv.DEFAULT_GRID))) == 6
    assert set(cv.DEFAULT_GRID) == {"max_features", "min_samples_leaf"}
    from src.baselines import RF_SPEC
    assert RF_SPEC["n_estimators"] == 200
    assert RF_SPEC["n_jobs"] == 1 and RF_SPEC["random_state"] == cv.SEED
