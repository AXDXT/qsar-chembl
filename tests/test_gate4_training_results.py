"""Audit real training-search caches without refitting or held-out scoring."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import shutil

import numpy as np
import pytest
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

from src import gate4_training as training
from src import gate4_combined as combined
from src import nested_cv as cv
from src.retrieve import _sha256


def _digest(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
    return hashlib.sha256(encoded.encode()).hexdigest()


def _metrics(predictions):
    observed = [p["observed_pIC50"] for p in predictions]
    predicted = [p["predicted_pIC50"] for p in predictions]
    return {"rmse": float(np.sqrt(mean_squared_error(observed, predicted))),
            "mae": mean_absolute_error(observed, predicted),
            "r2": r2_score(observed, predicted)}


def _combined_fit_plan(manifest):
    """Expected fitting/validation indices, in resumable execution order."""
    result = {}
    for number, outer in enumerate(manifest["outer_plan"]):
        for candidate in range(len(manifest["candidates"])):
            for fold_number, fold in enumerate(manifest["inner_plan"][number]):
                assert set(fold["train"] + fold["validation"]).isdisjoint(outer["validation"])
                result[f"outer{number}-inner-candidate{candidate}-fold{fold_number}"] = (
                    fold, f"outer{number}-inner-fold{fold_number}")
        result[f"outer{number}-evaluation"] = (outer, f"outer{number}-evaluation")
    for candidate in range(len(manifest["candidates"])):
        for number, fold in enumerate(manifest["final_plan"]):
            result[f"final-inner-candidate{candidate}-fold{number}"] = (
                fold, f"final-inner-fold{number}")
    return result


@pytest.fixture(scope="module")
def saved_controls():
    target = "CHEMBL284"
    inputs = (training.DATA_RAW / target, training.DATA_CURATED / target,
              training.RESULTS / "splits" / target, training.DATA_FEATURES / target)
    return inputs, training.load_training_inputs(*inputs)


@pytest.mark.parametrize("split", ["scaffold", "random"])
def test_logged_control_fits_match_grouped_plans_and_report(saved_controls, split):
    _, prepared = saved_controls
    directory = training.RESULTS / "gate4_training" / "CHEMBL284" / split / "descriptors"
    manifest = json.loads((directory / "manifest.json").read_text())
    identity = _digest(manifest)
    rows = prepared["splits"][split]["rows"]
    expected_rows = [
        {key: row[key] for key in ("inchikey", "scaffold", "pIC50_median", "partition")}
        for row in rows]
    assert manifest["rows"] == expected_rows
    assert manifest["source_checksums"] == prepared["splits"][split]["source_checksums"]
    for filename, digest in manifest["code_sha256"].items():
        assert digest == _sha256(Path(training.__file__).with_name(filename))
    assert manifest["folds"] == 5
    assert manifest["candidates"] == list(cv.ParameterGrid(cv.DEFAULT_GRID))
    assert manifest["estimator"] == cv.RandomForestRegressor(**cv.RF_SPEC).get_params()
    records = {}
    for path in directory.glob("*.json"):
        if path.name == "manifest.json" or path.name.startswith("._"):
            continue
        record = json.loads(path.read_text())
        assert record["name"] == path.stem
        assert record["identity"] == identity
        assert record["sha256"] == _digest(record["result"])
        records[path.stem] = record["result"]
    report = records.pop("report")
    assert report["identity_sha256"] == identity
    assert report["heldout_evaluations"] == 0
    assert report["grid"] == manifest["candidates"]
    expected_files, oof, selections, elapsed = set(), [], Counter(), []

    def check_fit(name, fold):
        expected_files.add(name)
        fit = records[name]
        train, validation = fold["train"], fold["validation"]
        assert set(train).isdisjoint(validation)
        assert {rows[i]["scaffold"] for i in train}.isdisjoint(
            rows[i]["scaffold"] for i in validation)
        predictions = fit["predictions"]
        assert [p["inchikey"] for p in predictions] == [rows[i]["inchikey"] for i in validation]
        for p, i in zip(predictions, validation, strict=True):
            assert p["scaffold"] == rows[i]["scaffold"]
            assert p["observed_pIC50"] == float(rows[i]["pIC50_median"])
        assert fit["metrics"] == pytest.approx(_metrics(predictions))
        assert fit["n_jobs"] == 4
        assert fit["fit_predict_seconds"] > 0
        elapsed.append(fit["fit_predict_seconds"])
        return fit

    def check_selection(selection, plan, prefix):
        assert [c["parameters"] for c in selection["candidates"]] == manifest["candidates"]
        for candidate, score in enumerate(selection["candidates"]):
            fits = [check_fit(f"{prefix}-candidate{candidate}-fold{number}", fold)
                    for number, fold in enumerate(plan)]
            assert score["fold_metrics"] == [fit["metrics"] for fit in fits]
            assert score["fit_predict_seconds"] == [fit["fit_predict_seconds"] for fit in fits]
            assert score["mean_validation_rmse"] == pytest.approx(
                np.mean([fit["metrics"]["rmse"] for fit in fits]))
        winner = min(range(len(selection["candidates"])),
                     key=lambda i: selection["candidates"][i]["mean_validation_rmse"])
        assert selection["selected_candidate"] == winner
        assert selection["selected_parameters"] == manifest["candidates"][winner]

    for number, outer in enumerate(manifest["outer_plan"]):
        assert sorted(outer["train"] + outer["validation"]) == list(range(len(rows)))
        inner_plan = manifest["inner_plan"][number]
        assert len(inner_plan) == 5
        assert Counter(i for fold in inner_plan for i in fold["validation"]) == Counter(
            outer["train"])
        for inner in inner_plan:
            assert sorted(inner["train"] + inner["validation"]) == sorted(outer["train"])
            assert set(inner["train"] + inner["validation"]).isdisjoint(outer["validation"])
        result = report["outer_folds"][number]
        check_selection(result["selection"], inner_plan, f"outer{number}-inner")
        selections[_digest(result["selection"]["selected_parameters"])] += 1
        evaluation = check_fit(f"outer{number}-evaluation", outer)
        assert result["metrics"] == evaluation["metrics"]
        oof.extend({**p, "outer_fold": number} for p in evaluation["predictions"])
    assert Counter(i for fold in manifest["final_plan"] for i in fold["validation"]) == Counter(
        range(len(rows)))
    for fold in manifest["final_plan"]:
        assert sorted(fold["train"] + fold["validation"]) == list(range(len(rows)))
    check_selection(report["final_selection"], manifest["final_plan"], "final-inner")
    assert set(records) == expected_files
    assert report["oof_predictions"] == oof
    assert Counter(p["inchikey"] for p in oof) == Counter(r["inchikey"] for r in rows)
    assert report["oof_metrics"] == pytest.approx(_metrics(oof))
    assert report["oof_bootstrap"] == cv.clustered_intervals(
        oof, repeats=manifest["bootstrap_repeats"])
    stability = report["selection_stability"]
    assert stability["all_outer_folds_agree"] == (len(selections) == 1)
    assert stability["distinct_selections"] == len(selections)
    assert {_digest(c["parameters"]): c["outer_folds"] for c in stability["counts"]} == selections
    for metric, spread in report["outer_fold_spread"].items():
        values = [r["metrics"][metric] for r in report["outer_folds"]]
        assert spread["values"] == values
        assert spread["mean"] == pytest.approx(np.mean(values))
        assert spread["sample_sd"] == pytest.approx(np.std(values, ddof=1))
    assert report["timing"]["cumulative_fit_predict_seconds"] == pytest.approx(sum(elapsed))
    assert report["timing"]["completion_invocation_seconds"] > 0
    assert report["timing"]["n_jobs"] == 4
    assert report["timing"]["search_workers"] == 1


def test_production_controls_resume_without_fitting(saved_controls, monkeypatch):
    inputs, _ = saved_controls

    def forbidden(*args, **kwargs):
        pytest.fail("A completed real-data descriptor control retrained")

    monkeypatch.setattr(cv, "_fit", forbidden)
    monkeypatch.setattr(cv.RandomForestRegressor, "fit", forbidden)
    reports = training.run_descriptor_control(
        *inputs, training.RESULTS / "gate4_training" / "CHEMBL284", n_jobs=1)
    assert set(reports) == {"scaffold", "random"}
    assert all(report["heldout_evaluations"] == 0 for report in reports.values())


@pytest.mark.parametrize("split", ["scaffold", "random"])
def test_production_combined_resume_reuses_completed_fits_before_first_missing(
        saved_controls, monkeypatch, tmp_path, split):
    """Probe resume read-only, stopping before an incomplete fit can do any work."""
    inputs, _ = saved_controls
    source = training.RESULTS / "gate4_training" / "CHEMBL284" / split
    cache_root = (tmp_path / "snapshot").resolve()
    # Another room participant may append fits. Work on immutable snapshots,
    # never instrument or write to the active production cache.
    for block in ("combined", "fingerprint_plan"):
        destination = cache_root / split / block
        destination.mkdir(parents=True)
        for path in list((source / block).glob("*.json")):
            if not path.name.startswith("._"):
                shutil.copyfile(path, destination / path.name)
    directory = cache_root / split / "combined"
    manifest = json.loads((directory / "manifest.json").read_text())
    cached = {p.stem for p in directory.glob("*.json")
              if p.stem not in {"manifest", "report"} and not p.name.startswith("._")}
    visited = []
    actual_get = cv._Cache.get

    class FirstMissingFit(Exception):
        pass

    def stop_before_missing(cache, name, compute):
        if cache.directory == directory and name != "report":
            if not (directory / f"{name}.json").exists():
                raise FirstMissingFit(name)
            visited.append(name)
        return actual_get(cache, name, compute)

    def forbidden(*args, **kwargs):
        pytest.fail("A completed real-data combined search repeated fitting or featurisation")

    monkeypatch.setattr(cv, "_fit", forbidden)
    monkeypatch.setattr(cv.RandomForestRegressor, "fit", forbidden)
    monkeypatch.setattr(combined, "_combined_values", forbidden)
    monkeypatch.setattr(cv._Cache, "get", stop_before_missing)

    def resume():
        return combined.run_combined_search(
            *inputs, cache_root, splits=(split,), n_jobs=1)

    if (directory / "report.json").exists():
        reports = resume()
        saved = json.loads((directory / "report.json").read_text())["result"]
        assert reports == {split: saved}
        assert saved["heldout_evaluations"] == 0
        assert visited == []
    else:
        with pytest.raises(FirstMissingFit):
            resume()
        prefix = []
        for name in _combined_fit_plan(manifest):
            if name not in cached:
                break
            prefix.append(name)
        assert visited == prefix
        assert len(visited) > 0


@pytest.mark.parametrize("split", ["scaffold", "random"])
def test_completed_combined_fits_are_valid_even_before_search_completion(saved_controls, split):
    """An interrupted search's work is checked, without pretending it is complete."""
    inputs, prepared = saved_controls
    directory = training.RESULTS / "gate4_training" / "CHEMBL284" / split / "combined"
    manifest = json.loads((directory / "manifest.json").read_text())
    plan_dir = directory.parent / "fingerprint_plan"
    plan_manifest = json.loads((plan_dir / "manifest.json").read_text())
    plan_record = json.loads((plan_dir / "report.json").read_text())
    plan = plan_record["result"]
    assert plan_record["identity"] == _digest(plan_manifest) == plan["identity_sha256"]
    assert plan_record["sha256"] == _digest(plan)
    assert plan["heldout_evaluations"] == 0
    control = json.loads((directory.parent / "descriptors/manifest.json").read_text())
    for key in ("outer_plan", "inner_plan", "final_plan"):
        assert manifest[key] == plan_manifest[key] == control[key]
    for key in ("seed", "bootstrap_repeats"):
        assert manifest[key] == control[key]
    rows = prepared["splits"][split]["rows"]
    structures = {r["inchikey"]: r["canonical_smiles"]
                  for r in training._read_csv(inputs[1] / "compounds.csv")}
    assert manifest["rows"] == [{**r, "canonical_smiles": structures[r["inchikey"]]} for r in rows]
    assert manifest["descriptor_sha256"] == hashlib.sha256(
        prepared["splits"][split]["descriptors"].tobytes(order="C")).hexdigest()
    assert all(r["partition"] == "train" for r in manifest["rows"])
    assert manifest["source_checksums"] == prepared["splits"][split]["source_checksums"]
    assert manifest["fingerprint_plan_sha256"] == _digest(plan)
    assert manifest["candidates"] == list(cv.ParameterGrid(cv.DEFAULT_GRID))
    assert manifest["estimator"] == cv.RandomForestRegressor(**cv.RF_SPEC).get_params()
    for filename, digest in manifest["code_sha256"].items():
        assert digest == _sha256(Path(training.__file__).with_name(filename))
    expected = _combined_fit_plan(manifest)
    checked = 0
    for path in directory.glob("*.json"):
        if path.stem in {"manifest", "report"} or path.name.startswith("._"):
            continue
        record = json.loads(path.read_text())
        assert record["name"] == path.stem
        assert record["identity"] == _digest(manifest)
        assert record["sha256"] == _digest(record["result"])
        fold, decision_name = expected[path.stem]
        fit = record["result"]
        assert set(fold["train"]).isdisjoint(fold["validation"])
        assert {rows[i]["scaffold"] for i in fold["train"]}.isdisjoint(
            rows[i]["scaffold"] for i in fold["validation"])
        assert [p["inchikey"] for p in fit["predictions"]] == [
            rows[i]["inchikey"] for i in fold["validation"]]
        for p, i in zip(fit["predictions"], fold["validation"], strict=True):
            assert p["scaffold"] == rows[i]["scaffold"]
            assert p["observed_pIC50"] == float(rows[i]["pIC50_median"])
        assert fit["metrics"] == pytest.approx(_metrics(fit["predictions"]))
        assert fit["fingerprint_decision"] == decision_name
        decision = plan["fit_decisions"][decision_name]
        assert decision["train_indices"] == fold["train"]
        assert fit["fpSize"] == decision["selected_fpSize"]
        assert fit["fit_predict_seconds"] > 0
        assert isinstance(fit["n_jobs"], int) and not isinstance(fit["n_jobs"], bool)
        assert fit["n_jobs"] >= 1
        checked += 1
    assert checked > 0


@pytest.mark.parametrize("split", ["scaffold", "random"])
def test_completed_combined_report_is_derived_from_all_cached_fits(saved_controls, split):
    """Rebuild selections and uncertainty independently of the search orchestrator."""
    _, prepared = saved_controls
    directory = training.RESULTS / "gate4_training" / "CHEMBL284" / split / "combined"
    manifest = json.loads((directory / "manifest.json").read_text())
    records = {}
    for path in directory.glob("*.json"):
        if path.name == "manifest.json" or path.name.startswith("._"):
            continue
        wrapper = json.loads(path.read_text())
        assert wrapper["name"] == path.stem
        assert wrapper["identity"] == _digest(manifest)
        assert wrapper["sha256"] == _digest(wrapper["result"])
        records[path.stem] = wrapper["result"]
    expected_fits = set(_combined_fit_plan(manifest))
    if "report" not in records and set(records) < expected_fits:
        pytest.skip(f"{split} combined report pending ({len(records)}/{len(expected_fits)} fits)")
    report = records.pop("report")
    assert set(records) == expected_fits
    assert report["identity_sha256"] == _digest(manifest)
    assert report["heldout_evaluations"] == 0
    assert report["feature_block"] == "descriptors+Morgan"
    assert report["grid"] == manifest["candidates"]
    plan = json.loads((directory.parent / "fingerprint_plan/report.json").read_text())["result"]
    assert report["fingerprint_plan_identity"] == plan["identity_sha256"]
    assert report["fingerprint_plan_sha256"] == _digest(plan)
    assert report["final_fpSize"] == plan["fit_decisions"]["final-fit"]["selected_fpSize"]

    def check_selection(selection, fold_plan, prefix):
        assert len(selection["candidates"]) == len(manifest["candidates"])
        scores = []
        for candidate, parameters in enumerate(manifest["candidates"]):
            fits = [records[f"{prefix}-candidate{candidate}-fold{number}"]
                    for number in range(len(fold_plan))]
            score = selection["candidates"][candidate]
            assert score["parameters"] == parameters
            metrics = [_metrics(fit["predictions"]) for fit in fits]
            for logged, actual in zip(score["fold_metrics"], metrics, strict=True):
                assert logged == pytest.approx(actual)
            mean_rmse = float(np.mean([m["rmse"] for m in metrics]))
            assert score["mean_validation_rmse"] == pytest.approx(mean_rmse)
            assert score["fold_fpSizes"] == [f["fpSize"] for f in fits]
            assert score["fit_predict_seconds"] == [f["fit_predict_seconds"] for f in fits]
            scores.append(mean_rmse)
        winner = min(range(len(scores)), key=lambda i: scores[i])
        assert selection["selected_candidate"] == winner
        assert selection["selected_parameters"] == manifest["candidates"][winner]

    assert len(report["outer_folds"]) == len(manifest["outer_plan"]) == 5
    oof, selections = [], Counter()
    for number, result in enumerate(report["outer_folds"]):
        assert result["outer_fold"] == number
        check_selection(result["selection"], manifest["inner_plan"][number],
                        f"outer{number}-inner")
        selections[_digest(result["selection"]["selected_parameters"])] += 1
        evaluation = records[f"outer{number}-evaluation"]
        assert result["metrics"] == evaluation["metrics"]
        assert result["fpSize"] == evaluation["fpSize"]
        assert result["fit_predict_seconds"] == evaluation["fit_predict_seconds"]
        oof.extend({**p, "outer_fold": number} for p in evaluation["predictions"])
    check_selection(report["final_selection"], manifest["final_plan"], "final-inner")
    assert report["oof_predictions"] == oof
    rows = prepared["splits"][split]["rows"]
    assert Counter(p["inchikey"] for p in oof) == Counter(r["inchikey"] for r in rows)
    assert report["oof_metrics"] == pytest.approx(_metrics(oof))
    assert report["oof_bootstrap"] == cv.clustered_intervals(
        oof, repeats=manifest["bootstrap_repeats"])
    stability = report["selection_stability"]
    assert stability["all_outer_folds_agree"] == (len(selections) == 1)
    assert stability["distinct_selections"] == len(selections)
    assert {_digest(c["parameters"]): c["outer_folds"] for c in stability["counts"]} == selections
    for metric, spread in report["outer_fold_spread"].items():
        values = [r["metrics"][metric] for r in report["outer_folds"]]
        assert spread["values"] == values
        assert spread["mean"] == pytest.approx(np.mean(values))
        assert spread["sample_sd"] == pytest.approx(np.std(values, ddof=1))
    assert report["timing"]["cumulative_fit_predict_seconds"] == pytest.approx(
        sum(r["fit_predict_seconds"] for r in records.values()))
    assert report["timing"]["completion_invocation_seconds"] > 0
    assert report["timing"]["search_workers"] == 1
    assert report["timing"]["n_jobs"] >= 1
