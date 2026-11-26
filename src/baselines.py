"""Gate 3 fixed, untuned mean and descriptor-only forest baselines.

Run after src.features: .venv/bin/python -m src.baselines CHEMBL284
Reads the exact saved Gate 2 assignments. No split generation, search or fitted
preprocessing occurs here. Fingerprints are cached for Gate 4, never fitted here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import sklearn
from sklearn.dummy import DummyRegressor
from sklearn.ensemble import RandomForestRegressor

from src.config import DATA_CURATED, DATA_FEATURES, DATA_RAW, RESULTS, SEED
from src.features import DESCRIPTOR_NAMES, _read_csv, load_features
from src.leakage_probe import BOOTSTRAP_REPEATS, clustered_intervals, point_metrics
from src.retrieve import MANIFEST, _sha256, verify_raw
from src.splits import _write_csv, assert_split_integrity

# Declared before the first run, never chosen by scores. A conventional bagged,
# fully grown forest establishes an intentionally simple descriptor-only floor.
# More trees reduce Monte Carlo variation; this count is fixed, not searched.
# One worker makes floating-point prediction accumulation deterministic.
RF_SPEC = dict(n_estimators=200, criterion="squared_error", max_depth=None,
               min_samples_split=2, min_samples_leaf=1, min_weight_fraction_leaf=0.0,
               max_features=1.0, max_leaf_nodes=None, min_impurity_decrease=0.0,
               bootstrap=True, oob_score=False, n_jobs=1, random_state=SEED,
               verbose=0, warm_start=False, ccp_alpha=0.0, max_samples=None)
MODEL_NAMES = {"mean": "Training mean", "descriptor_rf": "Descriptor RF"}


def make_models() -> dict:
    return {"mean": DummyRegressor(strategy="mean"),
            "descriptor_rf": RandomForestRegressor(**RF_SPEC)}


def predict_partition(rows: list[dict], descriptors: np.ndarray, assignment: list[str],
                      scaffold_disjoint: bool = False) -> dict[str, list[dict]]:
    """Fit on TRAIN ONLY, then label train resubstitution and test predictions explicitly.

    No scaler, imputer, variance filter or other fitted transform is used. Cached
    descriptor values are stateless per molecule; non-finite input is refused.
    """
    assert_split_integrity(rows, assignment, scaffold_disjoint)
    descriptors = np.asarray(descriptors, dtype=float)
    if descriptors.shape != (len(rows), len(DESCRIPTOR_NAMES)):
        raise ValueError("Descriptor shape must match the compounds and declared columns")
    y = np.array([float(row["pIC50_median"]) for row in rows])
    if not np.isfinite(descriptors).all() or not np.isfinite(y).all():
        raise ValueError("Baselines require finite descriptors and labels; no imputation")
    train = np.array([side == "train" for side in assignment])
    predictions = {}
    for name, model in make_models().items():
        model.fit(descriptors[train], y[train])
        values = model.predict(descriptors)
        predictions[name] = [
            {"inchikey": row["inchikey"], "scaffold": row["scaffold"],
             "partition": side, "observed_pIC50": float(observed),
             "predicted_pIC50": float(predicted)}
            for row, side, observed, predicted in zip(rows, assignment, y, values, strict=True)]
    return predictions


def run_baselines(raw_dir: Path, curated_dir: Path, split_dir: Path,
                  feature_dir: Path, out_dir: Path) -> dict:
    raw_dir, curated_dir, split_dir, feature_dir, out_dir = map(
        Path, (raw_dir, curated_dir, split_dir, feature_dir, out_dir))
    dest = out_dir.resolve()
    for source in (raw_dir, curated_dir, split_dir, feature_dir, DATA_RAW):
        source = source.resolve()
        if dest == source or source in dest.parents or dest in source.parents:
            raise ValueError("Baseline outputs must be separate from all input trees")
    if not verify_raw(raw_dir):
        raise ValueError("Raw checksum mismatch; baselines refused")
    curated = json.loads((curated_dir / MANIFEST).read_text())
    audit = json.loads((split_dir / "audit.json").read_text())
    features = json.loads((feature_dir / MANIFEST).read_text())
    raw_hash = _sha256(raw_dir / MANIFEST)
    compound_hash = _sha256(curated_dir / "compounds.csv")
    audit_hash = _sha256(split_dir / "audit.json")
    if (curated["source"]["raw_manifest_sha256"] != raw_hash
            or audit["raw_manifest_sha256"] != raw_hash
            or audit["curated_compounds_sha256"] != compound_hash
            or features["source"]["raw_manifest_sha256"] != raw_hash
            or features["source"]["curated_compounds_sha256"] != compound_hash
            or features["source"]["split_audit_sha256"] != audit_hash
            or len({curated["target_chembl_id"], audit["target_chembl_id"],
                    features["target_chembl_id"]}) != 1):
        raise ValueError("Baseline inputs disagree on raw/curated/split/feature provenance")
    sources = [raw_dir / MANIFEST, curated_dir / MANIFEST,
               split_dir / "audit.json", feature_dir / MANIFEST]
    for directory, manifest in ((curated_dir, curated), (split_dir, audit), (feature_dir, features)):
        for name, meta in manifest["files"].items():
            path = directory / name
            if _sha256(path) != meta["sha256"]:
                raise ValueError(f"Input checksum mismatch: {name}")
            sources.append(path)
    before = {path: _sha256(path) for path in sources}
    keys, descriptors, _ = load_features(feature_dir)
    compounds = sorted(_read_csv(curated_dir / "compounds.csv"), key=lambda r: r["inchikey"])
    if (keys != [r["inchikey"] for r in compounds]
            or len(keys) != curated["curated_compounds"]
            or descriptors.shape != (len(keys), len(DESCRIPTOR_NAMES))):
        raise ValueError("Cached descriptor rows do not reconcile with curated compounds")
    report = {
        "target_chembl_id": curated["target_chembl_id"],
        "command": f".venv/bin/python -m src.baselines {curated['target_chembl_id']}",
        "purpose": "Gate 3 baseline floor only; fixed, untuned hyperparameters; no Gate 4 model",
        "seed": SEED, "sklearn_version": sklearn.__version__, "numpy_version": np.__version__,
        "rdkit_version": features["rdkit_version"],
        "source": {"raw_manifest_sha256": raw_hash, "curated_compounds_sha256": compound_hash,
                   "split_audit_sha256": audit_hash,
                   "feature_manifest_sha256": before[feature_dir / MANIFEST],
                   "assignment_sha256": {name: audit["files"][f"{name}_assignments.csv"]["sha256"]
                                         for name in ("random", "scaffold")}},
        "features": list(DESCRIPTOR_NAMES),
        "preprocessing": "Stateless per-molecule cached descriptors only. No scaling, imputation, "
                         "variance filtering or feature selection. Non-finite values abort. "
                         "Both estimators, including the mean, fit training rows only. "
                         "Equal compound weights; no weighting by repeated report count.",
        "models": {name: {"class": type(model).__name__, "parameters": model.get_params()}
                   for name, model in make_models().items()},
        "hyperparameter_rationale": "Declared before the first run; no search. The mean ignores "
            "chemistry and is the trivial floor. The forest uses 200 bootstrap trees to reduce "
            "Monte Carlo variation, all descriptors at each node, fully grown trees and unit "
            "minimum leaf size as a conventional untuned regression baseline. Its flexibility "
            f"makes train/test gaps visible. One worker and seed {SEED} fix the execution. These "
            "are engineering defaults, not evidence of optimality.",
        "uncertainty": f"95% percentile intervals from {BOOTSTRAP_REPEATS} whole-test-scaffold "
            f"bootstrap draws, compound-weighted metrics, seed {SEED}; same draws for both models. "
            "Fixed fitted models, training sets and splits: no retraining or split-seed uncertainty. "
            "Cross-core similarity can still couple groups. Intervals are descriptive, not a "
            "significance test. Degenerate R2 draws are excluded and counted in the log.",
        "training_uncertainty": "Not estimated: train scores are resubstitution diagnostics of "
            "the fitted estimator, not an estimate of generalisation; test intervals are reported.",
        "splits": {},
        "disclaimer": "Literature/patent-reported IC50; no efficacy, safety or suitability claim.",
    }
    artifacts = {}
    for split in ("random", "scaffold"):
        saved = _read_csv(split_dir / f"{split}_assignments.csv")
        if ([r["inchikey"] for r in saved] != keys
                or any(r["pIC50_median"] != c["pIC50_median"]
                       for r, c in zip(saved, compounds, strict=True))):
            raise ValueError(f"Saved {split} row order or labels disagree with curated features")
        assignment = [r["partition"] for r in saved]
        rows = [{**c, "scaffold": r["scaffold"]} for c, r in zip(compounds, saved, strict=True)]
        predictions = predict_partition(rows, descriptors, assignment, split == "scaffold")
        result = {"train_compounds": assignment.count("train"),
                  "test_compounds": assignment.count("test"),
                  "training_mean_pIC50": float(np.mean([float(r["pIC50_median"])
                       for r, side in zip(rows, assignment) if side == "train"])), "models": {}}
        for side in ("train", "test"):
            if result[f"{side}_compounds"] != audit["splits"][split][side]["compounds"]:
                raise ValueError(f"Saved {split} partition counts disagree with audit")
        for model, values in predictions.items():
            train = [p for p in values if p["partition"] == "train"]
            test = [p for p in values if p["partition"] == "test"]
            result["models"][model] = {"train": point_metrics(train), "test": point_metrics(test),
                                       "test_bootstrap": clustered_intervals(test)}
            artifacts[f"{split}_{model}_predictions.csv"] = values
        report["splits"][split] = result
    if not verify_raw(raw_dir) or any(_sha256(path) != digest for path, digest in before.items()):
        raise ValueError("Baseline inputs changed during evaluation")
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, values in artifacts.items():
        _write_csv(out_dir / name, values, list(values[0]))
    report["inputs_verified_before_and_after"] = True
    report["files"] = {name: {"sha256": _sha256(out_dir / name)} for name in artifacts}
    (out_dir / "baselines.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    write_markdown(report, out_dir / "baselines.md")
    return report


def metrics_table(report: dict) -> str:
    lines = ["| Split / baseline | Train / test n | Train R² | Train RMSE | Train MAE | "
             "Test R² (95% interval) | Test RMSE (95% interval) | Test MAE (95% interval) |",
             "|---|---:|---:|---:|---:|---:|---:|---:|"]
    for split, result in report["splits"].items():
        for model, score in result["models"].items():
            cells = [f"{split} / {MODEL_NAMES[model]}",
                     f"{result['train_compounds']:,} / {result['test_compounds']:,}"]
            cells += [f"{score['train'][metric]:.3f}" if score['train'][metric] is not None
                      else "undefined" for metric in ("r2", "rmse", "mae")]
            for metric in ("r2", "rmse", "mae"):
                value = score["test"][metric]
                ci = score["test_bootstrap"]["metrics"][metric]["ci95"]
                cells.append(f"{value:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]"
                             if value is not None and ci else "undefined")
            lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def write_markdown(report: dict, path: Path) -> None:
    path.write_text("\n".join([
        "# Gate 3 fixed baselines", "", f"Command: `{report['command']}`.", "",
        metrics_table(report), "", "RMSE and MAE are pIC50 log units; R² is dimensionless.", "",
        report["uncertainty"], "", report["training_uncertainty"], "",
        report["preprocessing"], "", report["hyperparameter_rationale"], "",
        "Full estimator parameters (including library defaults) and input/output checksums are "
        "in baselines.json. Each predictions CSV contains both partitions with explicit labels. "
        "Saved Gate 2 assignments are read directly, never regenerated by this module.", "",
        "The Gate 2 fingerprint nearest-neighbour probe remains a separate structural leakage "
        "diagnostic. These baselines supersede it as the floor Gate 4 must beat, not as the "
        "evidence for the Gate 2 evaluation contrast. No fingerprint model is fitted here.", "",
        report["disclaimer"], ""]), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    target = parser.parse_args(argv).target_id.strip().upper()
    report = run_baselines(DATA_RAW / target, DATA_CURATED / target,
                           RESULTS / "splits" / target, DATA_FEATURES / target,
                           RESULTS / "baselines" / target)
    print(metrics_table(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
