"""Gate 2 only: a fixed Tanimoto 1-nearest-neighbour leakage probe, not a baseline.

Run .venv/bin/python -m src.leakage_probe CHEMBL284 after src.splits.
No tuning, model selection, repeated splits or assay-row weighting is performed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import DataStructs, rdBase

from src.config import DATA_CURATED, DATA_RAW, RESULTS, SEED
from src.retrieve import MANIFEST, _sha256, verify_raw
from src.splits import (FINGERPRINT_SPEC, NEAR_DUPLICATE_THRESHOLD, _write_csv,
                        assert_split_integrity, describe_structures, make_partitions)

BOOTSTRAP_REPEATS = 2000
PROBE_SPEC = {
    "name": "Fixed Tanimoto 1-nearest-neighbour leakage probe (not a baseline)",
    "fingerprint": FINGERPRINT_SPEC,
    "prediction": "Training pIC50_median at maximum Tanimoto; average all exact maximum ties",
    "fit": "Store training fingerprints and compound-median labels only; no tuning",
}


def predict_partition(rows: list[dict], fingerprints: list, assignment: list[str],
                      scaffold_disjoint: bool = False) -> list[dict]:
    """Fit by retaining training rows; held-out labels never enter predictions."""
    assert_split_integrity(rows, assignment, scaffold_disjoint)
    if len(rows) != len(fingerprints):
        raise ValueError("Fingerprint count must match compound count")
    train = [i for i, side in enumerate(assignment) if side == "train"]
    train_fp = [fingerprints[i] for i in train]
    train_labels = np.array([float(rows[i]["pIC50_median"]) for i in train])
    train_scaffolds = {rows[i]["scaffold"] for i in train}
    predictions = []
    for i, side in enumerate(assignment):
        if side != "test":
            continue
        similarities = np.array(DataStructs.BulkTanimotoSimilarity(fingerprints[i], train_fp))
        maximum = float(similarities.max())
        ties = np.flatnonzero(similarities == maximum)
        predicted = float(train_labels[ties].mean())
        predictions.append({
            "inchikey": rows[i]["inchikey"], "scaffold": rows[i]["scaffold"],
            "observed_pIC50": float(rows[i]["pIC50_median"]), "predicted_pIC50": predicted,
            "nearest_train_inchikey": rows[train[int(ties[0])]]["inchikey"],
            "nearest_tanimoto": maximum, "maximum_similarity_ties": len(ties),
            "shared_train_scaffold": rows[i]["scaffold"] in train_scaffolds,
            "near_train_neighbour": maximum >= NEAR_DUPLICATE_THRESHOLD,
            "n_records": int(rows[i]["n_records"]),
        })
    return predictions


def point_metrics(predictions: list[dict]) -> dict:
    observed = np.array([p["observed_pIC50"] for p in predictions], dtype=float)
    predicted = np.array([p["predicted_pIC50"] for p in predictions], dtype=float)
    if not len(observed) or not np.isfinite(observed).all() or not np.isfinite(predicted).all():
        raise ValueError("Metrics require nonempty finite observations and predictions")
    error = predicted - observed
    denominator = float(np.square(observed - observed.mean()).sum())
    return {"rmse": float(np.sqrt(np.square(error).mean())),
            "mae": float(np.abs(error).mean()),
            "r2": 1 - float(np.square(error).sum()) / denominator if denominator > 0 else None}


def clustered_intervals(predictions: list[dict], repeats: int = BOOTSTRAP_REPEATS,
                        seed: int = SEED) -> dict:
    """Percentile intervals: resample whole held-out exact scaffolds with replacement.

    Sufficient statistics preserve compound-weighted RMSE/MAE/R2 within each draw.
    Empty (acyclic) scaffolds are one group, matching the split rule. Predictions and
    training data are fixed. Degenerate R2 draws are excluded and counted explicitly.
    """
    if repeats < 2:
        raise ValueError("Need at least two bootstrap draws")
    point_metrics(predictions)  # validate before forming sufficient statistics
    groups = sorted({p["scaffold"] for p in predictions})
    stats = []
    for group in groups:
        selected = [p for p in predictions if p["scaffold"] == group]
        y = np.array([p["observed_pIC50"] for p in selected])
        e = np.array([p["predicted_pIC50"] - p["observed_pIC50"] for p in selected])
        stats.append([len(y), y.sum(), np.square(y).sum(), np.square(e).sum(), np.abs(e).sum()])
    stats = np.array(stats)
    rng = np.random.default_rng(seed)
    samples = {name: [] for name in ("rmse", "mae", "r2")}
    for _ in range(repeats):
        n, sy, sy2, se2, sae = stats[rng.integers(len(groups), size=len(groups))].sum(axis=0)
        samples["rmse"].append(float(np.sqrt(se2 / n)))
        samples["mae"].append(float(sae / n))
        denominator = sy2 - sy * sy / n
        if denominator > np.finfo(float).eps * max(1.0, abs(sy2)):
            samples["r2"].append(float(1 - se2 / denominator))
    return {"scaffold_groups": len(groups), "draws": repeats, "seed": seed,
            "metrics": {name: {"ci95": np.quantile(values, [.025, .975]).tolist() if values else None,
                               "valid_draws": len(values)} for name, values in samples.items()}}


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def run_probe(raw_dir: Path, curated_dir: Path, split_dir: Path, out_dir: Path) -> dict:
    raw_dir, curated_dir, split_dir, out_dir = map(Path, (raw_dir, curated_dir, split_dir, out_dir))
    dest = out_dir.resolve()
    for source in (raw_dir.resolve(), curated_dir.resolve(), split_dir.resolve(), DATA_RAW.resolve()):
        if dest == source or source in dest.parents or dest in source.parents:
            raise ValueError("Probe outputs must be separate from raw, curated and split inputs")
    if not verify_raw(raw_dir):
        raise ValueError("Raw checksum mismatch; probe refused")
    manifest = json.loads((curated_dir / MANIFEST).read_text())
    audit = json.loads((split_dir / "audit.json").read_text())
    raw_hash, compound_hash = _sha256(raw_dir / MANIFEST), _sha256(curated_dir / "compounds.csv")
    if (manifest["source"]["raw_manifest_sha256"] != raw_hash
            or audit["raw_manifest_sha256"] != raw_hash
            or audit["curated_compounds_sha256"] != compound_hash):
        raise ValueError("Probe inputs disagree on raw/curated provenance")
    for directory, files in ((curated_dir, manifest["files"]), (split_dir, audit["files"])):
        for name, meta in files.items():
            if _sha256(directory / name) != meta["sha256"]:
                raise ValueError(f"Input checksum mismatch: {name}")
    # Snapshot all source files used by the probe; recheck before publishing results.
    sources = [raw_dir / MANIFEST, curated_dir / MANIFEST, split_dir / "audit.json"]
    sources += [curated_dir / name for name in manifest["files"]]
    sources += [split_dir / name for name in audit["files"]]
    before = {path: _sha256(path) for path in sources}
    compounds = _read_csv(curated_dir / "compounds.csv")
    if len(compounds) != manifest["curated_compounds"]:
        raise ValueError("Curated count disagrees with manifest")
    rows, fingerprints = describe_structures(compounds)
    partitions = make_partitions(rows)
    report = {
        "target_chembl_id": manifest["target_chembl_id"],
        "command": f".venv/bin/python -m src.leakage_probe {manifest['target_chembl_id']}",
        "purpose": "Gate 2 leakage demonstration only; not a Gate 3 baseline or selected model",
        "probe": PROBE_SPEC, "seed": SEED, "rdkit_version": rdBase.rdkitVersion,
        "raw_manifest_sha256": raw_hash, "curated_compounds_sha256": compound_hash,
        "split_audit_sha256": before[split_dir / "audit.json"],
        "uncertainty": "95% percentile bootstrap, 2000 resamples of whole test exact-scaffold "
            "groups, compound-weighted metrics; fixed training set, split and predictions. "
            "Not uncertainty across split seeds or retraining. Cross-core similarity can "
            "still couple groups. Intervals are descriptive, not a significance test of the gap.",
        "splits": {},
        "disclaimer": "Literature/patent-reported IC50; no efficacy, safety or suitability claim.",
    }
    artifacts = {}
    for name, assignment in partitions.items():
        saved = _read_csv(split_dir / f"{name}_assignments.csv")
        expected = [{"inchikey": r["inchikey"], "partition": side, "scaffold": r["scaffold"],
                     "pIC50_median": r["pIC50_median"]} for r, side in zip(rows, assignment)]
        if saved != expected:
            raise ValueError(f"Saved {name} assignment/labels disagree with deterministic split")
        predictions = predict_partition(rows, fingerprints, assignment, name == "scaffold")
        strata = {}
        for field in ("shared_train_scaffold", "near_train_neighbour"):
            for value in (True, False):
                selected = [p for p in predictions if p[field] == value]
                strata[f"{field}={value}"] = {"compounds": len(selected),
                    "metrics": point_metrics(selected) if selected else None}
        report["splits"][name] = {
            "train_compounds": assignment.count("train"), "test_compounds": len(predictions),
            "metrics": point_metrics(predictions), "bootstrap": clustered_intervals(predictions),
            "strata": strata,
            "training_replicate_rows_collapsed": sum(int(r["n_records"]) - 1
                for r, side in zip(rows, assignment) if side == "train"),
        }
        artifacts[f"{name}_predictions.csv"] = predictions
    sizes = Counter(r["scaffold"] for r in rows)
    report["largest_scaffold_groups"] = [
        {"scaffold": core, "compounds": count,
         "random_test_compounds": sum(r["scaffold"] == core and side == "test"
                                       for r, side in zip(rows, partitions["random"]))}
        for core, count in sorted(sizes.items(), key=lambda item: (-item[1], item[0]))[:5]]
    report["rmse_gap_scaffold_minus_random"] = (
        report["splits"]["scaffold"]["metrics"]["rmse"] - report["splits"]["random"]["metrics"]["rmse"])
    if not verify_raw(raw_dir) or any(_sha256(path) != digest for path, digest in before.items()):
        raise ValueError("Probe inputs changed during evaluation")
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, predictions in artifacts.items():
        _write_csv(out_dir / name, predictions, list(predictions[0]))
    report["inputs_verified_before_and_after"] = True
    report["files"] = {name: {"sha256": _sha256(out_dir / name)} for name in artifacts}
    (out_dir / "probe.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_markdown(report, audit, out_dir / "probe.md")
    return report


def write_markdown(report: dict, audit: dict, path: Path) -> None:
    lines = ["# Gate 2 fixed leakage probe (not a baseline)", "", f"Command: `{report['command']}`.",
             "", report["probe"]["name"], "", report["probe"]["prediction"] + ". " + report["probe"]["fit"] + ".",
             "Features: nonchiral Morgan radius 2, 2048 bits; identical for both saved splits. "
             "All compounds carry equal weight after median aggregation. No parameter or split search.", "",
             "| Metric (95% interval) | Random | Scaffold |", "|---|---:|---:|"]
    a, b = (report["splits"][name] for name in ("random", "scaffold"))
    for metric in ("rmse", "mae", "r2"):
        values = []
        for result in (a, b):
            point, ci = result["metrics"][metric], result["bootstrap"]["metrics"][metric]["ci95"]
            values.append(f"{point:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]" if point is not None and ci else "undefined")
        lines.append(f"| {metric.upper()} | {' | '.join(values)} |")
    lines += ["", "RMSE and MAE are pIC50 log units; R2 is dimensionless.", "", report["uncertainty"],
              f"Test bootstrap groups: {a['bootstrap']['scaffold_groups']} random and "
              f"{b['bootstrap']['scaffold_groups']} scaffold. Seed {report['seed']}.", "",
              "## What this dataset lets the probe reuse", ""]
    ra, sa = (audit["splits"][name] for name in ("random", "scaffold"))
    lines += [f"Random holds out {a['test_compounds']} compounds, but {ra['test_compounds_with_train_scaffold']} "
              f"already have a training core ({ra['shared_scaffolds']} shared cores); "
              f"{ra['test_compounds_with_near_duplicate']} have a training neighbour at Tanimoto >=0.90. "
              f"Scaffold holds out {b['test_compounds']} compounds with zero shared cores, yet retains "
              f"{sa['test_compounds_with_near_duplicate']} such close-neighbour compounds. "
              "A stored analogue's label can therefore interpolate a familiar series under random allocation. "
              "Exact InChIKey isolation alone does not prevent this structural information reuse.", "",
              "Largest exact-core groups (canonical scaffold strings and sizes are in probe.json): "
              + ", ".join(f"{g['compounds']} compounds ({g['random_test_compounds']} random-test)"
                          for g in report["largest_scaffold_groups"]) + ".", "",
              f"Observed scaffold-minus-random RMSE: {report['rmse_gap_scaffold_minus_random']:+.3f} log units. "
              "The contrast measures this probe's sensitivity to the evaluation boundary. It does not isolate "
              "a causal leakage penalty: the test compounds and activity distributions also differ. "
              "No held-out labels affected features, assignments, predictions or parameter choices.", "",
              "| Split / diagnostic subgroup | Test compounds | RMSE |", "|---|---:|---:|"]
    for name in ("random", "scaffold"):
        for label, stratum in report["splits"][name]["strata"].items():
            score = f"{stratum['metrics']['rmse']:.3f}" if stratum["metrics"] else "undefined"
            lines.append(f"| {name}: {label} | {stratum['compounds']} | {score} |")
    lines += ["", "Subgroups explain the boundary; they were not used to tune or exclude compounds. "
              "They differ in chemistry and size and are not controlled experiments. "
              "Close neighbours need not agree in measured activity; compare their errors above "
              "rather than assuming similarity guarantees reliable prediction.", "",
              "Over-reported compounds cannot leak identical rows here: one median label per standardised "
              "InChIKey enters each split, and exact overlap is asserted zero. Training-side repeat rows "
              f"collapsed before fitting: {a['training_replicate_rows_collapsed']} random / "
              f"{b['training_replicate_rows_collapsed']} scaffold. The probe never weights those rows. "
              "Scaffold holdout also removes shared connectivity keys; residual fingerprint similarity "
              "can reflect related cores, stereochemistry ignored by the fingerprint, or bit collisions.", "",
              "These intervals do not turn conditional replicate dispersion into a universal RMSE floor. "
              "The predicted targets are compound medians, whereas replicate SD describes individual "
              "reported measurements on a selected repeated subset. A lower RMSE alone cannot establish leakage.",
              "", report["disclaimer"], ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    target = parser.parse_args(argv).target_id.strip().upper()
    report = run_probe(DATA_RAW / target, DATA_CURATED / target,
                       RESULTS / "splits" / target, RESULTS / "leakage_probe" / target)
    for name, result in report["splits"].items():
        r2 = result["metrics"]["r2"]
        r2_text = f"{r2:.3f}" if r2 is not None else "undefined"
        print(f"{name} fixed leakage probe: RMSE {result['metrics']['rmse']:.3f}; "
              f"MAE {result['metrics']['mae']:.3f}; R2 {r2_text}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
