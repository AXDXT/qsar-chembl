"""Gate 4 training-only, scaffold-grouped nested forest search.

Call run_nested_cv separately for each saved split and feature block, passing ONLY
rows whose saved partition is 'train', aligned stateless features, and upstream
checksums. Scaffold labels come from Gate 2's src.splits.describe_structures via
the saved assignment CSVs; the empty scaffold remains one group. This module
never loads or evaluates the saved held-out sets and never chooses a bit length.

Every completed fit is atomically cached. Resume with identical arguments and the
same cache directory. Changed inputs/settings/code abort; use a NEW directory for
a clean rebuild. No cache deletion or retraining is implicit. Single writer only.
The final search selects parameters on all training rows, but does not fit or
persist the final model: that and the single held-out evaluation belong downstream.
"""
from __future__ import annotations

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import platform
import tempfile
from time import perf_counter

import joblib
import numpy as np
from rdkit import rdBase
import sklearn
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import GroupKFold, ParameterGrid

from src.baselines import RF_SPEC
from src.config import SEED
from src.leakage_probe import BOOTSTRAP_REPEATS, clustered_intervals, point_metrics
from src.retrieve import _sha256
from src.splits import assert_split_integrity

# Eight candidates, fixed before results. Feature subsampling spans sqrt(p), a
# modest fraction, a larger fraction and all features: sparse Morgan bits can
# otherwise crowd out the dense descriptor block. 1.0 with leaf 1 IS the Gate 3
# forest, so every block's search contains the untuned baseline and tuning can
# never be forced away from it. Leaf size tests modest regularisation against
# fully grown trees. Keep the baseline's 200 trees and seed fixed so the tuning
# control changes only the two declared axes. Worker count changes speed only.
DEFAULT_GRID = {"max_features": ["sqrt", 0.1, 0.33, 1.0], "min_samples_leaf": [1, 3]}
GRID_RATIONALE = {
    "max_features": "sqrt(p), 10%, 33% and 100% vary feature competition and tree diversity, "
                    "especially for dense descriptors alongside sparse Morgan bits. 100% is the "
                    "Gate 3 setting and, on the 12-column descriptor block, the only value that "
                    "lets a split consider every descriptor (sqrt and 33% both round to 3 there, "
                    "10% to 1; the duplicate resolves by grid order).",
    "min_samples_leaf": "1 versus 3 tests modest regularisation of the Gate 3 train/test gap.",
    "baseline_in_grid": "max_features=1.0 with min_samples_leaf=1 reproduces the Gate 3 forest, "
                        "so the tuned control's candidate set contains the untuned baseline.",
    "n_estimators": "Fixed at the Gate 3 count, not selected; more trees reduce Monte Carlo "
                    "variation. Keeping the count fixed isolates the two search axes.",
    "compute": "Default grid: eight candidates; fits run one after another, each forest using "
               "n_jobs workers. Completed fits are cached individually, including timing, so "
               "interruption loses at most one fit.",
}


def _json(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)


def _digest(value) -> str:
    return hashlib.sha256(_json(value).encode()).hexdigest()


def _atomic_json(path: Path, value: dict) -> None:
    """Never expose a partially written completed-fit record on interruption."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix=".pending-", dir=path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as stream:
            stream.write(_json(value) + "\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


class _Cache:
    def __init__(self, directory: Path, identity: dict):
        self.directory = Path(directory)
        self.identity = _digest(identity)
        manifest = self.directory / "manifest.json"
        if manifest.exists():
            if json.loads(manifest.read_text()) != identity:
                raise ValueError("Nested-CV cache inputs/settings/code changed; use a new cache directory")
        else:
            # Refuse orphaned fit records rather than adopting them into a new run.
            if self.directory.exists() and any(
                    not p.name.startswith("._") for p in self.directory.iterdir()):
                raise ValueError("Nested-CV cache has no manifest; use a new cache directory")
            _atomic_json(manifest, identity)

    def get(self, name: str, compute):
        path = self.directory / f"{name}.json"
        if path.exists():
            record = json.loads(path.read_text())
            if (record.get("identity") != self.identity or record.get("name") != name
                    or record.get("sha256") != _digest(record.get("result"))):
                raise ValueError(f"Nested-CV cache checksum/identity mismatch: {name}")
            return record["result"]
        result = compute()
        _atomic_json(path, {"identity": self.identity, "name": name,
                            "sha256": _digest(result), "result": result})
        return result


def _folds(rows: list[dict], indices: list[int], count: int) -> list[dict]:
    """Map local GroupKFold indices back to training-pool indices; audit each fold."""
    subset = [rows[i] for i in indices]
    groups = [r["scaffold"] for r in subset]
    if len(set(groups)) < count:
        raise ValueError(f"Need at least {count} training scaffold groups in every CV loop")
    result, seen = [], []
    for train, validation in GroupKFold(n_splits=count).split(indices, groups=groups):
        train, validation = list(map(int, train)), list(map(int, validation))
        if (len(set(train)) != len(train) or len(set(validation)) != len(validation)
                or set(train) & set(validation)
                or set(train) | set(validation) != set(range(len(indices)))):
            raise ValueError("CV fold must partition its training pool exactly once")
        validation_set = set(validation)
        assignment = ["test" if i in validation_set else "train" for i in range(len(indices))]
        assert_split_integrity(subset, assignment, scaffold_disjoint=True)
        result.append({"train": [indices[i] for i in train],
                       "validation": [indices[i] for i in validation]})
        seen.extend(validation)
    if len(result) != count or sorted(seen) != list(range(len(indices))):
        raise ValueError("CV validation folds must cover each training compound once")
    return result


def _fit(rows, x, y, fold, parameters, n_jobs=1) -> dict:
    # Per-tree seeds are drawn before dispatch, so workers change speed, not the trees;
    # prediction sums may differ only in floating-point order (tested to 1e-12).
    start = perf_counter()
    model = RandomForestRegressor(**{**RF_SPEC, **parameters, "n_jobs": n_jobs})
    model.fit(x[fold["train"]], y[fold["train"]])
    predicted = model.predict(x[fold["validation"]])
    predictions = [{"inchikey": rows[i]["inchikey"], "scaffold": rows[i]["scaffold"],
                    "observed_pIC50": float(y[i]), "predicted_pIC50": float(value)}
                   for i, value in zip(fold["validation"], predicted, strict=True)]
    return {"metrics": point_metrics(predictions), "predictions": predictions,
            "fit_predict_seconds": perf_counter() - start, "n_jobs": n_jobs}


def _search(rows, x, y, folds, candidates, cache, prefix, n_jobs) -> dict:
    scores = []
    for candidate, parameters in enumerate(candidates):
        results = []
        for number, fold in enumerate(folds):
            result = cache.get(f"{prefix}-candidate{candidate}-fold{number}",
                              lambda: _fit(rows, x, y, fold, parameters, n_jobs=n_jobs))
            results.append(result)
        scores.append({"parameters": parameters,
                       "mean_validation_rmse": float(np.mean([r["metrics"]["rmse"] for r in results])),
                       "fold_metrics": [r["metrics"] for r in results],
                       "fit_predict_seconds": [r["fit_predict_seconds"] for r in results]})
    # ParameterGrid order breaks exact ties deterministically, without outer scores.
    winner = min(range(len(scores)), key=lambda i: scores[i]["mean_validation_rmse"])
    return {"selected_parameters": candidates[winner], "selected_candidate": winner,
            "candidates": scores}


def run_nested_cv(rows: list[dict], values: np.ndarray, cache_dir: Path, *,
                  feature_spec: dict, source_checksums: dict[str, str],
                  grid: dict | None = None, folds: int = 5,
                  bootstrap_repeats: int = BOOTSTRAP_REPEATS, n_jobs: int = 1) -> dict:
    """Estimate generalisation and select final parameters WITHOUT held-out rows.

    rows must retain 'partition=train' from the saved assignments, plus inchikey,
    scaffold and pIC50_median. The caller verifies upstream provenance and subsets
    the saved assignments and aligned features BEFORE calling. Each feature block
    gets its own cache directory; feature_spec records its column order, fingerprint
    settings and any training-only bit-length decision. No fitted preprocessing.

    Both CV loops use GroupKFold, including within the random split's training pool:
    their estimand is transfer across exact scaffolds, not random interpolation.
    Selection minimises equally weighted mean inner-fold RMSE (not pooled RMSE).
    Outer predictions are strictly out of fold; full-training selection is separate
    and never uses outer scores. The bootstrap conditions on fixed OOF predictions;
    outer-fold spread also reflects refitting on overlapping training subsets.

    n_jobs is the per-forest worker count. It is deliberately NOT part of the cache
    identity, so a resumed run on a machine with different cores reuses earlier fits.
    """
    if not rows or any(r.get("partition") != "train" for r in rows):
        raise ValueError("Nested CV accepts saved training rows only; held-out rows are forbidden")
    if len({r["inchikey"] for r in rows}) != len(rows) or any(not r["inchikey"] for r in rows):
        raise ValueError("Training InChIKeys must be nonempty and unique")
    if any(not isinstance(r["scaffold"], str) for r in rows):
        raise ValueError("Use saved Gate 2 scaffold strings, including the empty scaffold")
    x = np.asarray(values, dtype=np.float64)
    y = np.array([float(r["pIC50_median"]) for r in rows])
    if (x.ndim != 2 or x.shape[0] != len(rows) or x.shape[1] < 1
            or not np.isfinite(x).all() or not np.isfinite(y).all()):
        raise ValueError("Require aligned, nonempty finite features and labels; no imputation")
    if not feature_spec or not source_checksums:
        raise ValueError("Feature specification and upstream checksums are required")
    if not isinstance(folds, int) or isinstance(folds, bool) or folds < 2:
        raise ValueError("Need at least two CV folds")
    if bootstrap_repeats < 2:
        raise ValueError("Need at least two bootstrap draws")
    if not isinstance(n_jobs, int) or isinstance(n_jobs, bool) or n_jobs < 1:
        raise ValueError("n_jobs must be a positive worker count, logged as given")
    grid = DEFAULT_GRID if grid is None else grid
    if set(grid) != {"max_features", "min_samples_leaf"}:
        raise ValueError("Search only max_features and min_samples_leaf; tree count is fixed")
    candidates = list(ParameterGrid(grid))
    for parameters in candidates:
        RandomForestRegressor(**{**RF_SPEC, **parameters})._validate_params()
    indices = list(range(len(rows)))
    outer = _folds(rows, indices, folds)
    inner = [_folds(rows, fold["train"], folds) for fold in outer]
    final = _folds(rows, indices, folds)
    identity = {
        "schema": 1, "rows": [{k: r[k] for k in ("inchikey", "scaffold", "pIC50_median", "partition")}
                               for r in rows],
        "features_sha256": hashlib.sha256(x.tobytes(order="C")).hexdigest(),
        "feature_shape": list(x.shape), "feature_spec": feature_spec,
        "source_checksums": source_checksums, "candidates": candidates,
        "estimator": RandomForestRegressor(**RF_SPEC).get_params(),
        "versions": {"python": platform.python_version(), "numpy": np.__version__,
                     "sklearn": sklearn.__version__, "rdkit": rdBase.rdkitVersion,
                     "joblib": joblib.__version__},
        "code_sha256": {name: _sha256(Path(__file__).with_name(name)) for name in
                        ("nested_cv.py", "baselines.py", "config.py", "features.py",
                         "splits.py", "leakage_probe.py", "retrieve.py")},
        "folds": folds, "outer_plan": outer, "inner_plan": inner, "final_plan": final,
        "bootstrap_repeats": bootstrap_repeats, "seed": SEED,
    }
    # Normalise JSON-compatible tuples etc. before comparing to a loaded manifest.
    identity = json.loads(_json(identity))
    cache = _Cache(cache_dir, identity)

    def complete():
        start = perf_counter()
        outer_results, predictions = [], []
        for number, fold in enumerate(outer):
            selection = _search(rows, x, y, inner[number], candidates, cache,
                                f"outer{number}-inner", n_jobs)
            evaluation = cache.get(f"outer{number}-evaluation", lambda: _fit(
                rows, x, y, fold, selection["selected_parameters"], n_jobs=n_jobs))
            predictions.extend({**r, "outer_fold": number} for r in evaluation["predictions"])
            outer_results.append({"outer_fold": number, "selection": selection,
                                  "metrics": evaluation["metrics"],
                                  "fit_predict_seconds": evaluation["fit_predict_seconds"]})
        final_selection = _search(rows, x, y, final, candidates, cache, "final-inner", n_jobs)
        counts = Counter(_json(r["selection"]["selected_parameters"]) for r in outer_results)
        spread = {}
        for metric in ("rmse", "mae", "r2"):
            valid = [r["metrics"][metric] for r in outer_results if r["metrics"][metric] is not None]
            spread[metric] = {"values": [r["metrics"][metric] for r in outer_results],
                              "mean": float(np.mean(valid)) if valid else None,
                              "sample_sd": float(np.std(valid, ddof=1)) if len(valid) > 1 else None}
        searches = [r["selection"] for r in outer_results] + [final_selection]
        fit_seconds = sum(sum(c["fit_predict_seconds"]) for s in searches for c in s["candidates"])
        fit_seconds += sum(r["fit_predict_seconds"] for r in outer_results)
        return {
            "identity_sha256": cache.identity, "grid": candidates, "grid_rationale": GRID_RATIONALE,
            "selection_rule": "Minimum unweighted mean inner-fold RMSE; ties use logged grid order.",
            "fold_grouping": "GroupKFold in both loops, using saved Gate 2 exact scaffold labels; "
                             "acyclic empty scaffolds stay together, even for random-split training.",
            "outer_folds": outer_results, "outer_fold_spread": spread,
            "oof_predictions": predictions, "oof_metrics": point_metrics(predictions),
            "oof_bootstrap": clustered_intervals(predictions, repeats=bootstrap_repeats),
            "uncertainty": "Whole-scaffold percentile bootstrap conditions on these fixed OOF "
                           "predictions and folds; excludes selection/refit/split-seed uncertainty. "
                           "Fold SD is descriptive, not a confidence interval: training folds overlap.",
            "selection_stability": {"all_outer_folds_agree": len(counts) == 1,
                "distinct_selections": len(counts), "counts": [
                    {"parameters": json.loads(key), "outer_folds": value}
                    for key, value in sorted(counts.items())]},
            "final_selection": final_selection,
            "timing": {"cumulative_fit_predict_seconds": fit_seconds,
                       "completion_invocation_seconds": perf_counter() - start,
                       "note": "Cumulative fit times include reused work; invocation time may be "
                               "a resumed run. Cached reports retain their original timings. "
                               "n_jobs is the worker count of the invocation that completed "
                               "the report; each cached fit records its own.",
                       "n_jobs": n_jobs, "search_workers": 1},
            "heldout_evaluations": 0,
            "disclaimer": "Literature/patent-reported IC50; no efficacy, safety or suitability claim.",
        }

    return cache.get("report", complete)
