"""Verified saved-data adapter and training-only Gate 4 descriptor control.

Run: .venv/bin/python -m src.gate4_training CHEMBL284 --n-jobs 4
No final estimator is fitted and no saved test prediction or metric is computed.
Fingerprint length selection and the combined model remain separate work.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from rdkit import rdBase

from src.config import DATA_CURATED, DATA_FEATURES, DATA_RAW, RESULTS
from src.features import DESCRIPTOR_NAMES, _read_csv, load_features
from src.nested_cv import run_nested_cv
from src.retrieve import MANIFEST, _sha256, verify_raw
from src.splits import assert_split_integrity


def load_training_inputs(raw_dir: Path, curated_dir: Path, split_dir: Path,
                         feature_dir: Path) -> dict:
    """Verify provenance/alignment, then expose ONLY saved training rows to search.

    Assignment labels are read to check agreement with curation, never to score
    held-out rows. Reuse Gate 2's saved scaffold strings, including the empty
    scaffold; do not regenerate assignments or derive a new scaffold definition.
    File checks cover both partitions, but no held-out metric is calculated.
    """
    raw_dir, curated_dir, split_dir, feature_dir = map(
        Path, (raw_dir, curated_dir, split_dir, feature_dir))
    if not verify_raw(raw_dir):
        raise ValueError("Raw checksum mismatch; Gate 4 training refused")
    paths = {"raw_manifest": raw_dir / MANIFEST,
             "curated_manifest": curated_dir / MANIFEST,
             "split_audit": split_dir / "audit.json",
             "feature_manifest": feature_dir / MANIFEST}
    before = {path: _sha256(path) for path in paths.values()}
    curated = json.loads(paths["curated_manifest"].read_text())
    audit = json.loads(paths["split_audit"].read_text())
    features = json.loads(paths["feature_manifest"].read_text())
    compound_hash = _sha256(curated_dir / "compounds.csv")
    raw_hash = before[paths["raw_manifest"]]
    if (curated["source"]["raw_manifest_sha256"] != raw_hash
            or audit["raw_manifest_sha256"] != raw_hash
            or audit["curated_compounds_sha256"] != compound_hash
            or features["source"]["raw_manifest_sha256"] != raw_hash
            or features["source"]["curated_compounds_sha256"] != compound_hash
            or features["source"]["split_audit_sha256"] != before[paths["split_audit"]]
            or len({curated["target_chembl_id"], audit["target_chembl_id"],
                    features["target_chembl_id"]}) != 1):
        raise ValueError("Gate 4 input provenance disagrees across raw/curated/split/features")
    if features["rdkit_version"] != rdBase.rdkitVersion:
        raise ValueError("Feature RDKit version changed; rebuild features before training")
    for directory, manifest in ((curated_dir, curated), (split_dir, audit),
                                (feature_dir, features)):
        for name, meta in manifest["files"].items():
            path = directory / name
            digest = _sha256(path)
            if digest != meta["sha256"]:
                raise ValueError(f"Input checksum mismatch: {name}")
            before[path] = digest

    keys, descriptors, bits = load_features(feature_dir)
    compounds = sorted(_read_csv(curated_dir / "compounds.csv"), key=lambda r: r["inchikey"])
    if (keys != [r["inchikey"] for r in compounds]
            or len(keys) != len(set(keys)) or not all(keys)
            or len(keys) != curated["curated_compounds"]
            or descriptors.shape != (len(keys), len(DESCRIPTOR_NAMES))
            or bits.shape != (len(keys), features["fingerprint_block"]["fpSize"])
            or not np.isfinite(descriptors).all()):
        raise ValueError("Cached feature rows/shape do not reconcile with curated compounds")
    # Curated manifests contain a run timestamp. Verify it above and below, but
    # key the search on actual data and deterministic upstream manifests so an
    # identical run_all rebuild does not invalidate expensive completed fits.
    source = {name: before[path] for name, path in paths.items() if name != "curated_manifest"}
    source["curated_compounds"] = compound_hash
    source["training_adapter_code"] = _sha256(Path(__file__))
    result = {"target_chembl_id": curated["target_chembl_id"], "splits": {}}
    scaffold_labels = None
    for split in ("random", "scaffold"):
        path = split_dir / f"{split}_assignments.csv"
        saved = _read_csv(path)
        if ([r["inchikey"] for r in saved] != keys
                or any(r["pIC50_median"] != c["pIC50_median"]
                       for r, c in zip(saved, compounds, strict=True))):
            raise ValueError(f"Saved {split} row order or labels disagree with curated features")
        labels = [r["scaffold"] for r in saved]
        if scaffold_labels is not None and labels != scaffold_labels:
            raise ValueError("Saved splits disagree on Gate 2 scaffold labels")
        scaffold_labels = labels
        assignment = [r["partition"] for r in saved]
        assert_split_integrity(saved, assignment, scaffold_disjoint=split == "scaffold")
        for side in ("train", "test"):
            if assignment.count(side) != audit["splits"][split][side]["compounds"]:
                raise ValueError(f"Saved {split} partition counts disagree with audit")
        selected = [i for i, side in enumerate(assignment) if side == "train"]
        result["splits"][split] = {
            "rows": [dict(saved[i]) for i in selected],
            "descriptors": descriptors[selected].copy(),
            "source_checksums": {**source, "assignments": before[path]},
        }
    if not verify_raw(raw_dir) or any(_sha256(path) != digest for path, digest in before.items()):
        raise ValueError("Gate 4 inputs changed during training-input preparation")
    return result


def run_descriptor_control(raw_dir: Path, curated_dir: Path, split_dir: Path,
                           feature_dir: Path, cache_dir: Path, *,
                           splits: tuple[str, ...] = ("scaffold", "random"),
                           n_jobs: int = 1) -> dict:
    """Resume the tuned descriptor control; cache one search per saved split.

    Choose a new cache directory for a clean rebuild. Never overwrite inputs or
    silently delete incompatible fits. The nested-CV report logs the full grid,
    fold selections, timing, worker count and training-only OOF uncertainty.
    """
    if not splits or len(set(splits)) != len(splits) or set(splits) - {"random", "scaffold"}:
        raise ValueError("Choose unique saved split names: scaffold and/or random")
    dest = Path(cache_dir).resolve()
    for source in (raw_dir, curated_dir, split_dir, feature_dir, DATA_RAW):
        source = Path(source).resolve()
        if dest == source or source in dest.parents or dest in source.parents:
            raise ValueError("Gate 4 cache must be separate from all input trees")
    inputs = load_training_inputs(raw_dir, curated_dir, split_dir, feature_dir)
    reports = {}
    for split in splits:
        training = inputs["splits"][split]
        reports[split] = run_nested_cv(
            training["rows"], training["descriptors"], dest / split / "descriptors",
            feature_spec={"block": "descriptors", "columns": list(DESCRIPTOR_NAMES),
                          "saved_split": split, "preprocessing": "None; stateless descriptors"},
            source_checksums=training["source_checksums"], n_jobs=n_jobs)
    return reports


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    parser.add_argument("--split", choices=("scaffold", "random", "both"), default="both")
    parser.add_argument("--n-jobs", type=int, default=1)
    parser.add_argument("--cache-dir", type=Path,
                        help="New directory for a clean rebuild; existing compatible fits resume")
    args = parser.parse_args(argv)
    target = args.target_id.strip().upper()
    selected = ("scaffold", "random") if args.split == "both" else (args.split,)
    reports = run_descriptor_control(
        DATA_RAW / target, DATA_CURATED / target, RESULTS / "splits" / target,
        DATA_FEATURES / target, args.cache_dir or RESULTS / "gate4_training" / target,
        splits=selected, n_jobs=args.n_jobs)
    for split, report in reports.items():
        print(f"{split} descriptor control complete: {report['identity_sha256']}; "
              "training-only nested CV; no held-out evaluation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
