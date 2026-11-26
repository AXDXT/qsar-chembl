"""Training-only, fold-local Morgan bit-length decisions for Gate 4.

Run: .venv/bin/python -m src.gate4_fingerprints CHEMBL284
This prepares a checksummed plan; it fits no estimator and scores no held-out set.
The combined search must use each fit's own decision, NOT the final-fit decision
throughout CV. Existing descriptor-control caches are not changed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import json
from pathlib import Path
from time import perf_counter

from rdkit import Chem, rdBase
import sklearn

from src.config import DATA_CURATED, DATA_FEATURES, DATA_RAW, RESULTS
from src.features import FOLDING_SENSITIVITY_BITS, _read_csv, folding_evidence
from src.gate4_training import load_training_inputs
from src.nested_cv import _Cache, _folds, _json
from src.retrieve import _sha256
from src.splits import FINGERPRINT_SPEC


SELECTION_RULE = (
    "Minimise total within-compound identifiers lost on folding across the fit's "
    "training structures; ties choose the smallest bit length. No labels, validation "
    "structures, validation scores or held-out rows enter a decision."
)
RATIONALE = (
    "Use the Gate 3 folding-sensitivity lengths as a fixed laptop-memory budget. "
    "Within that budget retain as many distinct unfolded Morgan identifiers per "
    "compound as possible. The tie rule avoids extra columns without measured "
    "collision reduction. This is a representation-fidelity criterion, not evidence "
    "that a longer fingerprint predicts better. Distinct sparse identifiers may "
    "themselves collide; this does not recover chemical information already lost "
    "there, or encode stereochemistry."
)


def _collision_evidence(rows: list[dict]) -> dict:
    """Stateless per-molecule measurements; no aggregate selection occurs here."""
    started = perf_counter()
    records = []
    for row in rows:
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None or Chem.MolToInchiKey(mol) != row["inchikey"]:
            raise ValueError(f"Curated SMILES/InChIKey mismatch: {row['inchikey']}")
        # Reuse Gate 3's exact generators/definition rather than approximating
        # collision loss by modulo arithmetic or inventing a new fingerprint.
        evidence = folding_evidence([mol])
        identifiers = evidence["distinct_substructure_identifiers"]
        records.append({"inchikey": row["inchikey"], "identifiers": identifiers,
                        "lost_by_bit_length": {
                            str(bits): identifiers - evidence["by_bit_length"][str(bits)]["occupied_bits"]
                            for bits in FOLDING_SENSITIVITY_BITS}})
    return {"rows": records, "measurement_seconds": perf_counter() - started}


def _decision(evidence: list[dict], indices: list[int]) -> dict:
    """Aggregate ONLY the indices used to fit this particular estimator."""
    if not indices or len(set(indices)) != len(indices):
        raise ValueError("Bit-length selection needs unique nonempty fitting indices")
    candidates = []
    for bits in FOLDING_SENSITIVITY_BITS:
        losses = [evidence[i]["lost_by_bit_length"][str(bits)] for i in indices]
        candidates.append({"fpSize": bits, "lost_identifiers": sum(losses),
                           "compounds_losing_identifiers": sum(value > 0 for value in losses)})
    winner = min(candidates, key=lambda c: (c["lost_identifiers"], c["fpSize"]))
    return {"train_indices": indices, "training_compounds": len(indices),
            "selected_fpSize": winner["fpSize"], "candidates": candidates}


def plan_fingerprints(rows: list[dict], cache_dir: Path, *,
                      source_checksums: dict[str, str], folds: int = 5) -> dict:
    """Plan every nested fit and the eventual final fit from saved training rows.

    Per-row collision measurements can be shared because they are stateless. Every
    aggregate bit-length decision excludes that fit's validation rows, including
    the inner validation fold. Labels are deliberately absent from this identity
    and from all selection functions (upstream file checksums still pin the data).
    The fold plans reuse src.nested_cv._folds and saved Gate 2 scaffold strings.
    """
    if not rows or any(row.get("partition") != "train" for row in rows):
        raise ValueError("Fingerprint selection accepts saved training rows only")
    if (len({r["inchikey"] for r in rows}) != len(rows)
            or any(not r["inchikey"] or not isinstance(r["scaffold"], str)
                   or not r.get("canonical_smiles") for r in rows)):
        raise ValueError("Require unique keys, saved scaffold strings and curated SMILES")
    if not isinstance(folds, int) or isinstance(folds, bool) or folds < 2:
        raise ValueError("Need at least two grouped folds")
    if not source_checksums:
        raise ValueError("Upstream checksums are required")
    indices = list(range(len(rows)))
    outer = _folds(rows, indices, folds)
    inner = [_folds(rows, fold["train"], folds) for fold in outer]
    final = _folds(rows, indices, folds)
    identity = json.loads(_json({
        "schema": 1,
        "rows": [{k: r[k] for k in ("inchikey", "canonical_smiles", "scaffold", "partition")}
                 for r in rows],
        "source_checksums": source_checksums,
        "fingerprint_spec": FINGERPRINT_SPEC,
        "candidate_bits": FOLDING_SENSITIVITY_BITS,
        "selection_rule": SELECTION_RULE, "rationale": RATIONALE,
        "versions": {"rdkit": rdBase.rdkitVersion, "sklearn": sklearn.__version__},
        "code_sha256": {name: _sha256(Path(__file__).with_name(name)) for name in
                        ("gate4_fingerprints.py", "gate4_training.py", "nested_cv.py",
                         "features.py", "splits.py", "retrieve.py", "config.py")},
        "folds": folds, "outer_plan": outer, "inner_plan": inner, "final_plan": final,
    }))
    cache = _Cache(cache_dir, identity)

    def complete():
        start = perf_counter()
        measured = cache.get("compound-collisions", lambda: _collision_evidence(rows))
        evidence = measured["rows"]
        if [r["inchikey"] for r in evidence] != [r["inchikey"] for r in rows]:
            raise ValueError("Collision evidence row alignment changed")
        decisions = {}
        for number, fold in enumerate(outer):
            for inner_number, inner_fold in enumerate(inner[number]):
                decisions[f"outer{number}-inner-fold{inner_number}"] = _decision(
                    evidence, inner_fold["train"])
            decisions[f"outer{number}-evaluation"] = _decision(evidence, fold["train"])
        for number, fold in enumerate(final):
            decisions[f"final-inner-fold{number}"] = _decision(evidence, fold["train"])
        decisions["final-fit"] = _decision(evidence, indices)
        counts = Counter(d["selected_fpSize"] for d in decisions.values())
        return {"identity_sha256": cache.identity, "selection_rule": SELECTION_RULE,
                "rationale": RATIONALE, "fit_decisions": decisions,
                "selected_bit_length_counts": {str(k): v for k, v in sorted(counts.items())},
                "measurement_seconds": measured["measurement_seconds"],
                "completion_invocation_seconds": perf_counter() - start,
                "heldout_evaluations": 0, "estimators_fitted": 0,
                "disclaimer": "Literature-reported IC50; no efficacy, safety or suitability claim."}

    return cache.get("report", complete)


def run_fingerprint_plans(raw_dir: Path, curated_dir: Path, split_dir: Path,
                          feature_dir: Path, cache_dir: Path, *,
                          splits: tuple[str, ...] = ("scaffold", "random")) -> dict:
    """Join checksum-verified curated structures onto ONLY saved training keys."""
    if not splits or len(set(splits)) != len(splits) or set(splits) - {"scaffold", "random"}:
        raise ValueError("Choose unique saved split names: scaffold and/or random")
    dest = Path(cache_dir).resolve()
    for source in (raw_dir, curated_dir, split_dir, feature_dir, DATA_RAW):
        source = Path(source).resolve()
        if dest == source or source in dest.parents or dest in source.parents:
            raise ValueError("Fingerprint plans must be separate from all input trees")
    inputs = load_training_inputs(raw_dir, curated_dir, split_dir, feature_dir)
    path = Path(curated_dir) / "compounds.csv"
    structures = {r["inchikey"]: r["canonical_smiles"] for r in _read_csv(path)}
    if any(_sha256(path) != inputs["splits"][s]["source_checksums"]["curated_compounds"]
           for s in splits):
        raise ValueError("Curated structures changed during fingerprint preparation")
    reports = {}
    for split in splits:
        training = inputs["splits"][split]
        rows = [{**r, "canonical_smiles": structures[r["inchikey"]]} for r in training["rows"]]
        reports[split] = plan_fingerprints(
            rows, dest / split / "fingerprint_plan",
            source_checksums=training["source_checksums"])
    return reports


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    parser.add_argument("--split", choices=("scaffold", "random", "both"), default="both")
    parser.add_argument("--cache-dir", type=Path)
    args = parser.parse_args(argv)
    target = args.target_id.strip().upper()
    reports = run_fingerprint_plans(
        DATA_RAW / target, DATA_CURATED / target, RESULTS / "splits" / target,
        DATA_FEATURES / target, args.cache_dir or RESULTS / "gate4_training" / target,
        splits=("scaffold", "random") if args.split == "both" else (args.split,))
    for split, report in reports.items():
        print(f"{split} fingerprint plan complete: {report['identity_sha256']}; "
              "training-only; no model or held-out evaluation")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
