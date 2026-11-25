"""Gate 2 split design and structural leakage audit; no model is fitted.

Run after curation: .venv/bin/python -m src.splits CHEMBL284
Only the checksum-verified curated compounds are used, never raw activity rows.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, DataStructs, rdBase
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

from src.config import DATA_CURATED, DATA_RAW, RESULTS, SEED
from src.retrieve import MANIFEST, _sha256, verify_raw

TEST_FRACTION = 0.2
NEAR_DUPLICATE_THRESHOLD = 0.90
FINGERPRINT_SPEC = {"type": "Morgan bit vector", "radius": 2, "fpSize": 2048,
                    "includeChirality": False}


def describe_structures(compounds: list[dict]) -> tuple[list[dict], list]:
    """Read already-standardised structures; fail rather than silently lose a row."""
    if len(compounds) < 2:
        raise ValueError("At least two curated compounds are needed")
    rows, fingerprints, seen = [], [], set()
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048,
                                                          includeChirality=False)
    for original in sorted(compounds, key=lambda r: r["inchikey"]):
        row = dict(original)
        key = row["inchikey"]
        if not key or key in seen:
            raise ValueError(f"Empty or duplicate curated InChIKey: {key!r}")
        seen.add(key)
        mol = Chem.MolFromSmiles(row["canonical_smiles"])
        if mol is None or Chem.MolToInchiKey(mol) != key:
            raise ValueError(f"Curated SMILES/InChIKey mismatch: {key}")
        if not np.isfinite(float(row["pIC50_median"])):
            raise ValueError(f"Nonfinite curated pIC50: {key}")
        # Stereo variants must not obtain separate scaffold groups. Atom/bond types
        # remain intact: this is the exact Murcko scaffold, not the generic skeleton.
        row["scaffold"] = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        row["connectivity_key"] = key.split("-")[0]
        rows.append(row)
        fingerprints.append(generator.GetFingerprint(mol))
    return rows, fingerprints


def make_partitions(rows: list[dict], seed: int = SEED,
                    test_fraction: float = TEST_FRACTION) -> dict[str, list[str]]:
    """Seeded, label-independent assignments aligned with the supplied row order.

    Shuffle scaffold groups, then include a complete group when it reduces the
    distance to the requested test count. No labels enter the allocation. The empty
    scaffold of acyclic molecules is ONE group. Singletons follow the same rule.
    """
    n = len(rows)
    if n < 2 or not 0 < test_fraction < 1:
        raise ValueError("Need >=2 rows and a test fraction strictly between 0 and 1")
    target = min(n - 1, max(1, int(round(n * test_fraction))))
    random_test = set(np.random.default_rng(seed).permutation(n)[:target].tolist())
    groups = defaultdict(list)
    for i, row in enumerate(rows):
        groups[row["scaffold"]].append(i)
    if len(groups) < 2:
        raise ValueError("Cannot hold out a scaffold with fewer than two scaffold groups")
    names = sorted(groups)
    order = np.random.default_rng(seed).permutation(len(names))
    scaffold_test = set()
    for position in order:
        indices = groups[names[position]]
        candidate = len(scaffold_test) + len(indices)
        if candidate < n and abs(candidate - target) < abs(len(scaffold_test) - target):
            scaffold_test.update(indices)
    if not scaffold_test:
        # All groups may exceed twice the target: keep the smallest complete group.
        smallest = min(names, key=lambda name: (len(groups[name]), name))
        scaffold_test.update(groups[smallest])
    return {name: ["test" if i in test else "train" for i in range(n)]
            for name, test in (("random", random_test), ("scaffold", scaffold_test))}


def assert_split_integrity(rows: list[dict], assignment: list[str],
                           scaffold_disjoint: bool = False) -> None:
    """Explicit failures remain active under python -O as well as in pytest."""
    if len(rows) != len(assignment) or set(assignment) != {"train", "test"}:
        raise AssertionError("Every compound must occur once, with nonempty train and test")
    train = [r for r, side in zip(rows, assignment) if side == "train"]
    test = [r for r, side in zip(rows, assignment) if side == "test"]
    overlap = {r["inchikey"] for r in train} & {r["inchikey"] for r in test}
    if overlap:
        raise AssertionError(f"InChIKey leakage across split: {sorted(overlap)}")
    if len({r["inchikey"] for r in rows}) != len(rows):
        raise AssertionError("Duplicate curated InChIKey within a partition")
    if scaffold_disjoint and ({r["scaffold"] for r in train} & {r["scaffold"] for r in test}):
        raise AssertionError("Scaffold leakage across scaffold split")


def _distribution(rows: list[dict]) -> dict:
    values = np.array([float(r["pIC50_median"]) for r in rows])
    return {"min": float(values.min()), "p05": float(np.quantile(values, .05)),
            "median": float(np.median(values)), "p95": float(np.quantile(values, .95)),
            "max": float(values.max()), "mean": float(values.mean())}


def audit_split(rows: list[dict], fingerprints: list, assignment: list[str],
                scaffold_disjoint: bool = False) -> tuple[dict, list[dict], list[dict]]:
    """All cross-boundary similarities, not a subsample; write only flagged pairs.

    Nonchiral fingerprints intentionally catch stereoisomers. The first InChIKey
    block is checked separately, since a bit-vector similarity is not an identity
    test. A near duplicate is an audit flag, not an assertion of identical chemistry.
    """
    assert_split_integrity(rows, assignment, scaffold_disjoint)
    if len(fingerprints) != len(rows):
        raise ValueError("Fingerprint count must match compound count")
    train_indices = [i for i, side in enumerate(assignment) if side == "train"]
    test_indices = [i for i, side in enumerate(assignment) if side == "test"]
    train_fp = [fingerprints[i] for i in train_indices]
    train_rows = [rows[i] for i in train_indices]
    test_rows = [rows[i] for i in test_indices]
    train_scaffolds = {r["scaffold"] for r in train_rows}
    test_scaffolds = {r["scaffold"] for r in test_rows}
    train_connectivity = {r["connectivity_key"] for r in train_rows}
    global_sizes = Counter(r["scaffold"] for r in rows)
    pairs, nearest = [], []
    connectivity_pairs = 0
    for i in test_indices:
        similarities = DataStructs.BulkTanimotoSimilarity(fingerprints[i], train_fp)
        best = int(np.argmax(similarities))
        nearest.append({"test_inchikey": rows[i]["inchikey"],
                        "nearest_train_inchikey": train_rows[best]["inchikey"],
                        "tanimoto": similarities[best]})
        for j, similarity in zip(train_indices, similarities):
            same_connectivity = rows[i]["connectivity_key"] == rows[j]["connectivity_key"]
            connectivity_pairs += int(same_connectivity)
            if similarity >= NEAR_DUPLICATE_THRESHOLD or same_connectivity:
                pairs.append({"test_inchikey": rows[i]["inchikey"],
                              "train_inchikey": rows[j]["inchikey"], "tanimoto": similarity,
                              "same_connectivity": same_connectivity,
                              "same_scaffold": rows[i]["scaffold"] == rows[j]["scaffold"]})
    maxima = np.array([r["tanimoto"] for r in nearest])
    summary = {
        "inchikey_overlap": 0,
        "shared_scaffolds": len(train_scaffolds & test_scaffolds),
        "test_compounds_with_train_scaffold": sum(r["scaffold"] in train_scaffolds for r in test_rows),
        "shared_connectivity_keys": len(train_connectivity & {r["connectivity_key"] for r in test_rows}),
        "same_connectivity_cross_boundary_pairs": connectivity_pairs,
        "near_duplicate_cross_boundary_pairs": sum(p["tanimoto"] >= NEAR_DUPLICATE_THRESHOLD for p in pairs),
        "test_compounds_with_near_duplicate": int((maxima >= NEAR_DUPLICATE_THRESHOLD).sum()),
        "test_nearest_train_tanimoto": {"min": float(maxima.min()),
                                        "median": float(np.median(maxima)),
                                        "p95": float(np.quantile(maxima, .95)),
                                        "max": float(maxima.max())},
    }
    for side, selected in (("train", train_rows), ("test", test_rows)):
        summary[side] = {
            "compounds": len(selected), "scaffolds": len({r["scaffold"] for r in selected}),
            "global_singleton_scaffolds": sum(global_sizes[r["scaffold"]] == 1 for r in selected),
            "acyclic_compounds": sum(r["scaffold"] == "" for r in selected),
            "pIC50": _distribution(selected),
        }
    low, high = summary["train"]["pIC50"]["min"], summary["train"]["pIC50"]["max"]
    summary["test_compounds_outside_train_pIC50_range"] = sum(
        not low <= float(r["pIC50_median"]) <= high for r in test_rows)
    return summary, pairs, nearest


def _write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def split_target(raw_dir: Path, curated_dir: Path, out_dir: Path) -> dict:
    raw_dir, curated_dir, out_dir = map(Path, (raw_dir, curated_dir, out_dir))
    # Reject writes into either source tree, including the immutable raw root.
    dest = out_dir.resolve()
    for source in (raw_dir.resolve(), curated_dir.resolve(), DATA_RAW.resolve()):
        if dest == source or source in dest.parents or dest in source.parents:
            raise ValueError("Split outputs must be separate from raw and curated inputs")
    if not verify_raw(raw_dir):
        raise ValueError("Raw checksum mismatch; split refused")
    manifest = json.loads((curated_dir / MANIFEST).read_text())
    raw_manifest_hash = _sha256(raw_dir / MANIFEST)
    if manifest["source"]["raw_manifest_sha256"] != raw_manifest_hash:
        raise ValueError("Curated dataset does not match the frozen raw manifest")
    for name, meta in manifest["files"].items():
        if _sha256(curated_dir / name) != meta["sha256"]:
            raise ValueError(f"Curated checksum mismatch: {name}")
    with (curated_dir / "compounds.csv").open(newline="", encoding="utf-8") as fh:
        compounds = list(csv.DictReader(fh))
    if len(compounds) != manifest["curated_compounds"]:
        raise ValueError("Curated compound count disagrees with manifest")
    rows, fingerprints = describe_structures(compounds)
    partitions = make_partitions(rows)
    report = {
        "target_chembl_id": manifest["target_chembl_id"], "seed": SEED,
        "command": f".venv/bin/python -m src.splits {manifest['target_chembl_id']}",
        "requested_test_fraction": TEST_FRACTION, "rdkit_version": rdBase.rdkitVersion,
        "raw_manifest_sha256": raw_manifest_hash,
        "curated_compounds_sha256": _sha256(curated_dir / "compounds.csv"),
        "fingerprint": FINGERPRINT_SPEC, "near_duplicate_threshold": NEAR_DUPLICATE_THRESHOLD,
        "scaffold_backend": "RDKit native MurckoDecompose",
        "scaffold_backend_source_sha256": _sha256(Path(MurckoScaffold.__file__)),
        "scaffold_rule": "Exact Bemis-Murcko atom/bond scaffold without stereochemistry; "
                         "all acyclic compounds share the empty scaffold. Seeded group order; "
                         "greedily approach 20% test without splitting groups. Singletons "
                         "follow the same allocation rule. No pIC50 used for assignment.",
        "similarity_cluster_decision": "Deferred: exact scaffold holdout directly tests "
            "analogue-series reuse, while the exhaustive similarity audit measures what it "
            "misses. A similarity-cluster split would answer the additional question of "
            "transfer beyond a chosen fingerprint threshold, which is outside this "
            "random-versus-scaffold leakage comparison. Scaffold holdout is not claimed "
            "to guarantee low similarity; residual pairs are retained and reported.",
        "splits": {},
        "disclaimer": "Literature/patent-reported IC50; no efficacy, safety or suitability claim.",
    }
    artifacts = {}
    for name, assignment in partitions.items():
        summary, pairs, nearest = audit_split(rows, fingerprints, assignment, name == "scaffold")
        report["splits"][name] = summary
        artifacts[f"{name}_assignments.csv"] = ([
            {"inchikey": r["inchikey"], "partition": side, "scaffold": r["scaffold"],
             "pIC50_median": r["pIC50_median"]} for r, side in zip(rows, assignment)],
            ["inchikey", "partition", "scaffold", "pIC50_median"])
        artifacts[f"{name}_near_pairs.csv"] = (pairs, ["test_inchikey", "train_inchikey",
                                                        "tanimoto", "same_connectivity", "same_scaffold"])
        artifacts[f"{name}_nearest_train.csv"] = (nearest, ["test_inchikey", "nearest_train_inchikey", "tanimoto"])
    if not verify_raw(raw_dir):
        raise ValueError("Raw files changed during splitting")
    if _sha256(curated_dir / "compounds.csv") != report["curated_compounds_sha256"]:
        raise ValueError("Curated compounds changed during splitting")
    out_dir.mkdir(parents=True, exist_ok=True)
    for name, (records, fields) in artifacts.items():
        _write_csv(out_dir / name, records, fields)
    report["raw_verified_before_and_after"] = True
    report["files"] = {name: {"sha256": _sha256(out_dir / name)} for name in artifacts}
    (out_dir / "audit.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    write_markdown(report, out_dir / "audit.md")
    return report


def write_markdown(report: dict, path: Path) -> None:
    lines = [f"# Structural leakage audit: {report['target_chembl_id']}", "",
             f"Command: `{report['command']}`; seed {report['seed']}. No model fitted.", "",
             report["scaffold_rule"], "", f"Backend: {report['scaffold_backend']}; "
             f"RDKit {report['rdkit_version']}. The wrapper source checksum is logged in audit.json.",
             "", report["similarity_cluster_decision"], "",
             "| Property | Random | Scaffold |", "|---|---:|---:|"]
    a, b = (report["splits"][name] for name in ("random", "scaffold"))
    for label, key in (("InChIKey overlap (asserted zero)", "inchikey_overlap"),
                       ("Shared scaffolds", "shared_scaffolds"),
                       ("Test compounds with a training scaffold", "test_compounds_with_train_scaffold"),
                       ("Shared connectivity keys", "shared_connectivity_keys"),
                       ("Cross-boundary pairs with Tanimoto >= 0.90", "near_duplicate_cross_boundary_pairs"),
                       ("Test compounds with Tanimoto >= 0.90 neighbour", "test_compounds_with_near_duplicate"),
                       ("Test compounds outside training pIC50 range", "test_compounds_outside_train_pIC50_range")):
        lines.append(f"| {label} | {a[key]} | {b[key]} |")
    for side in ("train", "test"):
        for key in ("compounds", "scaffolds", "global_singleton_scaffolds", "acyclic_compounds"):
            lines.append(f"| {side}: {key} | {a[side][key]} | {b[side][key]} |")
        for key in ("min", "p05", "median", "p95", "max"):
            lines.append(f"| {side}: pIC50 {key} | {a[side]['pIC50'][key]:.3f} | {b[side]['pIC50'][key]:.3f} |")
    lines += ["", "Similarity: nonchiral Morgan radius 2, 2048 bits; exhaustive train/test "
              "Tanimoto comparisons. The 0.90 threshold is a fixed audit flag, not proof of "
              "chemical identity. Each test compound's nearest training neighbour and every "
              "flagged pair (including shared InChIKey connectivity blocks) are saved in CSV.", "",
              "The random split can reuse shared cores and close analogues despite zero exact "
              "compound overlap. The scaffold split blocks exact cores but can retain similar "
              "ring substitutions or fingerprint collisions. Counts above measure those risks; "
              "the fixed probe and its uncertainty remain a separate Gate 2 deliverable.", "",
              "pIC50 summaries diagnose distribution shift; they never determine assignments. "
              "No split is retried to improve target coverage or a score.", "", report["disclaimer"], ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    args = parser.parse_args(argv)
    target = args.target_id.strip().upper()
    report = split_target(DATA_RAW / target, DATA_CURATED / target, RESULTS / "splits" / target)
    for name, audit in report["splits"].items():
        print(f"{name}: {audit['train']['compounds']} train / {audit['test']['compounds']} test; "
              f"{audit['inchikey_overlap']} exact overlap; "
              f"{audit['test_compounds_with_near_duplicate']} test compounds with >=0.90 neighbour")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
