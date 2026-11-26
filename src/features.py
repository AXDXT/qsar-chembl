"""Gate 3 featurisation: a deterministic descriptor and Morgan cache. No model is fitted.

Run after src.splits: .venv/bin/python -m src.features CHEMBL284

Both feature blocks are computed only from the standardised structures in
data/curated/<TARGET>/compounds.csv, never from raw SMILES, and never from any label.
Featurisation is stateless: every row is a pure function of one standardised molecule,
so there is nothing that could be fitted on the full dataset or on held-out rows.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import Crippen, Descriptors, rdFingerprintGenerator, rdMolDescriptors

from src.config import DATA_CURATED, DATA_FEATURES, DATA_RAW, RESULTS, SEED
from src.retrieve import MANIFEST, _sha256, verify_raw
from src.splits import FINGERPRINT_SPEC, _write_csv, describe_structures

# Fixed, declared, untuned. Twelve interpretable 2D physicochemical descriptors: the
# Lipinski rule-of-five terms, the Veber flexibility/polarity terms, and size,
# polarisability, ring and saturation counts. Every column is a monotone, nameable
# property so a Gate 4 chemical interpretation can discuss it without a surrogate
# model. No 3D descriptor is included: those need conformer generation, which adds a
# stochastic embedding step this gate has no reason to introduce.
DESCRIPTORS = (
    ("MolWt", Descriptors.MolWt,
     "Average molecular weight; the size term of Lipinski's rule of five."),
    ("HeavyAtomCount", rdMolDescriptors.CalcNumHeavyAtoms,
     "Size as an integer atom count, free of hydrogen-weight and isotope conventions."),
    ("MolLogP", Crippen.MolLogP,
     "Crippen atom-contribution octanol/water logP; the lipophilicity term of the rule of five."),
    ("MolMR", Crippen.MolMR,
     "Crippen molar refractivity; polarisability and bulk, the Ghose-filter companion to logP."),
    ("TPSA", rdMolDescriptors.CalcTPSA,
     "Topological polar surface area; the polarity term of Veber's permeability rules."),
    ("NumLipinskiHBD", rdMolDescriptors.CalcNumLipinskiHBD,
     "Lipinski hydrogen-bond donors (N-H and O-H count), the rule-of-five donor definition."),
    ("NumLipinskiHBA", rdMolDescriptors.CalcNumLipinskiHBA,
     "Lipinski hydrogen-bond acceptors (nitrogen and oxygen count), the rule-of-five definition."),
    ("NumRotatableBonds", rdMolDescriptors.CalcNumRotatableBonds,
     "Conformational flexibility; the Veber rotatable-bond term."),
    ("RingCount", rdMolDescriptors.CalcNumRings,
     "Number of SSSR rings; total ring content."),
    ("NumAromaticRings", rdMolDescriptors.CalcNumAromaticRings,
     "Aromatic ring count; aromatic versus saturated ring character."),
    ("FractionCSP3", rdMolDescriptors.CalcFractionCSP3,
     "Fraction of sp3-hybridised carbons; saturation and three-dimensional character."),
    ("NumHeteroatoms", rdMolDescriptors.CalcNumHeteroatoms,
     "Heteroatom count; polar functionality that TPSA alone does not separate."),
)
DESCRIPTOR_NAMES = tuple(name for name, _, _ in DESCRIPTORS)

# Alternative foldings are measured for the record only. The bit length is NOT selected
# by any score: 2,048 is inherited from the Gate 2 similarity audit and leakage probe so
# that the audited near-duplicate structure describes the model's actual input space.
FOLDING_SENSITIVITY_BITS = (1024, 2048, 4096)

DESCRIPTOR_BLOCK_NOTES = {
    "nan_inf_policy": "Every descriptor value is checked with numpy.isfinite. A non-finite "
        "value is written to featurisation_failures.csv and aborts the run. Nothing is "
        "imputed, clipped or dropped, so no imputation statistic exists that could have "
        "been fitted on the full dataset.",
    "scaling_policy": "None. No centring, scaling, variance filtering or feature selection "
        "is applied at Gate 3. The descriptor baseline is a tree ensemble, which is "
        "invariant to monotone per-column rescaling, so a scaler would add a fitted "
        "object with no benefit. Any transform a later gate adds must be fitted on "
        "training rows only.",
    "statelessness": "Provably stateless: each row is a function of one molecule alone. "
        "No row, label, split or dataset-level statistic enters another row's features. "
        "test_features.py checks this by featurising single molecules in isolation.",
    "omitted_aliphatic_ring_count": "NumAliphaticRings is deliberately omitted: RDKit "
        "partitions SSSR rings into aromatic and aliphatic, so it is exactly "
        "RingCount - NumAromaticRings. The identity is measured, not assumed; see "
        "collinearity.aliphatic_ring_identity_max_deviation.",
}


def descriptor_values(mol: Chem.Mol) -> list[float]:
    """Twelve descriptors for one standardised molecule, in declared column order."""
    return [float(function(mol)) for _, function, _ in DESCRIPTORS]


def fingerprint_matrix(fingerprints: list) -> np.ndarray:
    """Dense 0/1 uint8 matrix from the same bit vectors the Gate 2 audit compared."""
    bits = FINGERPRINT_SPEC["fpSize"]
    matrix = np.zeros((len(fingerprints), bits), dtype=np.uint8)
    for i, fingerprint in enumerate(fingerprints):
        if fingerprint.GetNumBits() != bits:
            raise ValueError(f"Fingerprint {i} has {fingerprint.GetNumBits()} bits, expected {bits}")
        for bit in fingerprint.GetOnBits():
            matrix[i, bit] = 1
    return matrix


def folding_evidence(molecules: list[Chem.Mol]) -> dict:
    """Measure how much radius-2 substructure detail each bit length can keep.

    Counts distinct unfolded Morgan identifiers against occupied bits after folding.
    This records the cost of the inherited 2,048-bit choice; it selects nothing.
    """
    radius, chirality = FINGERPRINT_SPEC["radius"], FINGERPRINT_SPEC["includeChirality"]
    unfolded = rdFingerprintGenerator.GetMorganGenerator(radius=radius, includeChirality=chirality)
    identifiers = [set(unfolded.GetSparseCountFingerprint(m).GetNonzeroElements()) for m in molecules]
    evidence = {
        "radius": radius,
        "include_chirality": chirality,
        "distinct_substructure_identifiers": len(set().union(*identifiers)),
        "mean_identifiers_per_compound": float(np.mean([len(s) for s in identifiers])),
        "by_bit_length": {},
    }
    for bits in FOLDING_SENSITIVITY_BITS:
        generator = rdFingerprintGenerator.GetMorganGenerator(
            radius=radius, fpSize=bits, includeChirality=chirality)
        on_bits = [set(generator.GetFingerprint(m).GetOnBits()) for m in molecules]
        evidence["by_bit_length"][str(bits)] = {
            "occupied_bits": len(set().union(*on_bits)),
            "mean_on_bits_per_compound": float(np.mean([len(s) for s in on_bits])),
            "compounds_losing_at_least_one_identifier": sum(
                len(ids) > len(bitset) for ids, bitset in zip(identifiers, on_bits)),
        }
    return evidence


def collinearity(matrix: np.ndarray, names: tuple[str, ...]) -> dict:
    """Largest absolute Pearson correlation between descriptor columns.

    Reported for the training rows of each saved split, never for test rows, and used
    for nothing: no column is dropped, combined or reweighted on the strength of it.
    """
    if matrix.ndim != 2 or matrix.shape[1] != len(names):
        raise ValueError("Collinearity needs one column per descriptor name")
    deviation = matrix - matrix.mean(axis=0)
    scale = np.sqrt(np.square(deviation).sum(axis=0))
    constant = [name for name, value in zip(names, scale) if value == 0]
    keep = [i for i, value in enumerate(scale) if value > 0]
    normalised = deviation[:, keep] / scale[keep]
    correlation = normalised.T @ normalised
    np.fill_diagonal(correlation, 0.0)
    flat = int(np.argmax(np.abs(correlation)))
    a, b = divmod(flat, len(keep))
    return {"rows": int(matrix.shape[0]), "constant_columns": constant,
            "max_abs_pearson_r": float(correlation[a, b]),
            "max_abs_pearson_pair": sorted((names[keep[a]], names[keep[b]]))}


def _read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def load_features(feature_dir: Path) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Read a cache back: InChIKeys in row order, descriptors, unpacked fingerprints."""
    feature_dir = Path(feature_dir)
    manifest = json.loads((feature_dir / MANIFEST).read_text())
    for name, meta in manifest["files"].items():
        if _sha256(feature_dir / name) != meta["sha256"]:
            raise ValueError(f"Feature cache checksum mismatch: {name}")
    rows = _read_csv(feature_dir / "descriptors.csv")
    names = tuple(manifest["descriptor_block"]["names"])
    if names != DESCRIPTOR_NAMES:
        raise ValueError("Cached descriptor columns differ from the declared block")
    descriptors = np.array([[float(row[name]) for name in names] for row in rows])
    packed = np.load(feature_dir / "fingerprint_bits.npy")
    bits = np.unpackbits(packed, axis=1)[:, :manifest["fingerprint_block"]["fpSize"]]
    if descriptors.shape[0] != bits.shape[0] != len(rows):
        raise ValueError("Cached feature blocks disagree on row count")
    return [row["inchikey"] for row in rows], descriptors, bits


def build_features(raw_dir: Path, curated_dir: Path, split_dir: Path, out_dir: Path) -> dict:
    raw_dir, curated_dir, split_dir, out_dir = map(Path, (raw_dir, curated_dir, split_dir, out_dir))
    dest = out_dir.resolve()
    for source in (raw_dir.resolve(), curated_dir.resolve(), split_dir.resolve(), DATA_RAW.resolve()):
        if dest == source or source in dest.parents or dest in source.parents:
            raise ValueError("Feature outputs must be separate from raw, curated and split inputs")
    if not verify_raw(raw_dir):
        raise ValueError("Raw checksum mismatch; featurisation refused")
    manifest = json.loads((curated_dir / MANIFEST).read_text())
    audit = json.loads((split_dir / "audit.json").read_text())
    raw_hash, compound_hash = _sha256(raw_dir / MANIFEST), _sha256(curated_dir / "compounds.csv")
    if (manifest["source"]["raw_manifest_sha256"] != raw_hash
            or audit["raw_manifest_sha256"] != raw_hash
            or audit["curated_compounds_sha256"] != compound_hash):
        raise ValueError("Feature inputs disagree on raw/curated provenance")
    for directory, files in ((curated_dir, manifest["files"]), (split_dir, audit["files"])):
        for name, meta in files.items():
            if _sha256(directory / name) != meta["sha256"]:
                raise ValueError(f"Input checksum mismatch: {name}")
    if audit["fingerprint"] != FINGERPRINT_SPEC:
        raise ValueError("Saved split audit used a different fingerprint specification")
    sources = [raw_dir / MANIFEST, curated_dir / MANIFEST, split_dir / "audit.json"]
    sources += [curated_dir / name for name in manifest["files"]]
    sources += [split_dir / name for name in audit["files"]]
    before = {path: _sha256(path) for path in sources}

    compounds = _read_csv(curated_dir / "compounds.csv")
    if len(compounds) != manifest["curated_compounds"]:
        raise ValueError("Curated count disagrees with manifest")
    # describe_structures re-reads the standardised SMILES, re-derives the InChIKey and
    # refuses any mismatch, so the features below cannot come from an unstandardised
    # structure. It also fixes the row order and builds the audited fingerprints.
    rows, fingerprints = describe_structures(compounds)
    molecules = [Chem.MolFromSmiles(row["canonical_smiles"]) for row in rows]

    failures, values = [], []
    for row, mol in zip(rows, molecules):
        try:
            if mol is None:
                raise ValueError("standardised SMILES did not parse")
            computed = descriptor_values(mol)
            nonfinite = [name for name, value in zip(DESCRIPTOR_NAMES, computed)
                         if not np.isfinite(value)]
            if nonfinite:
                raise ValueError(f"non-finite descriptors: {','.join(nonfinite)}")
        except Exception as error:  # reported, never silently dropped
            failures.append({"inchikey": row["inchikey"],
                             "canonical_smiles": row["canonical_smiles"],
                             "reason": f"{type(error).__name__}: {error}"})
            continue
        values.append(computed)
    out_dir.mkdir(parents=True, exist_ok=True)
    _write_csv(out_dir / "featurisation_failures.csv", failures,
               ["inchikey", "canonical_smiles", "reason"])
    if failures:
        raise ValueError(f"{len(failures)} compounds failed featurisation; see "
                         f"{out_dir / 'featurisation_failures.csv'}. Count must reconcile "
                         f"with the {manifest['curated_compounds']} curated compounds.")
    descriptors = np.array(values, dtype=np.float64)
    bits = fingerprint_matrix(fingerprints)
    if descriptors.shape != (len(rows), len(DESCRIPTOR_NAMES)) or bits.shape[0] != len(rows):
        raise ValueError("Feature matrices do not reconcile with the curated compound count")

    order = [row["inchikey"] for row in rows]
    alignment = {}
    for name in ("random", "scaffold"):
        saved = [r["inchikey"] for r in _read_csv(split_dir / f"{name}_assignments.csv")]
        if saved != order:
            raise ValueError(f"Saved {name} split row order differs from the featurised order; "
                             "the split must not be regenerated or reordered")
        alignment[name] = True
    train_collinearity = {}
    for name in ("random", "scaffold"):
        assignment = _read_csv(split_dir / f"{name}_assignments.csv")
        selected = [i for i, r in enumerate(assignment) if r["partition"] == "train"]
        train_collinearity[name] = collinearity(descriptors[selected], DESCRIPTOR_NAMES)

    aliphatic = np.array([float(rdMolDescriptors.CalcNumAliphaticRings(m)) for m in molecules])
    ring = descriptors[:, DESCRIPTOR_NAMES.index("RingCount")]
    aromatic = descriptors[:, DESCRIPTOR_NAMES.index("NumAromaticRings")]

    report = {
        "target_chembl_id": manifest["target_chembl_id"],
        "command": f".venv/bin/python -m src.features {manifest['target_chembl_id']}",
        "purpose": "Gate 3 feature cache. No model is fitted, scored or selected here.",
        "seed": SEED,
        "seed_role": "Recorded for provenance only. Featurisation draws no random numbers, "
                     "so the cache is byte-identical on every rerun.",
        "determinism": "No run timestamp is recorded, which is what lets a rerun be "
                       "byte-identical. Provenance is carried by the input checksums below.",
        "rdkit_version": rdBase.rdkitVersion,
        "numpy_version": np.__version__,
        "source": {
            "curated_dir": f"data/curated/{manifest['target_chembl_id']}",
            "structures": "canonical_smiles from compounds.csv: sanitised, largest fragment, "
                          "neutralised and InChIKey-verified by Gate 2. Raw SMILES are never read.",
            "raw_manifest_sha256": raw_hash,
            "curated_compounds_sha256": compound_hash,
            "split_audit_sha256": before[split_dir / "audit.json"],
        },
        "reconciliation": {
            "curated_compounds": manifest["curated_compounds"],
            "rows_featurised": int(descriptors.shape[0]),
            "featurisation_failures": len(failures),
            "non_finite_descriptor_values": 0,
            "compounds_dropped": manifest["curated_compounds"] - int(descriptors.shape[0]),
        },
        "row_order": "Ascending standardised full InChIKey, identical to the saved split "
                     "assignment CSVs, so the cached matrices join to them positionally.",
        "split_row_order_matches_saved_assignments": alignment,
        "descriptor_block": {
            "count": len(DESCRIPTOR_NAMES),
            "names": list(DESCRIPTOR_NAMES),
            "shape": list(descriptors.shape),
            "rationale": {name: why for name, _, why in DESCRIPTORS},
            "selection_rule": "Fixed and declared before any model is fitted. Chosen for "
                              "named, monotone interpretability, not by any score.",
            **DESCRIPTOR_BLOCK_NOTES,
        },
        "fingerprint_block": {
            **FINGERPRINT_SPEC,
            "shape": list(bits.shape),
            "storage": "numpy.packbits along the bit axis; load_features unpacks it.",
            "bit_length_justification": "2,048 bits is inherited unchanged from the Gate 2 "
                "similarity audit and leakage probe, so the audited near-duplicate and "
                "shared-core counts describe this model's actual input space. It was not "
                "chosen by any metric; folding_evidence records what it costs.",
            "folding_evidence": folding_evidence(molecules),
        },
        "collinearity": {
            "train_rows_only": train_collinearity,
            "note": "Computed on each split's training rows only, with no label involved, "
                    "and used for nothing. Reported so a Gate 4 interpretation cannot "
                    "mistake two near-duplicate columns for independent evidence.",
            "aliphatic_ring_identity_max_deviation": float(
                np.abs(ring - aromatic - aliphatic).max()),
        },
        "disclaimer": "Literature/patent-reported IC50; no efficacy, safety or suitability claim.",
    }
    if not verify_raw(raw_dir) or any(_sha256(path) != digest for path, digest in before.items()):
        raise ValueError("Feature inputs changed during featurisation")
    _write_csv(out_dir / "descriptors.csv",
               [{"inchikey": key, **dict(zip(DESCRIPTOR_NAMES, row))}
                for key, row in zip(order, descriptors.tolist())],
               ["inchikey", *DESCRIPTOR_NAMES])
    np.save(out_dir / "fingerprint_bits.npy", np.packbits(bits, axis=1), allow_pickle=False)
    report["inputs_verified_before_and_after"] = True
    report["files"] = {name: {"sha256": _sha256(out_dir / name)} for name in
                       ("descriptors.csv", "fingerprint_bits.npy", "featurisation_failures.csv")}
    (out_dir / MANIFEST).write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                                    encoding="utf-8")
    write_markdown(report, out_dir / "features.md")
    return report


def write_markdown(report: dict, path: Path) -> None:
    descriptor, fingerprint = report["descriptor_block"], report["fingerprint_block"]
    reconciliation = report["reconciliation"]
    lines = [f"# Gate 3 feature cache: {report['target_chembl_id']}", "",
             f"Command: `{report['command']}`. {report['purpose']}", "",
             report["determinism"], "", report["source"]["structures"], "",
             "| Reconciliation | Count |", "|---|---:|",
             f"| Curated compounds in | {reconciliation['curated_compounds']:,} |",
             f"| Rows featurised | {reconciliation['rows_featurised']:,} |",
             f"| Featurisation failures | {reconciliation['featurisation_failures']} |",
             f"| Non-finite descriptor values | {reconciliation['non_finite_descriptor_values']} |",
             f"| Compounds dropped | {reconciliation['compounds_dropped']} |", "",
             f"Every failure, if any, is listed in `featurisation_failures.csv`; the run aborts "
             f"rather than continuing with a short matrix. {report['row_order']}", "",
             f"## Descriptor block ({descriptor['count']} columns)", "",
             descriptor["selection_rule"], "", "| Descriptor | Why it is in the set |", "|---|---|"]
    lines += [f"| {name} | {why} |" for name, why in descriptor["rationale"].items()]
    lines += ["", descriptor["omitted_aliphatic_ring_count"],
              f"Measured maximum deviation: "
              f"{report['collinearity']['aliphatic_ring_identity_max_deviation']:.1f}.", "",
              f"NaN/inf: {descriptor['nan_inf_policy']}", "",
              f"Scaling: {descriptor['scaling_policy']}", "",
              f"Statelessness: {descriptor['statelessness']}", "",
              "| Split (training rows only) | Rows | Largest absolute correlation | Pair |",
              "|---|---:|---:|---|"]
    for name, result in report["collinearity"]["train_rows_only"].items():
        lines.append(f"| {name} | {result['rows']:,} | {result['max_abs_pearson_r']:.3f} | "
                     f"{' / '.join(result['max_abs_pearson_pair'])} |")
    lines += ["", report["collinearity"]["note"], "",
              f"## Fingerprint block ({fingerprint['fpSize']:,} bits)", "",
              f"{fingerprint['type']}, radius {fingerprint['radius']}, "
              f"includeChirality={fingerprint['includeChirality']}, via the current "
              f"rdFingerprintGenerator API. {fingerprint['bit_length_justification']}", "",
              "| Bit length | Occupied bits | Mean on-bits per compound | "
              "Compounds losing >=1 identifier |", "|---:|---:|---:|---:|"]
    evidence = fingerprint["folding_evidence"]
    for bits, result in evidence["by_bit_length"].items():
        lines.append(f"| {int(bits):,} | {result['occupied_bits']:,} | "
                     f"{result['mean_on_bits_per_compound']:.1f} | "
                     f"{result['compounds_losing_at_least_one_identifier']:,} |")
    lines += ["", f"The curated set contains {evidence['distinct_substructure_identifiers']:,} "
              f"distinct radius-<={evidence['radius']} substructure identifiers, "
              f"{evidence['mean_identifiers_per_compound']:.1f} per compound on average. "
              "Folding maps several identifiers onto one bit; the final column counts the "
              "compounds where that happens. The alternative lengths are measured for the "
              "record and were not scored, so nothing here selects a bit length.", "",
              "Neither block distinguishes stereoisomers: the fingerprint ignores chirality "
              "by Gate 2's choice, and no descriptor counts stereocentres. 286 curated "
              "compounds share a connectivity skeleton with another compound, so that "
              "limitation is measurable rather than hypothetical.", "",
              report["disclaimer"], ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    target = parser.parse_args(argv).target_id.strip().upper()
    report = build_features(DATA_RAW / target, DATA_CURATED / target,
                            RESULTS / "splits" / target, DATA_FEATURES / target)
    reconciliation = report["reconciliation"]
    print(f"featurised {reconciliation['rows_featurised']} of "
          f"{reconciliation['curated_compounds']} curated compounds; "
          f"{reconciliation['featurisation_failures']} failures; "
          f"descriptors {tuple(report['descriptor_block']['shape'])}, "
          f"fingerprints {tuple(report['fingerprint_block']['shape'])}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
