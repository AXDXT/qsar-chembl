"""Offline Gate 3 featurisation checks.

These target the mistakes that would actually invalidate the baselines: features taken
from raw instead of standardised structures, a statistic computed across the whole
dataset or on held-out rows, a silently regenerated or reordered split, a molecule
dropped without a report, and a cache that is not reproducible.
"""
import csv
import json

import numpy as np
import pytest
from rdkit import Chem

from src import curate, features, retrieve, splits
from src.retrieve import read_jsonl_gz
from test_splits import _curated_fixture


def _inputs(tmp_path):
    raw, curated = _curated_fixture(tmp_path)
    split_dir = tmp_path / "splits"
    splits.split_target(raw, curated, split_dir)
    return raw, curated, split_dir


def _read_csv(path):
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


CACHE = features.DATA_FEATURES / "CHEMBL284"
CURATED = features.DATA_CURATED / "CHEMBL284"
SPLIT_DIR = features.RESULTS / "splits" / "CHEMBL284"


def test_descriptor_block_is_declared_fixed_and_free_of_the_exact_redundant_column():
    assert features.DESCRIPTOR_NAMES == tuple(name for name, _, _ in features.DESCRIPTORS)
    assert len(set(features.DESCRIPTOR_NAMES)) == len(features.DESCRIPTOR_NAMES) == 12
    assert all(why.endswith(".") for _, _, why in features.DESCRIPTORS)
    assert "NumAliphaticRings" not in features.DESCRIPTOR_NAMES
    # The omission claim is an identity, so it must hold for every curated structure.
    manifest = json.loads((CACHE / retrieve.MANIFEST).read_text())
    assert manifest["collinearity"]["aliphatic_ring_identity_max_deviation"] == 0.0


def test_fingerprint_spec_is_the_audited_gate2_spec_not_a_redefinition():
    audit = json.loads((SPLIT_DIR / "audit.json").read_text())
    assert features.FINGERPRINT_SPEC is splits.FINGERPRINT_SPEC
    assert audit["fingerprint"] == splits.FINGERPRINT_SPEC
    cached = json.loads((CACHE / retrieve.MANIFEST).read_text())["fingerprint_block"]
    assert {key: cached[key] for key in splits.FINGERPRINT_SPEC} == splits.FINGERPRINT_SPEC


@pytest.mark.parametrize("smiles", ["CCO", "c1ccccc1C(=O)NC1CCCCC1", "N[C@@H](C)C(=O)O"])
def test_descriptors_are_stateless_and_identical_molecule_by_molecule(smiles):
    """A row must not depend on which other molecules are present in the batch."""
    alone = features.descriptor_values(Chem.MolFromSmiles(smiles))
    padding = ["C", "CCCCCCCCCCCC", "c1ccc2ccccc2c1"]
    batch = [features.descriptor_values(Chem.MolFromSmiles(s)) for s in [*padding, smiles]]
    assert batch[-1] == alone
    assert features.descriptor_values(Chem.MolFromSmiles(smiles)) == alone


def test_cached_dpp4_descriptors_recompute_from_the_standardised_structures():
    compounds = {row["inchikey"]: row["canonical_smiles"]
                 for row in _read_csv(CURATED / "compounds.csv")}
    cached = _read_csv(CACHE / "descriptors.csv")
    assert len(cached) == len(compounds) == 3965
    for row in cached[::200]:
        expected = features.descriptor_values(Chem.MolFromSmiles(compounds[row["inchikey"]]))
        assert [float(row[name]) for name in features.DESCRIPTOR_NAMES] == expected


def test_features_come_from_standardised_structures_not_raw_smiles():
    """A salt's raw SMILES gives different descriptors; the cache must hold the parent's."""
    records = _read_csv(CURATED / "records_curated.csv")
    salts = {r["molecule_chembl_id"]: r for r in records if int(r["fragments_before"]) > 1}
    assert salts, "fixture expects curated records that were desalted at Gate 2"
    raw_smiles = {row["molecule_chembl_id"]: row["canonical_smiles"]
                  for row in read_jsonl_gz(features.DATA_RAW / "CHEMBL284" / "activities.jsonl.gz")
                  if row["molecule_chembl_id"] in salts}
    cached = {row["inchikey"]: row for row in _read_csv(CACHE / "descriptors.csv")}
    compared = 0
    for molecule_id, record in list(salts.items())[:25]:
        raw, standardised = raw_smiles[molecule_id], record["canonical_smiles_standardised"]
        assert "." in raw, "a desalted record's raw SMILES should still be a mixture"
        raw_row = features.descriptor_values(Chem.MolFromSmiles(raw))
        standardised_row = features.descriptor_values(Chem.MolFromSmiles(standardised))
        assert raw_row != standardised_row
        assert [float(cached[record["inchikey"]][n]) for n in features.DESCRIPTOR_NAMES] \
            == standardised_row
        compared += 1
    assert compared >= 10


def test_cached_row_order_matches_both_saved_split_assignments_positionally():
    order = [row["inchikey"] for row in _read_csv(CACHE / "descriptors.csv")]
    assert order == sorted(order)
    for name in ("random", "scaffold"):
        saved = [row["inchikey"] for row in _read_csv(SPLIT_DIR / f"{name}_assignments.csv")]
        assert saved == order


def test_cached_collinearity_is_measured_on_training_rows_only():
    manifest = json.loads((CACHE / retrieve.MANIFEST).read_text())
    audit = json.loads((SPLIT_DIR / "audit.json").read_text())
    for name, result in manifest["collinearity"]["train_rows_only"].items():
        assert result["rows"] == audit["splits"][name]["train"]["compounds"]
        assert result["rows"] < manifest["reconciliation"]["rows_featurised"]


def test_cached_fingerprints_equal_the_gate2_bit_vectors_and_load_back_as_bits():
    keys, descriptors, bits = features.load_features(CACHE)
    assert descriptors.shape == (3965, 12) and bits.shape == (3965, 2048)
    assert set(np.unique(bits)) <= {0, 1}
    assert np.isfinite(descriptors).all()
    rows, fingerprints = splits.describe_structures(_read_csv(CURATED / "compounds.csv"))
    assert keys == [row["inchikey"] for row in rows]
    for i in range(0, len(rows), 400):
        assert sorted(np.flatnonzero(bits[i]).tolist()) == sorted(fingerprints[i].GetOnBits())


def test_load_features_refuses_a_tampered_cache(tmp_path):
    for name in ("descriptors.csv", "fingerprint_bits.npy", "featurisation_failures.csv",
                 retrieve.MANIFEST):
        (tmp_path / name).write_bytes((CACHE / name).read_bytes())
    assert features.load_features(tmp_path)[0]
    with (tmp_path / "descriptors.csv").open("a", encoding="utf-8") as fh:
        fh.write("tampered\n")
    with pytest.raises(ValueError, match="checksum mismatch"):
        features.load_features(tmp_path)


def test_collinearity_reports_the_extreme_pair_and_names_constant_columns():
    matrix = np.array([[1., 2., 5., 0.], [2., 4., 1., 0.], [3., 6.1, 9., 0.]])
    result = features.collinearity(matrix, ("a", "b", "c", "flat"))
    assert result["constant_columns"] == ["flat"]
    assert result["max_abs_pearson_pair"] == ["a", "b"]
    assert result["max_abs_pearson_r"] == pytest.approx(1.0, abs=1e-3)
    assert result["rows"] == 3
    negated = features.collinearity(matrix * np.array([1., -1., 1., 1.]), ("a", "b", "c", "flat"))
    assert negated["max_abs_pearson_r"] == pytest.approx(-result["max_abs_pearson_r"], abs=1e-9)


def test_folding_evidence_measures_collisions_without_choosing_a_bit_length():
    molecules = [Chem.MolFromSmiles(s) for s in
                 ("CC(=O)Oc1ccccc1C(=O)O", "c1ccc2ccccc2c1", "CCCCCCCCCCCCN")]
    evidence = features.folding_evidence(molecules)
    assert evidence["radius"] == splits.FINGERPRINT_SPEC["radius"]
    assert set(evidence["by_bit_length"]) == {"1024", "2048", "4096"}
    identifiers = evidence["distinct_substructure_identifiers"]
    occupancies = [result["occupied_bits"] for result in evidence["by_bit_length"].values()]
    assert occupancies == sorted(occupancies)  # longer vectors cannot occupy fewer bits
    assert all(occupied <= identifiers for occupied in occupancies)
    assert features.folding_evidence(molecules) == evidence


def test_end_to_end_cache_is_byte_reproducible_and_leaves_every_input_untouched(tmp_path):
    raw, curated, split_dir = _inputs(tmp_path)
    before = {p: p.read_bytes() for directory in (raw, curated, split_dir)
              for p in directory.iterdir() if not p.name.startswith("._")}
    out = tmp_path / "features"
    report = features.build_features(raw, curated, split_dir, out)
    first = {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert report == features.build_features(raw, curated, split_dir, out)
    assert first == {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert all(p.read_bytes() == contents for p, contents in before.items())
    for name, meta in report["files"].items():
        assert retrieve._sha256(out / name) == meta["sha256"]
    assert retrieve.MANIFEST not in report["files"]  # the manifest cannot hash itself
    assert report["reconciliation"]["featurisation_failures"] == 0
    assert report["reconciliation"]["compounds_dropped"] == 0
    assert report["reconciliation"]["rows_featurised"] == report["reconciliation"]["curated_compounds"]
    assert report["split_row_order_matches_saved_assignments"] == {"random": True, "scaffold": True}
    assert "utc" not in json.dumps(report).lower()  # a timestamp would break reproducibility


def test_a_nonfinite_descriptor_is_reported_and_aborts_instead_of_being_imputed(tmp_path, monkeypatch):
    raw, curated, split_dir = _inputs(tmp_path)
    broken = (("MolWt", lambda mol: float("nan"), "injected fault."), *features.DESCRIPTORS[1:])
    monkeypatch.setattr(features, "DESCRIPTORS", broken)
    out = tmp_path / "features"
    with pytest.raises(ValueError, match="failed featurisation"):
        features.build_features(raw, curated, split_dir, out)
    reported = _read_csv(out / "featurisation_failures.csv")
    assert len(reported) == len(_read_csv(curated / "compounds.csv"))
    assert all("non-finite descriptors: MolWt" in row["reason"] for row in reported)
    assert not (out / "descriptors.csv").exists() and not (out / retrieve.MANIFEST).exists()


@pytest.mark.parametrize("corrupt", ["raw", "curated", "split", "provenance",
                                     "fingerprint_spec", "split_order"])
def test_corrupt_or_regenerated_inputs_fail_before_a_cache_is_published(tmp_path, corrupt):
    raw, curated, split_dir = _inputs(tmp_path)
    if corrupt == "raw":
        (raw / retrieve.RAW_TARGET).write_text("{}")
    elif corrupt == "curated":
        with (curated / "compounds.csv").open("a") as fh:
            fh.write("corrupted\n")
    elif corrupt == "split":
        with (split_dir / "random_assignments.csv").open("a") as fh:
            fh.write("corrupted\n")
    elif corrupt == "split_order":
        # A re-randomised or reordered split keeps its checksum but breaks the join.
        path = split_dir / "random_assignments.csv"
        rows = _read_csv(path)
        rows.reverse()
        splits._write_csv(path, rows, list(rows[0]))
        audit = json.loads((split_dir / "audit.json").read_text())
        audit["files"][path.name]["sha256"] = retrieve._sha256(path)
        (split_dir / "audit.json").write_text(json.dumps(audit))
    else:
        audit = json.loads((split_dir / "audit.json").read_text())
        if corrupt == "provenance":
            audit["raw_manifest_sha256"] = "wrong"
        else:
            audit["fingerprint"] = {**splits.FINGERPRINT_SPEC, "fpSize": 1024}
        (split_dir / "audit.json").write_text(json.dumps(audit))
    with pytest.raises(ValueError, match="checksum|provenance|disagree|row order|specification"):
        features.build_features(raw, curated, split_dir, tmp_path / "features")
    assert not (tmp_path / "features" / "descriptors.csv").exists()


def test_cache_cannot_be_written_into_any_source_tree(tmp_path):
    raw, curated, split_dir = _inputs(tmp_path)
    for dest in (raw, raw / "derived", curated, split_dir, split_dir / "derived", tmp_path):
        with pytest.raises(ValueError, match="separate"):
            features.build_features(raw, curated, split_dir, dest)


def test_published_dpp4_cache_matches_its_own_manifest_and_the_frozen_inputs():
    manifest = json.loads((CACHE / retrieve.MANIFEST).read_text())
    for name, meta in manifest["files"].items():
        assert retrieve._sha256(CACHE / name) == meta["sha256"]
    assert manifest["source"]["curated_compounds_sha256"] == \
        retrieve._sha256(CURATED / "compounds.csv")
    assert manifest["source"]["split_audit_sha256"] == retrieve._sha256(SPLIT_DIR / "audit.json")
    assert manifest["source"]["raw_manifest_sha256"] == \
        retrieve._sha256(features.DATA_RAW / "CHEMBL284" / retrieve.MANIFEST)
    assert manifest["seed"] == curate.SEED == 42
    assert _read_csv(CACHE / "featurisation_failures.csv") == []
