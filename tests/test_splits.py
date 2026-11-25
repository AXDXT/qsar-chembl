"""Offline tests of split isolation, deterministic allocation and the leakage audit."""
import csv
import json

import pytest

from src import curate, retrieve, splits
from test_curate import _row
from test_retrieve import _pull


def _compounds(smiles):
    compounds = []
    for i, smi in enumerate(smiles):
        structure = curate.standardise_smiles(smi)
        assert structure["ok"]
        compounds.append({"inchikey": structure["inchikey"],
                          "canonical_smiles": structure["canonical_smiles"],
                          "pIC50_median": str(5 + i / 10)})
    return compounds


SMILES = ["c1ccccc1", "Cc1ccccc1", "Oc1ccccc1", "c1ccncc1", "Cc1ccncc1",
          "C1CCCCC1", "CC1CCCCC1", "c1ccoc1", "CCO", "CCCC"]


@pytest.mark.parametrize("smiles,expected", [
    ("CCO", ""), ("C", ""), ("C=C", ""),
    ("Cc1ccccc1", "c1ccccc1"),
    ("Oc1ccncc1", "c1ccncc1"),
    ("Cc1cc(Oc2nccc(CCC)c2)ccc1", "c1ccc(Oc2ccccn2)cc1"),
    ("c1ccccc1CCc1ccccc1", "c1ccc(CCc2ccccc2)cc1"),
    ("CC(=O)C1CCCCC1", "C1CCCCC1"),
    ("O=C1CCCCC1", "O=C1CCCCC1"),
    ("C=C1CCCCC1", "C=C1CCCCC1"),
    ("Cn1cccc1", "c1cc[nH]c1"),
    ("CS(=O)(=O)n1ccnc1", "c1c[nH]cn1"),
    ("CC1CC2CCC1C2", "C1CC2CCC1C2"),
    ("c1ccc2ccccc2c1", "c1ccc2ccccc2c1"),
    ("C[C@H]1CCCCO1", "C1CCOCC1"),
])
def test_native_scaffold_handles_rings_linkers_and_exocyclic_bonds(smiles, expected):
    from rdkit import Chem
    from rdkit.Chem.Scaffolds import MurckoScaffold
    canonical = Chem.MolToSmiles(Chem.MolFromSmiles(expected), isomericSmiles=False)
    assert MurckoScaffold.MurckoScaffoldSmiles(smiles, includeChirality=False) == canonical


def test_splits_are_reproducible_disjoint_and_compound_complete():
    rows, _ = splits.describe_structures(_compounds(SMILES))
    assignments = splits.make_partitions(rows)
    assert assignments == splits.make_partitions(rows)
    assert assignments != splits.make_partitions(rows, seed=splits.SEED + 1)
    for name, assignment in assignments.items():
        splits.assert_split_integrity(rows, assignment, scaffold_disjoint=name == "scaffold")
        assert len(assignment) == len(rows)
    assert assignments["random"].count("test") == 2


def test_structure_order_is_canonical_and_does_not_depend_on_input_order():
    compounds = _compounds(SMILES)
    rows, _ = splits.describe_structures(compounds)
    reversed_rows, _ = splits.describe_structures(list(reversed(compounds)))
    assert rows == reversed_rows
    assert splits.make_partitions(rows) == splits.make_partitions(reversed_rows)


def test_scaffold_allocation_does_not_consult_pic50():
    rows, _ = splits.describe_structures(_compounds(SMILES))
    without_labels = [{"scaffold": r["scaffold"]} for r in rows]
    assert splits.make_partitions(rows) == splits.make_partitions(without_labels)


def test_acyclic_compounds_and_stereoisomers_stay_in_whole_scaffold_groups():
    compounds = _compounds(["CCO", "CCCC", "C[C@H]1CCCCO1", "C[C@@H]1CCCCO1", "c1ccccc1"])
    rows, _ = splits.describe_structures(compounds)
    assignment = splits.make_partitions(rows)["scaffold"]
    for scaffold in {r["scaffold"] for r in rows}:
        assert len({side for r, side in zip(rows, assignment) if r["scaffold"] == scaffold}) == 1
    assert sum(r["scaffold"] == "" for r in rows) == 2
    assert len({r["scaffold"] for r in rows}) == 3


def test_impossible_single_scaffold_split_fails_instead_of_splitting_the_core():
    rows, _ = splits.describe_structures(_compounds(["c1ccccc1", "Cc1ccccc1"]))
    with pytest.raises(ValueError, match="two scaffold groups"):
        splits.make_partitions(rows)


def test_smallest_whole_group_fallback_keeps_both_partitions_nonempty():
    rows = [{"scaffold": "A"}] * 3 + [{"scaffold": "B"}] * 4
    assignment = splits.make_partitions(rows, test_fraction=.01)["scaffold"]
    assert assignment == ["test"] * 3 + ["train"] * 4


@pytest.mark.parametrize("assignment", [["train"], ["train", "train"], ["train", "other"]])
def test_missing_rows_empty_partitions_or_unknown_assignments_fail(assignment):
    rows = [{"inchikey": "A", "scaffold": "A"}, {"inchikey": "B", "scaffold": "B"}]
    with pytest.raises(AssertionError, match="Every compound"):
        splits.assert_split_integrity(rows, assignment)


def test_exact_identity_leakage_fails_loudly():
    rows = [{"inchikey": "SAME", "scaffold": "A"}] * 2
    with pytest.raises(AssertionError, match="InChIKey leakage"):
        splits.assert_split_integrity(rows, ["train", "test"])


def test_scaffold_leakage_fails_even_with_different_inchikeys():
    rows = [{"inchikey": "A", "scaffold": "CORE"}, {"inchikey": "B", "scaffold": "CORE"}]
    with pytest.raises(AssertionError, match="Scaffold leakage"):
        splits.assert_split_integrity(rows, ["train", "test"], scaffold_disjoint=True)


def test_stereoisomer_leakage_is_reported_despite_distinct_full_keys():
    rows, fp = splits.describe_structures(_compounds(["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"]))
    audit, pairs, nearest = splits.audit_split(rows, fp, ["train", "test"])
    assert audit["inchikey_overlap"] == 0
    assert audit["shared_connectivity_keys"] == 1
    assert audit["same_connectivity_cross_boundary_pairs"] == 1
    assert audit["test_compounds_with_near_duplicate"] == 1
    assert len(pairs) == len(nearest) == 1
    assert pairs[0]["same_connectivity"] is True
    assert nearest[0]["tanimoto"] == pytest.approx(1)


def test_shared_connectivity_is_reported_independently_of_similarity_threshold(monkeypatch):
    rows, fp = splits.describe_structures(_compounds(["C[C@H](N)C(=O)O", "C[C@@H](N)C(=O)O"]))
    monkeypatch.setattr(splits.DataStructs, "BulkTanimotoSimilarity", lambda *args: [.5])
    audit, pairs, _ = splits.audit_split(rows, fp, ["train", "test"])
    assert len(pairs) == 1 and pairs[0]["same_connectivity"]
    assert audit["near_duplicate_cross_boundary_pairs"] == 0


def test_duplicate_or_mismatched_curated_identity_is_refused():
    compounds = _compounds(["CCO", "CCCC"])
    with pytest.raises(ValueError, match="duplicate"):
        splits.describe_structures([compounds[0], compounds[0]])
    compounds[0]["canonical_smiles"] = "N"
    with pytest.raises(ValueError, match="mismatch"):
        splits.describe_structures(compounds)


def _curated_fixture(tmp_path):
    records = [_row(activity_id=i, canonical_smiles=smi, standard_value=str(10 + i))
               for i, smi in enumerate(SMILES)]
    manifest = _pull(tmp_path / "raw", records=records)
    curate.curate_target(manifest.parent, tmp_path / "curated", tmp_path / "curation")
    return manifest.parent, tmp_path / "curated"


def test_end_to_end_outputs_are_reproducible_checksummed_and_sources_unmodified(tmp_path):
    raw, curated = _curated_fixture(tmp_path)
    before = {p: p.read_bytes() for directory in (raw, curated) for p in directory.iterdir()
              if not p.name.startswith("._")}
    out = tmp_path / "splits"
    report = splits.split_target(raw, curated, out)
    output = {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert report == splits.split_target(raw, curated, out)
    assert output == {p.name: p.read_bytes() for p in out.iterdir() if not p.name.startswith("._")}
    assert all(p.read_bytes() == contents for p, contents in before.items())
    assert retrieve.verify_raw(raw)
    for name, meta in report["files"].items():
        assert retrieve._sha256(out / name) == meta["sha256"]
    for name in ("random", "scaffold"):
        with (out / f"{name}_assignments.csv").open() as fh:
            assignments = list(csv.DictReader(fh))
        assert len(assignments) == len(SMILES)
        assert report["splits"][name]["inchikey_overlap"] == 0
    assert report["splits"]["scaffold"]["shared_scaffolds"] == 0


@pytest.mark.parametrize("corrupt", ["raw", "curated", "provenance"])
def test_bad_source_checksums_or_provenance_are_refused_before_writing(tmp_path, corrupt):
    raw, curated = _curated_fixture(tmp_path)
    if corrupt == "raw":
        (raw / retrieve.RAW_TARGET).write_text("{}")
    elif corrupt == "curated":
        with (curated / "compounds.csv").open("a") as fh:
            fh.write("corrupted\n")
    else:
        manifest = json.loads((curated / retrieve.MANIFEST).read_text())
        manifest["source"]["raw_manifest_sha256"] = "wrong"
        (curated / retrieve.MANIFEST).write_text(json.dumps(manifest))
    with pytest.raises(ValueError, match="checksum|manifest"):
        splits.split_target(raw, curated, tmp_path / "out")
    assert not (tmp_path / "out").exists()


def test_output_cannot_overwrite_source_files(tmp_path):
    raw, curated = _curated_fixture(tmp_path)
    for output in (raw, raw / "derived", curated, curated / "derived", tmp_path):
        with pytest.raises(ValueError, match="separate"):
            splits.split_target(raw, curated, output)


def test_saved_dpp4_splits_have_zero_identity_overlap_and_complete_membership():
    """Fail loudly if the actual saved coursework splits acquire exact leakage."""
    out = splits.RESULTS / "splits" / "CHEMBL284"
    report = json.loads((out / "audit.json").read_text())
    with (splits.DATA_CURATED / "CHEMBL284" / "compounds.csv").open() as fh:
        expected_keys = {r["inchikey"] for r in csv.DictReader(fh)}
    for name in ("random", "scaffold"):
        with (out / f"{name}_assignments.csv").open() as fh:
            rows = list(csv.DictReader(fh))
        assert {r["inchikey"] for r in rows} == expected_keys
        splits.assert_split_integrity(rows, [r["partition"] for r in rows], name == "scaffold")
        for side in ("train", "test"):
            assert sum(r["partition"] == side for r in rows) == report["splits"][name][side]["compounds"]
    for filename, metadata in report["files"].items():
        assert retrieve._sha256(out / filename) == metadata["sha256"]
