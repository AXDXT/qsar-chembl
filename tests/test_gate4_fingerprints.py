"""Fold-local collision selection, provenance and resumable plan guards."""
import json

import pytest
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator

from src import gate4_fingerprints as fp
from src import nested_cv, splits
from test_baselines import _cached_inputs


def _toy():
    rows = []
    for i in range(20):
        smiles = "C" * (i + 1)
        rows.append({"inchikey": Chem.MolToInchiKey(Chem.MolFromSmiles(smiles)),
                     "canonical_smiles": smiles, "partition": "train",
                     "scaffold": "" if i < 2 else f"core{i // 2}"})
    return rows


def _run(path, rows=None, **kwargs):
    return fp.plan_fingerprints(_toy() if rows is None else rows, path,
                               source_checksums=kwargs.pop("source_checksums", {"data": "toy"}),
                               **kwargs)


def _audit(directory, rows):
    manifest = json.loads((directory / "manifest.json").read_text())
    record = json.loads((directory / "report.json").read_text())
    report = record["result"]
    assert record["identity"] == nested_cv._digest(manifest) == report["identity_sha256"]
    assert record["sha256"] == nested_cv._digest(report)
    evidence = json.loads((directory / "compound-collisions.json").read_text())["result"]
    assert [r["inchikey"] for r in evidence["rows"]] == [r["inchikey"] for r in rows]
    # Independently regenerate every compound's collision losses from RDKit.
    for row, result in zip(rows, evidence["rows"], strict=True):
        molecule = Chem.MolFromSmiles(row["canonical_smiles"])
        generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, includeChirality=False)
        identifiers = len(generator.GetSparseCountFingerprint(molecule).GetNonzeroElements())
        assert result["identifiers"] == identifiers
        for size, lost in result["lost_by_bit_length"].items():
            generator = rdFingerprintGenerator.GetMorganGenerator(
                radius=2, includeChirality=False, fpSize=int(size))
            assert lost == identifiers - generator.GetFingerprint(molecule).GetNumOnBits()

    expected = {"final-fit": list(range(len(rows)))}
    for number, outer in enumerate(manifest["outer_plan"]):
        expected[f"outer{number}-evaluation"] = outer["train"]
        for inner_number, inner in enumerate(manifest["inner_plan"][number]):
            expected[f"outer{number}-inner-fold{inner_number}"] = inner["train"]
            assert set(inner["train"]).isdisjoint(outer["validation"])
    for number, fold in enumerate(manifest["final_plan"]):
        expected[f"final-inner-fold{number}"] = fold["train"]
    assert set(report["fit_decisions"]) == set(expected)
    for name, decision in report["fit_decisions"].items():
        indices = expected[name]
        assert decision["train_indices"] == indices
        assert decision["training_compounds"] == len(indices)
        for candidate in decision["candidates"]:
            losses = [evidence["rows"][i]["lost_by_bit_length"][str(candidate["fpSize"])]
                      for i in indices]
            assert candidate["lost_identifiers"] == sum(losses)
            assert candidate["compounds_losing_identifiers"] == sum(x > 0 for x in losses)
        best = min(decision["candidates"], key=lambda c: (c["lost_identifiers"], c["fpSize"]))
        assert decision["selected_fpSize"] == best["fpSize"]
    assert report["estimators_fitted"] == report["heldout_evaluations"] == 0
    return report, manifest


def test_every_decision_is_rederived_on_exact_fitting_indices_without_labels(tmp_path):
    rows = _toy()  # Deliberately no pIC50 field at all.
    report = _run(tmp_path / "plan", rows)
    assert report == _audit(tmp_path / "plan", rows)[0]
    assert report["fit_decisions"]["final-fit"]["selected_fpSize"] == 1024  # No benefit on toy alkanes.
    labelled = [{**r, "pIC50_median": object()} for r in rows]
    assert report == _run(tmp_path / "plan", labelled)  # Unreadable labels cannot enter selection.


def test_decision_minimises_total_loss_and_uses_size_only_to_break_ties():
    evidence = [{"lost_by_bit_length": {"1024": 8, "2048": 1, "4096": 1}},
                {"lost_by_bit_length": {"1024": 0, "2048": 1, "4096": 1}}]
    # Fewer affected compounds is not the criterion: preserving more identifiers is.
    decision = fp._decision(evidence, [0, 1])
    assert decision["selected_fpSize"] == 2048
    assert decision["candidates"][0]["compounds_losing_identifiers"] == 1
    assert decision["candidates"][1]["compounds_losing_identifiers"] == 2
    assert fp._decision(evidence, [1])["selected_fpSize"] == 1024
    for indices in ([], [0, 0]):
        with pytest.raises(ValueError, match="unique nonempty"):
            fp._decision(evidence, indices)


def test_validation_collision_changes_cannot_influence_its_outer_or_inner_choices(tmp_path, monkeypatch):
    rows = _toy()
    outer = nested_cv._folds(rows, list(range(len(rows))), 5)
    excluded = set(outer[0]["validation"])
    altered = [False]

    def synthetic_evidence(rows):
        return {"measurement_seconds": 0., "rows": [
            {"inchikey": r["inchikey"], "identifiers": 100,
             "lost_by_bit_length": {"1024": 99 if i in excluded and altered[0] else 0,
                                    "2048": 50 if i in excluded and altered[0] else 0,
                                    "4096": 0}} for i, r in enumerate(rows)]}

    monkeypatch.setattr(fp, "_collision_evidence", synthetic_evidence)
    original = _run(tmp_path / "original", rows)
    altered[0] = True
    changed = _run(tmp_path / "changed", rows)
    assert original["fit_decisions"]["final-fit"]["selected_fpSize"] == 1024
    assert changed["fit_decisions"]["final-fit"]["selected_fpSize"] == 4096
    names = [n for n in original["fit_decisions"] if n.startswith("outer0-")]
    assert len(names) == 6
    for name in names:
        assert original["fit_decisions"][name] == changed["fit_decisions"][name]


def test_completed_measurements_survive_interruption_and_resume_never_remeasures(tmp_path, monkeypatch):
    actual = fp._decision
    monkeypatch.setattr(fp, "_decision", lambda *a: (_ for _ in ()).throw(KeyboardInterrupt()))
    with pytest.raises(KeyboardInterrupt):
        _run(tmp_path / "plan")
    assert (tmp_path / "plan/compound-collisions.json").exists()
    assert not (tmp_path / "plan/report.json").exists()
    monkeypatch.setattr(fp, "_decision", actual)
    monkeypatch.setattr(fp, "_collision_evidence", lambda *a: pytest.fail("Remeasured"))
    report = _run(tmp_path / "plan")
    monkeypatch.setattr(fp, "_decision", lambda *a: pytest.fail("Repeated decisions"))
    assert report == _run(tmp_path / "plan")


@pytest.mark.parametrize("changed", ["structure", "source", "code", "version", "folds", "candidates"])
def test_changed_inputs_or_rules_cannot_silently_reuse_plan(tmp_path, monkeypatch, changed):
    _run(tmp_path / "plan")
    rows, kwargs = _toy(), {}
    if changed == "structure":
        rows[0]["canonical_smiles"] = "N"
    elif changed == "source":
        kwargs["source_checksums"] = {"data": "changed"}
    elif changed == "code":
        monkeypatch.setattr(fp, "_sha256", lambda path: "changed")
    elif changed == "version":
        monkeypatch.setattr(fp.rdBase, "rdkitVersion", "changed")
    elif changed == "folds":
        kwargs["folds"] = 2
    else:
        monkeypatch.setattr(fp, "FOLDING_SENSITIVITY_BITS", (2048, 4096))
    with pytest.raises(ValueError, match="inputs/settings/code changed"):
        _run(tmp_path / "plan", rows, **kwargs)


@pytest.mark.parametrize("which", ["report", "compound-collisions"])
def test_tampered_cache_is_refused(tmp_path, monkeypatch, which):
    if which == "compound-collisions":
        monkeypatch.setattr(fp, "_decision", lambda *a: (_ for _ in ()).throw(KeyboardInterrupt()))
        with pytest.raises(KeyboardInterrupt):
            _run(tmp_path / "plan")
    else:
        _run(tmp_path / "plan")
    path = tmp_path / "plan" / f"{which}.json"
    record = json.loads(path.read_text())
    record["sha256"] = "wrong"
    path.write_text(json.dumps(record))
    with pytest.raises(ValueError, match="checksum/identity"):
        _run(tmp_path / "plan")


@pytest.mark.parametrize("invalid", ["test", "partition", "duplicate", "smiles", "groups", "folds", "source"])
def test_invalid_inputs_are_refused_before_cache_creation(tmp_path, invalid):
    rows, kwargs = _toy(), {}
    if invalid == "test":
        rows[0]["partition"] = "test"
    elif invalid == "partition":
        del rows[0]["partition"]
    elif invalid == "duplicate":
        rows[0]["inchikey"] = rows[1]["inchikey"]
    elif invalid == "smiles":
        rows[0]["canonical_smiles"] = ""
    elif invalid == "groups":
        for row in rows:
            row["scaffold"] = "same"
    elif invalid == "folds":
        kwargs["folds"] = True
    else:
        kwargs["source_checksums"] = {}
    with pytest.raises(ValueError):
        _run(tmp_path / "plan", rows, **kwargs)
    assert not (tmp_path / "plan").exists()


def test_structure_key_mismatch_is_refused(tmp_path):
    rows = _toy()
    rows[0]["canonical_smiles"] = "N"
    with pytest.raises(ValueError, match="SMILES/InChIKey mismatch"):
        _run(tmp_path / "plan", rows)
    assert not (tmp_path / "plan/report.json").exists()


def test_adapter_passes_only_saved_training_structures_and_checks_postread_hash(tmp_path, monkeypatch):
    inputs = _cached_inputs(tmp_path)
    calls = []

    def spy(rows, cache_dir, **kwargs):
        calls.append(rows)
        assert all(r["partition"] == "train" for r in rows)
        assert all(Chem.MolToInchiKey(Chem.MolFromSmiles(r["canonical_smiles"])) == r["inchikey"]
                   for r in rows)
        return {"identity_sha256": "toy"}

    monkeypatch.setattr(fp, "plan_fingerprints", spy)
    monkeypatch.setattr(splits, "make_partitions", lambda *a: pytest.fail("Regenerated splits"))
    assert set(fp.run_fingerprint_plans(*inputs, tmp_path / "cache")) == {"scaffold", "random"}
    assert len(calls) == 2
    actual = fp._read_csv

    def changed(path):
        result = actual(path)
        with path.open("a") as stream:
            stream.write("\n")
        return result

    monkeypatch.setattr(fp, "_read_csv", changed)
    with pytest.raises(ValueError, match="changed during"):
        fp.run_fingerprint_plans(*inputs, tmp_path / "cache")
    assert len(calls) == 2


def test_adapter_rejects_input_overlap_and_bad_split_requests(tmp_path):
    inputs = _cached_inputs(tmp_path)
    for dest in (*inputs, inputs[0] / "plan", tmp_path):
        with pytest.raises(ValueError, match="separate"):
            fp.run_fingerprint_plans(*inputs, dest)
    for requested in ((), ("random", "random"), ("wrong",)):
        with pytest.raises(ValueError, match="unique saved split"):
            fp.run_fingerprint_plans(*inputs, tmp_path / "cache", splits=requested)


@pytest.mark.parametrize("split", ["scaffold", "random"])
def test_production_plan_matches_descriptor_folds_and_saved_train_keys(split, monkeypatch):
    target = "CHEMBL284"
    directory = fp.RESULTS / "gate4_training" / target / split
    manifest = json.loads((directory / "fingerprint_plan/manifest.json").read_text())
    rows = manifest["rows"]
    report, plan = _audit(directory / "fingerprint_plan", rows)
    control = json.loads((directory / "descriptors/manifest.json").read_text())
    for key in ("outer_plan", "inner_plan", "final_plan"):
        assert plan[key] == control[key]
    assert [r["inchikey"] for r in rows] == [r["inchikey"] for r in control["rows"]]
    monkeypatch.setattr(fp, "_collision_evidence", lambda *a: pytest.fail("Production remeasured"))
    resumed = fp.run_fingerprint_plans(
        fp.DATA_RAW / target, fp.DATA_CURATED / target, fp.RESULTS / "splits" / target,
        fp.DATA_FEATURES / target, fp.RESULTS / "gate4_training" / target, splits=(split,))
    assert resumed[split] == report
