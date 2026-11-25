"""Offline tests for src.curate: the cascade, the standardiser and the noise floor.

Every test runs on small hand-written records, so nothing here depends on the frozen pull
or on the network. The one test that touches data/raw only reads it.
"""
import copy
import csv
import json
import math

import pytest

from src import curate, retrieve
from test_retrieve import _pull


def _row(**kw):
    """A record that passes every filter unless a test overrides a field."""
    base = {
        "activity_id": 1, "assay_type": "B", "standard_type": "IC50",
        "target_organism": "Homo sapiens", "target_tax_id": 9606,
        "standard_units": "nM", "standard_value": "10.0", "standard_relation": "=",
        "canonical_smiles": "CCO", "molecule_chembl_id": "CHEMBL1",
        "parent_molecule_chembl_id": "CHEMBL1", "assay_chembl_id": "A1",
        "document_chembl_id": "D1", "src_id": 1,
    }
    base.update(kw)
    return base


# ----------------------------------------------------------------------------- pIC50

def test_pic50_matches_the_definition():
    assert curate.pic50_from_nm(1e9) == pytest.approx(0.0)      # 1 M
    assert curate.pic50_from_nm(1.0) == pytest.approx(9.0)      # 1 nM
    assert curate.pic50_from_nm(10.0) == pytest.approx(8.0)


# ------------------------------------------------------------------- unit harmonisation

@pytest.mark.parametrize("unit,value,expected_nm", [
    ("nM", "10", 10.0), ("uM", "10", 10_000.0), ("mM", "1", 1_000_000.0),
    ("M", "1", 1e9), ("pM", "1000", 1.0),
])
def test_molar_units_are_converted_to_nm(unit, value, expected_nm):
    log = curate.run_cascade([_row(standard_units=unit, standard_value=value)])
    assert log["records"][0]["value_nM"] == pytest.approx(expected_nm)


def test_mass_concentration_is_refused_and_reported_not_silently_dropped():
    log = curate.run_cascade([_row(standard_units="ug.mL-1", units="ug ml-1")])
    step = next(s for s in log["steps"] if s["step"] == "4_units_harmonised_to_nM")
    assert step["records_out"] == 0
    assert step["removed_breakdown"] == {"unit_not_convertible=ug.mL-1": 1}
    assert log["unit_refusals"][0]["standard_units"] == "ug.mL-1"
    assert "molar" in log["unit_refusals"][0]["reason"]


@pytest.mark.parametrize("value", [None, "", "not a number", "0", "-5"])
def test_values_that_cannot_become_a_pic50_are_refused_and_itemised(value):
    log = curate.run_cascade([_row(standard_value=value)])
    step = next(s for s in log["steps"] if s["step"] == "4_units_harmonised_to_nM")
    assert step["records_out"] == 0
    assert len(log["unit_refusals"]) == 1


# ------------------------------------------------------------------------- censoring

def test_censored_and_missing_relations_are_removed_and_counted_separately():
    rows = ([_row(activity_id=i, standard_relation="=") for i in range(3)]
            + [_row(activity_id=10, standard_relation=">"),
               _row(activity_id=11, standard_relation="<"),
               _row(activity_id=12, standard_relation=">>"),
               _row(activity_id=13, standard_relation=None)])
    log = curate.run_cascade(rows)
    step = next(s for s in log["steps"] if s["step"] == "5_censored_removed")
    assert step["records_in"] == 7 and step["records_out"] == 3
    assert step["censored_rows_removed"] == 3
    assert step["missing_relation_rows_removed"] == 1
    assert step["removed_breakdown"]["relation=<missing>"] == 1


def test_censoring_reconciliation_accounts_for_every_row_of_the_gate1_figure():
    """A value-less row counts against the Gate 1 relation rule but leaves earlier here."""
    rows = [_row(activity_id=1), _row(activity_id=2, standard_relation=">"),
            _row(activity_id=3, standard_relation="<"),
            # usable unit, no value and no relation: inside Gate 1's count, gone at step 4 here
            _row(activity_id=4, standard_relation=None, standard_value=None)]
    log = curate.run_cascade(rows)
    rec = curate.censoring_reconciliation(rows, log)

    assert rec["gate1_base_rows_recomputed"] == 4
    assert rec["gate1_rows_failing_exact_relation_recomputed"] == 3
    assert rec["gate2_step5_rows_removed"] == 2
    assert rec["difference"] == 1
    assert rec["difference_explained_by"]["rows"] == [4]
    assert rec["difference_explained_by"]["all_have_null_value"] is True
    assert all(rec["identity_checks"][k] is True for k in rec["identity_checks"]
               if isinstance(rec["identity_checks"][k], bool))


def test_censoring_reconciliation_is_exact_when_nothing_is_value_less():
    rows = [_row(activity_id=1), _row(activity_id=2, standard_relation=">")]
    log = curate.run_cascade(rows)
    rec = curate.censoring_reconciliation(rows, log)
    assert rec["difference"] == 0
    assert rec["difference_explained_by"]["count"] == 0
    assert rec["identity_checks"]["gate1_count_equals_step5_plus_difference"] is True
    assert rec["identity_checks"]["surviving_rows_agree_with_gate1"] is True


# -------------------------------------------------------------------- standardisation

def test_salt_is_reduced_to_its_largest_fragment():
    salt = curate.standardise_smiles("Cl.CCO")
    parent = curate.standardise_smiles("CCO")
    assert salt["ok"] and parent["ok"]
    assert salt["inchikey"] == parent["inchikey"]
    assert salt["fragments_before"] == 2 and parent["fragments_before"] == 1


def test_charged_form_is_neutralised_to_the_same_compound():
    anion = curate.standardise_smiles("CC(=O)[O-]")
    acid = curate.standardise_smiles("CC(=O)O")
    assert anion["ok"] and acid["ok"]
    assert anion["inchikey"] == acid["inchikey"]
    assert anion["formal_charge_before"] == -1


def test_standardisation_is_idempotent():
    once = curate.standardise_smiles("Cl.C[N+](C)(C)CC(=O)[O-]")
    twice = curate.standardise_smiles(once["canonical_smiles"])
    assert once["ok"] and twice["ok"]
    assert once["canonical_smiles"] == twice["canonical_smiles"]
    assert once["inchikey"] == twice["inchikey"]


def test_stereoisomers_stay_different_compounds():
    left = curate.standardise_smiles("C[C@H](N)C(=O)O")
    right = curate.standardise_smiles("C[C@@H](N)C(=O)O")
    assert left["inchikey"] != right["inchikey"]
    assert left["inchikey_skeleton"] == right["inchikey_skeleton"]


@pytest.mark.parametrize("smiles,reason", [
    (None, "missing_smiles"), ("", "missing_smiles"), ("   ", "missing_smiles"),
    ("not-a-smiles", "parse_failed"), ("c1ccccc1C(", "parse_failed"), ("C(C)(C)(C)(C)C", "sanitisation_failed"),
])
def test_structure_failures_are_named_never_silent(smiles, reason):
    result = curate.standardise_smiles(smiles)
    assert result["ok"] is False
    assert result["failure"].startswith(reason)


def test_failed_structures_appear_in_the_failure_log():
    log = curate.run_cascade([_row(activity_id=1, canonical_smiles="CCO"),
                              _row(activity_id=2, canonical_smiles=None)])
    step = next(s for s in log["steps"] if s["step"] == "6_structure_standardised")
    assert step["records_in"] == 2 and step["records_out"] == 1
    assert step["failed_rows"] == 1
    assert log["structure_failures"][0]["activity_id"] == 2
    assert "missing_smiles" in log["structure_failures"][0]["failure"]


# ----------------------------------------------------------------- duplicate resolution

def test_duplicates_collapse_to_the_median_and_record_the_spread():
    # IC50 1, 10 and 10000 nM -> pIC50 9, 8 and 5; median 8, spread 4.
    rows = [_row(activity_id=i, standard_value=v, assay_chembl_id=f"A{i}", document_chembl_id=f"D{i}")
            for i, v in enumerate(["1", "10", "10000"])]
    log = curate.run_cascade(rows)
    assert len(log["compounds"]) == 1
    compound = log["compounds"][0]
    assert compound["pIC50_median"] == pytest.approx(8.0)
    assert compound["pIC50_min"] == pytest.approx(5.0)
    assert compound["pIC50_max"] == pytest.approx(9.0)
    assert compound["pIC50_spread"] == pytest.approx(4.0)
    assert compound["n_records"] == 3
    assert compound["flag_spread_over_1_log"] is True
    assert compound["n_assays"] == 3 and compound["n_documents"] == 3


def test_tight_replicates_are_not_flagged():
    rows = [_row(activity_id=i, standard_value=v) for i, v in enumerate(["10", "12"])]
    compound = curate.run_cascade(rows)["compounds"][0]
    assert compound["pIC50_spread"] < curate.LARGE_SPREAD_LOG_UNITS
    assert compound["flag_spread_over_1_log"] is False


def test_salt_and_parent_become_one_compound_not_two():
    """The Gate 1 open question: a salt form must not count as a second compound."""
    rows = [_row(activity_id=1, molecule_chembl_id="CHEMBL_SALT",
                 parent_molecule_chembl_id="CHEMBL_PARENT", canonical_smiles="Cl.CCO", standard_value="10"),
            _row(activity_id=2, molecule_chembl_id="CHEMBL_PARENT",
                 parent_molecule_chembl_id="CHEMBL_PARENT", canonical_smiles="CCO", standard_value="10")]
    log = curate.run_cascade(rows)
    assert len(log["compounds"]) == 1
    compound = log["compounds"][0]
    assert compound["n_molecule_chembl_ids"] == 2
    assert compound["molecule_chembl_ids"] == "CHEMBL_PARENT;CHEMBL_SALT"
    decisions = curate.identity_decisions(log)
    assert decisions["compounds_merging_more_than_one_molecule_id"] == 1
    assert decisions["molecule_ids_absorbed_by_merging"] == 1


def test_stereoisomers_are_kept_apart_by_the_duplicate_rule():
    rows = [_row(activity_id=1, canonical_smiles="C[C@H](N)C(=O)O"),
            _row(activity_id=2, canonical_smiles="C[C@@H](N)C(=O)O")]
    log = curate.run_cascade(rows)
    assert len(log["compounds"]) == 2
    assert curate.identity_decisions(log)["stereoisomer_or_tautomer_pairs_sharing_a_skeleton"] == 1


# ------------------------------------------------------------------------ cascade shape

def test_each_filter_is_a_separate_step_in_the_required_order():
    log = curate.run_cascade([_row()])
    assert [s["step"] for s in log["steps"]] == [name for name, _ in curate.STEPS]
    for earlier, later in zip(log["steps"], log["steps"][1:]):
        assert later["records_in"] == earlier["records_out"]


def test_every_removed_record_is_accounted_for():
    rows = [_row(activity_id=1), _row(activity_id=2, assay_type="F"),
            _row(activity_id=3, standard_type="Ki"), _row(activity_id=4, target_organism="Rattus norvegicus"),
            _row(activity_id=5, standard_units="ug.mL-1"), _row(activity_id=6, standard_relation=">"),
            _row(activity_id=7, canonical_smiles="not-a-smiles")]
    log = curate.run_cascade(rows)
    for step in log["steps"][:-1]:  # the last step changes unit from records to compounds
        assert step["records_removed"] == sum(step["removed_breakdown"].values()), step["step"]
    assert log["steps"][-1]["records_out"] == 1


def test_input_records_are_never_mutated():
    rows = [_row(activity_id=1), _row(activity_id=2, standard_units="uM")]
    original = copy.deepcopy(rows)
    curate.run_cascade(rows)
    assert rows == original


def test_empty_input_invents_no_statistics():
    log = curate.run_cascade([])
    assert log["compounds"] == [] and log["records"] == []
    assert curate.replicate_noise([], [])["replicated_compounds"] == 0


# ------------------------------------------------------------------- target distribution

def test_assay_ceiling_check_detects_a_planted_ceiling_spike():
    """The ceiling check must be able to fire, not just report 'no spike' on this dataset."""
    ceilinged = [{"pIC50_median": 4.0} for _ in range(40)]
    spread = [{"pIC50_median": 6.0 + i / 100} for i in range(40)]
    summary = curate.target_variable_summary(ceilinged + spread)
    ceiling = summary["pathology_checks"]["assay_ceiling"]

    assert ceiling["largest_repeated_value_below_weak_cutoff"] == 4.0
    assert ceiling["compounds_sharing_it"] == 40
    assert ceiling["compounds_below_weak_cutoff"] == 40
    # The planted spike is the densest repeated value in the whole set, which is exactly the
    # signature the real dataset does NOT show.
    assert ceiling["largest_repeated_value_anywhere"] == 4.0


def test_assay_ceiling_check_reports_no_spike_when_the_weak_tail_is_spread_out():
    values = [{"pIC50_median": 3.0 + i / 50} for i in range(60)]
    ceiling = curate.target_variable_summary(values)["pathology_checks"]["assay_ceiling"]
    assert ceiling["compounds_sharing_it"] == 1
    assert "not proof" in ceiling["limit_of_this_check"]


def test_target_summary_does_not_blame_extremes_on_solubility():
    """Solubility is not a field in this pull, so it cannot be offered as the explanation."""
    summary = curate.target_variable_summary([{"pIC50_median": v} for v in (0.5, 6.0, 9.0)])
    extremes = summary["pathology_checks"]["implausible_extremes"]["finding"]
    assert "past the solubility" not in extremes
    assert "no solubility or top-concentration field" in extremes


# -------------------------------------------------------------------------- noise floor

def test_repeat_dispersion_is_measured_not_assumed():
    # One group differing by 0.3 log and one equal-valued group.
    rows = ([_row(activity_id=1, standard_value="10", canonical_smiles="CCO"),
             _row(activity_id=2, standard_value="20", canonical_smiles="CCO"),
             _row(activity_id=3, standard_value="10", canonical_smiles="CCC"),
             _row(activity_id=4, standard_value="10", canonical_smiles="CCC")])
    log = curate.run_cascade(rows)
    noise = curate.replicate_noise(log["compounds"], log["records"])
    assert noise["replicated_compounds"] == 2
    assert noise["groups_spread_below_0_01"] == 1
    assert noise["groups_spread_at_least_0_01"] == 1
    # The only differing group spans log10(20/10) = 0.301, so its SD is that over sqrt(2).
    assert noise["noise_floor"]["recommended_for_gate4"] == pytest.approx(math.log10(2) / math.sqrt(2))
    # Both sensitivity subsets contain the same group in this example.
    assert (noise["noise_floor"]["ceiling_including_irreconcilable_groups"]
            >= noise["noise_floor"]["recommended_for_gate4"])
    # Equal-valued groups reduce the all-groups descriptive SD in this example.
    assert noise["noise_floor"]["all_groups_figure_and_why_it_is_not_used"]["value"] < (
        noise["noise_floor"]["recommended_for_gate4"])


def test_conditional_dispersion_is_unchanged_by_additional_equal_reports():
    """Selection changes the estimand; the conditional statistic excludes equal groups."""
    differing = [_row(activity_id=1, standard_value="10", canonical_smiles="CCO"),
                 _row(activity_id=2, standard_value="20", canonical_smiles="CCO")]
    log_small = curate.run_cascade(differing)
    padded = differing + [_row(activity_id=100 + i, standard_value="10",
                               canonical_smiles=f"CC{'C' * (i + 1)}") for i in range(6)]
    log_padded = curate.run_cascade(padded)
    a = curate.replicate_noise(log_small["compounds"], log_small["records"])
    b = curate.replicate_noise(log_padded["compounds"], log_padded["records"])
    assert a["noise_floor"]["recommended_for_gate4"] == pytest.approx(
        b["noise_floor"]["recommended_for_gate4"])


def test_repeat_subset_and_histogram_boundaries_do_not_claim_exact_equality():
    widths = [0, 0.005, 0.01, 1, 3, 4]
    records = [{"inchikey": str(i), "pIC50": value}
               for i, width in enumerate(widths) for value in (0, width)]
    compounds = [{"n_records": 2, "pIC50_spread": width} for width in widths]
    noise = curate.replicate_noise(compounds, records)

    assert noise["groups_exactly_equal"] == 1
    assert noise["groups_spread_below_0_01"] == 2
    assert noise["groups_spread_at_least_0_01"] == 4
    assert noise["degrees_of_freedom_differing_under_1_log"] == 2
    assert noise["spread_histogram_log_units"] == {
        "0-0.01": 2, "0.01-0.3": 1, "0.3-1": 0, "1-2": 1, "2-3": 0, ">=3": 2,
    }
    assert noise["pooled_within_compound_sd_all_groups"] == pytest.approx(
        math.sqrt(sum(w ** 2 for w in widths) / 12))
    assert noise["noise_floor"]["recommended_for_gate4"] == pytest.approx(
        math.sqrt((0.01 ** 2 + 1) / 4))
    assert noise["noise_floor"]["status"].startswith("pending_human_decision")


# ------------------------------------------------------------------- end to end on fakes

def test_curate_target_writes_a_manifest_and_leaves_raw_untouched(tmp_path):
    records = [_row(activity_id=1, standard_value="10", canonical_smiles="CCO"),
               _row(activity_id=2, standard_value="30", canonical_smiles="CCO"),
               _row(activity_id=3, standard_value="100", canonical_smiles="c1ccccc1"),
               _row(activity_id=4, standard_relation=">", canonical_smiles="CCCC")]
    manifest = _pull(tmp_path / "raw", records=records)
    before = {p.name: p.read_bytes() for p in manifest.parent.iterdir() if not p.name.startswith("._")}

    log = curate.curate_target(manifest.parent, tmp_path / "curated", tmp_path / "results")

    assert before == {p.name: p.read_bytes() for p in manifest.parent.iterdir()
                      if not p.name.startswith("._")}
    assert retrieve.verify_raw(manifest.parent)
    assert len(log["compounds"]) == 2

    written = json.loads((tmp_path / "curated" / retrieve.MANIFEST).read_text())
    assert written["curated_compounds"] == 2
    assert written["curated_records"] == 3
    assert written["seed"] == curate.SEED
    assert written["source"]["raw_verified_before_and_after"] is True
    assert [s["step"] for s in written["cascade"]] == [name for name, _ in curate.STEPS]
    for name, meta in written["files"].items():
        assert meta["sha256"] == retrieve._sha256(tmp_path / "curated" / name)

    with (tmp_path / "curated" / "compounds.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2
    assert {"inchikey", "canonical_smiles", "pIC50_median", "n_records", "pIC50_spread"} <= set(rows[0])

    cascade = (tmp_path / "results" / "cascade.md").read_text()
    assert "5_censored_removed" in cascade
    for name in ("cascade.csv", "cascade_figure.csv", "curation_log.json", "replicate_noise.json",
                 "replicate_spreads.csv", "pic50_distribution.png", "cascade_funnel.png"):
        assert (tmp_path / "results" / name).exists(), name


def test_curation_refuses_a_corrupt_raw_pull(tmp_path):
    manifest = _pull(tmp_path / "raw")
    (manifest.parent / retrieve.RAW_TARGET).write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        curate.curate_target(manifest.parent, tmp_path / "curated", tmp_path / "results")
    assert not (tmp_path / "curated").exists()
