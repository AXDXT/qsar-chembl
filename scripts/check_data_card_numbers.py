"""Check DATA_CARD curation figures and the Gate 2 sign-off tables against the logs.

Usage (repo root): .venv/bin/python scripts/check_data_card_numbers.py
Reads the logged curation run, the curated CSVs and the frozen raw activities; exits
nonzero on any mismatch. Curation checks use the documented constants; sign-off table
checks read the Markdown and compare its cells directly with the logs. It writes nothing.
"""
from __future__ import annotations

import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
LOG = ROOT / "results/curation/CHEMBL284/curation_log.json"
CURATED = ROOT / "data/curated/CHEMBL284"
FEATURES = ROOT / "data/features/CHEMBL284"
RAW = ROOT / "data/raw/CHEMBL284"


def _table_rows(text: str) -> dict[str, list[list[str]]]:
    """Index every Markdown table row by its first cell, keeping duplicates visible."""
    rows: dict[str, list[list[str]]] = {}
    for line in text.splitlines():
        if line.startswith("|"):
            cells = [re.sub(r"\*\*", "", cell.strip()) for cell in line.strip("|").split("|")]
            rows.setdefault(cells[0], []).append(cells[1:])
    return rows


def gate3_feature_checks(text: str, features: dict) -> list[tuple[str, bool]]:
    """Re-derive every Gate 3 featurisation figure quoted in prose from the logged cache."""
    rows = _table_rows(text)
    reconciliation = features["reconciliation"]
    descriptor, fingerprint = features["descriptor_block"], features["fingerprint_block"]
    evidence = fingerprint["folding_evidence"]
    expected = {
        "Curated compounds in": [f"{reconciliation['curated_compounds']:,}"],
        "Rows featurised": [f"{reconciliation['rows_featurised']:,}"],
        "Featurisation failures": [str(reconciliation["featurisation_failures"])],
        "Non-finite descriptor values": [str(reconciliation["non_finite_descriptor_values"])],
        "Compounds dropped": [str(reconciliation["compounds_dropped"])],
        "Descriptor columns": [str(descriptor["count"])],
        "Morgan bits (radius 2, nonchiral)": [f"{fingerprint['fpSize']:,}"],
        "Distinct radius-≤2 substructure identifiers":
            [f"{evidence['distinct_substructure_identifiers']:,}"],
    }
    for bits in ("1024", "2048", "4096"):
        result = evidence["by_bit_length"][bits]
        expected[f"Occupied bits at {int(bits):,}"] = [f"{result['occupied_bits']:,}"]
        expected[f"Compounds losing ≥1 identifier at {int(bits):,}"] = [
            f"{result['compounds_losing_at_least_one_identifier']:,}"]
    for split, result in features["collinearity"]["train_rows_only"].items():
        expected[f"Largest absolute r among descriptors, {split} training rows"] = [
            f"{result['max_abs_pearson_r']:.3f} ({' / '.join(result['max_abs_pearson_pair'])})"]
    return [(f"Gate 3 feature table: {label}", rows.get(label) == [values])
            for label, values in expected.items()]


def signoff_table_checks(text: str, curation: dict, audit: dict, probe: dict) -> list[tuple[str, bool]]:
    """Compare the handoff's displayed tables with the run, including rounded CIs."""
    rows = _table_rows(text)
    expected = {}
    labels = ("Binding assays", "IC50 only", "Human target", "Finite positive molar units to nM",
              "Exact relation only", "Standardised structures", "Median per standardised full InChIKey")
    for label, step in zip(labels, curation["steps"], strict=True):
        values = [f"{step[key]:,}" for key in ("records_in", "records_out", "records_removed")]
        if label == labels[-1]:
            values = [values[0] + " records", values[1] + " compounds", values[2] + " rows collapsed"]
        expected[label] = values
    for label, field in (("Train / test compounds", "compounds"),
                         ("Train / test exact scaffolds", "scaffolds")):
        expected[label] = [f"{audit['splits'][split]['train'][field]:,} / "
                           f"{audit['splits'][split]['test'][field]:,}"
                           for split in ("random", "scaffold")]
    for label, field in (("Full InChIKey overlap, asserted zero", "inchikey_overlap"),
                         ("Shared exact scaffolds", "shared_scaffolds"),
                         ("Test compounds with Tanimoto ≥0.90 training neighbour", "test_compounds_with_near_duplicate")):
        expected[label] = [str(audit["splits"][split][field]) for split in ("random", "scaffold")]
    for label, metric in (("RMSE (95% descriptive interval), log units", "rmse"),
                          ("MAE (95% descriptive interval), log units", "mae"),
                          ("R² (95% descriptive interval)", "r2")):
        values = []
        for split in ("random", "scaffold"):
            result = probe["splits"][split]
            low, high = result["bootstrap"]["metrics"][metric]["ci95"]
            values.append(f"{result['metrics'][metric]:.3f} [{low:.3f}, {high:.3f}]")
        expected[label] = values
    # One row per label: missing and duplicate rows fail as well as wrong numbers.
    return [(f"sign-off table: {label}", rows.get(label) == [values])
            for label, values in expected.items()]


def main() -> int:
    d = json.loads(LOG.read_text())
    s = {x["step"]: x for x in d["steps"]}
    i = d["identity_decisions"]
    r = d["censoring_reconciliation"]
    t = d["target_variable"]
    n = d["replicate_noise"]
    checks: list[tuple[str, bool]] = []

    def ck(name: str, cond: object) -> None:
        checks.append((name, bool(cond)))

    ck("cascade counts", [(x["records_in"], x["records_out"], x["records_removed"]) for x in d["steps"]]
       == [(8397, 8174, 223), (8174, 5781, 2393), (5781, 5781, 0), (5781, 5650, 131),
           (5650, 4859, 791), (4859, 4855, 4), (4855, 3965, 890)])
    ck("28 other activity types", s["2_single_activity_type"]["removed_activity_types"] == 28)
    ck("unit refusals 78/48/5",
       sorted(s["4_units_harmonised_to_nM"]["removed_breakdown"].values()) == [5, 48, 78])
    ck("censored breakdown 732/54/4/1",
       sorted(s["5_censored_removed"]["removed_breakdown"].values()) == [1, 4, 54, 732])
    ck("855 salt or mixture rows", s["6_structure_standardised"]["salt_or_mixture_rows"] == 855)
    ck("2 rows neutralised", s["6_structure_standardised"]["charged_rows_neutralised"] == 2)
    ck("476 replicated compounds", s["7_duplicates_resolved"]["compounds_with_replicates"] == 476)
    ck("108 flagged compounds",
       s["7_duplicates_resolved"]["compounds_flagged_spread_over_1_log"] == 108)
    ck("4000 molecule / 3965 parent ids",
       (i["distinct_molecule_ids_in_curated_records"],
        i["distinct_parent_ids_in_curated_records"]) == (4000, 3965))
    ck("32 merged compounds absorbing 35 ids",
       (i["compounds_merging_more_than_one_molecule_id"],
        i["molecule_ids_absorbed_by_merging"]) == (32, 35))
    ck("286 additional compounds within shared connectivity groups",
       i["additional_compounds_sharing_connectivity_key"] == 286
       and i["stereoisomer_or_tautomer_pairs_sharing_a_skeleton"] == 286)
    ck("224 shared connectivity groups and 382 unordered pairs",
       (i["connectivity_key_groups_with_multiple_compounds"],
        i["compound_pairs_sharing_connectivity_key"]) == (224, 382))
    ck("src rows 4333 / 522", sorted(i["rows_by_src_id"].items()) == [("1", 4333), ("37", 522)])
    ck("validity comments 133 / 53",
       i["data_validity_comments"]["Outside typical range"] == 133
       and i["data_validity_comments"]["Potential transcription error"] == 53)
    ck("busiest compound 71 records", i["busiest_compound"]["n_records"] == 71)
    ck("busiest compound median 7.74", round(i["busiest_compound"]["pIC50_median"], 2) == 7.74)
    ck("busiest compound spread 4.14", round(i["busiest_compound"]["pIC50_spread"], 2) == 4.14)
    ck("censoring reconciliation 796 / 791 / 5",
       (r["gate1_rows_failing_exact_relation_recomputed"],
        r["gate2_step5_censored_rows_removed"], r["difference"]) == (796, 791, 5))
    ck("target summary statistics",
       [round(t[k], 2) for k in ("min", "p05", "median", "mean", "stdev", "p95", "max")]
       == [-0.04, 4.42, 7.17, 7.03, 1.33, 8.95, 10.10])
    ck("80 exact-integer pIC50 compounds", sum(t["exact_integer_pIC50_values"].values()) == 80)
    ck("integer pile-ups 8.0=33, 9.0=23",
       t["exact_integer_pIC50_values"]["8.0"] == 33 and t["exact_integer_pIC50_values"]["9.0"] == 23)
    ck("extremes 2 below pIC50 2, 108 below 4",
       (t["below_2_implausible"], t["below_4"]) == (2, 108))
    ceiling = t["pathology_checks"]["assay_ceiling"]
    ck("no ceiling spike: 6 of 359 share the densest weak value",
       (ceiling["compounds_sharing_it"], ceiling["compounds_below_weak_cutoff"]) == (6, 359))
    ck("densest repeated value anywhere 8.398 shared by 43",
       (round(ceiling["largest_repeated_value_anywhere"], 3),
        ceiling["compounds_sharing_largest_anywhere"]) == (8.398, 43))
    ck("extremes are not blamed on solubility",
       "past the solubility" not in t["pathology_checks"]["implausible_extremes"]["finding"])
    ck("noise 476 compounds / 1366 reports / 890 dof",
       (n["replicated_compounds"], n["replicate_measurements"],
        n["degrees_of_freedom"]) == (476, 1366, 890))
    ck("pooled SD all groups 0.828",
       round(n["pooled_within_compound_sd_all_groups"], 3) == 0.828)
    ck("pooled SD differing groups 1.011",
       round(n["pooled_within_compound_sd_differing_groups"], 3) == 1.011)
    ck("pooled SD differing, spread<=1 is 0.310",
       round(n["pooled_within_compound_sd_differing_under_1_log"], 3) == 0.310)
    ck("noise dof 597 / 208",
       (n["degrees_of_freedom_differing_groups"],
        n["degrees_of_freedom_differing_under_1_log"]) == (597, 208))
    ck("254 groups with spread >=0.01", n["groups_spread_at_least_0_01"] == 254)
    ck("222 groups with spread <0.01", n["groups_spread_below_0_01"] == 222)
    ck("56 groups spread >=3 log units", n["spread_histogram_log_units"][">=3"] == 56)
    ck("199 exactly equal groups", n["groups_exactly_equal"] == 199)
    ck("max spread 8.32", round(n["spread_max"], 2) == 8.32)

    # Figures quoted in the card that are computed from the files rather than read from the log.
    import sys
    sys.path.insert(0, str(ROOT))
    from src.retrieve import read_jsonl_gz  # noqa: E402  (after sys.path setup)

    compounds = pd.read_csv(CURATED / "compounds.csv")
    records = pd.read_csv(CURATED / "records_curated.csv")
    connectivity_counts = compounds.groupby("inchikey_skeleton").size()
    ck("connectivity group counts independently match curated CSV",
       (int((connectivity_counts > 1).sum()), int((connectivity_counts - 1).sum()),
        int((connectivity_counts * (connectivity_counts - 1) // 2).sum())) == (224, 286, 382))
    busiest_id, raw_rows = Counter(
        row["molecule_chembl_id"] for row in read_jsonl_gz(RAW / "activities.jsonl.gz")
    ).most_common(1)[0]
    surviving = records.loc[records.molecule_chembl_id == busiest_id]
    ck("busiest raw molecule is CHEMBL1422 with 96 rows", (busiest_id, raw_rows) == ("CHEMBL1422", 96))
    ck("63 of its rows survive to one compound",
       (len(surviving), surviving.inchikey.nunique()) == (63, 1))
    ck("3489 of 3965 compounds measured once",
       (int((compounds.n_records == 1).sum()), len(compounds)) == (3489, 3965))
    groups = records.groupby("inchikey").pIC50.agg(["count", "min", "max"])
    groups = groups.loc[groups["count"] > 1]
    spread = groups["max"] - groups["min"]
    ck("146 differing groups with spread <=1 log unit",
       int(((spread >= 0.01) & (spread <= 1)).sum()) == 146)

    ck("23 near-equal groups have nonzero spread", int(((spread > 0) & (spread < 0.01)).sum()) == 23)
    ck("29 groups at 3 and 27 above 3", (int((spread == 3).sum()), int((spread > 3).sum())) == (29, 27))

    checks.extend(signoff_table_checks(
        (ROOT / "GATE2_SIGNOFF.md").read_text(), d,
        json.loads((ROOT / "results/splits/CHEMBL284/audit.json").read_text()),
        json.loads((ROOT / "results/leakage_probe/CHEMBL284/probe.json").read_text())))
    checks.extend(gate3_feature_checks(
        (ROOT / "DATA_CARD.md").read_text(),
        json.loads((FEATURES / "manifest.json").read_text())))

    for name, ok in checks:
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    failed = [name for name, ok in checks if not ok]
    print(f"{len(checks) - len(failed)}/{len(checks)} data-card and sign-off checks match the logged run")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
