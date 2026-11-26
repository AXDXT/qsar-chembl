"""Check DATA_CARD and Gate 2–4 result tables against the logged runs.

Usage (repo root): .venv/bin/python scripts/check_data_card_numbers.py
Reads the logged curation run, the curated CSVs and the frozen raw activities; exits
nonzero on any mismatch. Curation checks use the documented constants; sign-off table
checks read the Markdown and compare its cells directly with the logs. It writes nothing.
"""
from __future__ import annotations

import hashlib
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


def gate3_baseline_checks(text: str, report: dict) -> list[tuple[str, bool]]:
    """Compare each baseline's displayed train/test scores, CIs and fixed choices to logs."""
    rows = _table_rows(text)
    expected = {}
    for split in ("random", "scaffold"):
        result = report["splits"][split]
        for model, label in (("mean", "Training mean"), ("descriptor_rf", "Descriptor RF")):
            score = result["models"][model]
            cells = [f"{result['train_compounds']:,} / {result['test_compounds']:,}"]
            cells += [f"{score['train'][metric]:.3f}" if score['train'][metric] is not None
                      else "undefined" for metric in ("r2", "rmse", "mae")]
            for metric in ("r2", "rmse", "mae"):
                value = score["test"][metric]
                ci = score["test_bootstrap"]["metrics"][metric]["ci95"]
                cells.append(f"{value:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]"
                             if value is not None and ci else "undefined")
            expected[f"{split} / {label}"] = cells
    params = report["models"]["descriptor_rf"]["parameters"]
    for label, key in (("Forest trees", "n_estimators"), ("Forest max_depth", "max_depth"),
                       ("Forest max_features", "max_features"),
                       ("Forest min_samples_split", "min_samples_split"),
                       ("Forest min_samples_leaf", "min_samples_leaf"),
                       ("Forest bootstrap", "bootstrap"), ("Estimator seed", "random_state"),
                       ("Forest workers", "n_jobs"), ("Forest criterion", "criterion")):
        expected[label] = [str(params[key])]
    draws = {score["test_bootstrap"]["draws"] for split in report["splits"].values()
             for score in split["models"].values()}
    expected["Test bootstrap draws"] = [str(next(iter(draws)))] if len(draws) == 1 else []
    return [(f"Gate 3 baseline table: {label}", rows.get(label) == [values])
            for label, values in expected.items()]


def gate3_contrast_checks(text: str, report: dict) -> list[tuple[str, bool]]:
    """Re-derive the quoted split and floor contrasts from the logged RMSEs.

    These are differences and ratios of logged point estimates, computed at full
    precision rather than from the rounded table cells, so a hand-written figure
    that was never derived from this run fails here.
    """
    rows = _table_rows(text)
    rmse = {(split, model, side): report["splits"][split]["models"][model][side]["rmse"]
            for split in ("random", "scaffold") for model in ("mean", "descriptor_rf")
            for side in ("train", "test")}
    expected = {
        "Scaffold minus random test RMSE, Descriptor RF":
            f"{rmse['scaffold', 'descriptor_rf', 'test'] - rmse['random', 'descriptor_rf', 'test']:+.3f}",
        "Scaffold minus random test RMSE, Training mean":
            f"{rmse['scaffold', 'mean', 'test'] - rmse['random', 'mean', 'test']:+.3f}",
    }
    for split in ("random", "scaffold"):
        expected[f"Descriptor RF test RMSE reduction vs mean, {split}"] = (
            f"{1 - rmse[split, 'descriptor_rf', 'test'] / rmse[split, 'mean', 'test']:.1%}")
        expected[f"Descriptor RF test/train RMSE ratio, {split}"] = (
            f"{rmse[split, 'descriptor_rf', 'test'] / rmse[split, 'descriptor_rf', 'train']:.2f}")
    return [(f"Gate 3 contrast: {label}", rows.get(label) == [[value]])
            for label, value in expected.items()]


def load_gate4_controls(directory: Path) -> dict:
    """Read completed controls, verifying the search identity and report checksum."""
    return _load_gate4_searches(directory, "descriptors")


def load_gate4_combined(directory: Path) -> dict:
    """Return each completed combined search while the other split may still be training."""
    return _load_gate4_searches(directory, "combined")


def _load_gate4_searches(directory: Path, block: str) -> dict:
    def digest(value):
        encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), allow_nan=False)
        return hashlib.sha256(encoded.encode()).hexdigest()

    runs = {}
    for split in ("scaffold", "random"):
        path = directory / split / block
        manifest = json.loads((path / "manifest.json").read_text())
        if block == "combined" and not (path / "report.json").exists():
            continue
        record = json.loads((path / "report.json").read_text())
        report = record["result"]
        identity = digest(manifest)
        if (record["name"] != "report" or record["identity"] != identity
                or record["sha256"] != digest(report)
                or report["identity_sha256"] != identity):
            raise ValueError(f"Gate 4 {split}/{block} checksum/identity mismatch")
        if (any(row["partition"] != "train" for row in manifest["rows"])
                or report["heldout_evaluations"] != 0):
            raise ValueError("Gate 4 search summary requires training-only results")
        runs[split] = (manifest, report)
    return runs


def gate4_control_checks(text: str, runs: dict) -> list[tuple[str, bool]]:
    """Compare training-only OOF scores, selections and compute figures with logs."""
    return _gate4_training_checks(text, runs, combined=False)


def gate4_combined_checks(text: str, runs: dict) -> list[tuple[str, bool]]:
    """Guard the combined search table, including fit-local representation choices."""
    if set(runs) != {"scaffold", "random"}:
        summary_rows = [label for label in _table_rows(text)
                        if label.startswith("Combined ") or "Tuned descriptors+Morgan (OOF)" in label]
        return [("Gate 4 combined score table awaits both completed reports", not summary_rows)]
    return _gate4_training_checks(text, runs, combined=True)


def _gate4_training_checks(text: str, runs: dict, *, combined: bool) -> list[tuple[str, bool]]:
    rows, expected = _table_rows(text), {}
    prefix = "Combined" if combined else "Control"
    model_label = "Tuned descriptors+Morgan" if combined else "Tuned descriptor control"

    def parameters(choice):
        return f"{choice['max_features']} / {choice['min_samples_leaf']}"

    protocol = {}
    for split in ("scaffold", "random"):
        manifest, report = runs[split]
        cells = [f"{len(manifest['rows']):,}"]
        for metric in ("r2", "rmse", "mae"):
            value = report["oof_metrics"][metric]
            ci = report["oof_bootstrap"]["metrics"][metric]["ci95"]
            cells.append(f"{value:.3f} [{ci[0]:.3f}, {ci[1]:.3f}]"
                         if value is not None and ci else "undefined")
        cells.append(f"{report['outer_fold_spread']['rmse']['sample_sd']:.3f}")
        expected[f"{split} / {model_label} (OOF)"] = cells
        protocol[split] = {
            "Control outer folds": str(len(manifest["outer_plan"])),
            "Control inner folds per outer fold": " / ".join(
                str(len(plan)) for plan in manifest["inner_plan"]),
            "Control final-selection folds": str(len(manifest["final_plan"])),
            "Control grid candidates": str(len(manifest["candidates"])),
            "Control forest trees": str(manifest["estimator"]["n_estimators"]),
            ("Control final-fit feature columns" if combined else "Control feature columns"):
                str(len(manifest["descriptor_columns"]) + report["final_fpSize"]
                if combined else manifest["feature_shape"][1]),
            "Control seed": str(manifest["seed"]),
            "Control bootstrap draws": str(report["oof_bootstrap"]["draws"]),
            "Control bootstrap scaffold groups": str(report["oof_bootstrap"]["scaffold_groups"]),
            "Control valid draws (R² / RMSE / MAE)": " / ".join(
                str(report["oof_bootstrap"]["metrics"][metric]["valid_draws"])
                for metric in ("r2", "rmse", "mae")),
            "Control full-training selection (max_features / leaf)":
                parameters(report["final_selection"]["selected_parameters"]),
            "Control distinct outer selections": str(report["selection_stability"]["distinct_selections"]),
            "Control all outer selections agree": str(report["selection_stability"]["all_outer_folds_agree"]),
            "Control cumulative fit/predict time (s)":
                f"{report['timing']['cumulative_fit_predict_seconds']:.2f}",
            "Control completion invocation wall time (s)":
                f"{report['timing']['completion_invocation_seconds']:.2f}",
            "Control report n_jobs": str(report["timing"]["n_jobs"]),
            "Control search workers": str(report["timing"]["search_workers"]),
            "Control held-out evaluations": str(report["heldout_evaluations"]),
        }
        for fold in report["outer_folds"]:
            label = f"Control outer fold {fold['outer_fold'] + 1} selection (max_features / leaf)"
            protocol[split][label] = parameters(fold["selection"]["selected_parameters"])
        if combined:
            # Every estimator's own bit length, not just the final fit's choice.
            selections = [f["selection"] for f in report["outer_folds"]] + [report["final_selection"]]
            bits = {b for s in selections for c in s["candidates"] for b in c["fold_fpSizes"]}
            bits.update(f["fpSize"] for f in report["outer_folds"])
            bits.add(report["final_fpSize"])
            protocol[split]["Control fit-local bit lengths"] = " / ".join(map(str, sorted(bits)))
            protocol[split]["Control final-fit Morgan bits"] = str(report["final_fpSize"])
        protocol[split] = {key.replace("Control", prefix, 1): value
                           for key, value in protocol[split].items()}
    for label in protocol["scaffold"]:
        expected[label] = [protocol[split][label] for split in ("scaffold", "random")]
    return [(f"Gate 4 {prefix.lower()} table: {label}", rows.get(label) == [values])
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
    baselines = json.loads((ROOT / "results/baselines/CHEMBL284/baselines.json").read_text())
    checks.extend(gate3_baseline_checks((ROOT / "DATA_CARD.md").read_text(), baselines))
    # The Gate 3 handoff repeats the same table, so it is held to the same logs.
    signoff = (ROOT / "GATE3_SIGNOFF.md").read_text()
    checks.extend((f"handoff {name}", ok)
                  for name, ok in gate3_baseline_checks(signoff, baselines))
    checks.extend(gate3_contrast_checks(signoff, baselines))
    checks.extend(gate4_control_checks(
        (ROOT / "DATA_CARD.md").read_text(),
        load_gate4_controls(ROOT / "results/gate4_training/CHEMBL284")))
    checks.extend(gate4_combined_checks(
        (ROOT / "DATA_CARD.md").read_text(),
        load_gate4_combined(ROOT / "results/gate4_training/CHEMBL284")))

    for name, ok in checks:
        print(f"{'OK  ' if ok else 'FAIL'} {name}")
    failed = [name for name, ok in checks if not ok]
    print(f"{len(checks) - len(failed)}/{len(checks)} data-card and sign-off checks match the logged run")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
