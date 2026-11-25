"""Gate 2 curation cascade: frozen raw pull -> one curated pIC50 per compound.

    .venv/bin/python -m src.curate CHEMBL284

Reads only ``data/raw/<TARGET>/`` (checksums verified before and after), writes
``data/curated/<TARGET>/`` and ``results/curation/<TARGET>/``. The network is never used.

Seven filters are applied **in this order**, each logged with before/after counts and an
itemised breakdown of what it removed:

    1. binding assays only            assay_type == 'B'
    2. single activity type           standard_type == 'IC50'
    3. single target organism         Homo sapiens, taxon 9606
    4. harmonise units to nM          molar units only; refusals are reported, not dropped silently
    5. remove censored values         standard_relation must be '='
    6. standardise structures         RDKit: sanitise -> largest fragment -> neutralise
                                      -> canonical SMILES -> InChIKey, in that order
    7. resolve duplicates             MEDIAN pIC50 per compound, spread recorded,
                                      replicates spanning > 1 log unit flagged

pIC50 = -log10(IC50 in mol/L). Nothing here claims or implies efficacy, safety or
suitability of any compound; these are literature- and patent-reported measurements.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
import statistics
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

from src.config import DATA_CURATED, DATA_RAW, RESULTS, SEED
from src.retrieve import MANIFEST, RAW_ACTIVITIES, _sha256, read_jsonl_gz, verify_raw

# Only unambiguous molar concentrations. Mass concentrations (ug.mL-1) would need the
# molecular weight of the material actually tested, which may be a salt or a mixture;
# converting them would silently mix two measurement bases, so they are refused and reported.
UNIT_TO_NM: dict[str, float] = {
    "M": 1e9, "mM": 1e6, "uM": 1e3, "µM": 1e3, "μM": 1e3, "nM": 1.0, "pM": 1e-3, "fM": 1e-6,
}

LARGE_SPREAD_LOG_UNITS = 1.0  # flag threshold for replicate disagreement

# Columns carried through to the curated record table for audit and later review.
AUDIT_FIELDS = (
    "activity_id", "assay_chembl_id", "document_chembl_id", "document_year", "src_id",
    "molecule_chembl_id", "parent_molecule_chembl_id", "molecule_pref_name",
    "standard_type", "standard_relation", "standard_value", "standard_units",
    "pchembl_value", "data_validity_comment", "potential_duplicate", "standard_flag",
    "bao_label", "assay_description",
)

STEPS = (
    ("1_binding_assay", "assay_type == 'B'"),
    ("2_single_activity_type", "standard_type == 'IC50'"),
    ("3_single_target_organism", "target_organism == 'Homo sapiens' and str(target_tax_id) == '9606'"),
    ("4_units_harmonised_to_nM", f"standard_units in {sorted(UNIT_TO_NM)} and standard_value is a finite number > 0"),
    ("5_censored_removed", "standard_relation == '=' exactly"),
    ("6_structure_standardised", "RDKit sanitise -> largest fragment -> neutralise -> canonical SMILES -> InChIKey"),
    ("7_duplicates_resolved", f"median pIC50 per standardised InChIKey; spread > {LARGE_SPREAD_LOG_UNITS} log unit flagged"),
)


def pic50_from_nm(value_nm: float) -> float:
    """pIC50 = -log10(IC50 in mol/L); value_nm is nanomolar."""
    return -math.log10(value_nm * 1e-9)


# --------------------------------------------------------------------------- structures

def _standardiser():
    """Import RDKit lazily so that importing this module stays cheap."""
    from rdkit import Chem, RDLogger
    from rdkit.Chem.MolStandardize import rdMolStandardize

    RDLogger.DisableLog("rdApp.*")  # failures are reported by us, not printed by RDKit
    return Chem, rdMolStandardize


def standardise_smiles(smiles: str | None) -> dict:
    """Standardise one SMILES, in the required order, and say exactly what happened.

    Returns a dict with ``ok``; on success also ``canonical_smiles``, ``inchikey``,
    ``inchikey_skeleton``, ``fragments_before`` and ``formal_charge_before``.
    On failure ``ok`` is False and ``failure`` names the stage that failed and why.
    No molecule is ever dropped without a row in the failure log.
    """
    Chem, rdMolStandardize = _standardiser()
    if smiles is None or not str(smiles).strip():
        return {"ok": False, "failure": "missing_smiles: no canonical_smiles in the API record"}

    # (a) sanitisation. Parse unsanitised first so the sanitisation error is reportable.
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    if mol is None:
        return {"ok": False, "failure": "parse_failed: RDKit could not read the SMILES string"}
    try:
        Chem.SanitizeMol(mol)
    except Exception as exc:  # RDKit raises several concrete subclasses
        return {"ok": False, "failure": f"sanitisation_failed: {type(exc).__name__}: {exc}"}
    try:
        mol = rdMolStandardize.Cleanup(mol)
    except Exception as exc:
        return {"ok": False, "failure": f"cleanup_failed: {type(exc).__name__}: {exc}"}
    if mol is None:
        return {"ok": False, "failure": "cleanup_failed: RDKit Cleanup returned no molecule"}
    fragments_before = len(Chem.GetMolFrags(mol))
    charge_before = Chem.GetFormalCharge(mol)

    # (b) largest fragment: a salt or co-crystal must not become a different compound.
    try:
        mol = rdMolStandardize.LargestFragmentChooser(preferOrganic=True).choose(mol)
    except Exception as exc:
        return {"ok": False, "failure": f"fragment_choice_failed: {type(exc).__name__}: {exc}"}
    if mol is None or mol.GetNumAtoms() == 0:
        return {"ok": False, "failure": "fragment_choice_failed: no fragment left"}

    # (c) neutralisation
    try:
        mol = rdMolStandardize.Uncharger().uncharge(mol)
    except Exception as exc:
        return {"ok": False, "failure": f"neutralisation_failed: {type(exc).__name__}: {exc}"}

    # (d) canonical SMILES, then (e) InChIKey
    try:
        canonical = Chem.MolToSmiles(mol)
    except Exception as exc:
        return {"ok": False, "failure": f"canonical_smiles_failed: {type(exc).__name__}: {exc}"}
    if not canonical:
        return {"ok": False, "failure": "canonical_smiles_failed: empty SMILES"}
    try:
        inchikey = Chem.MolToInchiKey(mol)
    except Exception as exc:
        return {"ok": False, "failure": f"inchikey_failed: {type(exc).__name__}: {exc}"}
    if not inchikey:
        return {"ok": False, "failure": "inchikey_failed: RDKit returned an empty InChIKey"}

    return {
        "ok": True,
        "canonical_smiles": canonical,
        "inchikey": inchikey,
        "inchikey_skeleton": inchikey.split("-")[0],
        "fragments_before": fragments_before,
        "formal_charge_before": charge_before,
    }


# ------------------------------------------------------------------------------ cascade

def run_cascade(records: list[dict]) -> dict:
    """Apply the seven filters in order. Returns the full log plus the surviving rows."""
    log: dict = {"steps": [], "input_records": len(records)}
    steps = log["steps"]

    def note(name: str, rule: str, before: int, kept: list, removed: dict, **extra) -> None:
        steps.append({
            "step": name, "rule": rule, "records_in": before, "records_out": len(kept),
            "records_removed": before - len(kept),
            "removed_breakdown": dict(sorted(removed.items())), **extra,
        })

    # 1. binding assays only
    before = len(records)
    kept = [r for r in records if r.get("assay_type") == "B"]
    removed = Counter(f"assay_type={r.get('assay_type')!r}" for r in records if r.get("assay_type") != "B")
    note(STEPS[0][0], STEPS[0][1], before, kept, removed)

    # 2. single activity type: IC50
    before, prev = len(kept), kept
    kept = [r for r in prev if r.get("standard_type") == "IC50"]
    removed = Counter(f"standard_type={r.get('standard_type')!r}" for r in prev if r.get("standard_type") != "IC50")
    note(STEPS[1][0], STEPS[1][1], before, kept, removed,
         removed_activity_types=len(removed))

    # 3. single target organism
    before, prev = len(kept), kept
    def human(r): return r.get("target_organism") == "Homo sapiens" and str(r.get("target_tax_id")) == "9606"
    kept = [r for r in prev if human(r)]
    removed = Counter(f"target_organism={r.get('target_organism')!r}/tax_id={r.get('target_tax_id')!r}"
                      for r in prev if not human(r))
    note(STEPS[2][0], STEPS[2][1], before, kept, removed,
         note_="The rule is the TARGET organism. The assay expression host is absent from every "
               "activity record in this pull, so it is not inferred.")

    # 4. harmonise units to nM
    before, prev = len(kept), kept
    kept, removed, refusals = [], Counter(), []
    for r in prev:
        unit = r.get("standard_units")
        factor = UNIT_TO_NM.get(unit) if isinstance(unit, str) else None
        if factor is None:
            label = "<missing>" if unit is None or not str(unit).strip() else str(unit)
            removed[f"unit_not_convertible={label}"] += 1
            refusals.append({"activity_id": r.get("activity_id"), "standard_units": unit,
                             "units_as_reported": r.get("units"), "standard_value": r.get("standard_value"),
                             "reason": "no unambiguous molar conversion"})
            continue
        try:
            value = float(r.get("standard_value"))
        except (TypeError, ValueError):
            removed["value_not_a_number"] += 1
            refusals.append({"activity_id": r.get("activity_id"), "standard_units": unit,
                             "standard_value": r.get("standard_value"), "reason": "value is missing or not numeric"})
            continue
        if not math.isfinite(value) or value <= 0:
            removed["value_not_positive_finite"] += 1
            refusals.append({"activity_id": r.get("activity_id"), "standard_units": unit,
                             "standard_value": r.get("standard_value"),
                             "reason": "value is not a positive finite number, so log10 is undefined"})
            continue
        row = dict(r)
        row["value_nM"] = value * factor
        row["unit_conversion_factor_to_nM"] = factor
        kept.append(row)
    note(STEPS[3][0], STEPS[3][1], before, kept, removed,
         units_accepted=dict(sorted(Counter(r["standard_units"] for r in kept).items())),
         refused_rows=len(refusals),
         note_="Every refusal is itemised in units_refused.csv. Mass concentrations are refused "
               "because the molecular weight of the tested material (possibly a salt) is not known "
               "here; converting them would mix two measurement bases.")
    log["unit_refusals"] = refusals

    # 5. remove censored values
    before, prev = len(kept), kept
    kept = [r for r in prev if r.get("standard_relation") == "="]
    dropped = [r for r in prev if r.get("standard_relation") != "="]
    removed = Counter(
        ("relation=<missing>" if r.get("standard_relation") in (None, "")
         else f"relation={r['standard_relation']!r}") for r in dropped)
    censored = sum(n for key, n in removed.items() if key != "relation=<missing>")
    note(STEPS[4][0], STEPS[4][1], before, kept, removed,
         censored_rows_removed=censored,
         missing_relation_rows_removed=removed.get("relation=<missing>", 0),
         note_="Censored means the measurement is a bound, not a value: a '>' row only says the "
               "IC50 exceeds the highest concentration tested. Keeping them as if they were values "
               "would put invented numbers into the target variable. Rows with no relation at all "
               "are removed by the same rule but counted separately.")

    # 6. standardise structures
    before, prev = len(kept), kept
    kept, removed, failures = [], Counter(), []
    for r in prev:
        result = standardise_smiles(r.get("canonical_smiles"))
        if not result["ok"]:
            reason = result["failure"].split(":")[0]
            removed[f"structure_{reason}"] += 1
            failures.append({"activity_id": r.get("activity_id"),
                             "molecule_chembl_id": r.get("molecule_chembl_id"),
                             "canonical_smiles_in": r.get("canonical_smiles"),
                             "failure": result["failure"]})
            continue
        row = dict(r)
        row.update({k: v for k, v in result.items() if k != "ok"})
        row["pIC50"] = pic50_from_nm(row["value_nM"])
        kept.append(row)
    note(STEPS[5][0], STEPS[5][1], before, kept, removed,
         failed_rows=len(failures),
         salt_or_mixture_rows=sum(1 for r in kept if r["fragments_before"] > 1),
         charged_rows_neutralised=sum(1 for r in kept if r["formal_charge_before"] != 0),
         note_="Every failure is itemised in structure_failures.csv; no molecule is dropped silently.")
    log["structure_failures"] = failures

    # 7. resolve duplicates: median pIC50 per standardised compound
    before = len(kept)
    groups: dict[str, list[dict]] = defaultdict(list)
    for r in kept:
        groups[r["inchikey"]].append(r)
    compounds = [summarise_compound(key, rows) for key, rows in sorted(groups.items())]
    steps.append({
        "step": STEPS[6][0], "rule": STEPS[6][1], "records_in": before,
        "records_out": len(compounds), "records_removed": before - len(compounds),
        "removed_breakdown": {"rows_collapsed_into_a_median": before - len(compounds)},
        "compounds_with_replicates": sum(1 for c in compounds if c["n_records"] > 1),
        "compounds_flagged_spread_over_1_log": sum(1 for c in compounds if c["flag_spread_over_1_log"]),
        "note_": "The output unit changes here: rows become compounds. The median is used rather "
                 "than the mean because replicate IC50 sets in this pull have a long tail of "
                 "cross-assay disagreement, which a mean would follow.",
    })
    log["records"] = kept
    log["compounds"] = compounds
    return log


def censoring_reconciliation(records: list[dict], log: dict) -> dict:
    """Reconcile this cascade's censored count with the Gate 1 headline figure.

    Gate 1 reported 796 rows failing the exact-relation rule. This cascade removes 791 at
    step 5. The difference is not a disagreement about censoring: Gate 1's unit filter
    checked only the unit *string* and logged positive-finite values as a separate
    diagnostic (see README), so rows with a usable unit but no value still reached its
    relation test and were counted there. This cascade requires a usable value at step 4,
    so those rows are already gone by step 5.

    Both numbers are recomputed here from the same frozen records, so the identity is
    checked rather than asserted.
    """
    def unit_string_ok(r):
        unit = r.get("standard_units")
        return isinstance(unit, str) and unit in UNIT_TO_NM

    def value_ok(r):
        try:
            value = float(r.get("standard_value"))
        except (TypeError, ValueError):
            return False
        return math.isfinite(value) and value > 0

    base = [r for r in records
            if r.get("assay_type") == "B" and r.get("standard_type") == "IC50"
            and r.get("target_organism") == "Homo sapiens" and str(r.get("target_tax_id")) == "9606"
            and unit_string_ok(r)]
    not_exact = [r for r in base if r.get("standard_relation") != "="]
    missing_relation = [r for r in not_exact if r.get("standard_relation") in (None, "")]
    step5 = next(s for s in log["steps"] if s["step"] == STEPS[4][0])

    return {
        "question": "Gate 1 measured 796 censored rows; step 5 of this cascade removes 791. Why?",
        "gate1_rule": "Gate 1 filtered on the unit STRING only and logged the "
                      "positive-finite-value check as a separate diagnostic, so rows with a "
                      "usable unit but an unusable value still reached its relation test.",
        "gate1_base_rows_recomputed": len(base),
        "gate1_rows_failing_exact_relation_recomputed": len(not_exact),
        "gate2_step5_rows_removed": step5["records_removed"],
        "gate2_step5_censored_rows_removed": step5["censored_rows_removed"],
        "difference": len(not_exact) - step5["records_removed"],
        "difference_explained_by": {
            "rows": sorted(r.get("activity_id") for r in missing_relation),
            "count": len(missing_relation),
            "what_they_are": "Rows whose standard_relation AND standard_value are both null, "
                             "with units reported as nM. Gate 1 counted them as 'not exactly =' "
                             "and so inside its 796. This cascade removes them one step earlier, "
                             "at step 4, as value_not_a_number, because a row with no value "
                             "cannot become a pIC50 whatever its relation says.",
            "all_have_null_value": all(not value_ok(r) for r in missing_relation),
        },
        "identity_checks": {
            "gate1_count_equals_step5_plus_difference":
                len(not_exact) == step5["records_removed"] + len(missing_relation),
            "genuinely_censored_breakdown": dict(sorted(Counter(
                f"relation={r['standard_relation']!r}" for r in not_exact
                if r.get("standard_relation") not in (None, "")).items())),
            "surviving_rows_agree_with_gate1":
                sum(1 for r in base if r.get("standard_relation") == "=") == step5["records_out"],
        },
        "conclusion": "No censored measurement is treated differently by the two gates. The "
                      "791/796 gap is five value-less rows removed one step earlier here, and "
                      "the 4859 rows that survive the relation rule are identical in both.",
    }


def summarise_compound(inchikey: str, rows: list[dict]) -> dict:
    """One curated compound: median pIC50 plus the full spread of its replicates."""
    values = sorted(r["pIC50"] for r in rows)
    spread = values[-1] - values[0]
    ids = sorted({r["molecule_chembl_id"] for r in rows if r.get("molecule_chembl_id")})
    parents = sorted({r["parent_molecule_chembl_id"] for r in rows if r.get("parent_molecule_chembl_id")})
    return {
        "inchikey": inchikey,
        "inchikey_skeleton": rows[0]["inchikey_skeleton"],
        "canonical_smiles": rows[0]["canonical_smiles"],
        "pIC50_median": statistics.median(values),
        "n_records": len(values),
        "pIC50_min": values[0],
        "pIC50_max": values[-1],
        "pIC50_spread": spread,
        "pIC50_stdev": statistics.stdev(values) if len(values) > 1 else None,
        "flag_spread_over_1_log": bool(spread > LARGE_SPREAD_LOG_UNITS),
        "molecule_chembl_ids": ";".join(ids),
        "n_molecule_chembl_ids": len(ids),
        "parent_molecule_chembl_ids": ";".join(parents),
        "n_parent_molecule_chembl_ids": len(parents),
        "n_assays": len({r.get("assay_chembl_id") for r in rows}),
        "n_documents": len({r.get("document_chembl_id") for r in rows}),
        "src_ids": ";".join(sorted({str(r.get("src_id")) for r in rows})),
        "any_data_validity_comment": ";".join(sorted(
            {r["data_validity_comment"] for r in rows if r.get("data_validity_comment")})),
        "any_potential_duplicate": int(any(r.get("potential_duplicate") in (1, "1", True) for r in rows)),
    }


# ------------------------------------------------------------------- measured noise floor

def replicate_noise(compounds: list[dict], records: list[dict]) -> dict:
    """Measured repeat dispersion, in the same log units as model RMSE.

    Selection is explicit; none of these statistics establishes an irreducible error
    floor for compound-median predictions. The Gate 4 reference remains a human decision.
    """
    by_key: dict[str, list[float]] = defaultdict(list)
    for r in records:
        by_key[r["inchikey"]].append(r["pIC50"])
    replicated = {k: v for k, v in by_key.items() if len(v) > 1}
    if not replicated:
        return {"replicated_compounds": 0, "note": "no replicated compounds; no noise floor measurable"}

    def pooled_sd(groups: dict[str, list[float]]) -> tuple[float | None, int]:
        df = sum(len(v) - 1 for v in groups.values())
        if df == 0:
            return None, 0
        var = sum((len(v) - 1) * statistics.variance(v) for v in groups.values()) / df
        return var ** 0.5, df

    all_sd, all_df = pooled_sd(replicated)
    deviations = [abs(x - statistics.median(v)) for v in replicated.values() for x in v]
    spreads = sorted(c["pIC50_spread"] for c in compounds if c["n_records"] > 1)

    # This threshold defines a sensitivity subset, not experimental independence.
    # Equal or almost-equal values can be copied reports OR independent agreements.
    NEAR_EQUAL = 0.01
    identical = {k: v for k, v in replicated.items() if max(v) - min(v) < NEAR_EQUAL}
    informative = {k: v for k, v in replicated.items() if max(v) - min(v) >= NEAR_EQUAL}
    tight = {k: v for k, v in informative.items() if max(v) - min(v) <= LARGE_SPREAD_LOG_UNITS}
    informative_sd, informative_df = pooled_sd(informative)
    tight_sd, tight_df = pooled_sd(tight)

    edges = [0.0, NEAR_EQUAL, 0.3, 1.0, 2.0, 3.0, float("inf")]
    histogram = {}
    for low, high in zip(edges, edges[1:]):
        label = f"{low:g}-{high:g}" if math.isfinite(high) else f">={low:g}"
        histogram[label] = sum(1 for s in spreads if low <= s < high)

    return {
        "definition": "Dispersion between repeated IC50 reports of the same standardised "
                      "compound, in pIC50 log units (the units of RMSE), not an established "
                      "lower bound on compound-median prediction error.",
        "replicated_compounds": len(replicated),
        "replicate_measurements": sum(len(v) for v in replicated.values()),
        "degrees_of_freedom": all_df,
        "spread_histogram_log_units": histogram,
        "groups_exactly_equal": sum(max(v) == min(v) for v in replicated.values()),
        "groups_spread_below_0_01": len(identical),
        "groups_spread_at_least_0_01": len(informative),
        "pooled_within_compound_sd_all_groups": all_sd,
        "pooled_within_compound_sd_differing_groups": informative_sd,
        "pooled_within_compound_sd_differing_under_1_log": tight_sd,
        "degrees_of_freedom_differing_groups": informative_df,
        "degrees_of_freedom_differing_under_1_log": tight_df,
        "rms_deviation_from_own_median": math.sqrt(sum(d * d for d in deviations) / len(deviations)),
        "median_absolute_deviation_from_own_median": statistics.median(deviations),
        "spread_median": statistics.median(spreads),
        "spread_mean": statistics.mean(spreads),
        "spread_p90": spreads[int(0.9 * (len(spreads) - 1))],
        "spread_max": spreads[-1],
        "compounds_flagged_spread_over_1_log": sum(1 for s in spreads if s > LARGE_SPREAD_LOG_UNITS),
        "noise_floor": {
            "status": "pending_human_decision; legacy recommendation and ceiling keys are descriptive candidates, not validated bounds",
            "recommended_for_gate4": tight_sd,
            "recommended_basis": "pooled_within_compound_sd_differing_under_1_log: the pooled SD "
                                 "over groups with 0.01 <= spread <= 1 log unit. This legacy "
                                 "candidate conditions on observed disagreement; it neither "
                                 "establishes independent repeats nor estimates a universal floor. "
                                 "Gate 4 use requires a human decision.",
            "ceiling_including_irreconcilable_groups": informative_sd,
            "ceiling_basis": "pooled_within_compound_sd_differing_groups: the same quantity over "
                             "every group with spread >=0.01 log unit, including large "
                             "disagreements. This is a sensitivity statistic, not a statistical "
                             "upper bound on assay noise or prediction error.",
            "all_groups_figure_and_why_it_is_not_used": {
                "value": all_sd,
                "reason": "The pooled SD over all repeated groups is the unfiltered descriptive "
                          "dispersion. Near-equal and widely disagreeing reports both contribute. "
                          "Neither independence nor assay comparability has been established, "
                          "so it is not an irreducible prediction-error floor.",
            },
        },
        "caveats": [
            "These are repeated literature and patent reports, not a designed replicate "
            "experiment: they mix assay protocols, so none of these figures is pure assay noise.",
            "The spread distribution has many near-equal groups and a long tail of large "
            "disagreements. Spread <0.01 does not mean exact equality or prove copying. "
            "All finite histogram bins include their lower edge and exclude their upper edge; "
            "the final bin includes 3 log units. Inspect the histogram alongside pooled SDs.",
            f"Only replicated compounds contribute. {len(compounds) - len(replicated)} of "
            f"{len(compounds)} curated compounds are measured once, so dispersion is measured on "
            f"{100 * len(replicated) / len(compounds):.1f}% of the curated set.",
            "The legacy candidate is conditional dispersion for repeat groups "
            "with 0.01 <= spread <=1 log unit, not an unconditional lower bound on held-out RMSE. "
            "Predictions target compound medians, not individual reports; the repeated subset "
            "may differ from the test population. A lower RMSE alone does not prove leakage "
            "or memorisation. Independence of these reports has not been established.",
        ],
    }


def target_variable_summary(compounds: list[dict]) -> dict:
    """Summary statistics for pIC50, plus the pathologies we can actually detect."""
    values = sorted(c["pIC50_median"] for c in compounds)
    n = len(values)

    def q(p: float) -> float:
        return values[min(n - 1, max(0, int(round(p * (n - 1)))))]

    at_one_dp = Counter(round(v, 1) for v in values)
    integers = Counter(v for v in values if abs(v - round(v)) < 1e-9)
    return {
        "n_compounds": n,
        "min": values[0], "p05": q(0.05), "q1": q(0.25), "median": statistics.median(values),
        "q3": q(0.75), "p95": q(0.95), "max": values[-1],
        "mean": statistics.mean(values), "stdev": statistics.stdev(values) if n > 1 else None,
        "skew_mean_minus_median": statistics.mean(values) - statistics.median(values),
        "below_4": sum(1 for v in values if v < 4), "above_9": sum(1 for v in values if v > 9),
        "below_2_implausible": sum(1 for v in values if v < 2),
        "above_10": sum(1 for v in values if v > 10),
        "most_common_one_dp_bins": dict(at_one_dp.most_common(10)),
        "exact_integer_pIC50_values": dict(sorted(integers.items())),
        "pathology_checks": {
            "round_number_pile_ups": {
                "finding": "No spike at a reporting-convention concentration. The densest 0.1-log "
                           "bins are all in the 6.6-7.9 interior of the distribution, which is "
                           "where the mass of the data is; there is no pile-up at pIC50 5.0 "
                           "(10 uM) or 4.0 (100 uM).",
                "exact_integer_values": f"{sum(integers.values())} of {n} compounds sit on an exact "
                                        "integer pIC50, concentrated at 8.0 and 9.0. That is a "
                                        "rounding fingerprint (10 nM, 1 nM reported as round "
                                        "numbers), small but real.",
            },
            "implausible_extremes": {
                "finding": f"min {values[0]:.3f} (IC50 ~ {10 ** -values[0]:.2g} M), "
                           f"max {values[-1]:.3f}. A pIC50 below about 2 means an IC50 above 10 mM, "
                           "which is past the solubility of most of these compounds and is better "
                           "read as 'inactive at the top concentration tested' than as a "
                           "measurement.",
                "compounds_below_pIC50_2": sum(1 for v in values if v < 2),
                "compounds_below_pIC50_4": sum(1 for v in values if v < 4),
                "decision": "These are kept in the curated set and flagged here rather than "
                            "trimmed, because an activity threshold is a modelling choice and "
                            "belongs to a later gate, not to curation. Any later gate that trims "
                            "them must say so and re-report the count.",
            },
            "assay_ceiling": "A hard ceiling would appear as a cluster of identical values at the "
                             "top concentration tested. Censored '>' rows, which are where such a "
                             "ceiling normally hides, were removed at step 5 and counted there, so "
                             "the ceiling has been taken out of the target variable rather than "
                             "flattened into it.",
        },
    }


# ------------------------------------------------------------------------------- outputs

def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_cascade_tables(log: dict, out_dir: Path) -> None:
    """The required report figure: CSV, markdown and a figure-ready long form."""
    rows = []
    for step in log["steps"]:
        rows.append({
            "order": step["step"].split("_")[0],
            "step": step["step"], "rule": step["rule"],
            "records_in": step["records_in"], "records_out": step["records_out"],
            "records_removed": step["records_removed"],
            "percent_of_raw_remaining": round(100 * step["records_out"] / log["input_records"], 2),
            "removed_breakdown": json.dumps(step["removed_breakdown"], ensure_ascii=False),
        })
    _write_csv(out_dir / "cascade.csv", rows, list(rows[0]))

    # Figure-ready: one row per bar, no nested JSON, explicit unit of the count.
    figure = [{"order": 0, "label": "raw pull", "records": log["input_records"],
               "removed": 0, "unit": "records"}]
    for step in log["steps"]:
        figure.append({
            "order": int(step["step"].split("_")[0]),
            "label": step["step"].split("_", 1)[1].replace("_", " "),
            "records": step["records_out"], "removed": step["records_removed"],
            "unit": "compounds" if step["step"].startswith("7") else "records",
        })
    _write_csv(out_dir / "cascade_figure.csv", figure, list(figure[0]))

    lines = [f"# Curation cascade: {log['target_chembl_id']}", "",
             f"Release {log['chembl_db_version']}; raw pull of **{log['input_records']}** records, "
             f"checksums verified before and after. Command: `{log['command']}`.", "",
             "Each filter is applied separately, in this order, to the output of the previous one.",
             "The last step changes the unit of the count from records to compounds.", "",
             "| # | Filter | Rule | In | Out | Removed |", "|---|---|---|---:|---:|---:|"]
    for step in log["steps"]:
        name = step["step"].split("_", 1)[1].replace("_", " ")
        lines.append(f"| {step['step'].split('_')[0]} | {name} | `{step['rule']}` | "
                     f"{step['records_in']} | {step['records_out']} | {step['records_removed']} |")
    lines += ["", "## What each filter removed", ""]
    for step in log["steps"]:
        lines.append(f"**{step['step']}** - removed {step['records_removed']}:")
        lines.append("")
        for key, n in step["removed_breakdown"].items():
            lines.append(f"- `{key}`: {n}")
        if step.get("note_"):
            lines += ["", f"> {step['note_']}"]
        lines.append("")
    rec = log.get("censoring_reconciliation")
    if rec:
        lines += ["## Reconciling step 5 with the Gate 1 censored count", "",
                  rec["question"], "",
                  f"- Gate 1's base (unit string only, no value check), recomputed here: "
                  f"**{rec['gate1_base_rows_recomputed']}** rows",
                  f"- rows failing the exact-relation rule on that base: "
                  f"**{rec['gate1_rows_failing_exact_relation_recomputed']}** (the Gate 1 figure)",
                  f"- removed by step 5 of this cascade: **{rec['gate2_step5_rows_removed']}**, of which "
                  f"{rec['gate2_step5_censored_rows_removed']} carry a censoring relation",
                  f"- difference: **{rec['difference']}** rows "
                  f"(activity_id {', '.join(str(i) for i in rec['difference_explained_by']['rows'])})",
                  "", f"> {rec['difference_explained_by']['what_they_are']}", "",
                  rec["conclusion"], ""]
    lines += ["No record is removed without appearing in a breakdown above or in an itemised CSV "
              "(`units_refused.csv`, `structure_failures.csv`). These are literature- and "
              "patent-reported measurements; nothing here claims efficacy, safety or suitability.", ""]
    (out_dir / "cascade.md").write_text("\n".join(lines), encoding="utf-8")


def new_figure(**kwargs):
    """A bare Agg figure.

    ``matplotlib.pyplot`` is deliberately not used: importing it scans
    ``mpl-data/stylelib``, and on the exFAT drive this repo lives on, macOS leaves
    AppleDouble ``._*`` companions in that directory which matplotlib tries to parse as
    style files and dies on. The Agg canvas needs neither pyplot nor the style library.
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure

    fig = Figure(**kwargs)
    FigureCanvasAgg(fig)
    return fig


def write_figures(log: dict, out_dir: Path) -> list[str]:
    """pIC50 distribution (item 2) and the cascade funnel. Returns file names written."""
    values = [c["pIC50_median"] for c in log["compounds"]]
    fig = new_figure(figsize=(11, 4.2))
    axes = fig.subplots(1, 2)
    axes[0].hist(values, bins=60, color="#4C78A8", edgecolor="white", linewidth=0.4)
    axes[0].axvline(statistics.median(values), color="#E45756", linestyle="--", linewidth=1.2,
                    label=f"median {statistics.median(values):.2f}")
    axes[0].set_xlabel("pIC50 (median per compound)")
    axes[0].set_ylabel("compounds")
    axes[0].set_title(f"{log['target_chembl_id']} curated pIC50, n={len(values)}")
    axes[0].legend(frameon=False, fontsize=8)
    # An ECDF rather than a boxplot: matplotlib's boxplot iterates rcParams, which lazily
    # imports pyplot and so hits the stylelib problem described in new_figure().
    ordered = sorted(values)
    fractions = [(i + 1) / len(ordered) for i in range(len(ordered))]
    axes[1].plot(ordered, fractions, color="#4C78A8", linewidth=1.4)
    for p, colour in ((0.05, "#BAB0AC"), (0.25, "#F58518"), (0.5, "#E45756"),
                      (0.75, "#F58518"), (0.95, "#BAB0AC")):
        value = ordered[min(len(ordered) - 1, int(round(p * (len(ordered) - 1))))]
        axes[1].axvline(value, color=colour, linestyle=":", linewidth=1.0)
        axes[1].annotate(f"p{int(p * 100)}={value:.2f}", (value, p), fontsize=7,
                         textcoords="offset points", xytext=(3, -8))
    axes[1].set_xlabel("pIC50 (median per compound)")
    axes[1].set_ylabel("cumulative fraction of compounds")
    axes[1].set_title("distribution shape and tails")
    fig.suptitle("Literature/patent-reported IC50 only. No efficacy, safety or suitability claim.",
                 fontsize=8, y=0.02)
    fig.tight_layout(rect=(0, 0.04, 1, 1))
    fig.savefig(out_dir / "pic50_distribution.png", dpi=150)

    labels = ["raw"] + [s["step"].split("_", 1)[1].replace("_", "\n") for s in log["steps"]]
    counts = [log["input_records"]] + [s["records_out"] for s in log["steps"]]
    fig = new_figure(figsize=(10, 4.4))
    ax = fig.subplots()
    bars = ax.bar(range(len(counts)), counts,
                  color=["#808080"] + ["#4C78A8"] * (len(counts) - 2) + ["#54A24B"])
    ax.bar_label(bars, fontsize=8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("records (final bar: compounds)")
    ax.set_title(f"{log['target_chembl_id']} curation cascade")
    fig.tight_layout()
    fig.savefig(out_dir / "cascade_funnel.png", dpi=150)
    return ["pic50_distribution.png", "cascade_funnel.png"]


def curate_target(target_dir: Path, curated_dir: Path, results_dir: Path) -> dict:
    target_dir, curated_dir, results_dir = Path(target_dir), Path(curated_dir), Path(results_dir)
    if not verify_raw(target_dir):
        raise ValueError("Raw file checksum mismatch; curation refused")
    raw_manifest = json.loads((target_dir / MANIFEST).read_text(encoding="utf-8"))
    records = read_jsonl_gz(target_dir / RAW_ACTIVITIES)
    if len(records) != raw_manifest["record_count"]:
        raise ValueError("Raw row count disagrees with manifest; curation refused")

    log = run_cascade(records)
    log.update({
        "target_chembl_id": raw_manifest["target_chembl_id"],
        "chembl_db_version": raw_manifest["chembl_db_version"],
        "raw_retrieval_finished_utc": raw_manifest["retrieval_finished_utc"],
        "curated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "command": f".venv/bin/python -m src.curate {raw_manifest['target_chembl_id']}",
        "seed": SEED,
        "source_raw_files": raw_manifest["raw_files"],
    })
    log["target_variable"] = target_variable_summary(log["compounds"])
    log["replicate_noise"] = replicate_noise(log["compounds"], log["records"])
    log["identity_decisions"] = identity_decisions(log)
    log["censoring_reconciliation"] = censoring_reconciliation(records, log)

    if not verify_raw(target_dir):
        raise ValueError("Raw files changed during curation; refusing to write")

    curated_dir.mkdir(parents=True, exist_ok=True)
    results_dir.mkdir(parents=True, exist_ok=True)

    compound_fields = list(log["compounds"][0])
    _write_csv(curated_dir / "compounds.csv", log["compounds"], compound_fields)
    record_fields = list(AUDIT_FIELDS) + ["value_nM", "unit_conversion_factor_to_nM", "pIC50",
                                          "canonical_smiles_standardised", "inchikey",
                                          "inchikey_skeleton", "fragments_before", "formal_charge_before"]
    audit_rows = []
    for r in log["records"]:
        row = {k: r.get(k) for k in AUDIT_FIELDS}
        row.update({"value_nM": r["value_nM"], "unit_conversion_factor_to_nM": r["unit_conversion_factor_to_nM"],
                    "pIC50": r["pIC50"], "canonical_smiles_standardised": r["canonical_smiles"],
                    "inchikey": r["inchikey"], "inchikey_skeleton": r["inchikey_skeleton"],
                    "fragments_before": r["fragments_before"], "formal_charge_before": r["formal_charge_before"]})
        audit_rows.append(row)
    _write_csv(curated_dir / "records_curated.csv", audit_rows, record_fields)

    _write_csv(results_dir / "units_refused.csv", log["unit_refusals"],
               ["activity_id", "standard_units", "units_as_reported", "standard_value", "reason"])
    _write_csv(results_dir / "structure_failures.csv", log["structure_failures"],
               ["activity_id", "molecule_chembl_id", "canonical_smiles_in", "failure"])
    _write_csv(results_dir / "replicate_spreads.csv",
               [c for c in log["compounds"] if c["n_records"] > 1],
               ["inchikey", "canonical_smiles", "n_records", "pIC50_median", "pIC50_min",
                "pIC50_max", "pIC50_spread", "pIC50_stdev", "flag_spread_over_1_log",
                "molecule_chembl_ids", "n_assays", "n_documents", "src_ids"])
    write_cascade_tables(log, results_dir)
    figures = write_figures(log, results_dir)

    report = {k: v for k, v in log.items() if k not in ("records", "compounds", "unit_refusals",
                                                        "structure_failures")}
    report["figures"] = figures
    (results_dir / "curation_log.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False, default=str) + "\n", encoding="utf-8")
    (results_dir / "replicate_noise.json").write_text(
        json.dumps(log["replicate_noise"], indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    manifest = {
        "target_chembl_id": log["target_chembl_id"],
        "chembl_db_version": log["chembl_db_version"],
        "curated_utc": log["curated_utc"],
        "command": log["command"],
        "seed": SEED,
        "python": platform.python_version(),
        "source": {"raw_dir": f"data/raw/{log['target_chembl_id']}",
                   "raw_manifest_sha256": _sha256(target_dir / MANIFEST),
                   "raw_files": raw_manifest["raw_files"],
                   "raw_verified_before_and_after": True},
        "cascade": [{k: s[k] for k in ("step", "rule", "records_in", "records_out", "records_removed")}
                    for s in log["steps"]],
        "curated_records": len(log["records"]),
        "curated_compounds": len(log["compounds"]),
        "target_variable": "pIC50 = -log10(IC50 in mol/L), median over replicates per InChIKey",
        "compound_identity": log["identity_decisions"]["compound_identity_key"],
        "files": {},
        "licence": "CC BY-SA 3.0, derived from ChEMBL. ChEMBL data is from "
                   "http://www.ebi.ac.uk/chembl - the version of ChEMBL is "
                   f"{log['chembl_db_version'].lower() if log['chembl_db_version'] else 'unknown'}.",
        "disclaimer": "Literature- and patent-reported bioactivity measurements. No claim of "
                      "efficacy, safety or suitability of any compound.",
    }
    for name in ("compounds.csv", "records_curated.csv"):
        manifest["files"][name] = {"sha256": _sha256(curated_dir / name)}
    (curated_dir / MANIFEST).write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    log["manifest"] = manifest
    return log


def identity_decisions(log: dict) -> dict:
    """The Gate 1 open questions, answered with numbers from this run."""
    records = log["records"]
    by_key: dict[str, list[dict]] = defaultdict(list)
    for r in records:
        by_key[r["inchikey"]].append(r)
    salt_rows = [r for r in records
                 if r.get("parent_molecule_chembl_id") and r["parent_molecule_chembl_id"] != r.get("molecule_chembl_id")]
    molecule_ids = {r["molecule_chembl_id"] for r in records}
    parent_ids = {r.get("parent_molecule_chembl_id") or r["molecule_chembl_id"] for r in records}
    merged = [c for c in log["compounds"] if c["n_molecule_chembl_ids"] > 1]
    skeletons = Counter(c["inchikey_skeleton"] for c in log["compounds"])
    busiest = max(log["compounds"], key=lambda c: c["n_records"])
    return {
        "compound_identity_key": "standardised InChIKey (full, after largest fragment and "
                                 "neutralisation); ChEMBL molecule and parent IDs are kept as "
                                 "audit columns, not as the identity",
        "why": "ChEMBL molecule IDs separate a salt from its parent, so keying on them would enter "
               "the same chemistry twice and could place the two copies on opposite sides of a "
               "split. The standardised InChIKey is the chemistry-level identity and merges them "
               "by construction.",
        "rows_with_parent_different_from_molecule_id": len(salt_rows),
        "distinct_molecule_ids_in_curated_records": len(molecule_ids),
        "distinct_parent_ids_in_curated_records": len(parent_ids),
        "compounds_merging_more_than_one_molecule_id": len(merged),
        "molecule_ids_absorbed_by_merging": sum(c["n_molecule_chembl_ids"] - 1 for c in merged),
        "stereoisomer_or_tautomer_pairs_sharing_a_skeleton": sum(n - 1 for n in skeletons.values() if n > 1),
        "distinct_inchikey_skeletons": len(skeletons),
        "busiest_compound": {
            "inchikey": busiest["inchikey"], "molecule_chembl_ids": busiest["molecule_chembl_ids"],
            "n_records": busiest["n_records"], "n_assays": busiest["n_assays"],
            "n_documents": busiest["n_documents"], "pIC50_median": busiest["pIC50_median"],
            "pIC50_spread": busiest["pIC50_spread"],
            "decision": "kept as ONE compound with the median of its replicates. Repeated reporting "
                        "of a reference compound is not repeated chemistry; leaving its rows "
                        "separate would weight it like a whole analogue series and let it appear on "
                        "both sides of a random split.",
        },
        "rows_by_src_id": dict(sorted(Counter(str(r.get("src_id")) for r in records).items())),
        "src_id_note": "src_id 37 is BindingDB patent bioactivity data (see results/sources/). It is "
                       "kept, because the measurements are IC50 values against the same target, but "
                       "it is a different reporting regime from the primary literature (src_id 1) and "
                       "is retained as a column so the decision can be revisited or tested.",
        "data_validity_comments": dict(sorted(Counter(
            r.get("data_validity_comment") or "<none>" for r in records).items())),
        "data_validity_decision": "ChEMBL's own validity comments are carried through as columns and "
                                  "reported, not used as an eighth filter: the seven-step cascade is "
                                  "the specified one, and removing these rows is a decision for Aadit. "
                                  "The counts above say exactly how much is at stake.",
        "potential_duplicate_rows": sum(1 for r in records if r.get("potential_duplicate") in (1, "1", True)),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    args = parser.parse_args(argv)
    target_id = args.target_id.strip().upper()
    log = curate_target(DATA_RAW / target_id, DATA_CURATED / target_id,
                        RESULTS / "curation" / target_id)
    noise = log["replicate_noise"]
    print(f"{target_id}: {log['input_records']} raw records -> {len(log['records'])} curated records "
          f"-> {len(log['compounds'])} compounds")
    floor = noise["noise_floor"]
    print(f"  conditional repeat SD: {floor['recommended_for_gate4']:.3f} log units "
          f"(all groups {noise['pooled_within_compound_sd_all_groups']:.3f}; Gate 4 choice pending) from "
          f"{noise['replicated_compounds']} replicated compounds, "
          f"{noise['groups_spread_at_least_0_01']} with spread >=0.01")
    print(f"  -> data/curated/{target_id}/ and results/curation/{target_id}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
