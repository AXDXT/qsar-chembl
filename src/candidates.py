"""Measure the Gate 1 target shortlist without curating any structures.

    .venv/bin/python -m src.candidates
    .venv/bin/python -m src.candidates CHEMBL220 CHEMBL203

Each target gets a dated, auditable count log; CSV and Markdown summarize that run.
The activity query is unfiltered except for target id, with a field projection to
avoid downloading structures. No source data or curated dataset is written.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import csv
from datetime import datetime, timezone
import importlib.metadata
import json
import math
from pathlib import Path

from src.config import RESULTS
from src.retrieve import (
    PAGE_SIZE, STATUS_URL, _configure_client, fetch_status, fetch_target, query_description,
)

CANDIDATES = ("CHEMBL220", "CHEMBL203", "CHEMBL4822", "CHEMBL284")
FIELDS = (
    "activity_id", "assay_type", "standard_type", "standard_units",
    "standard_relation", "standard_value", "target_organism", "target_tax_id",
    "molecule_chembl_id", "assay_chembl_id", "document_chembl_id",
)
# Only unambiguous molar concentration units. Mass/activity units need information
# we deliberately do not derive in Gate 1. Values are never converted or saved.
NM_FACTORS = {"M": 1e9, "mM": 1e6, "uM": 1e3, "µM": 1e3, "μM": 1e3,
              "nM": 1.0, "pM": 1e-3, "fM": 1e-6}
STAGES = (
    "all_records", "binding_assay", "ic50", "human_target",
    "molar_units", "exact_relation", "positive_finite_value",
)
FILTER_POLICY = {
    "binding_assay": "assay_type == 'B'",
    "ic50": "standard_type == 'IC50'",
    "human_target": "target_organism == 'Homo sapiens' and str(target_tax_id) == '9606'",
    "molar_units": f"standard_units in {list(NM_FACTORS)} (molar concentrations convertible to nM)",
    "exact_relation": "standard_relation == '='; censored, approximate and missing relations excluded",
    "positive_finite_value": "Diagnostic only, after the headline count: standard_value parses as a finite number > 0",
    "distinct_compounds": "unique nonempty molecule_chembl_id among surviving rows; no structure processing",
    "organism_interpretation": "Target organism, not assay expression host; candidates must have one human protein component",
    "limitations": "Count-only approximation of Gate 2: no salts, stereochemistry, duplicates, conflicting replicates, assay quality or structure checks",
    "unit_limitations": "Conservative whitelist of plain molar units; uncommon scaled unit strings are excluded pending interpretation in Gate 2",
}


def fetch_count_pages(query, total: int, target_id: str) -> list[dict]:
    """Read independent client-query slices, with at most four requests in flight.

    The caller orders by activity_id, so concatenating offset pages is stable.
    Each slice has its own cloned client query/session; no iterator is shared.
    """
    slices = [query[start:min(start + PAGE_SIZE, total)] for start in range(0, total, PAGE_SIZE)]

    def read_page(page):
        return [dict(rec) for rec in page]

    records = []
    with ThreadPoolExecutor(max_workers=4) as executor:
        for page in executor.map(read_page, slices):
            records.extend(page)
            print(f"{target_id}: {len(records)}/{total}", flush=True)
    return records


def count_records(records: list[dict]) -> dict:
    """Count sequential filters on unmodified records, without making curated data."""
    stages = dict.fromkeys(STAGES, 0)
    compounds: set[str] = set()
    assays: set[str] = set()
    documents: set[str] = set()
    missing_compound = 0
    for rec in records:
        stages["all_records"] += 1
        checks = (
            rec.get("assay_type") == "B",
            rec.get("standard_type") == "IC50",
            rec.get("target_organism") == "Homo sapiens" and str(rec.get("target_tax_id")) == "9606",
            rec.get("standard_units") in NM_FACTORS,
            rec.get("standard_relation") == "=",
        )
        for stage, passes in zip(STAGES[1:-1], checks):
            if not passes:
                break
            stages[stage] += 1
        else:
            if rec.get("molecule_chembl_id"):
                compounds.add(rec["molecule_chembl_id"])
            else:
                missing_compound += 1
            if rec.get("assay_chembl_id"):
                assays.add(rec["assay_chembl_id"])
            if rec.get("document_chembl_id"):
                documents.add(rec["document_chembl_id"])
            try:
                value = float(rec.get("standard_value"))
            except (TypeError, ValueError, OverflowError):
                continue
            if math.isfinite(value) and value > 0:
                stages["positive_finite_value"] += 1
    return {
        "cascade": stages,
        "filtered_records": stages["exact_relation"],
        "distinct_compounds": len(compounds),
        "survivors_missing_compound_id": missing_compound,
        "distinct_surviving_assays": len(assays),
        "distinct_surviving_documents": len(documents),
        "raw_distributions": {
            field: dict(sorted(Counter(
                str(rec[field]) if rec.get(field) is not None else "<missing>"
                for rec in records
            ).items()))
            for field in ("assay_type", "standard_type", "standard_units", "standard_relation", "target_organism")
        },
    }


def measure_target(target_id: str) -> dict:
    started = datetime.now(timezone.utc).isoformat(timespec="seconds")
    before = fetch_status()
    target = fetch_target(target_id)
    if not (target.get("target_type") == "SINGLE PROTEIN"
            and target.get("organism") == "Homo sapiens"
            and str(target.get("tax_id")) == "9606"
            and len(target.get("target_components", [])) == 1):
        raise ValueError(f"{target_id} is not a single-component human protein target")
    _configure_client()
    from chembl_webresource_client.new_client import new_client

    query = new_client.activity.filter(target_chembl_id=target_id).only(list(FIELDS)).order_by("activity_id")
    total = len(query)
    records = fetch_count_pages(query, total, target_id)
    after = fetch_status()
    if before.get("chembl_db_version") != after.get("chembl_db_version"):
        raise RuntimeError("ChEMBL release changed during counting; rerun")
    if len(records) != total:
        raise RuntimeError(f"Incomplete count: {len(records)} rows fetched, {total} expected")
    if len({rec["activity_id"] for rec in records}) != total:
        raise RuntimeError("Duplicate activity ids across pages; rerun before trusting counts")
    description = query_description(target_id)
    description["client_call"] += f".only({list(FIELDS)!r}).order_by('activity_id')"
    description["fields_requested"] = list(FIELDS)
    description["equivalent_url"] += "&only=" + ",".join(FIELDS) + "&order_by=activity_id"
    description["max_concurrent_pages"] = 4
    description["note"] = "All target activities fetched with a field projection; filters applied locally for counts only."
    return {
        "target_chembl_id": target_id,
        "target": target,
        "started_utc": started,
        "finished_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "status_endpoint": STATUS_URL,
        "status_before": before,
        "status_after": after,
        "client_version": importlib.metadata.version("chembl_webresource_client"),
        "query": description,
        "api_total_count": total,
        "filter_policy": FILTER_POLICY,
        "counts": count_records(records),
        "command": f".venv/bin/python -m src.candidates {target_id}",
    }


def write_summary(logs: list[dict], out_dir: Path, paths: list[Path]) -> None:
    rows = []
    for log, path in zip(logs, paths):
        counts = log["counts"]
        rows.append({
            "target_chembl_id": log["target_chembl_id"],
            "target_name": log["target"]["pref_name"],
            "chembl_release": log["status_before"]["chembl_db_version"],
            "measured_utc": log["finished_utc"],
            "total_records": log["api_total_count"],
            "approx_filtered_records": counts["filtered_records"],
            "approx_distinct_compounds": counts["distinct_compounds"],
            "surviving_assays": counts["distinct_surviving_assays"],
            "surviving_documents": counts["distinct_surviving_documents"],
            "run_log": path.name,
        })
    with (out_dir / "candidate_counts.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = [
        "# Gate 1 candidate counts", "",
        "Measured from ChEMBL; each linked log records the release, UTC times, target metadata, query and full count cascade.", "",
        "**Count-only approximation of the Gate 2 cascade**, not curated compound counts. "
        "Filters: binding assays, standard IC50, human target organism, molar units convertible to nM, exact '=' relation. "
        "A positive-finite-value diagnostic is logged separately and does not change the headline counts. "
        "The unit whitelist covers plain molar units; uncommon scaled strings are excluded pending interpretation. "
        "Distinct compounds means nonempty ChEMBL molecule IDs; no structures or replicates are curated.", "",
        "| Target | Protein | Total records | Approx. filtered records | Approx. distinct compounds | Log |",
        "|---|---|---:|---:|---:|---|",
    ]
    for row in rows:
        lines.append(f"| {row['target_chembl_id']} | {row['target_name']} | {row['total_records']} | "
                     f"{row['approx_filtered_records']} | {row['approx_distinct_compounds']} | "
                     f"[{row['chembl_release']}, {row['measured_utc']}]({row['run_log']}) |")
    lines += ["", "## Observed cautions", ""]
    for log in logs:
        counts = log["counts"]
        stages = counts["cascade"]
        relations_lost = stages["molar_units"] - stages["exact_relation"]
        lines.append(
            f"- **{log['target_chembl_id']}**: {counts['distinct_surviving_assays']} assays across "
            f"{counts['distinct_surviving_documents']} documents survive; assay conditions and replicate agreement "
            f"still need review. {relations_lost} rows drop at the relation step (including censored, approximate "
            f"or missing relations). Raw target-organism counts: "
            f"{json.dumps(counts['raw_distributions']['target_organism'], ensure_ascii=False)}. "
            f"{counts['filtered_records'] - counts['distinct_compounds']} surviving rows exceed the unique molecule "
            f"count; these are not independent compounds."
        )
    lines += ["", "Target choice remains provisional until Aadit checks class uniqueness and obtains instructor approval. "
              "A raw profile and the final recommendation are still due; none of these counts establishes compound "
              "efficacy, safety or suitability.", ""]
    (out_dir / "candidate_counts.md").write_text("\n".join(lines), encoding="utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_ids", nargs="*", default=CANDIDATES)
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "candidate_counts")
    args = parser.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    logs, paths = [], []
    for target_id in args.target_ids:
        target_id = target_id.strip().upper()
        log = measure_target(target_id)
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        path = args.out_dir / f"{target_id}_{stamp}.json"
        path.write_text(json.dumps(log, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        logs.append(log)
        paths.append(path)
        print(f"{target_id}: {log['counts']['distinct_compounds']} distinct compounds -> {path}", flush=True)
    write_summary(logs, args.out_dir, paths)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
