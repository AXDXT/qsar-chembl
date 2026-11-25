"""Describe a frozen raw pull without filtering, editing or curating it.

    .venv/bin/python -m src.profile CHEMBL284

Outputs JSON/Markdown summaries and per-molecule row counts under results/profiles/.
Repeated rows are observations of repeated reporting, not independent experiments.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
from statistics import median

from src.candidates import FILTER_POLICY, count_records
from src.config import DATA_RAW, RESULTS
from src.retrieve import MANIFEST, RAW_ACTIVITIES, _sha256, read_jsonl_gz, verify_raw

DISTRIBUTION_FIELDS = (
    "assay_type", "standard_type", "standard_units", "standard_relation",
    "target_organism", "assay_organism", "src_id",
)


def is_missing(value) -> bool:
    """Count absent/null/blank strings as missing; retain zero and False."""
    return value is None or (isinstance(value, str) and not value.strip())


def summarize_records(records: list[dict]) -> tuple[dict, Counter]:
    """Pure count-only summary; records and nested API values are never modified."""
    distributions = {
        field: dict(sorted(Counter(
            "<missing>" if is_missing(rec.get(field)) else str(rec[field])
            for rec in records
        ).items()))
        for field in DISTRIBUTION_FIELDS
    }
    missing = {}
    for field in sorted({field for rec in records for field in rec} | set(DISTRIBUTION_FIELDS)):
        absent = sum(field not in rec for rec in records)
        null = sum(field in rec and rec[field] is None for rec in records)
        blank = sum(isinstance(rec.get(field), str) and not rec[field].strip() for rec in records)
        missing[field] = {"absent": absent, "null": null, "blank": blank,
                          "total_missing": absent + null + blank}
    compounds = Counter(rec["molecule_chembl_id"] for rec in records
                        if not is_missing(rec.get("molecule_chembl_id")))
    documents = Counter(rec["document_chembl_id"] for rec in records
                        if not is_missing(rec.get("document_chembl_id")))
    counts = list(compounds.values())
    return {
        "record_count": len(records),
        "distributions": distributions,
        "missing_values_by_field": missing,
        "missing_definition": "Absent key, null, or blank string; empty containers are retained as API values.",
        "compound_row_counts": {
            "distinct_molecule_ids": len(compounds),
            "rows_missing_molecule_id": sum(is_missing(rec.get("molecule_chembl_id")) for rec in records),
            "molecules_with_multiple_rows": sum(n > 1 for n in counts),
            "rows_beyond_first_per_molecule": sum(n - 1 for n in counts),
            "minimum": min(counts) if counts else None,
            "median": median(counts) if counts else None,
            "maximum": max(counts) if counts else None,
            "histogram_rows_per_molecule": dict(sorted(Counter(counts).items())),
            "interpretation": "Raw rows per molecule ID across all endpoints and assays; not experimentally verified replicates.",
        },
        "documents": {
            "distinct_ids": len(documents),
            "rows_missing_id": sum(is_missing(rec.get("document_chembl_id")) for rec in records),
            "top_ten_by_rows": dict(sorted(documents.items(), key=lambda pair: (-pair[1], pair[0]))[:10]),
        },
        "count_only_gate2_approximation": count_records(records),
        "filter_policy": FILTER_POLICY,
    }, compounds


def write_markdown(profile: dict, path: Path) -> None:
    lines = [f"# Raw profile: {profile['target_chembl_id']}", "",
             f"Release: {profile['chembl_db_version']}. Raw retrieval finished: {profile['retrieval_finished_utc']}.",
             f"Command: `{profile['command']}`. Source hashes verified before and after profiling.", "",
             f"Raw records: **{profile['record_count']}**. No rows are removed or changed.", ""]
    for field, distribution in profile["distributions"].items():
        lines += [f"## {field}", "", "| API value | Rows |", "|---|---:|"]
        lines += [f"| `{value.replace('|', '&#124;')}` | {n} |" for value, n in distribution.items()]
        lines += [""]
    lines += ["## Repeated reporting per compound", "", profile["compound_row_counts"]["interpretation"], "",
              "Exact per-molecule counts are in `compound_row_counts.csv`.", "", "```json",
              json.dumps(profile["compound_row_counts"], indent=2), "```", "",
              "## Missing values", "", profile["missing_definition"], "",
              "| Field | Absent | Null | Blank | Total missing |", "|---|---:|---:|---:|---:|"]
    for field, counts in profile["missing_values_by_field"].items():
        lines.append(f"| {field} | {counts['absent']} | {counts['null']} | {counts['blank']} | {counts['total_missing']} |")
    lines += ["", "## Count-only Gate 2 approximation", "", "```json",
              json.dumps(profile["count_only_gate2_approximation"], indent=2), "```", "",
              "These counts precede structure checks, salt handling, replicate aggregation and assay-quality review. "
              "Literature-reported bioactivity is not a random sample of chemistry. "
              "No efficacy, safety or suitability claim follows from these observations.", ""]
    path.write_text("\n".join(lines), encoding="utf-8")


def profile_target(target_dir: Path, out_dir: Path) -> dict:
    target_dir, out_dir = Path(target_dir).resolve(), Path(out_dir).resolve()
    # Reject any overlap, including output to a parent of the frozen pull.
    if (out_dir == target_dir or out_dir in target_dir.parents
            or target_dir in out_dir.parents or out_dir == DATA_RAW.resolve()
            or DATA_RAW.resolve() in out_dir.parents):
        raise ValueError("Profile output must be separate from data/raw and the source directory")
    if not verify_raw(target_dir):
        raise ValueError("Raw file checksum mismatch; profile refused")
    manifest = json.loads((target_dir / MANIFEST).read_text(encoding="utf-8"))
    records = read_jsonl_gz(target_dir / RAW_ACTIVITIES)
    if len(records) != manifest["record_count"]:
        raise ValueError("Raw row count disagrees with manifest")
    summary, compounds = summarize_records(records)
    profile = {
        "target_chembl_id": manifest["target_chembl_id"],
        "chembl_db_version": manifest["chembl_db_version"],
        "retrieval_finished_utc": manifest["retrieval_finished_utc"],
        "profiled_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "source_manifest_sha256": _sha256(target_dir / MANIFEST),
        "source_raw_files": manifest["raw_files"],
        "command": f".venv/bin/python -m src.profile {manifest['target_chembl_id']}",
        **summary,
    }
    if not verify_raw(target_dir):
        raise ValueError("Raw files changed during profiling")
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "profile.json").write_text(json.dumps(profile, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    with (out_dir / "compound_row_counts.csv").open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow(["molecule_chembl_id", "raw_row_count"])
        writer.writerows(sorted(compounds.items()))
    write_markdown(profile, out_dir / "profile.md")
    return profile


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id")
    args = parser.parse_args(argv)
    target_id = args.target_id.strip().upper()
    out_dir = RESULTS / "profiles" / target_id
    profile = profile_target(DATA_RAW / target_id, out_dir)
    print(f"{target_id}: profiled {profile['record_count']} raw rows -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
