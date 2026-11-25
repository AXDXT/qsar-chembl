"""Resolve the ChEMBL ``src_id`` values that appear in a frozen raw pull.

    .venv/bin/python -m src.sources CHEMBL284

Reads the src_id distribution from ``data/raw/<TARGET>/`` (offline), then fetches the
description of each one from the ChEMBL ``source`` endpoint and saves the unmodified API
records, with row counts, to ``results/sources/<TARGET>/``.

This answers the Gate 1 open question "what is src_id 37?". It is **metadata only**: the
curation path in ``src/curate.py`` never calls this module and never needs the network.
"""
from __future__ import annotations

import argparse
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path

import requests

from src.config import CHEMBL_API_BASE, DATA_RAW, RESULTS
from src.retrieve import MANIFEST, RAW_ACTIVITIES, read_jsonl_gz, verify_raw

SOURCE_URL = f"{CHEMBL_API_BASE}/source.json"


def fetch_source(src_id: int, timeout: float = 60.0) -> dict:
    """Return the API's unmodified source record for one src_id."""
    resp = requests.get(SOURCE_URL, params={"src_id": src_id}, timeout=timeout)
    resp.raise_for_status()
    sources = resp.json().get("sources", [])
    if len(sources) != 1:
        raise ValueError(f"src_id {src_id}: expected one source record, got {len(sources)}")
    return sources[0]


def src_id_counts(target_dir: Path) -> dict[str, int]:
    """Offline: rows per src_id in the frozen pull, with checksums verified."""
    target_dir = Path(target_dir)
    if not verify_raw(target_dir):
        raise ValueError("Raw file checksum mismatch; source lookup refused")
    records = read_jsonl_gz(target_dir / RAW_ACTIVITIES)
    counts = Counter(str(rec.get("src_id")) for rec in records)
    return dict(sorted(counts.items(), key=lambda pair: int(pair[0])))


def describe_sources(target_dir: Path, out_dir: Path, fetch=fetch_source) -> dict:
    target_dir, out_dir = Path(target_dir), Path(out_dir)
    counts = src_id_counts(target_dir)
    manifest = json.loads((target_dir / MANIFEST).read_text(encoding="utf-8"))
    log = {
        "target_chembl_id": manifest["target_chembl_id"],
        "chembl_db_version": manifest["chembl_db_version"],
        "checked_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "endpoint": SOURCE_URL,
        "command": f".venv/bin/python -m src.sources {manifest['target_chembl_id']}",
        "rows_per_src_id": counts,
        "sources": {src_id: fetch(int(src_id)) for src_id in counts},
        "note": "Row counts are over the whole raw pull, before any curation filter. "
                "Source descriptions are the API's own text, unmodified.",
    }
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "sources.json").write_text(
        json.dumps(log, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    lines = [f"# Data sources in the {log['target_chembl_id']} raw pull", "",
             f"Release {log['chembl_db_version']}; checked {log['checked_utc']} with `{log['command']}`.",
             "Row counts are over all raw rows, before curation.", "",
             "| src_id | Short name | Description | Raw rows |", "|---|---|---|---:|"]
    for src_id, n in counts.items():
        src = log["sources"][src_id]
        lines.append(f"| {src_id} | {src.get('src_short_name')} | {src.get('src_description')} | {n} |")
    lines += ["", "Full API text, including each source's comment, is in `sources.json`.", ""]
    (out_dir / "sources.md").write_text("\n".join(lines), encoding="utf-8")
    return log


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("target_id", nargs="?", default="CHEMBL284")
    args = parser.parse_args(argv)
    target_id = args.target_id.strip().upper()
    out_dir = RESULTS / "sources" / target_id
    log = describe_sources(DATA_RAW / target_id, out_dir)
    print(f"{target_id}: {len(log['sources'])} src_id values described -> {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
