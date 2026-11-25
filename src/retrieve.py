"""Pull every ChEMBL bioactivity record for one target and freeze it in data/raw/.

Usage (one command per target):

    .venv/bin/python -m src.retrieve CHEMBL220
    .venv/bin/python -m src.retrieve CHEMBL220 --force   # overwrite an existing pull

Writes, under data/raw/<TARGET_ID>/:
    activities.jsonl.gz  one API activity record per line, exactly as returned
    target.json          the API target record (organism, type, components)
    manifest.json        release, retrieval time (UTC), query, counts, client version, sha256

Nothing in data/raw/ is ever edited after it is written. Downstream code reads it
and writes elsewhere (data/curated/, results/).
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import importlib.metadata
import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable

import requests

from src.config import CHEMBL_API_BASE, DATA_RAW

STATUS_URL = f"{CHEMBL_API_BASE}/status.json"
PAGE_SIZE = 1000  # API maximum; the client default of 20 means ~50x more requests
RAW_ACTIVITIES = "activities.jsonl.gz"
RAW_TARGET = "target.json"
MANIFEST = "manifest.json"

REQUIRED_MANIFEST_KEYS = (
    "target_chembl_id",
    "chembl_db_version",
    "chembl_release_date",
    "retrieval_started_utc",
    "retrieval_finished_utc",
    "query",
    "record_count",
    "api_total_count",
    "client",
    "raw_files",
)


def _configure_client() -> None:
    """Disable the client's on-disk cache (it could serve a previous release) and page at the max."""
    from chembl_webresource_client.settings import Settings

    settings = Settings.Instance()
    settings.CACHING = False
    settings.MAX_LIMIT = PAGE_SIZE
    settings.TIMEOUT = 60.0


def fetch_status() -> dict:
    resp = requests.get(STATUS_URL, timeout=60)
    resp.raise_for_status()
    return resp.json()


def fetch_target(target_id: str) -> dict:
    _configure_client()
    from chembl_webresource_client.new_client import new_client

    return dict(new_client.target.get(target_id))


def fetch_activities(target_id: str) -> tuple[list[dict], int]:
    """Return (records, total_count reported by the API for the same query)."""
    _configure_client()
    from chembl_webresource_client.new_client import new_client

    qs = new_client.activity.filter(target_chembl_id=target_id)
    total = len(qs)
    records = [dict(r) for r in qs]
    return records, total


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_jsonl_gz(records: Iterable[dict], path: Path) -> int:
    """Write records one per line; mtime=0 so identical content gives an identical file hash."""
    n = 0
    with open(path, "wb") as raw, gzip.GzipFile(fileobj=raw, mode="wb", mtime=0) as gz:
        for rec in records:
            gz.write((json.dumps(rec, ensure_ascii=False) + "\n").encode("utf-8"))
            n += 1
    return n


def read_jsonl_gz(path: Path) -> list[dict]:
    with gzip.open(path, "rt", encoding="utf-8") as fh:
        return [json.loads(line) for line in fh]


def query_description(target_id: str) -> dict:
    return {
        "resource": "activity",
        "client_call": f"new_client.activity.filter(target_chembl_id='{target_id}')",
        "filters": {"target_chembl_id": target_id},
        "equivalent_url": f"{CHEMBL_API_BASE}/activity.json?target_chembl_id={target_id}",
        "page_size": PAGE_SIZE,
        "client_cache": "disabled",
        "note": "No assay-type, activity-type, organism, unit or relation filter is applied at "
                "retrieval; all filtering happens downstream on a copy.",
    }


def pull_target(
    target_id: str,
    out_root: Path = DATA_RAW,
    force: bool = False,
    status_fn: Callable[[], dict] = fetch_status,
    target_fn: Callable[[str], dict] = fetch_target,
    activities_fn: Callable[[str], tuple[list[dict], int]] = fetch_activities,
) -> Path:
    """Pull and freeze one target. Returns the path of the manifest."""
    target_id = target_id.strip().upper()
    out_dir = Path(out_root) / target_id
    manifest_path = out_dir / MANIFEST
    if manifest_path.exists() and not force:
        raise FileExistsError(
            f"{manifest_path} exists; data/raw is immutable. Re-run with --force to replace the pull."
        )
    out_dir.mkdir(parents=True, exist_ok=True)

    started = datetime.now(timezone.utc)
    status_before = status_fn()
    target = target_fn(target_id)
    records, api_total = activities_fn(target_id)
    status_after = status_fn()
    finished = datetime.now(timezone.utc)

    if status_before.get("chembl_db_version") != status_after.get("chembl_db_version"):
        raise RuntimeError("ChEMBL release changed during the pull; re-run.")
    if len(records) != api_total:
        raise RuntimeError(f"Pulled {len(records)} records but API reported {api_total}; re-run.")

    act_path = out_dir / RAW_ACTIVITIES
    tgt_path = out_dir / RAW_TARGET
    n = write_jsonl_gz(records, act_path)
    tgt_path.write_text(json.dumps(target, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    manifest = {
        "target_chembl_id": target_id,
        "target_pref_name": target.get("pref_name"),
        "target_organism": target.get("organism"),
        "target_type": target.get("target_type"),
        "chembl_db_version": status_before.get("chembl_db_version"),
        "chembl_release_date": status_before.get("chembl_release_date"),
        "status_endpoint": STATUS_URL,
        "status_response": status_before,
        "retrieval_started_utc": started.isoformat(timespec="seconds"),
        "retrieval_finished_utc": finished.isoformat(timespec="seconds"),
        "query": query_description(target_id),
        "record_count": n,
        "api_total_count": api_total,
        "client": {
            "name": "chembl_webresource_client",
            "version": importlib.metadata.version("chembl_webresource_client"),
            "python": platform.python_version(),
        },
        "raw_files": {
            RAW_ACTIVITIES: {"sha256": _sha256(act_path), "format": "gzip JSON lines, one API record per line"},
            RAW_TARGET: {"sha256": _sha256(tgt_path), "format": "JSON, API target record"},
        },
        "command": f"python -m src.retrieve {target_id}" + (" --force" if force else ""),
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return manifest_path


def verify_raw(target_dir: Path) -> bool:
    """True iff every raw file still matches the sha256 recorded in its manifest."""
    manifest = json.loads((Path(target_dir) / MANIFEST).read_text(encoding="utf-8"))
    return all(
        _sha256(Path(target_dir) / name) == meta["sha256"] for name, meta in manifest["raw_files"].items()
    )


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("target_id", help="ChEMBL target id, e.g. CHEMBL220")
    ap.add_argument("--force", action="store_true", help="replace an existing pull for this target")
    args = ap.parse_args(argv)
    path = pull_target(args.target_id, force=args.force)
    m = json.loads(path.read_text(encoding="utf-8"))
    print(f"{m['target_chembl_id']} ({m['target_pref_name']}): {m['record_count']} records, "
          f"{m['chembl_db_version']} -> {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
