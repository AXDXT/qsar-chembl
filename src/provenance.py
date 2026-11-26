"""Fetch ChEMBL's own licence, attribution and release files, unmodified, with a manifest.

    .venv/bin/python -m src.provenance

Writes, under results/licence/:
    README, LICENSE, REQUIRED.ATTRIBUTION   byte-for-byte copies from the ChEMBL FTP "latest" dir
    status.json                             the API status record fetched in the same run
    manifest.json                           source URL, UTC fetch time and sha256 per file, plus
                                            whether the FTP README names the API's current release

README.md and DATA_CARD.md quote licence and citation text from these files only.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable

import requests

from src.config import RESULTS
from src.retrieve import STATUS_URL

FTP_LATEST = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest"
TERMS_FILES = ("README", "LICENSE", "REQUIRED.ATTRIBUTION")
MANIFEST = "manifest.json"


def _get(url: str) -> bytes:
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    return resp.content


def readme_release(readme: bytes) -> str | None:
    """Release string from the FTP README header, e.g. 'chembl_37'."""
    match = re.search(rb"Release:\s*(chembl_\d+)", readme, re.IGNORECASE)
    return match.group(1).decode() if match else None


def fetch_terms(out_dir: Path, get: Callable[[str], bytes] = _get) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    files = {}
    for name in TERMS_FILES:
        url = f"{FTP_LATEST}/{name}"
        body = get(url)
        (out_dir / name).write_bytes(body)
        files[name] = {
            "url": url,
            "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "sha256": hashlib.sha256(body).hexdigest(),
            "bytes": len(body),
        }
    status_body = get(STATUS_URL)
    (out_dir / "status.json").write_bytes(status_body)
    files["status.json"] = {
        "url": STATUS_URL,
        "fetched_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "sha256": hashlib.sha256(status_body).hexdigest(),
        "bytes": len(status_body),
    }
    api_version = json.loads(status_body).get("chembl_db_version")
    ftp_release = readme_release((out_dir / "README").read_bytes())
    manifest = {
        "command": ".venv/bin/python -m src.provenance",
        "files": files,
        "api_chembl_db_version": api_version,
        "ftp_readme_release": ftp_release,
        "release_matches": bool(api_version and ftp_release
                                and api_version.lower() == ftp_release.lower()),
    }
    (out_dir / MANIFEST).write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return manifest


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=RESULTS / "licence")
    args = parser.parse_args(argv)
    manifest = fetch_terms(args.out_dir)
    print(json.dumps({k: manifest[k] for k in
                      ("api_chembl_db_version", "ftp_readme_release", "release_matches")}))
    if not manifest["release_matches"]:
        print("WARNING: FTP 'latest' release differs from the API release; check before quoting.")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
