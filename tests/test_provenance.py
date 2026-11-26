import hashlib
import json

import pytest

from src import provenance

README = b"* Title: ChEMBL FTP Directory\n* Release:      chembl_37\n"


def fake_get(status_version):
    def get(url):
        if url.endswith("status.json"):
            return json.dumps({"chembl_db_version": status_version}).encode()
        if url.endswith("/README"):
            return README
        return b"body of " + url.encode()
    return get


def test_readme_release_parsed():
    assert provenance.readme_release(README) == "chembl_37"
    assert provenance.readme_release(b"no header") is None


def test_terms_saved_unmodified_with_manifest(tmp_path):
    manifest = provenance.fetch_terms(tmp_path, get=fake_get("ChEMBL_37"))
    assert manifest["release_matches"] is True
    for name in (*provenance.TERMS_FILES, "status.json"):
        entry = manifest["files"][name]
        body = (tmp_path / name).read_bytes()
        assert hashlib.sha256(body).hexdigest() == entry["sha256"]
        assert entry["url"].startswith("https://") and entry["fetched_utc"]
    assert json.loads((tmp_path / provenance.MANIFEST).read_text()) == manifest


def test_release_mismatch_flagged(tmp_path):
    assert provenance.fetch_terms(tmp_path, get=fake_get("ChEMBL_38"))["release_matches"] is False


@pytest.mark.live
def test_live_terms(tmp_path):
    manifest = provenance.fetch_terms(tmp_path)
    assert manifest["release_matches"]
    assert b"Attribution-ShareAlike" in (tmp_path / "LICENSE").read_bytes()
