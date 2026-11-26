"""Offline tests for src.sources: the API call is replaced by a fake."""
import json

import pytest

from src import sources, retrieve
from test_retrieve import _pull

FAKE_SOURCES = {
    1: {"src_id": 1, "src_short_name": "LITERATURE", "src_description": "Scientific Literature",
        "src_comment": "medchem journals", "src_url": None},
    37: {"src_id": 37, "src_short_name": "BINDINGDB",
         "src_description": "BindingDB Patent Bioactivity Data",
         "src_comment": "granted US patents", "src_url": "https://www.bindingdb.org/"},
}


def _records():
    return [{"activity_id": 1, "src_id": 1, "molecule_chembl_id": "CHEMBL1"},
            {"activity_id": 2, "src_id": 37, "molecule_chembl_id": "CHEMBL2"},
            {"activity_id": 3, "src_id": 37, "molecule_chembl_id": "CHEMBL3"}]


def test_counts_come_from_the_frozen_pull_and_are_ordered(tmp_path):
    manifest = _pull(tmp_path / "raw", records=_records())
    assert sources.src_id_counts(manifest.parent) == {"1": 1, "37": 2}


def test_descriptions_are_saved_with_row_counts_and_raw_is_untouched(tmp_path):
    manifest = _pull(tmp_path / "raw", records=_records())
    before = {p.name: p.read_bytes() for p in manifest.parent.iterdir() if not p.name.startswith("._")}

    log = sources.describe_sources(manifest.parent, tmp_path / "out",
                                   fetch=lambda src_id: dict(FAKE_SOURCES[src_id]))

    assert before == {p.name: p.read_bytes() for p in manifest.parent.iterdir()
                      if not p.name.startswith("._")}
    assert retrieve.verify_raw(manifest.parent)
    assert log["rows_per_src_id"] == {"1": 1, "37": 2}
    saved = json.loads((tmp_path / "out" / "sources.json").read_text())
    assert saved["sources"]["37"]["src_short_name"] == "BINDINGDB"
    table = (tmp_path / "out" / "sources.md").read_text()
    assert "BindingDB Patent Bioactivity Data" in table
    assert "| 37 | BINDINGDB" in table


def test_corrupt_raw_pull_is_rejected(tmp_path):
    manifest = _pull(tmp_path / "raw", records=_records())
    (manifest.parent / retrieve.RAW_TARGET).write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        sources.src_id_counts(manifest.parent)


def test_ambiguous_api_response_is_refused(monkeypatch):
    class FakeResponse:
        def raise_for_status(self): pass
        def json(self): return {"sources": []}

    monkeypatch.setattr(sources.requests, "get", lambda *a, **k: FakeResponse())
    with pytest.raises(ValueError, match="one source record"):
        sources.fetch_source(99)


@pytest.mark.live
def test_live_src_id_37_is_bindingdb():
    """The Gate 1 open question, checked against the API itself."""
    source = sources.fetch_source(37)
    assert source["src_short_name"] == "BINDINGDB"
