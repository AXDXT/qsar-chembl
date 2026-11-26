"""Offline tests for src.retrieve: the network is replaced by small fakes."""
import json

import pytest

from src import retrieve

FAKE_STATUS = {"chembl_db_version": "ChEMBL_TEST", "chembl_release_date": "2000-01-01", "status": "UP"}
FAKE_TARGET = {"target_chembl_id": "CHEMBL0", "pref_name": "Fake", "organism": "Homo sapiens",
               "target_type": "SINGLE PROTEIN"}
FAKE_RECORDS = [
    {"activity_id": 1, "standard_type": "IC50", "standard_relation": "=", "standard_value": "10.0",
     "standard_units": "nM", "assay_type": "B"},
    {"activity_id": 2, "standard_type": "Ki", "standard_relation": ">", "standard_value": None,
     "standard_units": None, "assay_type": "F"},
]


def _pull(tmp_path, records=FAKE_RECORDS, total=None, **kw):
    return retrieve.pull_target(
        "chembl0",
        out_root=tmp_path,
        status_fn=lambda: dict(FAKE_STATUS),
        target_fn=lambda tid: dict(FAKE_TARGET),
        activities_fn=lambda tid: ([dict(r) for r in records], len(records) if total is None else total),
        **kw,
    )


def test_manifest_written_and_complete(tmp_path):
    path = _pull(tmp_path)
    m = json.loads(path.read_text())
    for key in retrieve.REQUIRED_MANIFEST_KEYS:
        assert m.get(key) not in (None, "", {}), key
    assert m["target_chembl_id"] == "CHEMBL0"
    assert m["chembl_db_version"] == "ChEMBL_TEST"
    assert m["record_count"] == m["api_total_count"] == 2
    assert m["query"]["filters"] == {"target_chembl_id": "CHEMBL0"}
    assert m["retrieval_started_utc"].endswith("+00:00")
    assert m["client"]["name"] == "chembl_webresource_client" and m["client"]["version"]


def test_raw_records_round_trip_unmodified(tmp_path):
    path = _pull(tmp_path)
    assert retrieve.read_jsonl_gz(path.parent / retrieve.RAW_ACTIVITIES) == FAKE_RECORDS
    assert retrieve.verify_raw(path.parent)


def test_verify_raw_detects_mutation(tmp_path):
    path = _pull(tmp_path)
    (path.parent / retrieve.RAW_TARGET).write_text("{}")
    assert not retrieve.verify_raw(path.parent)


def test_refuses_to_overwrite_without_force(tmp_path):
    _pull(tmp_path)
    with pytest.raises(FileExistsError):
        _pull(tmp_path)
    _pull(tmp_path, force=True)


def test_incomplete_pull_is_rejected(tmp_path):
    with pytest.raises(RuntimeError):
        _pull(tmp_path, total=3)
    assert not (tmp_path / "CHEMBL0" / retrieve.MANIFEST).exists()


def test_identical_content_gives_identical_hash(tmp_path):
    a = json.loads(_pull(tmp_path / "a").read_text())
    b = json.loads(_pull(tmp_path / "b").read_text())
    assert a["raw_files"] == b["raw_files"]


@pytest.mark.live
def test_live_status_endpoint():
    status = retrieve.fetch_status()
    assert status["chembl_db_version"].startswith("ChEMBL_")
