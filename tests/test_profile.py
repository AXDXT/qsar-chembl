import copy
import csv
import json

import pytest

from src import profile, retrieve
from test_retrieve import _pull


def test_distributions_missing_values_and_repeated_reporting():
    records = [
        {"molecule_chembl_id": "CHEMBL1", "standard_units": "nM", "document_chembl_id": "D1", "value": 0},
        {"molecule_chembl_id": "CHEMBL1", "standard_units": None, "document_chembl_id": "D1", "value": " "},
        {"molecule_chembl_id": "CHEMBL2", "standard_units": "", "value": None},
        {"molecule_chembl_id": None},
    ]
    original = copy.deepcopy(records)
    summary, compounds = profile.summarize_records(records)
    assert records == original
    assert compounds == {"CHEMBL1": 2, "CHEMBL2": 1}
    assert summary["distributions"]["standard_units"] == {"<missing>": 3, "nM": 1}
    assert summary["missing_values_by_field"]["value"] == {"absent": 1, "null": 1, "blank": 1, "total_missing": 3}
    repeated = summary["compound_row_counts"]
    assert repeated["molecules_with_multiple_rows"] == 1
    assert repeated["rows_missing_molecule_id"] == 1
    assert repeated["rows_beyond_first_per_molecule"] == 1
    assert repeated["median"] == 1.5
    assert repeated["histogram_rows_per_molecule"] == {1: 1, 2: 1}
    assert summary["documents"] == {"distinct_ids": 1, "rows_missing_id": 2, "top_ten_by_rows": {"D1": 2}}


def test_profile_preserves_every_source_byte(tmp_path):
    manifest = _pull(tmp_path / "raw")
    before = {p.name: p.read_bytes() for p in manifest.parent.iterdir() if not p.name.startswith("._")}
    out_dir = tmp_path / "results"
    result = profile.profile_target(manifest.parent, out_dir)
    assert before == {p.name: p.read_bytes() for p in manifest.parent.iterdir() if not p.name.startswith("._")}
    assert retrieve.verify_raw(manifest.parent)
    assert result["record_count"] == 2
    assert result["source_manifest_sha256"] == retrieve._sha256(manifest)
    assert json.loads((out_dir / "profile.json").read_text())["record_count"] == 2
    assert "Raw records: **2**" in (out_dir / "profile.md").read_text()
    with (out_dir / "compound_row_counts.csv").open() as fh:
        assert list(csv.reader(fh)) == [["molecule_chembl_id", "raw_row_count"]]


def test_corrupt_source_is_rejected(tmp_path):
    manifest = _pull(tmp_path / "raw")
    (manifest.parent / retrieve.RAW_TARGET).write_text("{}")
    with pytest.raises(ValueError, match="checksum"):
        profile.profile_target(manifest.parent, tmp_path / "results")
    assert not (tmp_path / "results").exists()


@pytest.mark.parametrize("relative", [".", "nested", ".."])
def test_output_cannot_overlap_source(tmp_path, relative):
    manifest = _pull(tmp_path / "raw")
    with pytest.raises(ValueError, match="separate"):
        profile.profile_target(manifest.parent, manifest.parent / relative)


def test_empty_pull_has_no_invented_replicate_statistics():
    summary, compounds = profile.summarize_records([])
    assert not compounds
    assert summary["record_count"] == 0
    assert summary["compound_row_counts"]["median"] is None
    assert summary["compound_row_counts"]["maximum"] is None
