"""Offline checks of the count-only Gate 1 cascade, not structure curation."""
from copy import deepcopy
import csv
import json
import sys
from types import SimpleNamespace

import pytest

from src.candidates import CANDIDATES, NM_FACTORS, count_records, fetch_count_pages, main


def record(**changes):
    return dict({"assay_type": "B", "standard_type": "IC50", "target_organism": "Homo sapiens",
                 "target_tax_id": "9606", "standard_units": "nM", "standard_relation": "=",
                 "standard_value": "12.5", "molecule_chembl_id": "CHEMBL1"}, **changes)


@pytest.mark.parametrize("field,value", [
    ("assay_type", "F"), ("standard_type", "Ki"), ("target_organism", "Mus musculus"),
    ("target_tax_id", "10090"), ("standard_units", "ug.mL-1"),
    ("standard_relation", ">"), ("standard_relation", "<"), ("standard_relation", ">="),
    ("standard_relation", "<="), ("standard_relation", "~"), ("standard_relation", None),
])
def test_excludes_rows_failing_each_filter(field, value):
    result = count_records([record(**{field: value})])
    assert result["filtered_records"] == result["distinct_compounds"] == 0


@pytest.mark.parametrize("unit", NM_FACTORS)
def test_accepts_molar_units_without_changing_values(unit):
    records = [record(standard_units=unit)]
    original = deepcopy(records)
    assert count_records(records)["distinct_compounds"] == 1
    assert records == original


@pytest.mark.parametrize("value", [None, "NaN", "inf", "0", "-1", "invalid"])
def test_numeric_diagnostic_does_not_add_a_headline_filter(value):
    result = count_records([record(standard_value=value)])
    assert result["filtered_records"] == result["distinct_compounds"] == 1
    assert result["cascade"]["positive_finite_value"] == 0


def test_sequential_counts_and_unique_molecule_ids():
    records = [record(), record(), record(molecule_chembl_id="CHEMBL2"),
               record(molecule_chembl_id=None), record(assay_type="F"), record(standard_type="Ki"),
               record(standard_relation=">"), record(standard_value=None)]
    original = deepcopy(records)
    result = count_records(records)
    assert list(result["cascade"].values()) == [8, 7, 6, 6, 6, 5, 4]
    assert result["filtered_records"] == 5
    assert result["distinct_compounds"] == 2
    assert result["survivors_missing_compound_id"] == 1
    assert records == original


def test_empty_input_has_zero_counts():
    assert count_records([])["distinct_compounds"] == 0
    assert set(count_records([])["cascade"].values()) == {0}


def test_paging_includes_final_partial_page_in_order(monkeypatch):
    monkeypatch.setattr("src.candidates.PAGE_SIZE", 2)
    records = [{"activity_id": n} for n in range(5)]
    assert fetch_count_pages(records, len(records), "TEST") == records
    assert fetch_count_pages([], 0, "TEST") == []


def fake_api(monkeypatch, records, total=None, statuses=None):
    """Install a minimal client double without importing the live client module."""
    from src import candidates

    class Query:
        def only(self, fields):
            return self

        def order_by(self, field):
            assert field == "activity_id"
            return self

        def __len__(self):
            return len(records) if total is None else total

        def __getitem__(self, section):
            return records[section]

    query = Query()
    client = SimpleNamespace(activity=SimpleNamespace(filter=lambda **kwargs: query))
    monkeypatch.setitem(sys.modules, "chembl_webresource_client.new_client", SimpleNamespace(new_client=client))
    monkeypatch.setattr(candidates, "_configure_client", lambda: None)
    monkeypatch.setattr(candidates, "fetch_target", lambda tid: {
        "target_type": "SINGLE PROTEIN", "organism": "Homo sapiens", "tax_id": "9606",
        "target_components": [{}], "pref_name": "Fixture protein",
    })
    statuses = iter(statuses or [{"chembl_db_version": "ChEMBL_TEST"}] * 2)
    monkeypatch.setattr(candidates, "fetch_status", lambda: next(statuses))


@pytest.mark.parametrize("records,total,error", [
    ([record(activity_id=1)], 2, "Incomplete count"),
    ([record(activity_id=1), record(activity_id=1)], 2, "Duplicate activity ids"),
])
def test_measurement_rejects_incomplete_or_duplicate_pages(monkeypatch, records, total, error):
    from src.candidates import measure_target

    fake_api(monkeypatch, records, total=total)
    with pytest.raises(RuntimeError, match=error):
        measure_target("CHEMBL0")


def test_measurement_rejects_release_change(monkeypatch):
    from src.candidates import measure_target

    fake_api(monkeypatch, [record(activity_id=1)], statuses=[
        {"chembl_db_version": "ChEMBL_TEST_A"}, {"chembl_db_version": "ChEMBL_TEST_B"}])
    with pytest.raises(RuntimeError, match="release changed"):
        measure_target("CHEMBL0")


def test_measurement_logs_exact_query_and_verified_counts(monkeypatch):
    from src.candidates import measure_target

    fake_api(monkeypatch, [record(activity_id=1)])
    result = measure_target("CHEMBL0")
    assert result["api_total_count"] == result["counts"]["cascade"]["all_records"] == 1
    assert result["query"]["filters"] == {"target_chembl_id": "CHEMBL0"}
    assert "order_by('activity_id')" in result["query"]["client_call"]
    assert result["counts"]["distinct_compounds"] == 1


def test_cli_logs_queries_and_writes_measured_summary(monkeypatch, tmp_path):
    from src import candidates

    seen = []

    def measure(target_id):
        seen.append(target_id)
        return {"target_chembl_id": target_id, "target": {"pref_name": "Fixture protein"},
                "finished_utc": "2000-01-01T00:00:00+00:00", "api_total_count": 1,
                "status_before": {"chembl_db_version": "ChEMBL_TEST"},
                "query": {"filters": {"target_chembl_id": target_id}},
                "counts": count_records([record()])}

    monkeypatch.setattr(candidates, "measure_target", measure)
    assert main(["--out-dir", str(tmp_path)]) == 0
    assert seen == list(CANDIDATES)
    logs = list(tmp_path.glob("CHEMBL*.json"))
    assert len(logs) == len(CANDIDATES)
    assert json.loads(logs[0].read_text())["query"]["filters"]["target_chembl_id"] in CANDIDATES
    with (tmp_path / "candidate_counts.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == len(CANDIDATES)
    assert all(row["approx_distinct_compounds"] == "1" for row in rows)
    assert "Count-only approximation" in (tmp_path / "candidate_counts.md").read_text()
