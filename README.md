# Project 19 - QSAR Bioactivity Modelling (ChEMBL)

Second-year B.Tech CS (AI and ML) coursework, dataset-only. Status: **Gate 1 in progress**
(repository, retrieval, candidate-target counts). No curation, featurisation, splitting or
modelling yet - those are later gates.

This project models literature-reported bioactivity measurements. Nothing in it claims or
implies efficacy, safety or suitability of any compound.

## Layout

| Path | Contents |
|---|---|
| `src/config.py` | the single global seed (`SEED = 42`) and all paths |
| `src/retrieve.py` | target-parameterised ChEMBL pull -> `data/raw/<TARGET>/` + manifest |
| `src/candidates.py` | count-only comparison of four candidate human protein targets |
| `data/raw/` | **immutable**: unmodified API responses only, each with `manifest.json` (sha256 per file) |
| `data/curated/` | derived data (Gate 2 onwards), always regenerated from `data/raw/` by code |
| `results/` | logged outputs: counts, profiles, figures |
| `notebooks/` | exploration only; nothing a result depends on lives here |
| `tests/` | offline by default; `-m live` runs the tests that hit the ChEMBL API |

## Setup and run

```bash
.venv/bin/pip install -r requirements.txt   # Python 3.13; versions pinned from the working .venv
./run_all.sh                 # compile, offline tests, count shortlist and pull CHEMBL220 if missing
./run_all.sh <TARGET_ID>     # same for another ChEMBL target id
.venv/bin/python -m pytest -q -m live       # tests that need the live API
.venv/bin/python -m src.candidates          # live measured counts for all four candidates
.venv/bin/python -m src.candidates CHEMBL220 CHEMBL203  # compare selected targets
```

Pull one target directly: `.venv/bin/python -m src.retrieve CHEMBL220`. An existing pull is never
overwritten unless `--force` is given. `src.retrieve.verify_raw(path)` re-checks every raw file
against the sha256 in its manifest.

## Gate 1 target shortlist

The shortlist contains the specification example CHEMBL220 (acetylcholinesterase),
CHEMBL203 (EGFR), CHEMBL4822 (beta-secretase 1), and CHEMBL284 (DPP4). The counting
command verifies that each API target is a single-component human protein, then
fetches all its activities with a field projection. Dated JSON logs, a CSV table,
and a short summary are saved under `results/candidate_counts/`.
Each query orders by activity ID and reads independent pages with at most four
requests in flight; both the row total and unique activity IDs are checked before
a count log is accepted. `run_all.sh` reuses the committed candidate-count snapshot;
run `src.candidates` explicitly to measure a new dated snapshot.

The sequential headline filters are binding assays (`B`), standard `IC50`, human
target organism (taxon 9606), molar concentration units convertible to nM, and an
exact standard relation (`=`). A positive-finite-value diagnostic is logged
separately and does not change these headline counts. Distinct compounds are
unique nonempty ChEMBL molecule IDs among the surviving records. These are
**count-only approximations of Gate 2**, before structure checks, salt handling,
replicate aggregation, assay-quality review or conflict resolution. Nothing is
written to `data/curated/`, and no structures are standardised.

The human-organism rule refers to the target, not the assay's expression host.
Mass concentration units are excluded because converting them would require
additional compound information. The precise accepted units and exclusions are
recorded in each run log. Uncommon scaled unit strings are conservatively excluded
pending interpretation in Gate 2. The final target remains subject to class uniqueness
and instructor approval.

## Reproducibility conventions

- One seed, `src.config.SEED`; stochastic code calls `set_global_seed()` or passes `random_state=SEED`.
- The ChEMBL client's on-disk HTTP cache is disabled for retrieval so a stale release cannot be served.
- The manifest records the ChEMBL release (from the API status endpoint), UTC retrieval times,
  the exact client call and equivalent URL, the record count checked against the API's total,
  and the client and Python versions.

## Data source, licence and citation

To be filled from the EBI/ChEMBL pages themselves (Gate 1 item 3) - not yet verified.
