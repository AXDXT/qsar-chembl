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
| `data/raw/` | **immutable**: unmodified API responses only, each with `manifest.json` (sha256 per file) |
| `data/curated/` | derived data (Gate 2 onwards), always regenerated from `data/raw/` by code |
| `results/` | logged outputs: counts, profiles, figures |
| `notebooks/` | exploration only; nothing a result depends on lives here |
| `tests/` | offline by default; `-m live` runs the tests that hit the ChEMBL API |

## Setup and run

```bash
.venv/bin/pip install -r requirements.txt   # Python 3.13; versions pinned from the working .venv
./run_all.sh                 # compile, offline tests, pull CHEMBL220 if not already pulled
./run_all.sh <TARGET_ID>     # same for another ChEMBL target id
.venv/bin/python -m pytest -q -m live       # tests that need the live API
```

Pull one target directly: `.venv/bin/python -m src.retrieve CHEMBL220`. An existing pull is never
overwritten unless `--force` is given. `src.retrieve.verify_raw(path)` re-checks every raw file
against the sha256 in its manifest.

## Reproducibility conventions

- One seed, `src.config.SEED`; stochastic code calls `set_global_seed()` or passes `random_state=SEED`.
- The ChEMBL client's on-disk HTTP cache is disabled for retrieval so a stale release cannot be served.
- The manifest records the ChEMBL release (from the API status endpoint), UTC retrieval times,
  the exact client call and equivalent URL, the record count checked against the API's total,
  and the client and Python versions.

## Data source, licence and citation

To be filled from the EBI/ChEMBL pages themselves (Gate 1 item 3) - not yet verified.
