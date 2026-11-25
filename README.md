# Project 19 - QSAR Bioactivity Modelling (ChEMBL)

Second-year B.Tech CS (AI and ML) coursework, dataset-only. Status: **Gate 2 awaiting human
sign-off; Gate 3 fixed baselines recorded, review pending**. Human DPP4
(CHEMBL284) is locked by Aadit. Offline curation, the structural split audit and the fixed
leakage probe are implemented, with an executable visual notebook. Aadit's noise-reference
decision and final approval remain; see [the Gate 2 handoff](GATE2_SIGNOFF.md).
Gate 3 provides a deterministic, checksummed descriptor and Morgan fingerprint cache
built only from the standardised structures, plus training-mean and descriptor-only forest
baselines under both saved splits. The full train/test table and uncertainty are in
[DATA_CARD.md](DATA_CARD.md#gate-3-fixed-baselines) and
[`results/baselines/CHEMBL284/baselines.md`](results/baselines/CHEMBL284/baselines.md).
No Gate 4 model, hyperparameter search or combined feature model has been fitted.

This project models literature-reported bioactivity measurements. Nothing in it claims or
implies efficacy, safety or suitability of any compound.

## Layout

| Path | Contents |
|---|---|
| `src/config.py` | the single global seed (`SEED = 42`) and all paths |
| `src/retrieve.py` | target-parameterised ChEMBL pull -> `data/raw/<TARGET>/` + manifest |
| `src/candidates.py` | count-only comparison of four candidate human protein targets |
| `src/provenance.py` | saves ChEMBL's README, LICENSE, REQUIRED.ATTRIBUTION and API status to `results/licence/` |
| `src/profile.py` | read-only raw profile, missing values and row counts per molecule ID |
| `src/curate.py` | offline seven-step cascade, median compound pIC50, distributions and replicate dispersion |
| `src/splits.py` | reproducible compound random/scaffold splits and exhaustive structural leakage audit; no model fitting |
| `src/leakage_probe.py` | one fixed Tanimoto nearest-neighbour leakage probe, scaffold-group bootstrap uncertainty; not a Gate 3 baseline |
| `src/features.py` | Gate 3 featurisation: 12 declared RDKit descriptors and radius-2 Morgan bits from the standardised structures, cached with a manifest; no model fitted |
| `src/baselines.py` | fixed training-mean and descriptor-only Random Forest baselines; exact saved splits, train/test metrics and scaffold-group bootstrap test intervals |
| `DATA_CARD.md` | provenance, licence, measured composition and known biases |
| `data/raw/` | **immutable**: unmodified API responses only, each with `manifest.json` (sha256 per file) |
| `data/curated/` | derived data (Gate 2 onwards), always regenerated from `data/raw/` by code |
| `data/features/` | Gate 3 feature cache, regenerated from `data/curated/`; byte-identical on a rerun |
| `results/` | logged outputs: counts, profiles, figures |
| `notebooks/gate2_dpp4_curation_and_leakage.ipynb` | executed offline visual walkthrough; reruns the production modules from frozen raw data in temporary output directories |
| `scripts/run_gate2_notebook.py` | executes all notebook cells without a server/kernel listener and saves verified outputs |
| `scripts/check_data_card_numbers.py` | checks documented curation figures and compares the `GATE2_SIGNOFF.md`, Gate 3 feature and baseline tables directly with logged results |
| `tests/` | offline by default; `-m live` runs the tests that hit the ChEMBL API |

## Setup and run

```bash
.venv/bin/pip install -r requirements.txt   # Python 3.13; versions pinned from the working .venv
./run_all.sh                 # compile, offline tests, frozen DPP4 -> curation/splits/probe/features/baselines
./run_all.sh <TARGET_ID>     # same for another ChEMBL target id
.venv/bin/python -m pytest -q -m live       # tests that need the live API
.venv/bin/python -m src.candidates          # live measured counts for all four candidates
.venv/bin/python -m src.candidates CHEMBL220 CHEMBL203  # compare selected targets
.venv/bin/python -m src.profile CHEMBL284    # offline profile of a frozen pull
.venv/bin/python -m src.curate CHEMBL284     # offline cascade from checksum-verified raw data
.venv/bin/python -m src.splits CHEMBL284     # offline splits from checksum-verified curated compounds
.venv/bin/python -m src.leakage_probe CHEMBL284  # fixed probe on both saved splits, no tuning
.venv/bin/python -m src.features CHEMBL284   # Gate 3 feature cache from the standardised structures
.venv/bin/python -m src.baselines CHEMBL284  # fixed mean + descriptor RF on both exact saved splits
.venv/bin/python scripts/run_gate2_notebook.py  # execute/save the visual walkthrough offline
.venv/bin/python scripts/check_data_card_numbers.py  # curation figures and sign-off tables vs logs
```

The Gate 2 notebook can also be opened in Jupyter with the repo's `.venv` Python
kernel and run from top to bottom. Its final cell checks raw immutability, byte-identical
curated CSVs, split/similarity CSVs and predictions, and matching statistical summaries
against the committed results. The headless script above runs the same cells and fails
on any error; `--output /tmp/gate2.ipynb` saves a separate executed copy. No live API,
new dependency or additional model is used. Figures use the existing Agg helper to avoid
exFAT AppleDouble style-file errors.

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

Both build-room agents recommend **CHEMBL284 (DPP4)** provisionally: its measured
4,004 distinct molecule IDs after the headline count filters give room for later
curation, with a smaller raw pull than the other candidates. The recommendation
is based on dataset size and scope, not on any model score. The saved full pull
has 8,397 records; its raw profile independently reproduces the shortlist's
4,859 filtered rows and 4,004 molecule IDs. See `DATA_CARD.md` and
`results/profiles/CHEMBL284/` for composition and cautions.

The profiler verifies raw checksums before and after reading, writes only under
`results/`, and emits `profile.json`, `profile.md` and `compound_row_counts.csv`.
Counts of repeated rows do not establish experimentally independent replicates.
`run_all.sh <TARGET_ID>` retrieves and profiles another target without changing
the saved DPP4 pull or selecting a final coursework target.

## Reproducibility conventions

- One seed, `src.config.SEED`; stochastic code calls `set_global_seed()` or passes `random_state=SEED`.
- The ChEMBL client's on-disk HTTP cache is disabled for retrieval so a stale release cannot be served.
- The manifest records the ChEMBL release (from the API status endpoint), UTC retrieval times,
  the exact client call and equivalent URL, the record count checked against the API's total,
  and the client and Python versions.

## Data source, licence and citation

Checked 2026-11-25 (UTC) with `.venv/bin/python -m src.provenance`, which saves the files below
unmodified, with their sha256, under `results/licence/`. Everything in this section is quoted
or paraphrased from those files only.

- **Source:** ChEMBL, EMBL-EBI, <https://www.ebi.ac.uk/chembl/>, accessed through the REST API
  (`https://www.ebi.ac.uk/chembl/api/data`) via `chembl_webresource_client`.
- **Version:** `ChEMBL_37`, release date 2026-05-01, from the API status endpoint
  (`results/licence/status.json`). The FTP README header says "Release: chembl_37,
  Date: 01/05/2026". These two agree, and the manifest records that check (`release_matches: true`).
- **Licence:** Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0). Source:
  <https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/README> and the full licence text at
  <https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/LICENSE>.
- **Required attribution** (FTP README): the attribution must contain the URL of the resource and
  the release version, for example: *"ChEMBL data is from http://www.ebi.ac.uk/chembl - the version
  of ChEMBL is chembl_37."* If ChEMBL is built into other works, EBI asks that ChEMBL IDs be kept and
  the release number be shown clearly
  (<https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/REQUIRED.ATTRIBUTION>).
- **Citation EBI asks for** (REQUIRED.ATTRIBUTION, "primary current citation"):
  Mendez D, Gaulton A, Bento AP, Chambers J, De Veij M, Félix E, Magariños MP, Mosquera JF, Mutowo P,
  Nowotka M, Gordillo-Marañón M, Hunter F, Junco L, Mugumbate G, Rodriguez-Lopez M, Atkinson F,
  Bosc N, Radoux CJ, Segura-Cabrera A, Hersey A, Leach AR. ChEMBL: towards direct deposition of
  bioassay data. *Nucleic Acids Res.* 2019 47(D1):D930-D940. DOI: 10.1093/nar/gky1075

**ChEMBL data is from http://www.ebi.ac.uk/chembl - the version of ChEMBL is chembl_37.**

What this means for this repo (our reading, not legal advice): any ChEMBL-derived data we share
(`data/raw/`, `data/curated/`) is an adaptation under the ShareAlike term, so it must stay under
CC BY-SA 3.0 with the attribution above. ChEMBL IDs are kept throughout. The licence for this
repository's own *code* is still undecided (see "Needs decision" in the room file).

**Not verified:** the general EMBL-EBI terms of use for the web services, and whether the ChEMBL
web portal suggests a different or newer citation. Only the FTP files above were checked.
