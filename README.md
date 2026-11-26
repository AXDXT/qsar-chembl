# Project 19 - QSAR Bioactivity Modelling (ChEMBL)

Second-year B.Tech CS (AI and ML) coursework, dataset-only. Status: **Gate 2 awaiting human
sign-off; Gate 3 fixed baselines recorded, review pending**. Human DPP4
(CHEMBL284) is locked by Aadit. Offline curation, the structural split audit and the fixed
leakage probe are implemented, with an executable visual notebook. Aadit's noise-reference
decision and final approval remain; see [the Gate 2 handoff](GATE2_SIGNOFF.md).
Gate 3 provides a deterministic, checksummed descriptor and Morgan fingerprint cache
built only from the standardised structures, plus training-mean and descriptor-only forest
baselines under both saved splits. The full train/test table and uncertainty are in
[the Gate 3 handoff](GATE3_SIGNOFF.md),
[DATA_CARD.md](DATA_CARD.md#gate-3-fixed-baselines) and
[`results/baselines/CHEMBL284/baselines.md`](results/baselines/CHEMBL284/baselines.md).
Gate 4 now has a resumable scaffold-grouped nested-CV search and recorded
training-only searches for its tuned descriptor control under both saved splits.
The selected descriptor controls have been refitted and saved for both splits.
Both combined feature searches have partial production fit caches and remain in
progress. The held-out evaluation is pending.

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
| `src/nested_cv.py` | Gate 4 training-only grouped nested forest search, with incremental fit caches and outer-fold uncertainty |
| `src/gate4_training.py` | verifies saved inputs and runs the tuned descriptor control on training rows only; no held-out evaluation |
| `src/gate4_fingerprints.py` | records training-only, fold-local fingerprint bit-length decisions from collision loss; no model fitting |
| `src/gate4_combined.py` | training-only combined search consuming each fit's own fingerprint decision, with incremental fit caches |
| `src/gate4_refit.py` | refits the completed training search selection and verifies the saved forest on the same training rows |
| `src/gate4_model_io.py` | read-only loading of frozen forests, checking their completed search, current inputs, library versions and model checksum |
| `DATA_CARD.md` | provenance, licence, measured composition and known biases |
| `data/raw/` | **immutable**: unmodified API responses only, each with `manifest.json` (sha256 per file) |
| `data/curated/` | derived data (Gate 2 onwards), always regenerated from `data/raw/` by code |
| `data/features/` | Gate 3 feature cache, regenerated from `data/curated/`; byte-identical on a rerun |
| `results/` | logged outputs: counts, profiles, figures |
| `notebooks/gate2_dpp4_curation_and_leakage.ipynb` | executed offline visual walkthrough; reruns the production modules from frozen raw data in temporary output directories |
| `scripts/run_gate2_notebook.py` | executes all notebook cells without a server/kernel listener and saves verified outputs |
| `scripts/check_data_card_numbers.py` | checks documented curation figures, Gates 2–3 sign-off tables, Gate 3 features/baselines and Gate 4 training-only controls directly against logged results |
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
.venv/bin/python -m src.gate4_training CHEMBL284 --n-jobs 4  # training-only tuned descriptor control
.venv/bin/python -m src.gate4_fingerprints CHEMBL284  # training-only fold-local bit-length plans
.venv/bin/python -m src.gate4_combined CHEMBL284 --n-jobs 4  # resumable training-only combined search
.venv/bin/python -m src.gate4_combined CHEMBL284 --split scaffold --n-jobs 4  # resume scaffold only
.venv/bin/python -m src.gate4_combined CHEMBL284 --split random --n-jobs 4  # resume random only
.venv/bin/python scripts/run_gate2_notebook.py  # execute/save the visual walkthrough offline
.venv/bin/python scripts/check_data_card_numbers.py  # curation figures and sign-off tables vs logs
```

The Gate 4 control command is separate from `run_all.sh` while the combined-model
pipeline is being built. It defaults to both saved training pools, scaffold first;
`--split scaffold` or `--split random` runs just one. Both nested loops group on
the saved Gate 2 scaffold labels, including within random-split training. This
estimates transfer across scaffolds, even though the eventual random held-out
evaluation will describe interpolation. It computes no saved test prediction or
metric. Completed fits resume under `results/gate4_training/<TARGET>/`; a new
`--cache-dir /tmp/gate4-clean-build` starts a clean rebuild. Incompatible data,
settings, library versions or search/adapter code abort rather than silently reuse
the cache. A timestamp-only curation rerun retains the same search identity.
Each search directory's `report.json` wraps its results and checksum; its
`manifest.json` records data identities, feature specification and fold plans.
The recorded descriptor-control searches are under
[`results/gate4_training/CHEMBL284/`](results/gate4_training/CHEMBL284/).
Their outer-fold predictions cover only each saved training pool and estimate
transfer across scaffolds; they are not the final held-out results or a Gate 4
success-criteria assessment. `tests/test_gate4_training_results.py` checks every
completed fit against its fold plan, independently recomputes its metrics and
selections, checks reported uncertainty, and requires a rerun to reuse the caches
without fitting. Inner selections differ across outer folds; final parameters
are selected separately using grouped CV on the full saved training pool.
The [training-only control table](DATA_CARD.md#gate-4-training-only-descriptor-control)
reports OOF scores and intervals, fold selections, stability and compute timings.
The numeric documentation guard verifies it against checksummed cache reports;
mutation tests cover every displayed cell and interval endpoint.
Fingerprint plans are now recorded alongside the descriptor caches in each split's
`fingerprint_plan/` directory. The command above supports `--split` and
`--cache-dir` for the same resume/clean-rebuild workflow. It reuses the descriptor
search's grouped fold builder and Gate 3's Morgan folding measurements. Each
estimator's bit length minimises total within-compound identifier loss on its own
fitting structures, with smaller lengths breaking ties. No label or validation
structure enters that aggregate decision, including at the inner-fold level.
The fixed candidate lengths come from Gate 3's folding-sensitivity analysis; this
is a fidelity criterion within that memory budget, not a predictive-score claim.
Per-compound evidence, all fit-specific decisions, checksums and timings are logged.
The combined runner consumes these fit-specific decisions; using the final-fit
choice throughout CV would violate the intended boundary. The runner verifies the
plan again before search, checks every decision against that fit's training indices,
and records the applied bit length in each fit cache. It retains the control's grid,
selection rule and uncertainty helpers. After a split and feature block have a completed
search, `.venv/bin/python -m src.gate4_refit CHEMBL284 --split scaffold --block descriptors --n-jobs 4`
fits the full-training selection and saves the forest and checksum manifest under
`results/gate4_models/CHEMBL284/<split>/<block>/`. Both descriptor bundles have passed
reload checks on every training row. For a completed combined search, use `--block combined`.
Downstream code can use `from src.gate4_model_io import load_model`, then
`forest, metadata = load_model('CHEMBL284', 'scaffold', 'descriptors')`.
This verifies current input checksums, search code and library versions, the linked
completed search report, the binary checksum and the loaded forest's parameters
and feature count. It returns provenance including the seed, library versions and
input checksums without fitting, predicting or scoring. Missing or changed artifacts
raise an error. Keep each model bundle with its corresponding training-search
`manifest.json` and `report.json`; the loader also requires the verified input
snapshots. Load only this repository's trusted joblib files.
Partial production fits are saved in each split's `combined/` directory. The separate
split commands above were used to advance these searches; interruption preserves
completed fits. Complete combined OOF summaries, full-training parameter selections
and the final held-out evaluation remain pending.
`tests/test_gate4_training_results.py` independently audits every saved combined
fit's predictions, metrics, scaffold boundaries, fitting-row fingerprint decision,
input hashes and source identity. Its resume probes use cache snapshots and consume
cached fits without fitting or featurisation, stopping before the first incomplete fit.
The combined command
supports the same split, worker-count and clean-cache options as the control.
Scaffold-primary is Claude-and-GPT's working assumption pending Aadit's decision;
the random split is reported alongside it. The assessment split, repeat-dispersion
reference, code licence and final Gates 2 and 3 sign-off remain Aadit's decisions.

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
