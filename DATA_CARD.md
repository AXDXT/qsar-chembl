# Data card - ChEMBL bioactivity for one protein target

Status: **Gate 2 in progress; human DPP4 (CHEMBL284) locked by Aadit.**
The Gate 1 sections below describe the frozen raw pull. The Gate 2 split decisions are
recorded at the end with the fixed leakage probe. The executed walkthrough is
`notebooks/gate2_dpp4_curation_and_leakage.ipynb` (rerun with
`.venv/bin/python scripts/run_gate2_notebook.py`). The complete curation narrative,
noise-reference review and sign-off remain.

This dataset is literature-reported bioactivity measurements. It does not establish, and this
project does not claim, the efficacy, safety or suitability of any compound.

## Provenance

| Field | Value | Where it comes from |
|---|---|---|
| Source | ChEMBL, EMBL-EBI, <https://www.ebi.ac.uk/chembl/> | - |
| Access | REST API `https://www.ebi.ac.uk/chembl/api/data`, through `chembl_webresource_client` (version pinned in `requirements.txt`) | `src/retrieve.py` |
| Release | `ChEMBL_37`, released 2026-05-01 | API status endpoint, saved as `results/licence/status.json` on 2026-11-25 (UTC); every raw pull also records it in its own `manifest.json` |
| Retrieval date | 2026-11-25, 06:45:22 to 06:46:36 UTC | `data/raw/CHEMBL284/manifest.json` |
| Target | CHEMBL284, Dipeptidyl peptidase 4; Homo sapiens, taxon 9606; SINGLE PROTEIN, one component, accession P27487 | unmodified `data/raw/CHEMBL284/target.json`; locked by Aadit for Gate 2 |
| Query | `new_client.activity.filter(target_chembl_id='CHEMBL284')`, no other retrieval filters | manifest, client version 0.10.9; API total and downloaded count both 8,397 |
| Raw files | `activities.jsonl.gz` and `target.json`, with sha256 in the manifest | `.venv/bin/python -m src.retrieve CHEMBL284` |

## Licence and attribution

Checked 2026-11-25 (UTC) from the ChEMBL FTP release directory. Copies are saved, unmodified and
with their sha256, in `results/licence/` by `.venv/bin/python -m src.provenance`.

- **Licence:** Creative Commons Attribution-ShareAlike 3.0 Unported (CC BY-SA 3.0).
  <https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/LICENSE>, as described in
  <https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/README>.
- **Attribution:** "ChEMBL data is from http://www.ebi.ac.uk/chembl - the version of ChEMBL is
  chembl_37." EBI also asks that ChEMBL IDs be kept and the release number be shown clearly
  (<https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/REQUIRED.ATTRIBUTION>).
- **Citation:** Mendez D, Gaulton A, Bento AP, et al. ChEMBL: towards direct deposition of bioassay
  data. *Nucleic Acids Res.* 2019 47(D1):D930-D940. DOI: 10.1093/nar/gky1075. This is the "primary
  current citation" in REQUIRED.ATTRIBUTION, which also gives the full author list.
- **Derived data:** our reading, not legal advice, is that curated subsets are adaptations and must
  be shared under CC BY-SA 3.0 with the attribution above.
- **Not verified:** the general EMBL-EBI web-service terms of use, and any newer citation shown on
  the ChEMBL web portal.

## Composition

Measured by `.venv/bin/python -m src.profile CHEMBL284`; complete field distributions and missing
counts are in `results/profiles/CHEMBL284/profile.json` and `profile.md`. Source-file checksums
are verified before and after profiling. This is a count-only description; no records are curated.

| Property | Observed raw composition |
|---|---|
| Records / distinct molecule IDs | 8,397 / 5,752 |
| Assay type (API codes) | B: 8,174; A: 177; F: 46 |
| Activity types (selected) | IC50: 5,818; Inhibition: 734; Ki: 709; k_off: 325; kon: 325; full distribution in profile |
| Standard units (selected) | nM: 6,474; uM: 7; %: 880; s-1: 327; missing: 520; other unit strings in profile |
| Standard relations | `=`: 7,116; `>`: 839; `<`: 82; `>=`: 2; `<=`: 1; `>>`: 5; `~`: 2; missing: 350 |
| Target organisms | Homo sapiens: all 8,397; assay organism is absent from every activity response, so host organism is not inferred |
| Rows per molecule | minimum 1, median 1, maximum 96; 1,084 molecules have multiple rows; 2,645 rows beyond the first per molecule |
| Selected missing fields | standard_value: 1,003; pchembl_value: 3,006; canonical_smiles: 5; standard_relation: 350; standard_units: 520 |
| Identifiers and descriptions | no missing activity IDs, molecule IDs, assay IDs, document IDs, assay descriptions or target organisms |

Per-molecule counts are in `compound_row_counts.csv`. These are repeated reports across endpoints
and assays, not verified experimental replicates. An API null or absent key is missing; blank
strings are also counted as missing. Empty containers are retained, with no biological interpretation.

The count-only Gate 2 approximation independently reproduces the candidate-count snapshot:
8,397 total -> 8,174 binding -> 5,781 IC50 -> 5,781 human target -> 5,655 molar units ->
4,859 exact-relation records, representing **4,004 distinct molecule IDs**, 337 assays and 258
documents. All 4,859 pass the separate positive-finite-value diagnostic. No structures, salts,
conflicting measurements or assay quality have been checked; these are not final curated counts.

## Known biases

Literature-reported bioactivity is not a random sample of chemistry. The measured skews and
limitations visible here are:

- Binding-type rows dominate (8,174 of 8,397), but activity types and units vary. Even the headline
  subset spans 337 assays and 258 documents; assay conditions and endpoint comparability require
  later review. IC50 values must not be treated as interchangeable with Ki or kinetic endpoints.
- Censoring is substantial: 839 raw rows use `>` and 82 use `<`, with additional bound operators.
  The exact-relation rule drops 796 rows after binding/IC50/human/molar filtering. Removing bounds
  can change the observed activity distribution; no potency distribution is inferred at Gate 1.
- Repeated reporting is uneven (up to 96 rows for one molecule). Treating rows as independent
  compounds could overweight familiar compounds and leak their information across later splits.
  No replicate aggregation or split has been performed.
- Missing numeric values (1,003 rows) and structures (5 rows) need later quality review.
  Null validity comments do not prove a measurement is valid; a supplied SMILES does not prove
  a structure passes curation.
- All target records are human, which limits organism scope. Assay-host metadata is absent in
  this activity pull; target organism alone does not establish expression host or assay context.

Gate 2 and later work must examine these limitations before any modelling. Selecting this target
makes no claim about compound efficacy, safety or suitability.

## Gate 2 split design and structural audit

Command: `.venv/bin/python -m src.splits CHEMBL284`, also run by `./run_all.sh` after
curation. Input: the 3,965 standardised, aggregated compounds in `compounds.csv`, whose
checksum and link to the frozen raw manifest are checked before splitting. Neither raw
activity rows nor replicate counts weight the allocation. Full outputs, assignment CSVs,
nearest-neighbour CSVs and every flagged pair are in `results/splits/CHEMBL284/`.

- **Random:** seed 42, approximately 80/20, at compound level. This measures interpolation
  within the collected chemical series; it does not measure transfer to unseen series.
- **Scaffold:** exact Bemis–Murcko cores with atom/bond types retained, stereochemistry
  ignored. Groups are shuffled with seed 42, then whole groups are included if they
  bring the test count closer to 20%. This uses no activity labels. Singleton cores
  follow the same rule. All 16 acyclic compounds share the empty scaffold and stay
  together. Whole groups can prevent an exact 80/20 allocation; no retry is made.
- **Similarity clusters deferred:** this comparison asks whether reusing exact cores
  makes a random evaluation optimistic. A third split would test an additional,
  fingerprint-threshold-specific extrapolation question. The exhaustive similarity
  audit measures residual analogue overlap without claiming scaffold holdout eliminates it.
  A later evaluation requiring stronger similarity separation would need that third split.

| Observed property | Random | Scaffold |
|---|---:|---:|
| Train / test compounds | 3,172 / 793 | 3,163 / 802 |
| Full InChIKey overlap (asserted zero) | 0 | 0 |
| Shared scaffolds | 268 | 0 |
| Test compounds with training neighbour Tanimoto >=0.90 | 120 | 44 |

The similarity definition is fixed before any model: nonchiral Morgan bit fingerprints,
radius 2, 2,048 bits, with Tanimoto >=0.90 flagged. This catches stereoisomers and very
close analogues but can also flag fingerprint collisions. Shared first InChIKey blocks
(connectivity) are audited independently. Every test compound is compared with every
training compound. Near pairs are reported, not silently removed to improve an evaluation.
The random split shares a core with training for 582/793 test compounds. Exact-identity
isolation therefore does not establish independence of the chemical series. The scaffold
test still has 44 close-neighbour compounds: its result must be described as transfer to
unseen exact cores, with this residual similarity limitation.

The saved audit includes sizes, core counts, singleton handling, target quantiles and
test compounds outside the training pIC50 range. These diagnose distribution shift;
they never determine assignments.

### Fixed leakage probe, not a baseline

Command: `.venv/bin/python -m src.leakage_probe CHEMBL284`, also run by `./run_all.sh`.
The only fitted model is a fixed Tanimoto 1-nearest-neighbour probe on the same
nonchiral Morgan radius-2, 2048-bit fingerprints used by the structural audit.
Fitting stores training compound-median labels; predictions average the labels of
all neighbours tied at maximum similarity. No tuning, split retries, activity-based
feature selection or replicate weighting occurs. Saved assignment checksums,
structures, labels and deterministic allocation are verified before prediction;
zero InChIKey overlap and scaffold isolation are asserted again.

| Metric (95% descriptive interval) | Random | Scaffold |
|---|---:|---:|
| RMSE, log units | 0.970 [0.857, 1.088] | 1.129 [0.950, 1.308] |
| MAE, log units | 0.618 [0.558, 0.677] | 0.777 [0.654, 0.919] |
| R2 | 0.505 [0.381, 0.610] | 0.390 [0.105, 0.585] |

Intervals use 2,000 seed-42 bootstrap draws of whole held-out exact scaffold groups
(472 random / 341 scaffold), preserving compounds within groups. They condition on
these fixed training sets, predictions and splits. They do not include variability
from retraining or different split seeds; similarity across distinct cores can still
couple bootstrap groups. Intervals overlap and are not a significance test of the gap.
Full precision, input/output hashes, subgroup metrics and predictions are in
`results/leakage_probe/CHEMBL284/`.

The observed RMSE increase is **0.159 log units**. Under random allocation, the probe
can reuse a training core for 582/793 test compounds and a >=0.90-similarity neighbour
for 120/793. The largest five exact-core groups contain 133, 66, 60, 52 and 40
compounds; random places 26, 16, 10, 4 and 7 of them respectively in test, leaving
members of all five groups in training. This is the actual analogue-series reuse
that exact-identity checks alone miss. Random-test RMSE is 0.921 with a training
core versus 1.093 without one, and 0.819 with a close neighbour versus 0.994 without.
These diagnostic subsets differ in chemistry and are not controlled comparisons.

Scaffold allocation removes shared exact cores and shared connectivity keys but
retains 44 close-neighbour test compounds. Their RMSE is 1.161 versus 1.127 for the
other 758: closeness does not guarantee assay agreement. The overall gap is consistent
with easier series interpolation under random allocation, but also includes changed
chemistry and target distributions; it is not an isolated causal leakage penalty.
Over-reported compounds have already been collapsed to one median per standardised
InChIKey, so the probe cannot memorise duplicate rows across a boundary. The 723
random / 689 scaffold training-side repeat rows collapsed during curation carry no
extra fitting weight. This gate demonstrates evaluation sensitivity, not a model
selection result, and makes no efficacy, safety or suitability claim.

**RDKit pin decision:** use 2026.9.1. The 2026.9.2 macOS arm64 PyPI wheel raises
`RuntimeError: bad any cast` in native `MurckoDecompose` even for benzene, reproduced
locally and reported in [RDKit issue #9836](https://github.com/rdkit/rdkit/issues/9836).
The preceding wheel runs the native implementation. After changing the pin, the full
offline suite and pipeline pass, and both curated CSVs remain byte-identical to the
committed curation. No private reference algorithm is substituted for native Murcko.

This pin exists because of a **platform-specific packaging fault, not a bug in RDKit's
source**. On the issue thread the same machine works on the 2026.9.1 wheel and on the
conda-forge 2026.9.2 build, and Linux 2026.9.2 wheels are unaffected; a maintainer
suggested the PyPI builds may have switched to nanobind, and the report was redirected
to `kuelumbus/rdkit-pypi`. A reader on Linux, or on conda-forge, is therefore not
working around an RDKit defect and can relax the pin once the macOS arm64 wheel is
rebuilt — the scaffold algorithm itself is unchanged across these versions. Verified
against the primary source by Aadit; the local reproduction is gpt's.

**Open statistical review before sign-off:** the existing 0.310-log-unit noise figure
conditions on differing replicate groups with spread <=1 log unit. It describes that
selected subset, not an unconditional measurement floor for compound-median predictions.
Similar reported values alone do not establish independent measurements. A test RMSE
below this figure does not by itself prove leakage: the prediction target is an aggregated
median and the replicate subset may differ from the held-out population. The noise report's
caveat now explicitly states that a lower RMSE alone cannot prove leakage. Its legacy
`recommended_for_gate4` field and supporting interpretation still need review alongside
all-group dispersion before Gate 4 can cite it as a floor.
