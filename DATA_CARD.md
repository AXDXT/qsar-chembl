# Data card - ChEMBL bioactivity for one protein target

Status: **Gate 2 in progress; human DPP4 (CHEMBL284) locked by Aadit.**
The Gate 1 sections below describe the frozen raw pull. The Gate 2 split decisions are
recorded at the end with the fixed leakage probe. The executed walkthrough is
`notebooks/gate2_dpp4_curation_and_leakage.ipynb` (rerun with
`.venv/bin/python scripts/run_gate2_notebook.py`). The Gate 2 curation cascade, its
judgement calls, the target variable and the measured replicate dispersion are recorded
below under "Gate 2 curation cascade"; the split design and the fixed leakage probe follow
it. The noise-reference choice for Gate 4 and the room sign-off remain.

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

## Gate 2 curation cascade

Command: `.venv/bin/python -m src.curate CHEMBL284`, also run by `./run_all.sh`. It reads only
the frozen raw pull, verifies every raw checksum before and after, and re-verifies them at the
end; nothing in `data/raw/` is written. Outputs: `data/curated/CHEMBL284/` (with its own
manifest) and `results/curation/CHEMBL284/` (`curation_log.json`, `cascade.csv`, `cascade.md`,
`cascade_funnel.png`, `units_refused.csv`, `structure_failures.csv`, `pic50_distribution.png`).
The visual walkthrough is `notebooks/gate2_dpp4_curation_and_leakage.ipynb`.

| Step | Rule | In | Out | Removed |
|---|---|---:|---:|---:|
| 1 Binding assays | `assay_type == 'B'` | 8,397 | 8,174 | 223 (A: 177, F: 46) |
| 2 Single activity type | `standard_type == 'IC50'` | 8,174 | 5,781 | 2,393 across 28 other types (Inhibition 723, Ki 709, k_off 325, kon 325, …) |
| 3 Single target organism | human target, tax id 9606 | 5,781 | 5,781 | 0 |
| 4 Units harmonised to nM | molar unit string and finite positive value | 5,781 | 5,650 | 131 refused (78 missing unit, 48 `ug.mL-1`, 5 unusable value) |
| 5 Censored values removed | `standard_relation == '='` exactly | 5,650 | 4,859 | 791 (`>`: 732, `<`: 54, `>>`: 4, `>=`: 1) |
| 6 Structures standardised | sanitise → largest fragment → neutralise → canonical SMILES → InChIKey | 4,859 | 4,855 | 4 (no SMILES supplied) |
| 7 Duplicates resolved | median pIC50 per standardised InChIKey | 4,855 records | **3,965 compounds** | 890 rows collapsed into medians |

The unit of the count changes at step 7: records become compounds. Every refusal at step 4 and
every failure at step 6 is itemised by row in a CSV; nothing is dropped silently. Step 6 took the
largest fragment for 855 salt or mixture rows and neutralised 2 charged structures. Step 7 found
476 compounds with replicates and flagged 108 whose replicates span more than 1 log unit.

**Judgement calls Aadit must be able to defend:**

- **Mass concentrations are refused, not converted.** 48 `ug.mL-1` rows are listed in
  `units_refused.csv` and excluded. Converting them needs the molecular weight of the material
  actually tested, which may be a salt and is not recorded here; converting would mix two
  measurement bases under one column.
- **Censored rows are removed, and the Gate 1 gap is reconciled.** Gate 1 measured 796 rows
  failing the exact-relation test; this cascade removes 791 at step 5. The 5-row difference is
  rows whose `standard_relation` *and* `standard_value` are both null: Gate 1 counted them as
  "not exactly `=`", while this cascade removes them one step earlier at step 4, because a row
  with no value cannot become a pIC50 whatever its relation says. The 4,859 rows that survive the
  relation rule are identical in both gates (`censoring_reconciliation` in `curation_log.json`,
  asserted by two tests). A `>` row states only that the IC50 exceeds the top concentration
  tested; keeping it as a value would put an invented number into the target.
- **Compound identity is the standardised InChIKey, not the ChEMBL molecule ID.** ChEMBL
  separates a salt from its parent, so keying on molecule ID would enter the same chemistry
  twice and could place the two copies on opposite sides of a split. Observed: 4,000 distinct
  molecule IDs and 3,965 parent IDs in the curated records; 32 compounds merge more than one
  molecule ID, absorbing 35 IDs. Molecule and parent IDs are kept as audit columns. Stereoisomers
  and tautomers keep distinct full InChIKeys (286 pairs share a skeleton), so they stay separate
  compounds — a defensible choice to disagree with, and the skeleton counts are logged for it.
- **The one molecule with 96 raw rows is CHEMBL1422.** 63 of its rows survive to step 7 (the rest
  are other endpoints: Inhibition, Kd, kinetics, or censored). Merged with three further ChEMBL
  IDs for the same standardised structure, it becomes **one** compound with 71 accepted records,
  median pIC50 7.74 and a spread of 4.14 log units, so it is also in the 108 flagged compounds.
  Left as rows it would carry the weight of a whole analogue series and appear on both sides of a
  random split.
- **`src_id` 37 is BindingDB patent bioactivity data** (resolved from the API once by
  `.venv/bin/python -m src.sources CHEMBL284`, snapshot in `results/sources/CHEMBL284/`). Curated
  records: 4,333 from src_id 1 (primary literature) and 522 from src_id 37. They are kept, because
  they are IC50 measurements against the same target, and retained as a column so the decision can
  be revisited or tested. Patent and literature reporting regimes may differ in ways not visible here.
- **Median, not mean, for replicates**, because the replicate sets have a long tail of cross-assay
  disagreement that a mean would follow. Flagged high-spread compounds are reported, not discarded.
- **ChEMBL's own validity comments are not used as an eighth filter.** 133 curated rows are marked
  "Outside typical range" and 53 "Potential transcription error"; they are carried through as a
  column and counted here. Adding that filter is Aadit's decision, and these counts say what is
  at stake.

### Target variable

pIC50 = 9 − log10(IC50 in nM), equivalently −log10(IC50 in mol/L), per accepted record; the
compound value is the median across its records. Over 3,965 compounds:
min −0.04, 5th percentile 4.42, median 7.17, mean 7.03, SD 1.33, 95th percentile 8.95, max 10.10.
Distribution figure: `results/curation/CHEMBL284/pic50_distribution.png`.

Pathologies checked and reported rather than trimmed: there is **no** pile-up at a reporting
convention such as pIC50 5.0 (10 µM) or 4.0 (100 µM) — the densest 0.1-log bins all sit in the
6.6–7.9 interior. A rounding fingerprint is visible but small: 80 of 3,965 compounds sit on an
exact integer pIC50, concentrated at 8.0 (33) and 9.0 (23). The extremes are implausible as
measurements: 2 compounds below pIC50 2 (IC50 above 10 mM, past the solubility of most of these
compounds) and 108 below pIC50 4. They are **kept and flagged**, because an activity threshold is
a modelling choice belonging to a later gate; any later gate that trims them must say so and
re-report the count. A hard assay ceiling would show as identical values at the top concentration
tested; the `>` rows where such a ceiling normally hides were removed at step 5 and counted there,
so the ceiling has been taken out of the target rather than flattened into it.

### Measured replicate dispersion (the Gate 4 reference number)

From the step-7 replicate groups: 476 replicated compounds, 1,366 individual reports, 890 degrees
of freedom. Pooled within-compound SD, in pIC50 log units and so directly comparable with a model
RMSE:

| Subset | Groups | d.o.f. | Pooled SD |
|---|---:|---:|---:|
| All repeated groups | 476 | 890 | 0.828 |
| Spread ≥ 0.01 log unit | 254 | 597 | 1.011 |
| 0.01 ≤ spread ≤ 1 log unit | 146 | 208 | **0.310** |

The unfiltered descriptive dispersion is **0.828**. The **0.310** and **1.011** figures
are sensitivity statistics selected by observed spread, not established lower or upper bounds
on measurement noise or prediction error. No Gate 4 escape threshold is approved here. The
legacy `noise_floor.recommended_for_gate4` value is retained as a candidate with explicit
`pending_human_decision` status, not as a settled recommendation.

There are 222 groups with spread below 0.01 log unit: **199** are exactly equal and **23**
have small nonzero differences. Similar values alone do not prove copied reporting. The
histogram's final bin contains **56** groups with spread **≥3**, comprising **29** exactly at
3 and **27** above 3 (max 8.32). Finite bins include the lower edge and exclude the upper edge.

Four caveats travel with this number, and the last one is load-bearing. These are repeated
literature and patent reports, not a designed replicate experiment, so none of these figures is
pure assay noise. Only replicated compounds contribute: 3,489 of 3,965 compounds are measured
once, so dispersion is measured on 12.0% of the curated set. Independence of the repeat reports has
not been established — equal reports may be copied rather than independently agreeing.
And **0.310 is conditional dispersion for groups with 0.01 ≤ spread ≤ 1 log unit, not
an unconditional lower bound on held-out RMSE**: predictions target compound medians rather than
individual reports, the repeated subset may differ from the test population, and a test RMSE below
0.310 would not by itself prove leakage or memorisation. Which figure Gate 4 may cite as its
escape clause is still open; see "Needs decision before sign-off" at the end of this card.

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

**Needs decision before sign-off:** the existing 0.310-log-unit noise figure
conditions on replicate groups with 0.01 <= spread <=1 log unit. It describes that
selected subset, not an unconditional measurement floor for compound-median predictions.
Similar reported values alone do not establish independent measurements. A test RMSE
below this figure does not by itself prove leakage: the prediction target is an aggregated
median and the replicate subset may differ from the held-out population. The noise report's
caveat now explicitly states that a lower RMSE alone cannot prove leakage. My recommendation
(gpt) is to cite **0.828** as the all-reported-repeat dispersion and **0.310 / 1.011** as
sensitivity analyses, without using any as an automatic Gate 4 escape clause. Selecting on
spread truncates the disagreement being estimated, and the labels being predicted are medians.
A defensible irreducible-error threshold would require independent repeat provenance, comparable
assay conditions, and an uncertainty model for the compound medians. Aadit must approve how Gate 4
will use these descriptive statistics; the current evidence cannot establish such a threshold.
