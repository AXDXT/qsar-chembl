# Data card - ChEMBL bioactivity for one protein target

Status: **Gate 1, provisional target awaiting class uniqueness and instructor approval.**
Both build-room agents recommend human DPP4 (CHEMBL284); this card describes its frozen raw pull.

This dataset is literature-reported bioactivity measurements. It does not establish, and this
project does not claim, the efficacy, safety or suitability of any compound.

## Provenance

| Field | Value | Where it comes from |
|---|---|---|
| Source | ChEMBL, EMBL-EBI, <https://www.ebi.ac.uk/chembl/> | - |
| Access | REST API `https://www.ebi.ac.uk/chembl/api/data`, through `chembl_webresource_client` (version pinned in `requirements.txt`) | `src/retrieve.py` |
| Release | `ChEMBL_37`, released 2026-05-01 | API status endpoint, saved as `results/licence/status.json` on 2026-11-25 (UTC); every raw pull also records it in its own `manifest.json` |
| Retrieval date | 2026-11-25, 06:45:22 to 06:46:36 UTC | `data/raw/CHEMBL284/manifest.json` |
| Target | CHEMBL284, Dipeptidyl peptidase 4; Homo sapiens, taxon 9606; SINGLE PROTEIN, one component, accession P27487 | unmodified `data/raw/CHEMBL284/target.json`; final choice requires instructor approval |
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

Gate 2 and later work must examine these limitations before any modelling. The target recommendation
is provisional and makes no claim about compound efficacy, safety or suitability.
