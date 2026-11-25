# Data card - ChEMBL bioactivity for one protein target

Status: **draft, Gate 1.** Provenance and licence are filled in and verified. Composition and
known biases will be added from the raw-data profile of the recommended target (Gate 1 item 5).

This dataset is literature-reported bioactivity measurements. It does not establish, and this
project does not claim, the efficacy, safety or suitability of any compound.

## Provenance

| Field | Value | Where it comes from |
|---|---|---|
| Source | ChEMBL, EMBL-EBI, <https://www.ebi.ac.uk/chembl/> | - |
| Access | REST API `https://www.ebi.ac.uk/chembl/api/data`, through `chembl_webresource_client` (version pinned in `requirements.txt`) | `src/retrieve.py` |
| Release | `ChEMBL_37`, released 2026-05-01 | API status endpoint, saved as `results/licence/status.json` on 2026-11-25 (UTC); every raw pull also records it in its own `manifest.json` |
| Retrieval date | recorded per pull, in UTC, in `data/raw/<TARGET>/manifest.json` | `src/retrieve.py` |
| Target | to be fixed after instructor approval; shortlist and counts in `results/candidate_counts/` | `src/candidates.py` |

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

To be filled from the raw profile (Gate 1 item 5): record counts, assay types, activity types,
units, relation operators, organisms, replicates per compound, missing values.

## Known biases

To be filled from the raw profile (Gate 1 item 5). Literature-reported bioactivity is not a random
sample of chemistry. This section will name the specific skews visible in the pulled data.
