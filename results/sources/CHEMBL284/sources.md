# Data sources in the CHEMBL284 raw pull

Release ChEMBL_37; checked 2026-11-25T07:04:10+00:00 with `.venv/bin/python -m src.sources CHEMBL284`.
Row counts are over all raw rows, before curation.

| src_id | Short name | Description | Raw rows |
|---|---|---|---:|
| 1 | LITERATURE | Scientific Literature | 7109 |
| 37 | BINDINGDB | BindingDB Patent Bioactivity Data | 1281 |
| 38 | PATENT | SureChEMBL Patent Bioactivity Data | 2 |
| 54 | DONATED_PROBES | SGC Frankfurt - Donated Chemical Probes | 2 |
| 65 | LIT_EUBOPEN_CGL | EUbOPEN Chemogenomic Library Literature Data | 1 |
| 72 | LIT_CHEM_PROBES | Chemical Probe data from Scientific Literature | 2 |

Full API text, including each source's comment, is in `sources.json`.
