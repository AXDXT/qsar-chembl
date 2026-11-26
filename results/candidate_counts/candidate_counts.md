# Gate 1 candidate counts

Measured from ChEMBL; each linked log records the release, UTC times, target metadata, query and full count cascade.

**Count-only approximation of the Gate 2 cascade**, not curated compound counts. Filters: binding assays, standard IC50, human target organism, molar units convertible to nM, exact '=' relation. A positive-finite-value diagnostic is logged separately and does not change the headline counts. The unit whitelist covers plain molar units; uncommon scaled strings are excluded pending interpretation. Distinct compounds means nonempty ChEMBL molecule IDs; no structures or replicates are curated.

| Target | Protein | Total records | Approx. filtered records | Approx. distinct compounds | Log |
|---|---|---:|---:|---:|---|
| CHEMBL220 | Acetylcholinesterase | 20132 | 7646 | 6134 | [ChEMBL_37, 2026-11-25T06:33:00+00:00](CHEMBL220_20261125T063300727419Z.json) |
| CHEMBL203 | Epidermal growth factor receptor | 58847 | 18381 | 10524 | [ChEMBL_37, 2026-11-25T06:36:17+00:00](CHEMBL203_20261125T063617396212Z.json) |
| CHEMBL4822 | Beta-secretase 1 | 27804 | 12893 | 8272 | [ChEMBL_37, 2026-11-25T06:38:01+00:00](CHEMBL4822_20261125T063801514352Z.json) |
| CHEMBL284 | Dipeptidyl peptidase 4 | 8397 | 4859 | 4004 | [ChEMBL_37, 2026-11-25T06:38:49+00:00](CHEMBL284_20261125T063849168172Z.json) |

## Observed cautions

- **CHEMBL220**: 965 assays across 649 documents survive; assay conditions and replicate agreement still need review. 885 rows drop at the relation step (including censored, approximate or missing relations). Raw target-organism counts: {"Homo sapiens": 20132}. 1512 surviving rows exceed the unique molecule count; these are not independent compounds.
- **CHEMBL203**: 1910 assays across 1022 documents survive; assay conditions and replicate agreement still need review. 5234 rows drop at the relation step (including censored, approximate or missing relations). Raw target-organism counts: {"Homo sapiens": 58847}. 7857 surviving rows exceed the unique molecule count; these are not independent compounds.
- **CHEMBL4822**: 702 assays across 445 documents survive; assay conditions and replicate agreement still need review. 702 rows drop at the relation step (including censored, approximate or missing relations). Raw target-organism counts: {"Homo sapiens": 27804}. 4621 surviving rows exceed the unique molecule count; these are not independent compounds.
- **CHEMBL284**: 337 assays across 258 documents survive; assay conditions and replicate agreement still need review. 796 rows drop at the relation step (including censored, approximate or missing relations). Raw target-organism counts: {"Homo sapiens": 8397}. 855 surviving rows exceed the unique molecule count; these are not independent compounds.

Target choice remains provisional until Aadit checks class uniqueness and obtains instructor approval. A raw profile and the final recommendation are still due; none of these counts establishes compound efficacy, safety or suitability.
