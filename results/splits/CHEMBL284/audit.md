# Structural leakage audit: CHEMBL284

Command: `.venv/bin/python -m src.splits CHEMBL284`; seed 42. No model fitted.

Exact Bemis-Murcko atom/bond scaffold without stereochemistry; all acyclic compounds share the empty scaffold. Seeded group order; greedily approach 20% test without splitting groups. Singletons follow the same allocation rule. No pIC50 used for assignment.

Backend: RDKit native MurckoDecompose; RDKit 2026.09.1. The wrapper source checksum is logged in audit.json.

Deferred: exact scaffold holdout directly tests analogue-series reuse, while the exhaustive similarity audit measures what it misses. A similarity-cluster split would answer the additional question of transfer beyond a chosen fingerprint threshold, which is outside this random-versus-scaffold leakage comparison. Scaffold holdout is not claimed to guarantee low similarity; residual pairs are retained and reported.

| Property | Random | Scaffold |
|---|---:|---:|
| InChIKey overlap (asserted zero) | 0 | 0 |
| Shared scaffolds | 268 | 0 |
| Test compounds with a training scaffold | 582 | 0 |
| Shared connectivity keys | 75 | 0 |
| Cross-boundary pairs with Tanimoto >= 0.90 | 189 | 69 |
| Test compounds with Tanimoto >= 0.90 neighbour | 120 | 44 |
| Test compounds outside training pIC50 range | 0 | 0 |
| train: compounds | 3172 | 3163 |
| train: scaffolds | 1233 | 1096 |
| train: global_singleton_scaffolds | 751 | 714 |
| train: acyclic_compounds | 14 | 16 |
| train: pIC50 min | -0.041 | -0.041 |
| train: pIC50 p05 | 4.472 | 4.518 |
| train: pIC50 median | 7.174 | 7.194 |
| train: pIC50 p95 | 8.959 | 8.931 |
| train: pIC50 max | 10.097 | 10.097 |
| test: compounds | 793 | 802 |
| test: scaffolds | 472 | 341 |
| test: global_singleton_scaffolds | 199 | 236 |
| test: acyclic_compounds | 2 | 0 |
| test: pIC50 min | 2.417 | 3.000 |
| test: pIC50 p05 | 4.265 | 4.174 |
| test: pIC50 median | 7.191 | 7.073 |
| test: pIC50 p95 | 8.886 | 8.959 |
| test: pIC50 max | 10.000 | 9.921 |

Similarity: nonchiral Morgan radius 2, 2048 bits; exhaustive train/test Tanimoto comparisons. The 0.90 threshold is a fixed audit flag, not proof of chemical identity. Each test compound's nearest training neighbour and every flagged pair (including shared InChIKey connectivity blocks) are saved in CSV.

The random split can reuse shared cores and close analogues despite zero exact compound overlap. The scaffold split blocks exact cores but can retain similar ring substitutions or fingerprint collisions. Counts above measure those risks; the fixed probe and its uncertainty remain a separate Gate 2 deliverable.

pIC50 summaries diagnose distribution shift; they never determine assignments. No split is retried to improve target coverage or a score.

Literature/patent-reported IC50; no efficacy, safety or suitability claim.
