# Gate 2 fixed leakage probe (not a baseline)

Command: `.venv/bin/python -m src.leakage_probe CHEMBL284`.

Fixed Tanimoto 1-nearest-neighbour leakage probe (not a baseline)

Training pIC50_median at maximum Tanimoto; average all exact maximum ties. Store training fingerprints and compound-median labels only; no tuning.
Features: nonchiral Morgan radius 2, 2048 bits; identical for both saved splits. All compounds carry equal weight after median aggregation. No parameter or split search.

| Metric (95% interval) | Random | Scaffold |
|---|---:|---:|
| RMSE | 0.970 [0.857, 1.088] | 1.129 [0.950, 1.308] |
| MAE | 0.618 [0.558, 0.677] | 0.777 [0.654, 0.919] |
| R2 | 0.505 [0.381, 0.610] | 0.390 [0.105, 0.585] |

RMSE and MAE are pIC50 log units; R2 is dimensionless.

95% percentile bootstrap, 2000 resamples of whole test exact-scaffold groups, compound-weighted metrics; fixed training set, split and predictions. Not uncertainty across split seeds or retraining. Cross-core similarity can still couple groups. Intervals are descriptive, not a significance test of the gap.
Test bootstrap groups: 472 random and 341 scaffold. Seed 42.

## What this dataset lets the probe reuse

Random holds out 793 compounds, but 582 already have a training core (268 shared cores); 120 have a training neighbour at Tanimoto >=0.90. Scaffold holds out 802 compounds with zero shared cores, yet retains 44 such close-neighbour compounds. A stored analogue's label can therefore interpolate a familiar series under random allocation. Exact InChIKey isolation alone does not prevent this structural information reuse.

Largest exact-core groups (canonical scaffold strings and sizes are in probe.json): 133 compounds (26 random-test), 66 compounds (16 random-test), 60 compounds (10 random-test), 52 compounds (4 random-test), 40 compounds (7 random-test).

Observed scaffold-minus-random RMSE: +0.159 log units. The contrast measures this probe's sensitivity to the evaluation boundary. It does not isolate a causal leakage penalty: the test compounds and activity distributions also differ. No held-out labels affected features, assignments, predictions or parameter choices.

| Split / diagnostic subgroup | Test compounds | RMSE |
|---|---:|---:|
| random: shared_train_scaffold=True | 582 | 0.921 |
| random: shared_train_scaffold=False | 211 | 1.093 |
| random: near_train_neighbour=True | 120 | 0.819 |
| random: near_train_neighbour=False | 673 | 0.994 |
| scaffold: shared_train_scaffold=True | 0 | undefined |
| scaffold: shared_train_scaffold=False | 802 | 1.129 |
| scaffold: near_train_neighbour=True | 44 | 1.161 |
| scaffold: near_train_neighbour=False | 758 | 1.127 |

Subgroups explain the boundary; they were not used to tune or exclude compounds. They differ in chemistry and size and are not controlled experiments. Close neighbours need not agree in measured activity; compare their errors above rather than assuming similarity guarantees reliable prediction.

Over-reported compounds cannot leak identical rows here: one median label per standardised InChIKey enters each split, and exact overlap is asserted zero. Training-side repeat rows collapsed before fitting: 723 random / 689 scaffold. The probe never weights those rows. Scaffold holdout also removes shared connectivity keys; residual fingerprint similarity can reflect related cores, stereochemistry ignored by the fingerprint, or bit collisions.

These intervals do not turn conditional replicate dispersion into a universal RMSE floor. The predicted targets are compound medians, whereas replicate SD describes individual reported measurements on a selected repeated subset. A lower RMSE alone cannot establish leakage.

Literature/patent-reported IC50; no efficacy, safety or suitability claim.
