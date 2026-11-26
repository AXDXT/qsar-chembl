# Gate 3 feature cache: CHEMBL284

Command: `.venv/bin/python -m src.features CHEMBL284`. Gate 3 feature cache. No model is fitted, scored or selected here.

No run timestamp is recorded, which is what lets a rerun be byte-identical. Provenance is carried by the input checksums below.

canonical_smiles from compounds.csv: sanitised, largest fragment, neutralised and InChIKey-verified by Gate 2. Raw SMILES are never read.

| Reconciliation | Count |
|---|---:|
| Curated compounds in | 3,965 |
| Rows featurised | 3,965 |
| Featurisation failures | 0 |
| Non-finite descriptor values | 0 |
| Compounds dropped | 0 |

Every failure, if any, is listed in `featurisation_failures.csv`; the run aborts rather than continuing with a short matrix. Ascending standardised full InChIKey, identical to the saved split assignment CSVs, so the cached matrices join to them positionally.

## Descriptor block (12 columns)

Fixed and declared before any model is fitted. Chosen for named, monotone interpretability, not by any score.

| Descriptor | Why it is in the set |
|---|---|
| MolWt | Average molecular weight; the size term of Lipinski's rule of five. |
| HeavyAtomCount | Size as an integer atom count, free of hydrogen-weight and isotope conventions. |
| MolLogP | Crippen atom-contribution octanol/water logP; the lipophilicity term of the rule of five. |
| MolMR | Crippen molar refractivity; polarisability and bulk, the Ghose-filter companion to logP. |
| TPSA | Topological polar surface area; the polarity term of Veber's permeability rules. |
| NumLipinskiHBD | Lipinski hydrogen-bond donors (N-H and O-H count), the rule-of-five donor definition. |
| NumLipinskiHBA | Lipinski hydrogen-bond acceptors (nitrogen and oxygen count), the rule-of-five definition. |
| NumRotatableBonds | Conformational flexibility; the Veber rotatable-bond term. |
| RingCount | Number of SSSR rings; total ring content. |
| NumAromaticRings | Aromatic ring count; aromatic versus saturated ring character. |
| FractionCSP3 | Fraction of sp3-hybridised carbons; saturation and three-dimensional character. |
| NumHeteroatoms | Heteroatom count; polar functionality that TPSA alone does not separate. |

NumAliphaticRings is deliberately omitted: RDKit partitions SSSR rings into aromatic and aliphatic, so it is exactly RingCount - NumAromaticRings. The identity is measured, not assumed; see collinearity.aliphatic_ring_identity_max_deviation.
Measured maximum deviation: 0.0.

NaN/inf: Every descriptor value is checked with numpy.isfinite. A non-finite value is written to featurisation_failures.csv and aborts the run. Nothing is imputed, clipped or dropped, so no imputation statistic exists that could have been fitted on the full dataset.

Scaling: None. No centring, scaling, variance filtering or feature selection is applied at Gate 3. The descriptor baseline is a tree ensemble, which is invariant to monotone per-column rescaling, so a scaler would add a fitted object with no benefit. Any transform a later gate adds must be fitted on training rows only.

Statelessness: Provably stateless: each row is a function of one molecule alone. No row, label, split or dataset-level statistic enters another row's features. test_features.py checks this by featurising single molecules in isolation.

| Split (training rows only) | Rows | Largest absolute correlation | Pair |
|---|---:|---:|---|
| random | 3,172 | 0.980 | HeavyAtomCount / MolWt |
| scaffold | 3,163 | 0.982 | HeavyAtomCount / MolWt |

Computed on each split's training rows only, with no label involved, and used for nothing. Reported so a Gate 4 interpretation cannot mistake two near-duplicate columns for independent evidence.

## Fingerprint block (2,048 bits)

Morgan bit vector, radius 2, includeChirality=False, via the current rdFingerprintGenerator API. 2,048 bits is inherited unchanged from the Gate 2 similarity audit and leakage probe, so the audited near-duplicate and shared-core counts describe this model's actual input space. It was not chosen by any metric; folding_evidence records what it costs.

| Bit length | Occupied bits | Mean on-bits per compound | Compounds losing >=1 identifier |
|---:|---:|---:|---:|
| 1,024 | 1,024 | 51.3 | 2,523 |
| 2,048 | 2,028 | 51.8 | 1,854 |
| 4,096 | 3,692 | 52.1 | 905 |

The curated set contains 9,309 distinct radius-<=2 substructure identifiers, 52.4 per compound on average. Folding maps several identifiers onto one bit; the final column counts the compounds where that happens. The alternative lengths are measured for the record and were not scored, so nothing here selects a bit length.

Neither block distinguishes stereoisomers: the fingerprint ignores chirality by Gate 2's choice, and no descriptor counts stereocentres. 286 curated compounds share a connectivity skeleton with another compound, so that limitation is measurable rather than hypothetical.

Literature/patent-reported IC50; no efficacy, safety or suitability claim.
