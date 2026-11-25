# Raw profile: CHEMBL284

Release: ChEMBL_37. Raw retrieval finished: 2026-11-25T06:46:36+00:00.
Command: `.venv/bin/python -m src.profile CHEMBL284`. Source hashes verified before and after profiling.

Raw records: **8397**. No rows are removed or changed.

## assay_type

| API value | Rows |
|---|---:|
| `A` | 177 |
| `B` | 8174 |
| `F` | 46 |

## standard_type

| API value | Rows |
|---|---:|
| `% Ctrl` | 1 |
| `% Inhibition of Control Values` | 1 |
| `Activity` | 98 |
| `Drug metabolism` | 18 |
| `EC50` | 4 |
| `IC5` | 2 |
| `IC50` | 5818 |
| `INH` | 1 |
| `Inhibition` | 734 |
| `K` | 62 |
| `K inact` | 2 |
| `K obs / 1` | 8 |
| `Ka` | 1 |
| `Kapp` | 1 |
| `Kd` | 26 |
| `Kdiss` | 1 |
| `Ki` | 709 |
| `Km` | 23 |
| `Ratio IC50` | 12 |
| `Ratio Ki` | 5 |
| `Selectivity index` | 68 |
| `Stability` | 80 |
| `T1/2` | 40 |
| `TIME` | 21 |
| `TdeltaS` | 1 |
| `Vmax` | 4 |
| `Vmax/Km` | 4 |
| `deltaG` | 1 |
| `deltaH` | 1 |
| `k_off` | 325 |
| `kon` | 325 |

## standard_units

| API value | Rows |
|---|---:|
| `%` | 880 |
| `/min` | 8 |
| `10'-2/s` | 3 |
| `10'-3/s` | 4 |
| `10'-4/s` | 12 |
| `10'-5/s` | 1 |
| `10'3/M/s` | 1 |
| `10'4/M/s` | 1 |
| `10'5/M/s` | 4 |
| `10'6/M/min` | 6 |
| `10'6/M/s` | 9 |
| `10'7/M/s` | 2 |
| `10'8/M/min` | 1 |
| `10^-4/s` | 6 |
| `10^2kJ/mol` | 2 |
| `10^3kJ/mol` | 1 |
| `10^4/M/s` | 1 |
| `10^5/M` | 1 |
| `10^5/M/s` | 5 |
| `<missing>` | 520 |
| `M-1 s-1` | 8 |
| `hr` | 61 |
| `nM` | 6474 |
| `nmol/hr/microg` | 4 |
| `s-1` | 327 |
| `uM` | 7 |
| `ug.mL-1` | 48 |

## standard_relation

| API value | Rows |
|---|---:|
| `<` | 82 |
| `<=` | 1 |
| `<missing>` | 350 |
| `=` | 7116 |
| `>` | 839 |
| `>=` | 2 |
| `>>` | 5 |
| `~` | 2 |

## target_organism

| API value | Rows |
|---|---:|
| `Homo sapiens` | 8397 |

## assay_organism

| API value | Rows |
|---|---:|
| `<missing>` | 8397 |

## src_id

| API value | Rows |
|---|---:|
| `1` | 7109 |
| `37` | 1281 |
| `38` | 2 |
| `54` | 2 |
| `65` | 1 |
| `72` | 2 |

## Repeated reporting per compound

Raw rows per molecule ID across all endpoints and assays; not experimentally verified replicates.

Exact per-molecule counts are in `compound_row_counts.csv`.

```json
{
  "distinct_molecule_ids": 5752,
  "rows_missing_molecule_id": 0,
  "molecules_with_multiple_rows": 1084,
  "rows_beyond_first_per_molecule": 2645,
  "minimum": 1,
  "median": 1.0,
  "maximum": 96,
  "histogram_rows_per_molecule": {
    "1": 4668,
    "2": 594,
    "3": 276,
    "4": 38,
    "5": 15,
    "6": 65,
    "7": 51,
    "8": 9,
    "9": 18,
    "10": 3,
    "11": 3,
    "13": 2,
    "14": 1,
    "20": 1,
    "22": 1,
    "27": 1,
    "36": 1,
    "40": 1,
    "43": 1,
    "46": 1,
    "72": 1,
    "96": 1
  },
  "interpretation": "Raw rows per molecule ID across all endpoints and assays; not experimentally verified replicates."
}
```

## Missing values

Absent key, null, or blank string; empty containers are retained as API values.

| Field | Absent | Null | Blank | Total missing |
|---|---:|---:|---:|---:|
| action_type | 0 | 7791 | 0 | 7791 |
| activity_comment | 0 | 7047 | 0 | 7047 |
| activity_id | 0 | 0 | 0 | 0 |
| activity_properties | 0 | 0 | 0 | 0 |
| assay_chembl_id | 0 | 0 | 0 | 0 |
| assay_description | 0 | 0 | 0 | 0 |
| assay_organism | 8397 | 0 | 0 | 8397 |
| assay_type | 0 | 0 | 0 | 0 |
| assay_variant_accession | 0 | 8397 | 0 | 8397 |
| assay_variant_mutation | 0 | 8397 | 0 | 8397 |
| bao_endpoint | 0 | 0 | 0 | 0 |
| bao_format | 0 | 0 | 0 | 0 |
| bao_label | 0 | 0 | 0 | 0 |
| canonical_smiles | 0 | 5 | 0 | 5 |
| data_validity_comment | 0 | 8018 | 0 | 8018 |
| data_validity_description | 0 | 8018 | 0 | 8018 |
| document_chembl_id | 0 | 0 | 0 | 0 |
| document_journal | 0 | 1286 | 0 | 1286 |
| document_year | 0 | 0 | 0 | 0 |
| ligand_efficiency | 0 | 3043 | 0 | 3043 |
| modality | 0 | 8397 | 0 | 8397 |
| molecule_chembl_id | 0 | 0 | 0 | 0 |
| molecule_pref_name | 0 | 7908 | 0 | 7908 |
| parent_molecule_chembl_id | 0 | 0 | 0 | 0 |
| pchembl_value | 0 | 3006 | 0 | 3006 |
| potential_duplicate | 0 | 0 | 0 | 0 |
| qudt_units | 0 | 927 | 0 | 927 |
| record_id | 0 | 0 | 0 | 0 |
| relation | 0 | 350 | 0 | 350 |
| src_id | 0 | 0 | 0 | 0 |
| standard_flag | 0 | 0 | 0 | 0 |
| standard_relation | 0 | 350 | 0 | 350 |
| standard_text_value | 0 | 8115 | 0 | 8115 |
| standard_type | 0 | 0 | 0 | 0 |
| standard_units | 0 | 520 | 0 | 520 |
| standard_upper_value | 0 | 8397 | 0 | 8397 |
| standard_value | 0 | 1003 | 0 | 1003 |
| target_chembl_id | 0 | 0 | 0 | 0 |
| target_organism | 0 | 0 | 0 | 0 |
| target_pref_name | 0 | 0 | 0 | 0 |
| target_tax_id | 0 | 0 | 0 | 0 |
| text_value | 0 | 8115 | 0 | 8115 |
| toid | 0 | 8397 | 0 | 8397 |
| type | 0 | 0 | 0 | 0 |
| units | 0 | 440 | 0 | 440 |
| uo_units | 0 | 927 | 0 | 927 |
| upper_value | 0 | 8395 | 0 | 8395 |
| value | 0 | 1003 | 0 | 1003 |

## Count-only Gate 2 approximation

```json
{
  "cascade": {
    "all_records": 8397,
    "binding_assay": 8174,
    "ic50": 5781,
    "human_target": 5781,
    "molar_units": 5655,
    "exact_relation": 4859,
    "positive_finite_value": 4859
  },
  "filtered_records": 4859,
  "distinct_compounds": 4004,
  "survivors_missing_compound_id": 0,
  "distinct_surviving_assays": 337,
  "distinct_surviving_documents": 258,
  "raw_distributions": {
    "assay_type": {
      "A": 177,
      "B": 8174,
      "F": 46
    },
    "standard_type": {
      "% Ctrl": 1,
      "% Inhibition of Control Values": 1,
      "Activity": 98,
      "Drug metabolism": 18,
      "EC50": 4,
      "IC5": 2,
      "IC50": 5818,
      "INH": 1,
      "Inhibition": 734,
      "K": 62,
      "K inact": 2,
      "K obs / 1": 8,
      "Ka": 1,
      "Kapp": 1,
      "Kd": 26,
      "Kdiss": 1,
      "Ki": 709,
      "Km": 23,
      "Ratio IC50": 12,
      "Ratio Ki": 5,
      "Selectivity index": 68,
      "Stability": 80,
      "T1/2": 40,
      "TIME": 21,
      "TdeltaS": 1,
      "Vmax": 4,
      "Vmax/Km": 4,
      "deltaG": 1,
      "deltaH": 1,
      "k_off": 325,
      "kon": 325
    },
    "standard_units": {
      "%": 880,
      "/min": 8,
      "10'-2/s": 3,
      "10'-3/s": 4,
      "10'-4/s": 12,
      "10'-5/s": 1,
      "10'3/M/s": 1,
      "10'4/M/s": 1,
      "10'5/M/s": 4,
      "10'6/M/min": 6,
      "10'6/M/s": 9,
      "10'7/M/s": 2,
      "10'8/M/min": 1,
      "10^-4/s": 6,
      "10^2kJ/mol": 2,
      "10^3kJ/mol": 1,
      "10^4/M/s": 1,
      "10^5/M": 1,
      "10^5/M/s": 5,
      "<missing>": 520,
      "M-1 s-1": 8,
      "hr": 61,
      "nM": 6474,
      "nmol/hr/microg": 4,
      "s-1": 327,
      "uM": 7,
      "ug.mL-1": 48
    },
    "standard_relation": {
      "<": 82,
      "<=": 1,
      "<missing>": 350,
      "=": 7116,
      ">": 839,
      ">=": 2,
      ">>": 5,
      "~": 2
    },
    "target_organism": {
      "Homo sapiens": 8397
    }
  }
}
```

These counts precede structure checks, salt handling, replicate aggregation and assay-quality review. Literature-reported bioactivity is not a random sample of chemistry. No efficacy, safety or suitability claim follows from these observations.
