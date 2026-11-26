# Curation cascade: CHEMBL284

Release ChEMBL_37; raw pull of **8397** records, checksums verified before and after. Command: `.venv/bin/python -m src.curate CHEMBL284`.

Each filter is applied separately, in this order, to the output of the previous one.
The last step changes the unit of the count from records to compounds.

| # | Filter | Rule | In | Out | Removed |
|---|---|---|---:|---:|---:|
| 1 | binding assay | `assay_type == 'B'` | 8397 | 8174 | 223 |
| 2 | single activity type | `standard_type == 'IC50'` | 8174 | 5781 | 2393 |
| 3 | single target organism | `target_organism == 'Homo sapiens' and str(target_tax_id) == '9606'` | 5781 | 5781 | 0 |
| 4 | units harmonised to nM | `standard_units in ['M', 'fM', 'mM', 'nM', 'pM', 'uM', 'µM', 'μM'] and standard_value is a finite number > 0` | 5781 | 5650 | 131 |
| 5 | censored removed | `standard_relation == '=' exactly` | 5650 | 4859 | 791 |
| 6 | structure standardised | `RDKit sanitise -> largest fragment -> neutralise -> canonical SMILES -> InChIKey` | 4859 | 4855 | 4 |
| 7 | duplicates resolved | `median pIC50 per standardised InChIKey; spread > 1.0 log unit flagged` | 4855 | 3965 | 890 |

## What each filter removed

**1_binding_assay** - removed 223:

- `assay_type='A'`: 177
- `assay_type='F'`: 46

**2_single_activity_type** - removed 2393:

- `standard_type='% Ctrl'`: 1
- `standard_type='% Inhibition of Control Values'`: 1
- `standard_type='Activity'`: 48
- `standard_type='EC50'`: 4
- `standard_type='IC5'`: 2
- `standard_type='INH'`: 1
- `standard_type='Inhibition'`: 723
- `standard_type='K inact'`: 2
- `standard_type='K obs / 1'`: 8
- `standard_type='K'`: 62
- `standard_type='Ka'`: 1
- `standard_type='Kapp'`: 1
- `standard_type='Kd'`: 26
- `standard_type='Kdiss'`: 1
- `standard_type='Ki'`: 709
- `standard_type='Km'`: 23
- `standard_type='Ratio IC50'`: 12
- `standard_type='Ratio Ki'`: 5
- `standard_type='Selectivity index'`: 68
- `standard_type='T1/2'`: 14
- `standard_type='TIME'`: 20
- `standard_type='TdeltaS'`: 1
- `standard_type='Vmax'`: 4
- `standard_type='Vmax/Km'`: 4
- `standard_type='deltaG'`: 1
- `standard_type='deltaH'`: 1
- `standard_type='k_off'`: 325
- `standard_type='kon'`: 325

**3_single_target_organism** - removed 0:


> The rule is the TARGET organism. The assay expression host is absent from every activity record in this pull, so it is not inferred.

**4_units_harmonised_to_nM** - removed 131:

- `unit_not_convertible=<missing>`: 78
- `unit_not_convertible=ug.mL-1`: 48
- `value_not_a_number`: 5

> Every refusal is itemised in units_refused.csv. Mass concentrations are refused because the molecular weight of the tested material (possibly a salt) is not known here; converting them would mix two measurement bases.

**5_censored_removed** - removed 791:

- `relation='<'`: 54
- `relation='>'`: 732
- `relation='>='`: 1
- `relation='>>'`: 4

> Censored means the measurement is a bound, not a value: a '>' row only says the IC50 exceeds the highest concentration tested. Keeping them as if they were values would put invented numbers into the target variable. Rows with no relation at all are removed by the same rule but counted separately.

**6_structure_standardised** - removed 4:

- `structure_missing_smiles`: 4

> Every failure is itemised in structure_failures.csv; no molecule is dropped silently.

**7_duplicates_resolved** - removed 890:

- `rows_collapsed_into_a_median`: 890

> The output unit changes here: rows become compounds. The median is used rather than the mean because replicate IC50 sets in this pull have a long tail of cross-assay disagreement, which a mean would follow.

## Reconciling step 5 with the Gate 1 censored count

Gate 1 measured 796 censored rows; step 5 of this cascade removes 791. Why?

- Gate 1's base (unit string only, no value check), recomputed here: **5655** rows
- rows failing the exact-relation rule on that base: **796** (the Gate 1 figure)
- removed by step 5 of this cascade: **791**, of which 791 carry a censoring relation
- difference: **5** rows (activity_id 744390, 1603666, 1603715, 1603720, 1603754)

> Rows whose standard_relation AND standard_value are both null, with units reported as nM. Gate 1 counted them as 'not exactly =' and so inside its 796. This cascade removes them one step earlier, at step 4, as value_not_a_number, because a row with no value cannot become a pIC50 whatever its relation says.

No censored measurement is treated differently by the two gates. The 791/796 gap is five value-less rows removed one step earlier here, and the 4859 rows that survive the relation rule are identical in both.

No record is removed without appearing in a breakdown above or in an itemised CSV (`units_refused.csv`, `structure_failures.csv`). These are literature- and patent-reported measurements; nothing here claims efficacy, safety or suitability.
