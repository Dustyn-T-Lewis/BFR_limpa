# 04_singscore_pheno_associations

Asks whether a pathway's movement tracks the phenotype, and whether it tracks it differently under
BFR than under high load.

| | |
|---|---|
| **Script** | `a_script/04_singscore_pheno_associations.R` |
| **Reads** | `singscore.rds`, `gene_sets.rds`, `proteins.rds`, `phenotype.csv`, `metadata.csv` |
| **Writes** | `c_data/pheno_associations.rds`, `04_singscore_pheno_associations.xlsx`, four figures |

## The design it leans on

Every participant contributed one BFR leg and one HLRT leg — 33 of 33, asserted at the top of the
script. So the differential question is asked **within participant**: for each person, the
BFR-minus-HLRT difference in score change against the same difference in outcome. One number per
participant, 32 complete pairs.

That removes every between-participant confound and roughly doubles the power against comparing
two independent groups of legs: at a true r of 0.4, 0.63 against 0.37.

**S06's left leg has T1 only** — never acquired, not a reinjection — and that leg is BFR, which is
why one arm contributes 32 legs and the other 33. The script asserts this rather than discovering
it again: 131 samples, 33 participants, 66 legs, treatment agreeing 66/66 between `metadata.csv`
and `phenotype.csv`, zero MS legs absent from the phenotype.

## Three tables, increasing strength of claim

| Sheet | What it asks | n |
|---|---|---:|
| `pooled` | does this pathway track adaptation at all, ignoring treatment | 65 legs |
| `by_treatment` | the same correlation inside each arm, side by side | 32 / 33 legs |
| `differential` | **the study question**, within participant | 32 pairs |

`by_treatment` is descriptive. Two separate confidence intervals are not a test of their
difference, which is what `differential` exists to provide.

Four outcomes, all leg-level change scores: `vl_csa`, `vl_echo`, `rf_csa`, `rf_echo`. Vastus
lateralis is the biopsied muscle and the primary; rectus femoris was not biopsied, so a set
tracking VL but not RF is more convincing than one tracking both.

## Read signal_check before any row

2,732 sets x 4 outcomes is 10,928 tests. Nominal p is the screening statistic and BH is reported
beside it, not used as a gate. `signal_check` gives, at each threshold, the observed hit count
against the count expected by chance.

**As measured: at p < 0.01 the differential analysis returns 12 to 25 hits against 27 expected,
every ratio below 1.0.** The screen is indistinguishable from noise. No row in `differential` is a
finding, and the top of that table is what a null table looks like.

With 32 pairs the study can resolve a paired correlation of about 0.5. This is a null result at
that resolution, not evidence that no difference exists.
