# 02_Differential_Expression / 01_Design

Builds the design matrix and the five contrasts, and measures the assumption they rest on.
Separate from the fit so a design problem surfaces in seconds.

| | |
|---|---|
| **Script** | `a_script/01_design.qmd` |
| **Reads** | `01_Preprocess/02_Quantification/c_data/proteins.rds` |
| **Writes** | `c_data/design.rds`, `c_data/correlation_strata.csv` |

## The design

```r
model.matrix(~ 0 + group + participant, data = targets)
```

131 rows by 36 columns: four group means, one per treatment-and-timepoint cell, plus 32 participant
dummies. Full rank, 95 residual degrees of freedom.

Participant is fixed, not random, because every comparison happens inside one person. With the term
in the design, treatments are only compared within the same person, which is what makes pre-to-post
paired. A random participant effect is for designs that also compare between people.

One participant contributes three cells rather than four, so the group columns are
participant-adjusted cell means rather than raw ones. The design is full rank regardless.

## Two assertions before anything is fitted

The design must be full rank, or the contrasts are not estimable. And every column name must survive
`make.names()`, because `makeContrasts()` parses its arguments as R code and a level named like
`2E-T1` would be read as subtraction.

## The correlation diagnostic

Samples share a participant, and within that they share a leg. The participant term removes the
first exactly; the second stays in the residual. If it were large, the within-leg comparisons would
be tested too conservatively and the between-leg ones too liberally.

`correlation_strata.csv` reports both. On this data the within-leg value is 0.012, so the fixed term
is sufficient and `02_Differential` calls `dpcDE()` without `block =`.

Two caveats. Both correlations are estimated on the bare expression matrix, so neither carries the
precision weights or sample weights the real fit uses; they describe the matrix, not the model. And
nothing automates the decision. If that number ever rose, the fix is `block = leg_id` in the next
notebook, but there is no branch, `leg_id` is not saved, and nothing would warn you. Read the row.
