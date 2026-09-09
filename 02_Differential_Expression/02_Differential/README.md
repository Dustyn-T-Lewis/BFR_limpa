# 02_Differential_Expression / 02_Differential

Normalises the protein matrix, fits the model, writes the results.

| | |
|---|---|
| **Script** | `a_script/02_differential.qmd` |
| **Reads** | `proteins.rds`, `01_Design/c_data/design.rds` |
| **Writes** | `c_data/protein_contrasts_long.csv`, `contrast_summary.csv`, `fit.rds` |
| **Beside it** | `a_script/02b_normalization_report.qmd`, run when the input matrix changes |

## The sequence

```r
proteins$E <- normalizeBetweenArrays(proteins$E, method = "cyclicloess")
fit <- dpcDE(proteins, design, sample.weights = TRUE, plot = TRUE)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, robust = TRUE)
```

limpa places normalisation here, between `dpcQuant()` and `dpcDE()`, and treats it as optional.
`cyclicloess` is used because most samples disagree with the cohort more at one end of the
abundance range than the other. A per-sample shift cannot fix that. The report beside this script
shows the MA plots the choice rests on, and states its cost.

`dpcDE()` reads the standard errors from `proteins.rds` and turns them into precision weights.
`sample.weights` estimates one quality weight per sample, because biopsies differed in how much
blood they carried. Contrasts are applied before moderation, which is why `eBayes` comes last.

One assertion runs first: that the design rows are the same samples in the same order as the matrix
columns. limma checks only that the counts match.

## Multiple testing

Benjamini-Hochberg within each contrast, never pooled. The five share participants and the
interaction is a difference of two of the others. `topTable()` is called with `sort.by = "p"`
because its default sorts on the B statistic, which orders differently.

## Reading a null contrast

Three numbers are reported per contrast. `propTrueNull()` estimates the null proportion from the
p-value distribution. Median |t| is 0.674 under a calibrated null, and `p_below_0.2` is 0.200.

Compare the between-leg contrasts against `BFR_vs_HLRT_at_T1`, not against those theoretical
values. That contrast is a true null by construction, so whatever it returns is what a null looks
like in this design. It returns median |t| near 0.54, and the interaction returns near 0.52.

## Effect size

`treat()` tests against a fold-change floor rather than against zero and controls an error rate,
unlike filtering on logFC afterwards. Three floors are reported so the counts read as a property of
the threshold rather than of the data.

`pi_score` is `P.Value^abs(logFC)`, the score of Xiao et al. (2014, PMID 22321699) written
inverted. It ranks within a contrast and selects nothing. It rewards large fold changes, and here
those come from the proteins with fewest detected precursors, so the observation count is printed
beside the ranking.

## The control

`BFR_vs_HLRT_at_T1` compares two untrained legs of the same person and should find nothing. Its row
count is asserted, because a name matching no contrast would otherwise look like a clean control.
The hit count is printed rather than asserted: if two untrained legs ever differ, that belongs in
the report, not in a halted render.

Results are written before the check runs, so a failure leaves the evidence on disk.
