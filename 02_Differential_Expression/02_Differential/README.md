# 02_Differential_Expression / 02_Differential

Fits the model and writes the results.

| | |
|---|---|
| **Script** | `a_script/02_differential.qmd` |
| **Reads** | `proteins.rds`, `01_Design/c_data/design.rds` |
| **Writes** | `c_data/fit.rds`, `c_data/02_differential.xlsx` |

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

## The sequence

```r
fit <- dpcDE(proteins, design, sample.weights = TRUE, plot = TRUE)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit, robust = TRUE)
```

`proteins.rds` arrives normalised. `dpcDE()` reads its standard errors and turns them into
precision weights; `lmFit()` would discard them. `sample.weights` estimates one quality weight per
sample, because biopsies differed in how much blood they carried. Contrasts are applied before
moderation, which is why `eBayes` comes last.

The weights are summarised per cell. They downweight, they do not exclude, and what matters is
that the downweighting falls evenly, since uneven weighting would act like an undeclared covariate.

One assertion runs first: that the design rows are the same samples in the same order as the matrix
columns. limma checks only that the counts match.

## What comes out

`02_differential.xlsx` holds two sheets. `DEP_matrix` is one row per protein with `logFC`,
`P.Value` and `adj.P.Val` for each of the five contrasts. `contrast_summary` counts hits per
contrast. `fit.rds` is the complete fitted model, so the long per-contrast form is not written.

## Multiple testing

Benjamini-Hochberg within each contrast, never pooled. The five share participants and the
interaction is a difference of two of the others. `topTable()` is called with `sort.by = "none"`
so all five contrasts come back in the protein order of the matrix and line up row for row;
its default would sort on the B statistic.

## Reading a null contrast

Three numbers are reported per contrast: `propTrueNull()`, median |t|, and the share of p-values
below 0.2. Under a calibrated null those sit at 1, 0.674 and 0.200. `propTrueNull()` saturates at
its ceiling, which is why the other two are reported beside it. All three read the same p-value
histogram, which is plotted next to them: a spike at zero is signal, a slope toward 1 is a
conservative test.

Compare the between-leg contrasts against `BFR_vs_HLRT_at_T1`, not against the theoretical values.
That contrast is a true null by construction, so whatever it returns is what a null looks like in
this design, and it has run conservative here. The honest reading of an empty interaction is that
no between-modality difference is detectable under a test that is under-calling, which is weaker
than saying the two training modes produce equivalent proteomes.

The notebook also prints the effect size the interaction had 80% power to detect. An empty
contrast says nothing about effects below that size, so quote the two together.

## Effect size

`treat()` tests against a fold-change floor rather than against zero and controls an error rate,
unlike filtering on logFC afterwards. One floor is fixed in advance at 1.15-fold, so the count
reads as a property of that threshold rather than of the data.

`pi_score` is `P.Value^abs(logFC)`, the score of Xiao et al. (2014, PMID 22321699) written
inverted. It ranks within a contrast and selects nothing. It rewards large fold changes, and here
those come from the proteins with fewest detected precursors, so the observation count is printed
beside the ranking.

It is applied only to the two training contrasts. Ranking the interaction would print a table of
its largest fold changes even though `propTrueNull` sits at its ceiling there, and such a table
gets read as a result list whatever the caption says.

## The control

`BFR_vs_HLRT_at_T1` compares two untrained legs of the same person and should find nothing. Its row
count is asserted, because a name matching no contrast would otherwise look like a clean control.
The hit count is printed rather than asserted: if two untrained legs ever differ, that belongs in
the report, not in a halted render. Results are written before the check runs, so a failure leaves
the evidence on disk.
