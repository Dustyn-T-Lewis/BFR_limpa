# 02_Differential_Expression / 02_Differential

**Planned.** This folder holds this README.

Fits the model and writes the results: one row per protein per comparison, a summary count per
comparison, and the fitted model.

## The sequence

Fit with limpa's function, apply the contrasts, then moderate the variance estimates. limpa does
the first step, limma the rest.

The fit reads the standard errors and turns them into weights, so a protein built mostly from
missing peptides counts for less. One quality weight per sample is estimated too, because
biopsies differed in how much blood they carried.

## Multiple testing

Adjusted within each comparison, never pooled across the five.

A secondary ranking score may be reported alongside, for ordering supplementary tables only. It
controls no error rate, so any count quoted against it must name its own threshold.

## The control

The before-training comparison between legs should find nothing. The stage asserts that, and
writes the tables to disk first, so a failure leaves the evidence rather than losing it.

If two untrained legs of one person differ, the cause is upstream. That failure is the finding.
Do not relax the check to get a clean render.
