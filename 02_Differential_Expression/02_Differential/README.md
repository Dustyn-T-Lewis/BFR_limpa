# 02_Differential_Expression / 02_Differential

**Planned.** This folder holds this README.

Fits the model and writes the results.

| | |
|---|---|
| **Will read** | the protein table, and the design from `01_Design` |
| **Will write** | one row per protein per comparison, a summary count per comparison, and the fitted model |
| **Read by** | `03_Pathway_Enrichment`, `04_Network`, `05_Figures` |

## The sequence

Fit with limpa's own function, apply the contrasts, then moderate the variance estimates. limpa
does the first step and limma does the rest.

The fit reads the standard errors from the protein table and turns them into weights, so a protein
built mostly from missing peptides counts for less than one measured directly. Nothing else has
access to that information. Handing the abundances to an ordinary linear model would treat every
protein as equally reliable.

One quality weight per sample is estimated as well. Biopsies differed in how much blood they
carried, and removing the contaminant rows afterwards does not undo what that did to the sample
underneath.

## Multiple testing

Adjusted within each comparison, never pooled across the five. They share participants and the
interaction is a difference of two of the others, so a pooled adjustment would not control
anything that can be stated.

A secondary ranking score may be reported alongside, as an ordering for supplementary tables only.
It controls no error rate, so any count quoted against it has to name its own threshold and cannot
borrow the language of false discovery rates.

## The control

The before-training comparison between legs should find nothing. The stage will assert that, and
write the result tables to disk before the assertion runs, so a failure leaves the evidence on
disk instead of losing it with the render.

If two untrained legs of one person differ, the cause is upstream: annotation, contamination, or a
sample mix-up. That failure is the finding. Do not relax the check to get a clean render.
