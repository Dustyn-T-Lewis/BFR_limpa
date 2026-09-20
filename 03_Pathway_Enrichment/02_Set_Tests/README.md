# 03_Pathway_Enrichment / 02_Set_Tests

Tests whether a set moved, and scores where each sample sits on it. This is the stage whose
p-values you quote.

| | |
|---|---|
| **Script** | `a_script/02_set_tests.qmd` |
| **Reads** | `proteins.rds`, `fit.rds`, `design.rds`, `01_Gene_Sets/c_data/pathway_inputs.rds`, `00_Input/phenotype.csv` |
| **Writes** | `c_data/set_tests.rds`, `02_set_tests.xlsx`, `fry_results.csv` |

```sh
quarto render 03_Pathway_Enrichment/02_Set_Tests/a_script/02_set_tests.qmd --output-dir ../b_reports
```

## Two tools, two questions

`fry()` asks whether a set moved on a contrast. It is a rotation test, so members that correlate
by construction — mitochondria, ribosome — do not inflate it, and it takes the participant design
and limpa's precision weights from `fit$EList$weights`. One call per contrast: it accepts a single
contrast vector, never a matrix.

`singscore` asks where each sample sits on each set, giving 1,032 by 131 numbers. It carries no
p-value. `simpleScore()` returns a score and a dispersion; the permutation route handles one set
per call at a thousand re-scorings each, and answers a per-sample question that never sees the
contrast. The scores are the output worth having, and fry supplies the significance.

## What it found

fry calls 265 sets on `BFR_post_vs_pre` and 226 on `HLRT_post_vs_pre`, and **nothing** on the
interaction, on `BFR_vs_HLRT_at_T2`, or on the negative control. Both legs trained; neither
differs from the other. That empty interaction agrees with the protein-level result and is bounded
by the same detectable-effect calculation.

Read it against `01_Gene_Sets`, where fgsea called one set significant on the negative control and
59 on `BFR_vs_HLRT_at_T2`. Same sets, same matrix, different assumption about whether proteins are
exchangeable.

## Whether fry needs a block term

Participant is already fixed in the design, so what remains is the correlation between the two
samples of one leg. Measured on the matrix fry actually tests, it is 0.027, and fry runs without
`block`.

The set-score matrix gives 0.046 for the same quantity. Both are near zero and they are not the
same number, which is the whole point: a correlation belongs to the matrix it was estimated on.
This project once handed fry the singscore consensus instead of the protein-level one and 29 sets
crossed the FDR line on the strength of it. Nothing here branches on either value, so read the row.

## The signature, and why its control does not pass

A set built from a contrast and tested on that same contrast is circular and always significant.
So the signature is defined on `HLRT_post_vs_pre` — 107 proteins, 61 up and 46 down at BH < 0.05 —
and tested on `BFR_post_vs_pre`, which asks whether restriction reproduced the high-load response.

It did, emphatically: both halves come back at FDR below 0.001. Restriction moves the same
proteins the same way.

**The negative control is not clean, so that result cannot be quoted yet.** Tested on
`BFR_vs_HLRT_at_T1` — two untrained legs of one person — the up half returns p = 0.030, FDR =
0.060. It does not cross 0.05, but it is not the nothing a true null should give, and the FDR
there is adjusted over two tests, which is barely an adjustment at all.

Two readings, and this stage cannot separate them. The signature selects proteins for large
effects, and those may be the proteins with the widest baseline spread between legs for reasons
that have nothing to do with training. Or there is a real pre-training difference between legs in
exactly those proteins, which would be an upstream problem. Either way the interaction result
that follows — up half at FDR 0.031, direction reversed — rests on a control that is leaking, and
should not be read as the study's answer.

The fixed set results above do not have this problem: fry returns nothing on the control across
all 1,032 sets. The difference is that those sets were defined before any contrast was seen.

## Reading the phenotype table

`score_vs_csa` ranks sets by how their change across the two timepoints tracks the change in
cross-sectional area, one row per leg. It carries no multiplicity correction and implies none:
it ranks 1,032 sets against one outcome, and any single row would need its own test.
