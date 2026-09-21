# 03_Pathway_Enrichment / 03_Sample_Scores

Scores every set in every sample with ssGSEA, tests the scores on the same design, and reads the
score changes against muscle growth. The output is a 1,040 × 131 matrix: a pathway turned into a
variable that can be plotted against the phenotype rather than only reported with a p-value.

| | |
|---|---|
| **Script** | `a_script/03_sample_scores.qmd` |
| **Reads** | `01_Gene_Sets/c_data/pathway_inputs.rds`, `00_Input/phenotype.csv` |
| **Writes** | `c_data/set_scores.rds`, `score_contrast_results.csv`, `score_phenotype_correlations.csv` |

```sh
quarto render 03_Pathway_Enrichment/03_Sample_Scores/a_script/03_sample_scores.qmd --output-dir ../b_reports
```

## What the conversion costs

ssGSEA ranks within a sample, so limpa's precision weights do not travel: a protein rebuilt
largely from missing values carries the same rank as a well-measured one. Scores are also relative
to the cohort, so adding or removing samples changes them. Error control for the contrasts
therefore stays with the whole-set tests; these scores answer a different question.

## Calibration

The scores go into the fixed-participant design with `lmFit`, `contrasts.fit` and
`eBayes(trend = TRUE)`, and the negative control is read first.

| Contrast | Sets at FDR < 0.05 | Median p | p < 0.05 |
|---|---|---|---|
| `BFR_vs_HLRT_at_T1` (control) | 0 | 0.510 | 3.2% |
| `BFR_post_vs_pre` | 69 | 0.281 | 21.6% |
| `HLRT_post_vs_pre` | 29 | 0.310 | 18.9% |
| `BFR_vs_HLRT_at_T2` | 0 | 0.467 | 6.2% |
| `interaction` | 0 | 0.537 | 2.4% |

A median of 0.510 on a true null is the closest to flat of any method run in this stage. The
training calls repeat the set tests: oxidative phosphorylation rises in both legs by almost the
same amount (score difference 0.0142 in BFR, 0.0132 in HLRT).

## Against the phenotype

The unit is a leg: one score change (T2 − T1) against one change in vastus lateralis CSA,
correlated within each treatment separately so each person enters each correlation once. Only
sets flagged by the score-level test are correlated.

Nothing survives FDR in either arm (smallest FDR 0.129). HLRT shows a coherent cluster before
correction, all positive: tRNA aminoacylation (ρ = 0.45), amino-acid metabolism (ρ = 0.48),
branched-chain amino-acid metabolism (ρ = 0.37). BFR has nothing comparable (smallest raw p 0.08).
With 33 legs per arm, 80% power needs ρ of about 0.47, so this is underpowered rather than null.

A correlation here relates two changes measured on the same leg. It is not evidence that the
pathway caused the growth.
