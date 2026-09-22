# 03_Pathway_Enrichment / 91_Sample_Scores_ssGSEA

Scores every set in every sample with ssGSEA, tests the scores on the same design, and reads the
score changes against muscle growth. The output is a 1,040 × 131 matrix: a pathway turned into a
variable that can be plotted against the phenotype rather than only reported with a p-value.

| | |
|---|---|
| **Script** | `a_script/91_sample_scores_ssgsea.qmd` |
| **Reads** | `90_Set_Test_Calibration/c_data/tested_sets.rds`, `01_Preprocess/02_Quantification/c_data/proteins.rds`, `02_Differential_Expression/01_Design/c_data/design.rds`, `00_Input/phenotype.csv` |
| **Writes** | `c_data/set_scores.rds`, `score_contrast_results.csv`, `score_phenotype_correlations.csv`, `score_phenotype_correlations_hallmark_all.csv`, `score_phenotype_arm_difference_hallmark.csv` |

```sh
quarto render 03_Pathway_Enrichment/91_Sample_Scores_ssGSEA/a_script/91_sample_scores_ssgsea.qmd --output-dir ../b_reports
```

Render `90_Set_Test_Calibration` first: this stage scores exactly the sets that one tested.

## What the conversion costs

ssGSEA ranks within a sample, so limpa's precision weights do not travel: a protein rebuilt
largely from missing values carries the same rank as a well-measured one. Scores are also relative
to the cohort, so adding or removing samples changes them. Error control for the contrasts
therefore stays with the whole-set tests; these scores answer a different question.

The main pathway stages score samples with singscore instead, which is sample-independent. The two
are alternatives for the same job; `92_Method_Comparison` lists where they are compared.

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

A median of 0.510 on a true null is the closest to flat of any method run here. The training
calls repeat the set tests: oxidative phosphorylation rises in both legs by almost the same amount
(score difference 0.0142 in BFR, 0.0132 in HLRT).

## Against the phenotype

The unit is a leg: one score change (T2 − T1) against one change in vastus lateralis CSA,
correlated within each treatment separately so each person enters each correlation once.

**Filtered (85 sets that moved with training).** Nothing survives FDR in either arm (smallest FDR
0.129). HLRT shows a coherent cluster before correction, all positive: tRNA aminoacylation
(ρ = 0.45), amino-acid metabolism (ρ = 0.48), branched-chain amino-acid metabolism (ρ = 0.37). BFR
has nothing comparable (smallest raw p 0.08).

**Unfiltered (all 41 Hallmark sets).** The filter selects on the mean change, and a set can track
growth without moving on average, so the same correlation is also run on every Hallmark set.
Nothing survives FDR (smallest FDR 0.59 in HLRT, 0.97 in BFR); raw p < 0.05 in 3 HLRT and 1 BFR
sets, against 2.1 expected by chance.

**Between arms.** Being significant in one arm and not the other is not a difference, so the two
correlations are compared directly with Fisher's z. No Hallmark set differs between arms after
correction (smallest FDR 0.85).

With 33 legs per arm, 80% power needs ρ of about 0.47, so this is underpowered rather than null.
`00_Input/phenotype.csv` is being corrected for eight legs in an open pull request; this stage
reads that file directly, so a rerun after the merge picks the correction up.

A correlation here relates two changes measured on the same leg. It is not evidence that the
pathway caused the growth.
