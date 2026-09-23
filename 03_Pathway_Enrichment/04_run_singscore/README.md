# 04_run_singscore

Scores every sample on every set. One number per set per sample, no p-value.

| | |
|---|---|
| **Script** | `a_script/04_run_singscore.R` |
| **Reads** | `gene_sets.rds`, `proteins.rds` |
| **Writes** | `c_data/singscore.rds`, `04_run_singscore.xlsx`, `set_scores.csv`, two figures and a two-page figure PDF |

## Why singscore

Rank-based and sample-independent: a sample's score does not move when the cohort changes, which
is what a paired within-participant design needs. Measured here, dropping 51 of 131 samples moves
a singscore value by 0 and a GSVA value by up to 0.34 on a range of 1.5.

It carries no p-value and never sees the contrast. `01_run_fgsea_and_fry` supplies significance; this
supplies the 1,990 x 131 matrix that `05_classify_and_associate_sets` consumes.

## What the structure check says

Participant identity dominates the raw scores — PC1 holds 33.6% of the variance and participant
explains 0.64 of it. That is why the phenotype analysis works on the within-leg change and never
on the raw value.

## Two cohort figures

singscore ships `plotDispersion` and `plotRankDensity`, and both draw one signature at a time.
Reporting here is at cohort scale, so `b_reports/` holds plain plots over the two matrices
`multiScore()` already returns.

`01_score_dispersion.png` puts every set at its mean score against its mean dispersion, coloured
by collection. Dispersion is how far a set's member ranks spread inside one sample: low means the
members sit together, high means they scatter. Reactome scores highest and disperses least.

| Collection | Sets | Median score | Median dispersion |
|---|---:|---:|---:|
| Reactome | 430 | 0.078 | 854 |
| KEGG_Legacy | 94 | 0.025 | 932 |
| Hallmark | 41 | 0.026 | 1010 |
| GO:BP | 1,366 | -0.007 | 1006 |
| GO Slim | 59 | -0.004 | 1045 |

`02_score_distribution.png` is the score distribution per study group over all sets, which is
where a global shift between groups would show if one existed.

## Handoff

```r
ss <- readRDS("03_Pathway_Enrichment/04_run_singscore/c_data/singscore.rds")
ss$scores          # 1,990 sets x 131 samples, sets in rows
ss$score_summary   # range and spread, across sets and across samples
ss$structure_check # variance per component and the participant share
ss$collection_spread # median score and dispersion per collection
```
