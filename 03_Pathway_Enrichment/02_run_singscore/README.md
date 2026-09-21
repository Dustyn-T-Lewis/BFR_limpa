# 02_run_singscore

Scores every sample on every set. One number per set per sample, no p-value.

| | |
|---|---|
| **Script** | `a_script/02_run_singscore.R` |
| **Reads** | `gene_sets.rds`, `proteins.rds`, `config.yml` |
| **Writes** | `c_data/singscore.rds`, `02_run_singscore.xlsx`, `set_scores.csv` |

## Why singscore

Rank-based and sample-independent: a sample's score does not move when the cohort changes, which
is what a paired within-participant design needs. Measured here, dropping 51 of 131 samples moves
a singscore value by 0 and a GSVA value by up to 0.34 on a range of 1.5.

It carries no p-value and never sees the contrast. `01_run_fgsea` supplies significance; this
supplies the 2,732 x 131 matrix that `04_singscore_pheno_associations` consumes.

## What the structure check says

Participant identity dominates the raw scores — PC1 holds 33.7% of the variance and participant
explains 0.65 of it. That is why the phenotype analysis works on the within-leg change and never
on the raw value.

## Handoff

```r
ss <- readRDS("03_Pathway_Enrichment/02_run_singscore/c_data/singscore.rds")
ss$scores          # 2,732 sets x 131 samples, sets in rows
ss$score_summary   # range and spread, across sets and across samples
ss$structure_check # variance per component and the participant share
```
