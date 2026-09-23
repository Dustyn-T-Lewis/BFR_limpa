# 04_run_singscore

One singscore per set per sample.

| | |
|---|---|
| Reads | `gene_sets.rds`, `proteins.rds` |
| Writes | `singscore.rds`, `set_scores.csv`, `04_run_singscore.xlsx`, 2 figures |

singscore ranks proteins within each sample, so a score does not depend on the cohort. Measured
here, dropping 51 of 131 samples changed singscore values by 0 and GSVA values by up to 0.34.

Participant dominates raw scores (PC1 33.6%, of which participant explains 0.64), so later steps
use within-leg change.

Figures: mean score against mean dispersion per set, coloured by collection; score distribution by
study group. Reactome scores highest and disperses least.

```r
ss <- readRDS("03_Pathway_Enrichment/04_run_singscore/c_data/singscore.rds")
ss$scores            # 1,990 sets x 131 samples
ss$collection_spread # median score and dispersion per collection
ss$structure_check   # variance per component, participant share
```
