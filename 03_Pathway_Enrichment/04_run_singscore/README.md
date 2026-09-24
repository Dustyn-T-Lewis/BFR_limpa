# 04_run_singscore

One singscore per set per sample.

| | |
|---|---|
| Reads | `gene_sets.rds`, `00_build_gene_sets.xlsx`, `proteins.rds` |
| Writes | `set_scores.rds`, `04_run_singscore.xlsx`, `04_run_singscore_figures.pdf` (S12) |

singscore ranks proteins within each sample, so a score does not depend on the cohort. Measured
here, dropping 51 of 131 samples changed singscore values by 0 and GSVA values by up to 0.34.

Participant dominates raw scores (PC1 33.6%, of which participant explains 0.64), so later steps
use within-leg change.

S12: (A) mean score against mean dispersion per set, coloured by collection; (B) score
distribution by study group. Reactome scores highest and disperses least.

The score matrix is written twice: as a sheet for reading, and as `set_scores.rds`, which
`05_classify_and_associate_sets` reads. singscore values carry exact rank ties, and a trip
through Excel changes the last bit of some of them, which breaks the ties and moves 05's
Wilcoxon and Spearman p-values by up to 10%.
