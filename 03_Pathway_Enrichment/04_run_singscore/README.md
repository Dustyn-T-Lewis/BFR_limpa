# 03_Pathway_Enrichment / 04_run_singscore

Scores every sample on every set with singscore.

- Reads: `00_build_gene_sets/c_data/gene_sets.rds`, `00_build_gene_sets.xlsx` (`set_catalog`,
  `protein_gene_map`), `01_Preprocess/02_Quantification/c_data/proteins.rds`
- Writes: `c_data/set_scores.rds` (read by `05_classify_and_associate_sets`),
  `c_data/04_run_singscore.xlsx` (`overview`, `collection_spread`, `structure_check`,
  `set_scores`), `b_reports/04_run_singscore_figures.pdf` (S12)
- Run: `Rscript 03_Pathway_Enrichment/04_run_singscore/a_script/04_run_singscore.R`, 11
  seconds.

## A sample's score does not depend on the cohort

singscore ranks proteins within each sample. It carries no p-value and never sees a contrast.
Measured here, dropping 51 of 131 samples changed singscore values by 0 and GSVA values by up to
0.34, and the paired design needs that stability.

Participant dominates raw scores (PC1 33.6%, of which participant explains 0.64), so later steps
use within-leg change.

S12: (A) mean score against mean dispersion per set, coloured by collection; (B) score
distribution by study group. Reactome scores highest and disperses least.

## The matrix reaches 05 as rds, not through the workbook

The score matrix is written twice: as a sheet for reading, and as `set_scores.rds`, which
`05_classify_and_associate_sets` reads. singscore values carry exact rank ties, and a trip through
Excel changes the last bit of some of them, which breaks the ties and moves 05's Wilcoxon and
Spearman p-values by up to 10%.
