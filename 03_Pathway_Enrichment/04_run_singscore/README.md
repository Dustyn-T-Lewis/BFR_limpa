# 04_run_singscore

This step scores every sample on every set with singscore, from
`00_build_gene_sets/c_data/gene_sets.rds`, the catalogue and protein map in
`00_build_gene_sets.xlsx`, and `01_Preprocess/02_Quantification/c_data/proteins.rds`.

The score matrix is written twice: as a sheet in `c_data/04_run_singscore.xlsx`, beside
per-collection medians and the principal-component check, and as `c_data/set_scores.rds`, which
`05_set_classification` and `06_set_association` read. singscore values carry exact rank ties, and
Excel changes the last bit of some, which breaks ties and moves 05's Wilcoxon and Spearman p by up
to 10%. S11 is in `b_reports/04_run_singscore_figures.pdf`.

```sh
Rscript 03_Pathway_Enrichment/04_run_singscore/a_script/04_run_singscore.R
```

singscore ranks within each sample, so a score has no p-value and never sees a contrast. Dropping 51
of the 131 samples changed singscore values by 0 and GSVA values by up to 0.34, and the paired
design needs that stability.

Participant dominates the raw scores (PC1 explains 33.6% of variance, participant 0.64 of PC1), so
later steps use within-leg change.

S11 plots each set's mean score against mean dispersion by collection, and the score distribution
per study group. Reactome scores highest and disperses least.
