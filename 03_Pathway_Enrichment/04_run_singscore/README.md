# 04_run_singscore

This step scores every sample on every set with singscore. It reads the sets from
`00_build_gene_sets/c_data/gene_sets.rds`, the catalogue and protein map from
`00_build_gene_sets.xlsx`, and the protein table from
`01_Preprocess/02_Quantification/c_data/proteins.rds`.

The score matrix is written twice. `c_data/04_run_singscore.xlsx` has it as a sheet for reading,
with the median score and dispersion per collection and the principal-component check below.
`c_data/set_scores.rds` is the copy `05_classify_and_associate_sets` reads. singscore values carry
exact rank ties, and a trip through Excel changes the last bit of some of them, which breaks the
ties and moves 05's Wilcoxon and Spearman p-values by up to 10%. Figure S12 is in
`b_reports/04_run_singscore_figures.pdf`.

```sh
Rscript 03_Pathway_Enrichment/04_run_singscore/a_script/04_run_singscore.R
```

singscore ranks proteins within each sample, so a score carries no p-value and never sees a
contrast. Dropping 51 of the 131 samples changed singscore values by 0 and GSVA values by up to
0.34, and the paired design needs that stability.

Participant dominates the raw scores. The first principal component explains 33.6% of the
variance, and participant explains 0.64 of that component, so later steps use the within-leg
change.

S12 plots each set's mean score against its mean dispersion, coloured by collection, and the score
distribution in each study group. Reactome scores highest and disperses least.
