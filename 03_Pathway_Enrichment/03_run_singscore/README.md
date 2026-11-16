# 03 · Run singscore

Scores every sample on every tested set with singscore, which ranks within each sample, so a
score has no p-value and never sees a contrast.

Inputs: `00_build_gene_sets/c_data/gene_sets.rds`, the catalogue and protein map in
`00_build_gene_sets.xlsx`, `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds`.

```sh
quarto render 03_Pathway_Enrichment/03_run_singscore/a_script/03_run_singscore.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/03_run_singscore.html`: the walkthrough.
- `b_reports/03_run_singscore_figures.pdf`: S5, each set's mean score against its mean
  dispersion (the spread of member ranks within a sample) over all samples.
- `c_data/set_scores.rds`: the set-by-sample score matrix, read by `04_` and `05_`. It travels as
  rds because an Excel round trip breaks singscore's exact rank ties.
- `c_data/03_run_singscore.xlsx`, read by no later step:
  - `collection_spread`: one row per collection: median score and median dispersion.
  - `structure_check`: one row per principal component: variance explained, participant share.
  - `set_scores`: the score matrix, one row per set.

Participant dominates the raw scores (PC1 explains 34.0% of variance, participant 0.64 of PC1),
so later steps use within-leg change.
