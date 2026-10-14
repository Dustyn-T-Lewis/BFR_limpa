# 04 · Run singscore

Scores every sample on every tested set with singscore, which ranks within each sample, so a
score has no p-value and never sees a contrast.

Inputs: `00_build_gene_sets/c_data/gene_sets.rds`, the catalogue and protein map in
`00_build_gene_sets.xlsx`, `01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 03_Pathway_Enrichment/04_run_singscore/a_script/04_run_singscore.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/04_run_singscore.html`: the walkthrough.
- `b_reports/04_run_singscore_figures.pdf`: S11, mean score against mean dispersion per set, and
  the score distribution per study group.
- `c_data/set_scores.rds`: the set-by-sample score matrix, read by `05_set_classification` and
  `06_set_association`. It travels as rds because an Excel round trip breaks singscore's exact
  rank ties.
- `c_data/04_run_singscore.xlsx`: the same matrix for reading, per-collection medians and the
  principal-component check.

Participant dominates the raw scores (PC1 explains 33.6% of variance, participant 0.64 of PC1),
so later steps use within-leg change.
