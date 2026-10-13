# 02 · Differential expression

Fits the protein estimates with `dpcDE()` and tests the five contrasts, with sample quality weights,
robust empirical Bayes and BH adjustment within each contrast. When the fitted-slope checkpoint
exists, it refits on it as a sensitivity check.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`, `01_Design/c_data/design.rds`.

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_differential.html`: the walkthrough, with diagnostic plots.
- `c_data/fit.rds`: the fitted contrasts, read by stage 03.
- `c_data/02_differential.xlsx`: per-protein results (`DEP_matrix`), discovery counts, sample
  weights, p-value summaries and fold-change-threshold counts.
