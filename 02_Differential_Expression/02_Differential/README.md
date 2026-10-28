# 02 · Differential expression

Fits the protein estimates with `dpcDE()` and tests the five contrasts with
`eBayes(robust = TRUE)` and BH adjustment within each contrast. `dpcDE()` uses the standard errors
as precision weights and, as in limpa's default, adds no sample quality weights. When the
fitted-slope quantification run exists, the notebook refits on it as a sensitivity check: 85 and
177 training hits against the primary 86 and 161 (76 and 151 shared), and none between legs under
either slope.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`, `01_Design/c_data/design.rds`.

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_differential.html`: the walkthrough, with diagnostic plots.
- `c_data/fit.rds`: the fitted contrasts, read by `03_Pathway_Enrichment/01_run_fgsea_and_fry`
  and `05_Mitochondria/02_test_mito`.
- `c_data/02_differential.xlsx`:
  - `DEP_matrix`: one row per protein, with log2 fold change, p and BH-adjusted p per contrast;
  - `contrast_summary`: discoveries per contrast at FDR 0.05 and 0.10;
  - `null_summary`: per contrast, estimated share of true nulls, median |t| and share of p below
    0.2.
