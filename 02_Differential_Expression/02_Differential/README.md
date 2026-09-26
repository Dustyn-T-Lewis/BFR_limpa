# 02_Differential_Expression / 02_Differential

Fits the model on the five contrasts and writes the results.

- Reads: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
  `02_Differential_Expression/01_Design/c_data/design.rds`
- Writes: `c_data/fit.rds` (read by `01_run_fgsea_and_fry`), `c_data/02_differential.xlsx`
  (`overview`, `DEP_matrix` with `logFC`, `P.Value` and `adj.P.Val` per protein and contrast,
  `contrast_summary`), `b_reports/02_differential.html`
- Run: `quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports`,
  17 seconds.
