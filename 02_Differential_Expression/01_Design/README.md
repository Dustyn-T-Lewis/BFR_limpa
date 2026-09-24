# 02_Differential_Expression / 01_Design

Builds the design matrix and the five contrasts, and measures the within-leg correlation they
assume is small.

- Reads: `01_Preprocess/02_Quantification/c_data/proteins.rds`
- Writes: `c_data/design.rds` (read by `02_Differential` and `01_run_fgsea_and_fry`),
  `c_data/01_design.xlsx` (`overview`, `correlation_strata`), `b_reports/01_design.html`
- Run: `quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd --output-dir ../b_reports`,
  a few seconds.
