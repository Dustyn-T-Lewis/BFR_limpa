# 01 · Design

Builds the design matrix (treatment-by-timepoint groups plus a fixed effect per participant) and
the five contrasts, and checks the residual within-leg correlation.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_design.html`: the walkthrough.
- `c_data/design.rds`: design, contrasts and the baseline contrast's name, read by
  `02_Differential` and later stages.
- `c_data/01_design.xlsx`: overview, samples, design, contrasts and the correlation check.
