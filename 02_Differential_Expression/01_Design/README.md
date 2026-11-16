# 01 · Design

Builds the design matrix (treatment-by-timepoint groups, a fixed effect per participant and the
blood index as a covariate) and the five contrasts, and checks the residual within-leg correlation.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_design.html`: the report.
- `c_data/design.rds`: design, contrasts and the baseline contrast's name, for later steps.
- `c_data/01_design.xlsx`:
  - `samples`: sample identifiers and study metadata.
  - `design`: the design matrix, one row per sample.
  - `contrasts`: contrast weights, one column per comparison.
  - `correlation_strata`: the within-leg residual correlation, participant fixed.

Results: 131 samples, 37 coefficients, 94 residual degrees of freedom. The within-leg residual
correlation is 0.015, so `02_Differential` fits without a leg block. The estimate ignores the
standard errors the fit uses.
