# 02_Differential_Expression / 01_Design

This step builds the design matrix and the five contrasts from the sample targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds`. It also measures the within-leg
correlation, which the design assumes is small.

The design and contrasts are saved as `c_data/design.rds`, which `02_Differential` and
`01_run_fgsea_and_fry` read. `c_data/01_design.xlsx` holds the correlation estimate, and the
rendered notebook is `b_reports/01_design.html`.

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd --output-dir ../b_reports
```
