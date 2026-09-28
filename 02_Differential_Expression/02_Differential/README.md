# 02_Differential_Expression / 02_Differential

This step fits the model to `01_Preprocess/02_Quantification/c_data/proteins.rds` with the design
from `02_Differential_Expression/01_Design/c_data/design.rds`, and tests the five contrasts.

`01_run_fgsea_and_fry` reads the fit, `c_data/fit.rds`. `c_data/02_differential.xlsx` holds `logFC`,
`P.Value` and `adj.P.Val` per protein and contrast (`DEP_matrix`) and hits per contrast at FDR 0.05
and 0.10; `b_reports/02_differential.html` is the report.

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

The notebook ends by refitting on the fitted-slope run (0.488): 82 BFR and 114 HLRT training hits
against 60 and 107 at the preset 0.7, with 57 and 99 in both. The other three contrasts have none
under either slope.