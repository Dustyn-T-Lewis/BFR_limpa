# 02_Differential_Expression / 02_Differential

This step fits the model to `01_Preprocess/02_Quantification/c_data/proteins.rds` with the design
from `02_Differential_Expression/01_Design/c_data/design.rds`, and tests the five contrasts.

The fit is saved as `c_data/fit.rds`, which `01_run_fgsea_and_fry` reads.
`c_data/02_differential.xlsx` has `logFC`, `P.Value` and `adj.P.Val` for every protein and
contrast in its `DEP_matrix` sheet, and the hits per contrast at FDR 0.05 and 0.10. The rendered
notebook is `b_reports/02_differential.html`.

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

The notebook ends by refitting the model on the fitted-slope run (0.488). BFR training has 82 hits
there against 60 at the preset 0.7, and HLRT training 114 against 107, with 57 and 99 found under
both. The other three contrasts have none under either slope.
