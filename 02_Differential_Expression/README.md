# 02 · Differential Expression

Takes `proteins.rds` and asks what happened to the proteins. The design sits in its own step, so
a design problem surfaces before anything is fitted.

| Step | Runs | Writes |
|---|---|---|
| [`01_Design`](01_Design/README.md) | design matrix, five contrasts, within-leg correlation | `design.rds`, `01_design.xlsx` |
| [`02_Differential`](02_Differential/README.md) | `dpcDE()`, contrasts, `eBayes(robust = TRUE)`, `treat()` | `fit.rds`, `02_differential.xlsx` |
| [`03_Phenotype`](03_Phenotype/README.md) | every protein against the ultrasound outcomes | `03_phenotype.xlsx` |
