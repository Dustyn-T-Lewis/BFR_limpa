# 02 · Differential Expression

This stage takes `proteins.rds` and asks what happened to the proteins. The design sits in its own
step, so a design problem surfaces before anything is fitted.

[`01_Design`](01_Design/README.md) builds the design matrix and the five contrasts and estimates
the within-leg correlation. It writes `design.rds` and `01_design.xlsx`.
[`02_Differential`](02_Differential/README.md) fits `dpcDE()`, applies the contrasts with
`eBayes(robust = TRUE)` and `treat()`, and writes `fit.rds` and `02_differential.xlsx`.
[`03_Phenotype`](03_Phenotype/README.md) tests every protein against the ultrasound outcomes and
writes `03_phenotype.xlsx`.
