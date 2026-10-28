# 02 · Differential Expression

This stage fits the protein model to `proteins.rds`. The design sits in its own step, so a design
problem surfaces before anything is fitted.

[`01_Design`](01_Design/README.md) builds the design matrix and the five contrasts and estimates the
within-leg correlation. It writes `design.rds` and `01_design.xlsx`.
[`02_Differential`](02_Differential/README.md) fits `dpcDE()`, applies the contrasts with
`eBayes(robust = TRUE)`, and writes `fit.rds` and `02_differential.xlsx`.
[`03_Protein_Association`](03_Protein_Association/README.md) fits each protein's change against the
ultrasound outcomes and writes `association_t.rds`, `association_designs.rds`,
`03_protein_association.xlsx`, one PDF per outcome and a chance summary.
[`04_Protein_Classification`](04_Protein_Classification/README.md) measures how well each protein
separates a participant's BFR leg from their HLRT leg in three tasks and writes
`04_protein_classification.xlsx`, ROC PDFs for two of the tasks and a chance summary.

## What each step answers

| Step | Question | Answer |
|---|---|---|
| `02_Differential` | Which proteins change? | 86 with BFR training and 161 with HLRT at FDR 0.05; none between the arms or in the interaction. |
| `03_Protein_Association` | Does a protein's change track the phenotype? | No protein survives BH. |
| `04_Protein_Classification` | Does a protein separate BFR from HLRT legs? | At chance at baseline, after training and in change. |
