# 02 · Differential Expression

This stage takes `proteins.rds` and asks what happened to the proteins. The design sits in its own
step, so a design problem surfaces before anything is fitted.

[`01_Design`](01_Design/README.md) builds the design matrix and the five contrasts and estimates the
within-leg correlation. It writes `design.rds` and `01_design.xlsx`.
[`02_Differential`](02_Differential/README.md) fits `dpcDE()`, applies the contrasts with
`eBayes(robust = TRUE)` and `treat()`, and writes `fit.rds` and `02_differential.xlsx`.
[`03_Protein_Association`](03_Protein_Association/README.md) correlates every protein with the
ultrasound outcomes and writes `03_protein_association.xlsx` and a PDF of every nominal pair.
[`04_Protein_Classification`](04_Protein_Classification/README.md) measures how well each protein
separates the study groups and writes one ROC PDF per task. The first two are Quarto notebooks,
since they hold the limpa steps; the last two are plain R scripts.

## What each step answers

| Step | Question | Answer |
|---|---|---|
| `02_Differential` | Which proteins change? | 60 with BFR training and 107 with HLRT at FDR 0.05; none between the arms or in the interaction. |
| `03_Protein_Association` | Does a protein's change track the phenotype? | No protein survives BH. |
| `04_Protein_Classification` | Does a protein separate the groups? | 53 on BFR training and 220 on HLRT after BH; BFR against HLRT sits at chance. |
