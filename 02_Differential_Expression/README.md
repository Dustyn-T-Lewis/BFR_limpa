# 02 · Differential Expression

This stage fits the protein model to `proteins.rds`. The design sits in its own step, so a design
problem surfaces before anything is fitted.

[`01_Design`](01_Design/README.md) builds the design matrix (group, participant and the blood index)
and the five contrasts and estimates the within-leg correlation. It writes `design.rds` and `01_design.xlsx`.
[`02_Differential`](02_Differential/README.md) fits `dpcDE()`, applies the contrasts with
`eBayes(robust = TRUE)`, and writes `fit.rds` and `02_differential.xlsx`.
[`03_Protein_Association`](03_Protein_Association/README.md) fits each protein's change against the
ultrasound outcomes and writes `association_t.rds`, `association_designs.rds`,
`03_protein_association.xlsx`, one PDF per outcome and a chance summary.
[`04_Protein_Classification`](04_Protein_Classification/README.md) measures how well each protein
separates a participant's BFR leg from their HLRT leg in three tasks and writes
`04_protein_classification.xlsx`, ROC PDFs for two of the tasks and a chance summary.
[`05_Response_Magnitude`](05_Response_Magnitude/README.md) asks whether the BFR response is smaller
overall than the HLRT response, with three tests, and writes `05_response_magnitude.xlsx` and one
figure.

## What each step answers

| Step | Question | Answer |
|---|---|---|
| `02_Differential` | Which proteins change? | 95 with BFR training and 218 with HLRT at FDR 0.05; none between the arms or in the interaction. |
| `03_Protein_Association` | Does a protein's change track the phenotype? | No protein survives BH. |
| `04_Protein_Classification` | Does a protein separate BFR from HLRT legs? | At chance at baseline, after training and in change. |
| `05_Response_Magnitude` | Is one arm's response smaller overall? | BFR about 0.77 of HLRT by signature score (p 0.031); cross-fitted magnitude and fry lean the same way without passing. |
