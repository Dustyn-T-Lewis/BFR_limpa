# 02 · Differential Expression

This stage fits the protein model to `proteins.rds`. The design sits in its own step, so a design
problem surfaces before anything is fitted.

[`01_Design`](01_Design/README.md) builds the design matrix (group, participant and the blood index)
and the five contrasts and estimates the within-leg correlation. It writes `design.rds` and
`01_design.xlsx`. [`02_Differential`](02_Differential/README.md) fits `dpcDE()`, applies the
contrasts with `eBayes(robust = TRUE)`, and writes `fit.rds` and `02_differential.xlsx`.
[`03_Protein_Association`](03_Protein_Association/README.md) fits each protein's change against the
ultrasound outcomes and writes `association_t.rds`, `association_designs.rds`,
`03_protein_association.xlsx`, one PDF per outcome and a chance summary. Whether one arm's response
is smaller overall is answered in `08_Plasticity`.

## What each step answers

| Step | Question | Answer |
|---|---|---|
| `02_Differential` | Which proteins change? | 95 with BFR training and 218 with HLRT at FDR 0.05, 51 in both, all moving the same way. Interaction leads at nominal p (BGN, PC, RPL36, MSN higher with HLRT); none passes BH between the arms or in the interaction. |
| `03_Protein_Association` | Does a protein's change track the phenotype? | Closest lead AGMAT with RF CSA (pooled, q 0.090); no protein survives BH. |
