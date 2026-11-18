# 02 · Differential Expression

Fits the protein models to `proteins.rds`. The design has its own step, so a design problem
surfaces before anything is fitted.

1. [`01_Design`](01_Design/README.md): how should samples be compared? It builds the design matrix
   (group, participant, blood index) and the five contrasts, and checks the within-leg correlation.
2. [`02_Differential`](02_Differential/README.md): which proteins change with training in each
   arm, differ between arms, or respond differently to the two arms? It fits `dpcDE()` and tests
   the contrasts.
3. [`03_Protein_Association`](03_Protein_Association/README.md): does a protein's training change
   track the same leg's change in muscle size, echo intensity or volume load? It fits slopes with
   `dpcDE()`.

Whether one arm's response is smaller overall is asked in `08_Plasticity`.
