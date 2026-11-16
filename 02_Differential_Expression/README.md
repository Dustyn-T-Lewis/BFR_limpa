# 02 · Differential Expression

Fits the protein models to `proteins.rds`. The design has its own step, so a design problem
surfaces before anything is fitted.

1. [`01_Design`](01_Design/README.md): design matrix (group, participant, blood index) and the five
   contrasts.
2. [`02_Differential`](02_Differential/README.md): `dpcDE()` fit and contrast tests.
3. [`03_Protein_Association`](03_Protein_Association/README.md): slopes of protein change on the
   ultrasound outcomes and volume load.

Results: training changes 95 proteins with BFR and 218 with HLRT at FDR 0.05, 51 in both, all
moving the same way. No protein passes BH between the arms or in the interaction. No protein's
change tracks a CSA or echo change after BH; with volume load, CA14 and TOMM40 fall and HSPA12A
rises across BFR legs (q 0.02), and HSPA12A differs between arms (q 0.034). Whether one arm's
response is smaller overall is answered in `08_Plasticity`.
