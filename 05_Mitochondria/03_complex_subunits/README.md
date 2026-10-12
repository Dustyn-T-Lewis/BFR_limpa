# 03_complex_subunits

This step shows every measured subunit and assembly factor of complexes I to V with its training
change in each leg, so a reader can see which proteins carry a complex. It reads `DEP_matrix` from
`02_Differential_Expression/02_Differential/c_data/02_differential.xlsx`, the membership sheet of
`01_build_mito_sets` and the complex tests of `02_test_mito`. It writes
`c_data/03_complex_subunits.xlsx` and `b_reports/03_complex_subunits.pdf`.

```sh
Rscript 05_Mitochondria/03_complex_subunits/a_script/03_complex_subunits.R
```

The idea comes from Sam's PR #26. His version tested BFR against HLRT with a paired t-test across
the subunits of each complex, which treats correlated proteins as independent replicates and makes
p too small. Here the bars are descriptive and each panel title carries fry's interaction FDR on
that complex's subunits from `02_test_mito`, the set-level test that allows for correlation.
Assembly factors have no set test, so their panels carry counts only.

The table holds 88 subunits and 43 assembly factors with their logFC and BH FDR in each leg and
for the interaction. No subunit's interaction reaches FDR 0.05, and fry's interaction FDR is 0.99
for every complex.
