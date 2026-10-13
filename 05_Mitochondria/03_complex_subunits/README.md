# 03_complex_subunits

This step shows every measured subunit and assembly factor of complexes I to V with its training
change in each leg, so a reader can see which proteins carry a complex. It reads `DEP_matrix` from
`02_Differential_Expression/02_Differential/c_data/02_differential.xlsx`, the membership sheet of
`01_build_mito_sets` and the complex tests of `02_test_mito`. It writes
`c_data/03_complex_subunits.xlsx` and `b_reports/03_complex_subunits.pdf`.

```sh
Rscript 05_Mitochondria/03_complex_subunits/a_script/03_complex_subunits.R
```

The PDF has three pages. The first is Sam's summary from PR #26: each complex's mean subunit change
in each leg, with standard error and every subunit as a point. The second and third draw every
subunit and every assembly factor by leg, ordered within complex by BFR minus HLRT.

Sam's version tested BFR against HLRT with a paired t-test across the subunits of each complex,
which treats correlated proteins as independent replicates and makes p too small (complex I came
out below 0.001). Here the bars are descriptive and each complex is labelled with fry's interaction
FDR from `02_test_mito`, the set-level test that allows for that correlation. Assembly factors have
no set test, so their panels carry counts only.

The table holds 88 subunits and 43 assembly factors with their logFC and BH FDR in each leg and
for the interaction. No subunit's interaction reaches FDR 0.05, and fry's interaction FDR is 0.99
for every complex.
