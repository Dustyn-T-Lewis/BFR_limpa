# 03 · Complex subunits

Shows every measured subunit and assembly factor of complexes I to V with its training change in
each leg, so a reader can see which proteins carry a complex. Descriptive only; the set-level test
is fry in `02_test_mito`.

Inputs: `DEP_matrix` from `02_Differential_Expression/02_Differential/c_data/02_differential.xlsx`,
the `membership` sheet of `01_build_mito_sets` and the `complex_tests` sheet of `02_test_mito`.

```sh
quarto render 05_Mitochondria/03_complex_subunits/a_script/03_complex_subunits.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/03_complex_subunits.html`: the walkthrough.
- `b_reports/03_complex_subunits.pdf`: page 1 stacks Sam's two summaries from PR #26, complex
  subunits and then assembly factors with "Other" (CYCS, HCCS and the respirasome assembly genes
  tied to no complex): mean change per group and leg as bars with the standard error across the
  proteins, proteins as points, axis labels giving measured of all MitoCarta members, and fry's
  interaction FDR above each complex's subunits. The standard error describes spread among
  correlated proteins and is not a test. Pages 2 and 3 show every subunit and every assembly
  factor by leg, ordered within complex by BFR minus HLRT.
- `c_data/03_complex_subunits.xlsx`, one sheet:
  - `subunits`: one row per complex, role and measured protein, with training logFC and BH FDR per
    leg and for the interaction (BFR minus HLRT) from the stage 02 fit: subunits, assembly factors
    and the "Other" proteins. TMEM70 is listed under both CI and CV.

No protein's interaction reaches FDR 0.05, and fry's interaction FDR is 0.93 for every complex.
