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
- `b_reports/03_complex_subunits.pdf`: Sam's summary from PR #26 (mean subunit change per complex
  and leg), then every subunit and every assembly factor by leg.
- `c_data/03_complex_subunits.xlsx`, one sheet:
  - `subunits`: one row per complex, role and measured protein, with training logFC and BH FDR per
    leg and for the interaction (BFR minus HLRT) from the stage 02 fit. 88 subunits and 42 assembly
    factors; TMEM70 is listed under both CI and CV, so the sheet has 131 rows.

No protein's interaction reaches FDR 0.05, and fry's interaction FDR is 0.99 for every complex.
