# 03 · Complex subunits

Shows every measured subunit and assembly factor of complexes I to V with its training change in
each leg, so a reader can see which proteins carry a complex (Sam's PR #26). Descriptive only; the
set-level test is fry in `02_test_mito`.

Inputs: the `protein_results` sheet of
`02_Differential_Expression/02_Differential/c_data/02_differential.xlsx`, the `membership` sheet of
`01_build_mito_sets` and the `complex_tests` sheet of `02_test_mito`.

```sh
quarto render 05_Mitochondria/03_complex_subunits/a_script/03_complex_subunits.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/03_complex_subunits.html`: the walkthrough.
- `b_reports/03_complex_subunits.pdf`: page 1, mean change per complex and leg with the standard
  error across proteins, for subunits (labelled with fry's interaction FDR) and for assembly
  factors with "Other" (CYCS, HCCS and respirasome assembly genes tied to no complex); pages 2
  and 3, every subunit and every assembly factor by leg, ordered by BFR minus HLRT.
- `c_data/03_complex_subunits.xlsx`:
  - `subunits`: per complex, role and measured protein, training logFC and BH FDR per leg and for
    the interaction from the stage 02 fit (TMEM70 is listed under both CI and CV).

The subunits of every complex rise on average in both legs, by 0.10 to 0.27 log2 with BFR and 0.12
to 0.19 with HLRT. Complex I rises a little less with BFR (0.10 against 0.15) and complex II more
(0.27 against 0.12). The largest single differences are in assembly factors and a few subunits:
CMC1 (BFR +0.70, HLRT −0.02), ATPAF2 and MT-CO1 rise more with BFR, while TMEM126A (−0.03, +0.40),
NDUFAF5 and COX8A rise more with HLRT. No protein's interaction reaches FDR 0.05, and fry's
interaction FDR is 0.73 for every complex; the standard errors describe spread among correlated
proteins and are not tests.
