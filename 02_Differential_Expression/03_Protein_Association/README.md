# 03 · Protein association

Tests whether each protein's pre-to-post change tracks the same leg's change in muscle size and
echo intensity. The association is a slope in the limpa model (`dpcDE()`, leg and arm terms,
`eBayes(robust = TRUE)`), fitted over all legs (pooled) and within participant (differential), with
BH within analysis and outcome. Its moderated t ranks proteins for gene set tests.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds` (models) and `proteins_adjusted.rds`
(paired changes and figures), `00_Input/phenotype.csv`.

```sh
quarto render 02_Differential_Expression/03_Protein_Association/a_script/03_protein_association.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/03_protein_association.html`: the walkthrough, with previews and the chance figure.
- `b_reports/03_association_<outcome>.pdf`: one per outcome (`vl_csa`, `vl_echo`, `rf_csa`,
  `rf_echo`), every protein at nominal p, twelve to a page. Pooled pages plot each leg, with the
  within-arm Spearman r inset; differential pages plot the BFR minus HLRT differences, one point
  per participant. Lines are least-squares fits (within each arm on pooled pages); panel headers
  give the tested slope's t, p and BH q within analysis and outcome.
- `b_reports/03_chance.pdf`: nominal hits against chance. Bar length is nominal hits divided by
  the chance count; red bars exceed 1, and labels give observed of tested.
- `c_data/association_t.rds`: moderated t, one row per protein, one column per analysis and
  outcome, for gene set tests.
- `c_data/association_designs.rds`: the eight design matrices, read by
  `03_Pathway_Enrichment/05_set_association` for fry.
- `c_data/03_protein_association.xlsx`:
  - `summary`: one row per analysis and outcome; nominal and BH hits against chance;
  - `associations`: one row per protein, analysis and outcome; slope, t, p and FDR, sorted by FDR.
    Read by `04_Protein_Classification`.

The closest lead is AGMAT, which rises more in legs whose RF CSA grew more (pooled, t 4.5,
q 0.090). SMTNL1, which falls with training in both arms, falls further in legs that grew more in
RF CSA, both pooled (t −3.4) and within participant (t −4.3, q 0.40). TJP1 runs against both VL
and RF CSA in the pooled model (t −3.7). Nothing survives BH, and nominal counts sit near or below
the 152 expected by chance from 3,042 proteins, a count that swings widely because proteins move
together: shuffling vl_csa across legs within arm gave 63 to 308 nominal hits (middle 90% of 40
shuffles; before the blood covariate).

Nominal and BH hits per analysis and outcome:

| Analysis | Outcome | Nominal | Ratio to chance | BH < 0.05 |
|---|---|---:|---:|---:|
| pooled | vl_csa | 137 | 0.90 | 0 |
| pooled | vl_echo | 162 | 1.07 | 0 |
| pooled | rf_csa | 124 | 0.82 | 0 |
| pooled | rf_echo | 120 | 0.79 | 0 |
| differential | vl_csa | 67 | 0.44 | 0 |
| differential | vl_echo | 101 | 0.66 | 0 |
| differential | rf_csa | 129 | 0.85 | 0 |
| differential | rf_echo | 116 | 0.76 | 0 |
