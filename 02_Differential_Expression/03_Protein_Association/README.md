# 03 · Protein association

Correlates each protein's pre-to-post change with the same leg's change in muscle size and echo
intensity (Spearman), pooled over all legs and as BFR minus HLRT within participant, with BH within
analysis and outcome.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`, `00_Input/phenotype.csv`,
`02_Differential/c_data/02_differential.xlsx` (training hits).

```sh
quarto render 02_Differential_Expression/03_Protein_Association/a_script/03_protein_association.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/03_protein_association.html`: the walkthrough, with previews and the chance figure.
- `b_reports/03_association_<outcome>.pdf`: one per outcome (`vl_csa`, `vl_echo`, `rf_csa`,
  `rf_echo`), every protein at nominal p, twelve to a page, with the within-arm correlations inset.
- `b_reports/03_chance.pdf`: nominal hits against chance.
- `c_data/03_protein_association.xlsx`: hits against chance, training overlap, every correlation
  (`correlations`, read by `04_Protein_Classification`) and the descriptive within-arm correlations.

| Analysis | Outcome | n | Nominal | Ratio to chance | BH < 0.05 |
|---|---|---:|---:|---:|---:|
| pooled | vl_csa | 65 | 183 | 1.20 | 0 |
| pooled | vl_echo | 65 | 127 | 0.83 | 0 |
| pooled | rf_csa | 65 | 134 | 0.88 | 0 |
| pooled | rf_echo | 65 | 113 | 0.74 | 0 |
| differential | vl_csa | 32 | 146 | 0.96 | 0 |
| differential | vl_echo | 32 | 130 | 0.85 | 0 |
| differential | rf_csa | 32 | 168 | 1.10 | 0 |
| differential | rf_echo | 32 | 116 | 0.76 | 0 |

Chance is 152 nominal hits from 3,042 proteins. Nothing survives BH (best adjusted p 0.281). Of the
130 training hits, 14 reach nominal p against pooled vl_csa where 7.8 are expected (one-sided Fisher
p 0.022), and 12 against pooled rf_csa where 5.7 are expected (p 0.011); the other six sit at
chance.
