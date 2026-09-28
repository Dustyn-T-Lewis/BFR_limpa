# 02_Differential_Expression / 03_Protein_Association

This step correlates each protein's pre-to-post change with the same leg's change in muscle size and
quality. It reads `01_Preprocess/02_Quantification/c_data/proteins.rds`, `00_Input/phenotype.csv`
and the training hits in `02_Differential/c_data/02_differential.xlsx`, and writes Spearman r, p and
FDR per protein and outcome, plus the counts below, to `c_data/03_protein_association.xlsx`.
`b_reports/` holds one PDF per outcome (`03_association_vl_csa.pdf` and so on), drawing every
protein at nominal p twelve to a page with the within-arm correlations inset, and `03_chance.pdf`,
the hits against chance. The `protein_by_arm` sheet repeats each correlation inside BFR and inside
HLRT, for description only.

```sh
Rscript 02_Differential_Expression/03_Protein_Association/a_script/03_protein_association.R
```

Pooled correlations use all 65 legs. Differential ones use 32 participants, correlating the BFR
minus HLRT change in protein with the same difference in outcome.

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

Chance is 152 nominal hits from 3,042 proteins. Nothing survives BH (best adjusted p 0.281); with 32
pairs the smallest resolvable paired correlation is about 0.5, so this is a null at that resolution.

Of the 130 training hits in either arm, 14 reach nominal p against pooled vl_csa where 7.8 are
expected (one-sided Fisher p 0.022), and 12 against pooled rf_csa where 5.7 are expected (p 0.011).
The other six pairs sit at chance.

Protein classification by ROC is `04_Protein_Classification`; `02_Differential` gives the
model-based group test. Set-level association is `03_Pathway_Enrichment/06_set_association`.
