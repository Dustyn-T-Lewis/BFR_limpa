# 01 · Run fgsea and fry

Tests every gene set on every contrast with fgsea and fry, and marks the significant fgsea sets
that are not redundant with a stronger one.

Inputs: `00_build_gene_sets/c_data/gene_sets.rds` and `00_build_gene_sets.xlsx`,
`02_Differential_Expression/02_Differential/c_data/fit.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`,
`01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 03_Pathway_Enrichment/01_run_fgsea_and_fry/a_script/01_run_fgsea_and_fry.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_run_fgsea_and_fry.html`: the walkthrough.
- `b_reports/01_run_fgsea_and_fry_figures.pdf`: S1 (ten strongest collapse survivors per
  contrast, all collections), S2–S6 (the same within each collection) and S7 (significant sets
  before and after collapse).
- `c_data/01_run_fgsea_and_fry.xlsx`: `set_tests`, one row per set, contrast and method, and
  `protein_results`, every protein and contrast from the fit; `02_`, `03_` and `08_` read both,
  `07_` reads `set_tests`.

Significant sets at an adjusted p below 0.05:

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 376 | 88 | 656 |
| HLRT_Post-Pre | 281 | 68 | 589 |
| BFR_Post-HLRT_Post | 159 | 52 | 0 |
| Modality_x_Time_Interaction | 26 | 10 | 0 |
| BFR_Pre-HLRT_Pre *(baseline)* | 1 | 1 | 0 |

fry finds nothing on the between-leg contrasts, so the interaction's 10 fgsea pathways are leads
rather than findings. The baseline contrast's one hit, `REACTOME_STRIATED_MUSCLE_CONTRACTION`, is
the floor for every count.
