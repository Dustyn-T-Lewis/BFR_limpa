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
  contrast, all collections) and S2 (significant sets before and after collapse). In S1 point
  size is the number of measured genes in the set and colour the collection; a contrast with no
  survivor has no panel. Each collection's own view is the dotplot in `07_pathway_maps`. S2 counts sets at FDR 0.05 before and after `collapsePathways`, which makes a set
  the child of a stronger one when it shows no enrichment of its own (p > 0.05).
- `c_data/01_run_fgsea_and_fry.xlsx`, both sheets read by `02_`, `06_` and `07_`:
  - `set_summary`: one row per contrast: significant sets by fgsea, after collapse, and by fry.
  - `protein_summary`: one row per contrast: proteins up and down at FDR 0.05, smallest FDR.
  - `set_tests`: one row per set, contrast and method; `main` marks collapse survivors.
  - `protein_results`: one row per protein and contrast from the fit, with the pi-score.

Significant sets at an adjusted p below 0.05, BH over all 1,990 sets in the contrast:

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 376 | 88 | 656 |
| HLRT_Post-Pre | 281 | 68 | 589 |
| BFR_Post-HLRT_Post | 159 | 52 | 0 |
| Modality_x_Time_Interaction | 26 | 10 | 0 |
| BFR_Pre-HLRT_Pre *(baseline)* | 1 | 1 | 0 |

fry finds no set on either between-leg contrast, BFR_Post-HLRT_Post or the interaction, so their
fgsea counts are leads only. The baseline contrast compares the legs before training, when no
difference is expected, so its one fgsea hit, `REACTOME_STRIATED_MUSCLE_CONTRACTION`, shows what a
contrast with nothing to find returns.
