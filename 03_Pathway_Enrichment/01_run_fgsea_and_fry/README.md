# 01 · Run fgsea and fry

Tests every gene set on every contrast with fgsea and fry, and marks the significant fgsea sets
that are not redundant with a stronger one. fry uses the stage 02 design, which carries the blood
index as a covariate.

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
  survivor has no panel. Each collection's own view is the dotplot in `07_pathway_maps`. S2
  counts sets at FDR 0.05 before and after `collapsePathways`, which makes a set the child of a
  stronger one when it shows no enrichment of its own (p > 0.05).
- `c_data/01_run_fgsea_and_fry.xlsx`; `02_`, `06_` and `07_` read `set_tests` and
  `protein_results`:
  - `set_summary`: one row per contrast: significant sets by fgsea, after collapse, and by fry.
  - `protein_summary`: one row per contrast: proteins up and down at FDR 0.05, smallest FDR.
  - `set_tests`: one row per set, contrast and method; `main` marks collapse survivors.
  - `protein_results`: one row per protein and contrast from the fit, with the pi-score.

Significant sets at an adjusted p below 0.05, BH over all 1,990 sets in the contrast:

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 433 | 102 | 635 |
| HLRT_Post-Pre | 279 | 62 | 410 |
| BFR_Post-HLRT_Post | 154 | 53 | 0 |
| Modality_x_Time_Interaction | 38 | 25 | 0 |
| BFR_Pre-HLRT_Pre *(baseline)* | 70 | 29 | 0 |

Training raises mitochondrial pathways in both legs. The top sets are the same in each arm: aerobic
respiration and electron transport (NES 3.2 with BFR, 2.8 with HLRT), oxidative phosphorylation,
mitochondrial gene expression and mitochondrial protein degradation, all with BH q below 10⁻⁹. The
BFR response in these sets is as strong as the HLRT one.

Between the arms, the interaction points to structural remodelling that HLRT training drives and
BFR training does not. The sets lower in BFR than in HLRT are cell motility, cell adhesion,
cytoskeleton organisation, wound healing, epithelial–mesenchymal transition and Rho GTPase effectors
(fgsea q 0.0008 to 0.015). These are exploratory leads: fry does not confirm any between-leg set,
and the baseline contrast, where no difference is expected, also gives fgsea 70 sets, so fgsea calls
between legs run high in these data.
