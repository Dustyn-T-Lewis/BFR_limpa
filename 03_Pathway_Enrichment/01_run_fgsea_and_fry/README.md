# 01 · Run fgsea and fry

Tests every gene set on every contrast with fgsea and fry, and marks the significant fgsea sets
not redundant with a stronger one. fry uses the stage 02 design, blood index included.

Inputs: `00_build_gene_sets/c_data/gene_sets.rds` and `00_build_gene_sets.xlsx`,
`02_Differential_Expression/02_Differential/c_data/fit.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`,
`01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 03_Pathway_Enrichment/01_run_fgsea_and_fry/a_script/01_run_fgsea_and_fry.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/01_run_fgsea_and_fry.html`: the walkthrough.
- `b_reports/01_run_fgsea_and_fry_figures.pdf`: S1, the ten strongest collapse survivors per
  contrast (size: measured genes, colour: collection); S2, significant sets before and after
  `collapsePathways`, which makes a set the child of a stronger one when it shows no enrichment
  of its own (p > 0.05).
- `c_data/01_run_fgsea_and_fry.xlsx`; `02_`, `06_` and `07_` read `set_tests` and
  `protein_results`:
  - `set_summary`: one row per contrast: significant sets by fgsea, after collapse, and by fry.
  - `protein_summary`: one row per contrast: proteins up and down at FDR 0.05, smallest FDR.
  - `set_tests`: one row per set, contrast and method; `main` marks collapse survivors.
  - `protein_results`: one row per protein and contrast from the fit, with the pi-score.

Sets at adjusted p < 0.05, BH over all 1,990 sets in the contrast:

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 433 | 102 | 635 |
| HLRT_Post-Pre | 279 | 62 | 410 |
| BFR_Post-HLRT_Post | 154 | 53 | 0 |
| Modality_x_Time_Interaction | 38 | 25 | 0 |
| BFR_Pre-HLRT_Pre (baseline) | 70 | 29 | 0 |

Training raises mitochondrial pathways in both legs, with the same top sets in each arm: aerobic
respiration and electron transport (NES 3.2 with BFR, 2.8 with HLRT), oxidative phosphorylation,
mitochondrial gene expression and mitochondrial protein degradation, all at BH q below 10⁻⁹.

Between the arms, fgsea puts cell motility, cell adhesion, cytoskeleton organisation, wound
healing, epithelial-mesenchymal transition and Rho GTPase effector sets lower in BFR than HLRT
(interaction q 0.0008 to 0.015). fry confirms no between-leg set, and the baseline contrast,
where no difference is expected, also gives 70 fgsea sets, so these are leads only.
