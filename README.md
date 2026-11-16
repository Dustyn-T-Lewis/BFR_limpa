# BFR Proteomics

DIA mass-spectrometry proteomics from a unilateral resistance training trial: each participant
trained one leg with blood-flow restriction (BFR) and the other with high load (HLRT). Both legs
were biopsied before and after training, giving 131 MS runs from 33 participants; one
post-training biopsy was not acquired. The primary comparison is the interaction: whether BFR
changes the muscle proteome differently from high load.

limpa models missing values instead of imputing them, giving every protein an estimate and a
standard error in every sample, and the standard error carries into the downstream statistics. So
nothing is filtered on missingness before quantification.

## Stages

| Stage | Contents | State |
|---|---|---|
| `00_Input/` | study data, no code | ready |
| `01_Preprocess/` | search output to a normalised protein table | ready |
| `02_Differential_Expression/` | design, five contrasts, protein association | ready |
| `03_Pathway_Enrichment/` | gene set tests, scores, classification, association and maps | ready |
| `04_Network/` | WGCNA co-expression modules and their phenotype associations | ready |
| `05_Mitochondria/` | mitochondrial content, complexes and pathways | ready |
| `06_Composition/` | slow fibre fraction (FibeRtypeR) and its training response | ready |
| `07_Multivariate/` | multilevel PCA, sPLS-DA and sPLS regression (mixOmics) | ready |
| `08_Plasticity/` | BFR against HLRT response as a whole | ready |
| `0X_Figures/` | manuscript panels, built last | planned |

## Layout

Each step holds `README.md`, `a_script/` (a Quarto notebook), `b_reports/` (the rendered HTML and
any figure PDF) and, where it writes tables, `c_data/` (one workbook, plus an `.rds` where a later
step needs the R object). Each step's README lists its workbook sheets. Data passes through disk,
so any step reruns from a fresh session. "Planned" means only a README exists.

## Running the pipeline

`00_Input/report.parquet` is too large for git. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet

quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports

d=02_Differential_Expression
quarto render $d/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render $d/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
quarto render $d/03_Protein_Association/a_script/03_protein_association.qmd \
  --output-dir ../b_reports

for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea 03_run_singscore \
         04_set_classification 05_set_association 06_leading_edge 07_pathway_maps; do
  quarto render 03_Pathway_Enrichment/$s/a_script/$s.qmd --output-dir ../b_reports
done

for s in 01_modules 02_module_traits; do
  quarto render 04_Network/$s/a_script/$s.qmd --output-dir ../b_reports
done

for s in 01_build_mito_sets 02_test_mito 03_complex_subunits; do
  quarto render 05_Mitochondria/$s/a_script/$s.qmd --output-dir ../b_reports
done

quarto render 06_Composition/02_fibre_type/a_script/02_fibre_type.qmd --output-dir ../b_reports

for s in 01_pca 02_splsda 03_spls; do
  quarto render 07_Multivariate/$s/a_script/$s.qmd --output-dir ../b_reports
done

quarto render 08_Plasticity/01_plasticity/a_script/01_plasticity.qmd --output-dir ../b_reports
```

`dpcQuant()` takes about 100 minutes, so it runs separately from
`01_Preprocess/02_Quantification/a_script/02_quantify_run.R`; its output is committed, so rerun it
only when the precursor matrix changes.

## Dependencies

Install Quarto and the R packages `pacman`, `knitr` and `rmarkdown`. Every notebook attaches its
packages with `pacman::p_load()`, which installs any that are missing. Packages called only as
`pkg::fun()` are not attached, so install them beforehand.

- Stages 01 and 02 attach `here`, `limpa`, `limma`, `dplyr`, `tidyr`, `tibble`, `purrr`,
  `stringr`, `readr`, `readxl`, `writexl` and `ggplot2`, and call `nanoparquet` (read by
  `readDIANN()`).
- Stage 03 adds `patchwork`, `clusterProfiler`, `enrichplot` and `ggvenn`, and calls `fgsea`,
  `singscore`, `msigdbr`, `GO.db`, `GSEABase`, `AnnotationDbi`, `pROC` and `enrichVolcano` 2.0.0
  or later (GitHub, `Dustyn-T-Lewis/enrichVolcano`, not on CRAN).
- Stage 04 adds `WGCNA` and `ggrepel`, and calls `ggnewscale`.
- Stages 05 and 06 need nothing new.
- Stage 07 adds `mixOmics` (Bioconductor).
- Stage 08 calls `RRHO2` (GitHub, `RRHO2/RRHO2`).
