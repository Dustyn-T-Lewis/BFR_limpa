# BFR Proteomics

DIA mass-spectrometry proteomics from a unilateral resistance training trial: each participant
trained one leg with blood-flow restriction and the other with conventional high load. Both legs
were biopsied before and after training, giving 131 MS runs from 33 participants; one post-training
biopsy from one leg was not acquired.

The primary comparison is the interaction: whether blood-flow restriction alters the muscle
proteome differently than heavy load.

## Stages

| Stage | Contents | State |
|---|---|---|
| `00_Input/` | study data, no code | ready |
| `01_Preprocess/` | search output to a normalised protein table | ready |
| `02_Differential_Expression/` | model fitting, five contrasts, protein association and classification | ready |
| `03_Pathway_Enrichment/` | gene set tests, per-sample set scores, set classification and phenotype association | ready |
| `04_Network/` | protein co-expression modules (WGCNA) | ready |
| `05_Mitochondria/` | mitochondrial content, complexes and pathways | ready |
| `06_Composition/` | fibre type, and DE adjusted for composition | export only |
| `0X_Figures/` | manuscript panels, built last | planned |

Each sub-stage holds `README.md`, `a_script/` (a Quarto notebook), `b_reports/` (the rendered HTML
report, beside any figure PDFs) and, where the step writes tables, `c_data/` (one workbook, plus an
`.rds` where a later step needs the R object). Each step's README lists its workbook sheets. Data
passes through disk, so any sub-stage reruns from a fresh session, and no notebook caches.
"Planned" means only a README exists; "export only" means only the step before the FibeRtypeR web
app exists.

## Running the pipeline

`00_Input/report.parquet` is too large for git. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet

quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports

quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
quarto render 02_Differential_Expression/03_Protein_Association/a_script/03_protein_association.qmd       --output-dir ../b_reports
quarto render 02_Differential_Expression/04_Protein_Classification/a_script/04_protein_classification.qmd --output-dir ../b_reports

for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea 03_run_singscore \
         04_set_classification 05_set_association 06_leading_edge 07_pathway_maps; do
  quarto render 03_Pathway_Enrichment/$s/a_script/$s.qmd --output-dir ../b_reports
done

quarto render 04_Network/01_modules/a_script/01_modules.qmd --output-dir ../b_reports

for s in 01_build_mito_sets 02_test_mito 03_complex_subunits; do
  quarto render 05_Mitochondria/$s/a_script/$s.qmd --output-dir ../b_reports
done

quarto render 06_Composition/01_fibre_type_export/a_script/01_fibre_type_export.qmd --output-dir ../b_reports
```

The block takes about ten minutes. `dpcQuant()` takes about 100 minutes, so it runs separately from
`01_Preprocess/02_Quantification/a_script/02_quantify_run.R`; its output is committed, so rerun it
only when the precursor matrix changes.

## Approach

Half of a DIA matrix is missing, and low-abundance peptides go missing more often. limpa models that
instead of imputing, giving every protein an estimate and a standard error in every sample, and the
standard error carries into the downstream statistics. So nothing is filtered on missingness before
quantification, and the protein matrix never goes to a standard linear model.

## Dependencies

Install Quarto and the R packages `pacman`, `knitr` and `rmarkdown` first. Every notebook then
attaches its packages with `pacman::p_load()`, which installs any that are missing. Packages a
notebook only calls as `pkg::fun()` are not attached, so install them once beforehand.

- Stages 01 and 02 attach `here`, `limpa`, `limma`, `readr`, `readxl`, `writexl`,
  `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr` and `ggplot2`, and call `nanoparquet` (read by
  `readDIANN()`) and `pROC`.
- Stage 03 adds `patchwork`, `clusterProfiler`, `enrichplot` and `ggvenn`, and calls `fgsea`,
  `singscore`, `msigdbr`, `GO.db`, `GSEABase`, `AnnotationDbi` and `enrichVolcano` 2.0.0 or later
  (GitHub, `Dustyn-T-Lewis/enrichVolcano`, not on CRAN).
- Stage 04 adds `WGCNA`.
- Stages 05 and 06 need nothing new.
