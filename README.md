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
| `04_Network/` | proteins that move together | planned |
| `05_Mitochondria/` | mitochondrial content, complexes and pathways | ready |
| `06_Composition/` | fibre type, and DE adjusted for composition | export only |
| `0X_Figures/` | manuscript panels, built last | planned |

Each sub-stage holds `README.md`, `a_script/` (code), `b_reports/` (an HTML report for the Quarto
notebooks, figure PDFs for the R scripts) and `c_data/` (one workbook of tables, plus an `.rds`
where a later step needs the R object). The limpa steps (stage 01, `01_Design` and
`02_Differential`) are Quarto notebooks; every step after the fit is a plain R script. Data passes
through disk, so any sub-stage reruns from a fresh session, and no notebook caches. "Planned" means
only a README exists; "export only" means only the step before the FibeRtypeR web app exists.

## Running the pipeline

`00_Input/report.parquet` is too large for git. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet

quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports

quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
Rscript 02_Differential_Expression/03_Protein_Association/a_script/03_protein_association.R
Rscript 02_Differential_Expression/04_Protein_Classification/a_script/04_protein_classification.R

for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea \
         03_enrich_scatter_fgsea 04_run_singscore 05_set_classification 06_set_association \
         07_leading_edge 08_pathway_maps; do
  Rscript 03_Pathway_Enrichment/$s/a_script/$s.R
done

Rscript 05_Mitochondria/01_build_mito_sets/a_script/01_build_mito_sets.R
Rscript 05_Mitochondria/02_test_mito/a_script/02_test_mito.R
Rscript 05_Mitochondria/03_complex_subunits/a_script/03_complex_subunits.R
```

The block takes about five minutes. `dpcQuant()` takes about 100 minutes, so it runs separately from
`01_Preprocess/02_Quantification/a_script/02_quantify_run.R`; its output is committed, so rerun it
only when the precursor matrix changes.

## Approach

Half of a DIA matrix is missing, and low-abundance peptides go missing more often. limpa models that
instead of imputing, giving every protein an estimate and a standard error in every sample, and the
standard error carries into the downstream statistics. So nothing is filtered on missingness before
quantification, and the protein matrix never goes to a standard linear model.

## Dependencies

Every script loads its packages with `pacman::p_load()`, which installs any that are missing.
Stages 01 and 02 use `pacman`, `here`, `limpa`, `limma`, `nanoparquet` (read by `readDIANN()`),
`readr`, `writexl`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr` and `ggplot2`. Stage 03 adds
`readxl`, `fgsea`, `singscore`, `msigdbr`, `GO.db`, `GSEABase`, `AnnotationDbi`, `pROC`, `ggrepel`,
`patchwork`, `clusterProfiler`, `enrichplot`, `eulerr` and `enrichVolcano` 2.0.0 or later (GitHub,
`Dustyn-T-Lewis/enrichVolcano`, not on CRAN). Stages 05 and 06 need nothing new.
