# BFR Proteomics

DIA mass-spectrometry proteomics from a unilateral resistance training trial. Each participant
trained one leg with blood-flow restriction and the other with conventional high load. Biopsies
came from both legs before and after training: 33 participants and 131 MS runs, four samples per
participant except one participant whose post-training biopsy from one leg was not acquired.

The primary comparison is the interaction: whether blood-flow restriction alters the muscle
proteome differently than heavy load.

## Stages 01 to 03 run; 04 and 05 are planned

| Stage | Contents | State |
|---|---|---|
| `00_Input/` | study data, no code | ready |
| `01_Preprocess/` | search output to a normalised protein table | ready |
| `02_Differential_Expression/` | model fitting, five contrasts, protein against phenotype | ready |
| `03_Pathway_Enrichment/` | gene set tests, per-sample set scores, set classification and phenotype association | ready |
| `04_Network/` | proteins that move together | planned |
| `05_Figures/` | manuscript panels | planned |

## Every sub-stage runs on its own

Each sub-stage holds `README.md`, `a_script/` (code), `b_reports/` (an HTML report in stages 01
and 02, one figure PDF per step in stage 03) and `c_data/` (one workbook of tables, plus an `.rds`
where a later step needs the R object). Stages 01 and 02 are Quarto notebooks; stage 03 is plain R
scripts. Data passes through disk, so any sub-stage re-runs from a fresh session. No notebook uses
knitr caching. "Planned" means only a README exists.

## One block runs the whole pipeline

`00_Input/report.parquet` is too large for git. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet

quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports

quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
quarto render 02_Differential_Expression/03_Phenotype/a_script/03_phenotype.qmd       --output-dir ../b_reports

for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea \
         03_enrich_scatter_fgsea 04_run_singscore 05_classify_and_associate_sets; do
  Rscript 03_Pathway_Enrichment/$s/a_script/$s.R
done
```

The block above ran in 155 seconds. `dpcQuant()` takes about 100 minutes and runs separately, from
`01_Preprocess/02_Quantification/a_script/02_quantify_run.R`. Its output is committed, so re-run
it only when the precursor matrix changes.

## limpa treats a missing value as evidence of low abundance

Half of a DIA matrix is missing, and low-abundance peptides are missed more often than
abundant ones. limpa models that relationship instead of imputing a replacement. Every protein
gets an estimate in every sample with a standard error, and the standard error propagates into
the downstream statistics.

So nothing is filtered on missingness before quantification, and the protein matrix never goes
to a standard linear model.

## Dependencies are not pinned

Every script loads its packages with `pacman::p_load()`, which installs any that are missing.
Stages 01 and 02 use `pacman`, `here`, `limpa`, `limma`, `nanoparquet` (read by `readDIANN()`),
`readr`, `writexl`, `dplyr`, `tidyr`, `tibble`, `purrr`, `stringr` and `ggplot2`. Stage 03 adds
`readxl`, `fgsea`, `singscore`, `msigdbr`, `GO.db`, `GSEABase`, `AnnotationDbi`, `pROC`, `ggrepel`,
`patchwork` and `enrichVolcano` 2.0.0 or later (GitHub, `Dustyn-T-Lewis/enrichVolcano`, not on
CRAN).
