# BFR Proteomics

DIA mass-spectrometry proteomics from a unilateral resistance training trial. Each
participant trained one leg with blood-flow restriction and the other with conventional
high load, with biopsies taken from both legs before and after training. 33 participants
and 131 MS runs: four samples per participant, except one participant whose
post-training biopsy from one leg was not acquired.

The primary comparison is the interaction: whether blood-flow restriction alters the
muscle proteome differently than heavy load.

## Stages

| Stage | Contents | State |
|---|---|---|
| `00_Input/` | study data, no code | ready |
| `01_Preprocess/` | search output to a normalised protein table | ready |
| `02_Differential_Expression/` | model fitting, five contrasts, protein against phenotype | ready |
| `03_Pathway_Enrichment/` | gene set tests, per-sample set scores, set classification and phenotype association | ready |
| `04_Network/` | proteins that move together | planned |
| `05_Figures/` | manuscript panels | planned |

Each sub-stage holds `a_script/` (code), `b_reports/` (HTML reports or figures) and `c_data/`
(outputs). Stages 01 and 02 are Quarto notebooks; stage 03 is plain R scripts. Data passes
through disk, so any sub-stage re-runs on its own. "Planned" means only a README exists.

## Data

`00_Input/report.parquet` is too large for git. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet
```

## Running the pipeline

```sh
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

Everything above takes minutes. `dpcQuant()` takes about 100 and runs separately, from
`01_Preprocess/02_Quantification/a_script/02_quantify_run.R`. Its output is committed, so re-run
it only when the precursor matrix changes.

## Approach

Quantification uses limpa. Roughly half of a DIA matrix is missing, and missingness is
not random: low-abundance peptides are missed more often than abundant ones. limpa models
that relationship and treats a missing value as evidence of low abundance rather than
imputing a replacement. Every protein receives an estimate in every sample along with a
standard error, and that uncertainty propagates into the downstream statistics.

Two consequences for the workflow: do not filter on missingness before quantification, and
do not pass the protein matrix to a standard linear model.

## Dependencies

`limpa`, `limma`, `here`, `nanoparquet`, `writexl`, `readr`, `dplyr`, `tidyr`, `tibble`, `purrr`,
`stringr`, `ggplot2`. Stage 03 adds `fgsea`, `singscore`, `msigdbr`, `GO.db`,
`GSEABase`, `AnnotationDbi`, `pROC`, `ggrepel`, `patchwork`, `qpdf` and `enrichVolcano` (not on
CRAN). Versions are not pinned.
