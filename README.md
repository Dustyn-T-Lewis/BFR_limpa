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
| `02_Differential_Expression/` | model fitting, five contrasts | ready |
| `03_Pathway_Enrichment/` | frozen gene sets, protein plots, enrichment tests | first sub-stage ready, tests planned |
| `04_Network/` | data-driven protein modules | planned |
| `05_Figures/` | manuscript panels | planned |

Each sub-stage contains `a_script/` (code), `b_reports/` (rendered HTML), and `c_data/`
(outputs). Stages run in order and pass data through disk, so any stage can be re-run on
its own. "Planned" means the folder currently contains only a README.

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

quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd --output-dir ../b_reports
```

These renders take minutes. The one expensive step is `dpcQuant()` at roughly 100 minutes,
which runs separately from
`01_Preprocess/02_Quantification/a_script/02_quantify_run.R`. Its output is committed, so
the renders above load it from disk. Re-run that script only when the precursor matrix
changes.

## Approach

Quantification uses **limpa**. Roughly half of a DIA matrix is missing, and missingness is
not random: low-abundance peptides are missed more often than abundant ones. limpa models
that relationship and treats a missing value as evidence of low abundance rather than
imputing a replacement. Every protein receives an estimate in every sample along with a
standard error, and that uncertainty propagates into the downstream statistics.

Two consequences for the workflow: do not filter on missingness before quantification, and
do not pass the protein matrix to a standard linear model.

## Dependencies

`limpa`, `limma`, `here`, `nanoparquet`, `writexl`, `dplyr`, `stringr`, `purrr`, `readr`,
`tibble`, `tidyr`, `ggplot2`. `03_Pathway_Enrichment` additionally requires `msigdbr`,
`enrichVolcano`, `digest`, and `Matrix`, and checks for them at runtime. Versions are not
pinned.
