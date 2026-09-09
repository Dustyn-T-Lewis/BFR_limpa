# BFR Proteomics

DIA mass-spectrometry proteomics from a unilateral resistance training trial. Each participant
trained one leg with blood-flow restriction and the other with conventional high load, with a
biopsy from both legs before and after. 131 MS runs, 33 participants, four samples each.

The question is the interaction: does restriction change the muscle proteome differently from
heavy load?

## Stages

| Stage | Does | State |
|---|---|---|
| `00_Input/` | study data, runs nothing | ready |
| `01_Preprocess/` | search output to a protein table | ready |
| `02_Differential_Expression/` | normalise, fit, test five contrasts | ready |
| `03_Pathway_Enrichment/` | which processes moved | planned |
| `04_Network/` | protein groups the data defines | planned |
| `05_Figures/` | manuscript panels | planned |

Each sub-stage holds `a_script/` (code), `b_reports/` (rendered HTML), `c_data/` (its outputs).
Stages run in order and pass data through disk, so any one can re-run alone. "Planned" means
the folder holds only a README.

## Get the data

`00_Input/report.parquet` is too large for git. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet
```

## Run it

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports

quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

Two scripts sit beside the pipeline and run only when needed: `02b_quantify_slope07.qmd` refits
quantification at a fixed detection slope, and `02b_normalization_report.qmd` decides which
normalisation the protein matrix needs.

About two hours, nearly all of it quantification. That step is cached.

## Approach

We use **limpa**. Half a DIA matrix is empty, and not at random: faint peptides go missing more
often than abundant ones. limpa fits that relationship and treats a missing value as evidence
the protein was low, rather than filling in a guess. Every protein gets a value in every sample
plus a standard error, and that uncertainty carries into the statistics.

Two rules follow. Never filter on missing values before quantification, and never hand the
protein matrix to an ordinary linear model.

Packages: `limpa`, `limma`, `nanoparquet`, and `dplyr`/`stringr`/`purrr`/`readr`/`tibble`.
Versions are not pinned.
