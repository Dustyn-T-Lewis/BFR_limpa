# BFR Proteomics

DIA mass-spectrometry proteomics from a unilateral resistance training trial. Each participant
trained one leg with blood-flow restriction and the other with conventional high load, with a
biopsy from both legs before and after. 33 participants and 131 MS runs: four samples each,
except one participant whose post-training biopsy on one leg was never acquired.

The question is the interaction: does restriction change the muscle proteome differently from
heavy load?

## Stages

| Stage | Does | State |
|---|---|---|
| `00_Input/` | study data, runs nothing | ready |
| `01_Preprocess/` | search output to a normalised protein table | ready |
| `02_Differential_Expression/` | fit the model, test five contrasts | ready |
| `03_Pathway_Enrichment/` | frozen gene sets, then which processes moved | ready |
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

quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd   --output-dir ../b_reports
quarto render 03_Pathway_Enrichment/02_Set_Tests/a_script/02_set_tests.qmd   --output-dir ../b_reports
quarto render 03_Pathway_Enrichment/03_Volcanoes/a_script/03_volcanoes.qmd   --output-dir ../b_reports
```

Minutes, not hours. No notebook computes anything slow. The one expensive step is `dpcQuant()` at
about 100 minutes, and it lives in `01_Preprocess/02_Quantification/a_script/02_quantify_run.R`.
Its output is committed, so these renders load it. Run that script by hand only when the precursor
matrix changes.

## Approach

We use **limpa**. Half a DIA matrix is empty, and not at random: faint peptides go missing more
often than abundant ones. limpa fits that relationship and treats a missing value as evidence
the protein was low, rather than filling in a guess. Every protein gets a value in every sample
plus a standard error, and that uncertainty carries into the statistics.

Two rules follow. Never filter on missing values before quantification, and never hand the
protein matrix to an ordinary linear model.

Packages: `limpa`, `limma`, `here`, `nanoparquet`, `writexl`, and `dplyr`, `stringr`, `purrr`,
`readr`, `tibble`, `tidyr`, `ggplot2`. `03_Pathway_Enrichment` adds `fgsea`, `singscore`,
`msigdbr` and `enrichVolcano`. Versions are not pinned.
