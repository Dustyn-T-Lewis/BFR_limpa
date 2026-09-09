# BFR Proteomics

DIA mass-spectrometry proteomics from a unilateral resistance training trial. Each participant
trained one leg with blood-flow restriction and the other with conventional high load, with a
muscle biopsy from both legs before and after. We want to know whether restriction changes the
muscle proteome differently from heavy load.

131 MS runs, 33 participants, four samples each. Because both legs belong to one person, every
participant is their own control, which cancels diet, sleep and hormones without having to
model them.

## Goal

Turn the DIA-NN search output into a protein table we can trust, test the five comparisons the
study design supports, then ask what the results mean biologically.

The primary question is the interaction: what does restriction add, over and above training
the leg hard? Everything in the pipeline exists to answer that one honestly.

## Stages

| Stage | Does | State |
|---|---|---|
| `00_Input/` | holds the study data. Runs nothing. | ready |
| `01_Preprocess/` | search output to a protein table with uncertainty | ready |
| `02_Differential_Expression/` | fit the model, test five contrasts | planned |
| `03_Pathway_Enrichment/` | ask which biological processes moved | planned |
| `04_Network/` | find protein groups the data itself defines | planned |
| `05_Figures/` | the manuscript panels | planned |

Stages run in order. Each one reads the previous stage's output from disk and writes its own,
so any stage can be re-run alone. Nothing is passed in memory.

"Planned" means the folder holds this README and nothing else. Pilot code exists offline and is
not ready to be read or relied on.

## Folder layout

Every sub-stage carries the same three folders:

- `a_script/` the code, one notebook per step
- `b_reports/` the rendered HTML and the figures it drew
- `c_data/` the tables and matrices that step produced

Only `a_script/` is committed. Reports and data regenerate from code, so `.gitignore` excludes
them. Four people committing their own render of the same HTML would make every pull request a
conflict in a file nobody reads.

Numbers in folder and file names are the run order, and numbering restarts inside each
sub-stage.

## Approach

We use **limpa** rather than a summarise-then-impute pipeline. Half a DIA matrix is empty, and
not at random: faint peptides go missing more often than abundant ones. limpa fits that
relationship and treats a missing value as evidence the protein was low, instead of filling in
a guess. Every protein gets a value in every sample plus a standard error saying how much to
trust it, and that uncertainty travels into the statistics.

The practical consequence is a rule that holds throughout: **never filter on missing values
before quantification, and never hand the protein matrix to an ordinary linear model.** Both
throw away the information limpa exists to use.

## Packages

`limpa` does the quantification and the differential testing. `limma` supplies the linear
modelling underneath it. `nanoparquet` is needed to read the DIA-NN report even though no line
calls it directly. Data handling is `dplyr`, `stringr`, `purrr`, `readr` and `tibble`. Later
stages will add gene-set and network packages, named in their own READMEs.

Package versions are not pinned. That is a deliberate trade for a small project, and it means a
re-run months from now may not reproduce exactly.

## Running it

Run from the project root, where the `.here` marker makes paths resolve.

```
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

`--output-dir` is what puts reports in `b_reports/`. No notebook sets it itself.

Budget about two hours, nearly all of it the quantification step. It is cached, so re-rendering
to fix text costs seconds.

## Before you start

`00_Input/report.parquet` is too large for git, so it is attached to a release. From the repo
root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet
```

Raw instrument files are held at the collaborating site and will be deposited in MassIVE, with
the accession released on publication.
