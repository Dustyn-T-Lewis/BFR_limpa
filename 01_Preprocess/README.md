# 01 · Preprocess

Turns the DIA-NN report into a protein matrix. Two sub-stages, run in order. Each one reads
the previous one's `c_data` folder and writes its own. Nothing passes in memory.

```
00_Input/report.parquet + metadata.csv
  01_Filtering       repair the annotation, then remove rows  -> precursors_filtered.rds
  02_Quantification  fit the DPC, roll up, filter             -> proteins.rds
                                                                 |
                                              02_Differential_Expression
```

Run both from the repo root:

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

## Why the stages are shaped this way

This is the limpa vignette's sequence with nothing added: `readDIANN`, the two precursor
filters, `dpc`, `dpcQuant`, `filterByDetection`, then `dpcDE` in the next stage. The only
code here that is not a limpa call repairs what the search FASTA broke and removes proteins
skeletal muscle does not express. Both have to run before `dpc()` sees the matrix.

limpa reverses two habits from older proteomics pipelines. It does not impute, and it does
not want you to filter missing values out. A precursor nobody detected is evidence that it
was faint, `dpc()` fits that relationship, and `dpcQuant()` then uses it. Filter on
missingness and you throw away the information the method exists to use.

## There is no normalization stage

limpa ships no normalization function. Its vignette does not normalize. Its published
DIA-NN case study does not normalize. `readDIANN()` reads DIA-NN's `Precursor.Normalised`
column by default, which is where the normalization already happened: DIA-NN applies a
retention-time correction during the search, and this pipeline uses that output.

An earlier version of this stage estimated one offset per sample and subtracted it before
the DPC. We removed it on 2026-09-07. The measurement behind it was a 1.37 log2 spread in
column medians, but that number is the range across 131 samples, which two extreme samples
can set on their own, and it sat on top of a correction DIA-NN had already applied.

limpa's documentation mentions normalization once, in `filterByDetection()`'s help: "after
`dpcQuant` but before normalization or `dpcDE`". That would put it on the protein matrix,
not the precursor matrix, and limpa still gives you no function for it.

Anything that refers to `sample_offsets.csv`, `normalization_estimators.csv` or
`precursors_normalized.rds` describes the old pipeline.

## The DPC slope

Reported as fitted. We never substitute a preset value. The limpa FAQ puts the acceptable
range at 0.1 to 1.0 and calls 0.7 to 0.9 typical for DIA-NN searched with match-between-runs.
limpa's own DIA-NN case study fitted 0.59 and used it. The vignette says the slope comes out
low when peptides vary a lot, which is what a within-subject design with heavy missingness
looks like. A low slope also recovers less from missing values, so it errs toward finding
nothing rather than toward finding too much.

## Packages

`limpa` for everything it covers, and `dplyr`, `stringr`, `purrr`, `readr` and `tibble` for
the handling around it. `readDIANN()` needs `nanoparquet` to read a DIA-NN v2 report, so
install it even though no line here calls it. `limma` enters at
`02_Differential_Expression`, not here.
