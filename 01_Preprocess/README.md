# 01 · Preprocess

DIA-NN output to a protein table with uncertainty attached.

```
00_Input/report.parquet + metadata.csv
  01_Filtering       clean up the peptide table  -> precursors_filtered.rds
  02_Quantification  peptides to proteins        -> proteins.rds
```

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

## What comes out

`proteins.rds` holds three things per protein per sample: abundance, a standard error, and how
many peptides were actually detected. Later stages read this file. The standard error is the
point, so keep the `.rds` rather than exporting abundances alone.

## Two orderings that matter

The search-database repair runs **before** any filter, because the fault it fixes makes limpa's
filters delete real muscle protein.

Nothing is filtered on missing values **before** quantification, because those gaps are what the
detection curve is fitted to.

## No normalization step

limpa ships no normalization function and its own worked examples do not normalize, because the
DIA-NN column we read is already normalized. An earlier version of this stage added a second
correction. Removing it changed nothing, which was the answer.

Quantification takes about 100 minutes and 11 GB. It is cached and keyed to its input file.
