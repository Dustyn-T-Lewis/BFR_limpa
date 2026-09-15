# 01 · Preprocess

DIA-NN output to a protein table with uncertainty attached.

```
00_Input/report.parquet + metadata.csv
  01_Filtering       clean up the precursor table       -> precursors_filtered.rds
  02_Quantification  precursors to proteins, normalised -> proteins.rds
```

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

Run them in order: the second reads the first's `.rds`.

## What comes out

`proteins.rds` holds three things per protein per sample: abundance, a standard error, and how
many precursors were detected. Abundances are cyclic-loess normalised at the end of
`02_Quantification`; the standard errors are left on the unnormalised scale, which that stage's
README states as an approximation. Later stages read this file. The standard error is the point, so
keep the `.rds` rather than exporting abundances alone.

## Two orderings that matter

The search-database repair runs **before** any filter, because the fault it fixes makes limpa's
filters delete real muscle protein.

Nothing is filtered on missing values **before** quantification, because those gaps are what the
detection curve is fitted to.

## Cost

Both notebooks render in seconds. `dpcQuant()` costs about 100 minutes per call, so it does not
live in a notebook at all: `02_Quantification/a_script/dpc_quant.R` computes it, you run that by
hand, and `02_quantify.qmd` only loads the result from `c_data/dpcQuant/`. That output is
committed, so a fresh clone renders without computing anything. No stage uses knitr caching, so
no `_cache/` directory is ever created.
