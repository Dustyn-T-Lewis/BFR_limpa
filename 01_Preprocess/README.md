# 01 · Preprocess

Turns the DIA-NN search output into a protein table with uncertainty attached. Two sub-stages.

```
00_Input/report.parquet + metadata.csv
  01_Filtering       clean up the peptide table    -> precursors_filtered.rds
  02_Quantification  peptides to proteins          -> proteins.rds
                                                      |
                                       02_Differential_Expression
```

| Sub-stage | What it does |
|---|---|
| `01_Filtering` | reads the report, repairs a fault in the search database, removes peptides that cannot be assigned to one protein, removes blood and skin proteins |
| `02_Quantification` | learns how missingness relates to abundance, then rolls peptides up to proteins |

Run both from the repo root:

```
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd        --output-dir ../b_reports
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

## What comes out

`proteins.rds` is the file every later stage reads. It holds three things per protein per
sample: an abundance, a standard error, and a count of how many peptides were actually detected
behind it. The standard error is the point. A protein reconstructed mostly from missing peptides
still gets a number, but a wide one, and the statistics downstream weight it accordingly.

Keep the `.rds`. A CSV of the abundances alone drops the standard errors, which makes the
analysis ordinary again.

## Why the order is what it is

This follows limpa's own sequence with nothing added: read, remove ambiguous peptides, fit the
detection curve, quantify, filter on detection.

Two orderings are load-bearing. **The search-database repair happens before any filter**, because
the fault it fixes makes limpa's filters delete real muscle protein. And **nothing is filtered on
missing values before quantification**, because those gaps are the data the detection curve is
fitted to.

## No normalization step

There isn't one, and that is deliberate. limpa ships no normalization function and its own
worked examples do not normalize, because the DIA-NN column we read is already normalized during
the search. An earlier version of this stage added a second correction on top. Removing it
changed essentially nothing, which was the answer.

## Packages

`limpa` for the quantification, `nanoparquet` so it can read the report, and
`dplyr`/`stringr`/`purrr`/`readr`/`tibble` for the handling around it. `limma` is not needed
until the next stage.

## Cost

The quantification step takes about 100 minutes on a laptop and holds roughly 11 GB. It is
cached and keyed to its input file, so it re-runs when the filtering output changes and not
otherwise. Filtering itself takes a few minutes.
