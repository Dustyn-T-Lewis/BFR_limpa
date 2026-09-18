# 01_Preprocess / 01_Filtering

Reads the DIA-NN report and removes rows. Drops no sample, filters nothing on missing values.

| | |
|---|---|
| **Script** | `a_script/01_filter.qmd` |
| **Reads** | `00_Input/report.parquet`, `00_Input/metadata.csv` |
| **Writes** | `c_data/precursors_filtered.rds`, `c_data/01_filter.xlsx` |

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports
```

Only `precursors_filtered.rds` is read downstream. The workbook holds two sheets, `filter_log`
and `contaminants_removed`, recording what this stage did.

## What it does

1. Read the report with `readDIANN()`.
2. Rename columns from MS run names to our sample IDs, dropping one discarded injection. Every
   metadata run must match a report column, or the sample would silently become all-NA.
3. Repair the search database fault and drop reversed decoys.
4. Remove peptides mapping to several proteins, then groups holding several proteins.
5. Remove blood, plasma, skin and antibody proteins, reporting each cell's contaminant share
   first, since uneven contamination would land on a contrast.

## The repair

The contaminant database appended to the search FASTA repeats proteins the human FASTA already
lists, so creatine kinase and myoglobin arrive tagged as contaminants and look ambiguous. Left
alone, step 4 deletes them along with most of the muscle. So we drop only the entries that are
contaminant-only, then strip the tag from the rest. The order matters: once the tag is gone, a
contaminant-only entry looks real.

## Contaminants

147 gene symbols plus four family patterns, defined in the notebook. Symbols were screened against
single-cell muscle expression, not whole-tissue values, which are themselves blood-contaminated.

A group counts as contamination only when **every** gene name on it does. One protein can carry
several names, so testing the first would get it wrong both ways.

## The guard

One assertion: all 28 symbols in `MUST_KEEP` must survive. It exists in order to fail, and catches
both a broken repair and a contaminant rule that reaches too far.
