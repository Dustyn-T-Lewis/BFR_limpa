# 01_Preprocess / 01_Filtering

Reads the DIA-NN report and removes rows. Drops no sample, filters nothing on missing values.

| | |
|---|---|
| **Script** | `a_script/01_filter.qmd` |
| **Reads** | `00_Input/report.parquet`, `00_Input/metadata.csv` |
| **Writes** | `c_data/precursors_filtered.rds` and two CSVs recording what went |

## What it does

1. Read the report with `readDIANN()`.
2. Rename columns from MS run names to our sample IDs, dropping one discarded injection.
3. Check every retained sample clears half the median depth. Removes nothing.
4. Repair the search database fault.
5. Remove peptides mapping to several proteins, then groups holding several proteins.
6. Remove blood, plasma and skin proteins, recording each sample's contaminant share by timepoint
   first, since uneven contamination between visits would land on the two time contrasts.
7. Assert 28 named proteins survived, and list what left.

## The repair

The search FASTA listed every protein twice, once tagged as a contaminant, which makes creatine
kinase and myoglobin look ambiguous. Left alone, step 5 deletes them along with most of the
muscle. So we drop only the entries that are contaminant-only, then strip the tag from the rest.
The order matters: once the tag is gone, a contaminant-only entry looks real.

## Contaminants

147 gene symbols plus four family patterns, in the notebook. Symbols were checked against
single-cell muscle expression, not whole-tissue values, which are themselves blood-contaminated.

A group counts as contamination only when **every** gene name on it does. One protein can carry
several names, so testing the first would get it wrong both ways.

## The guard

One assertion: all 28 symbols in `MUST_KEEP` must survive. It exists in order to fail, and
catches both a broken repair and a contaminant rule that reaches too far.
