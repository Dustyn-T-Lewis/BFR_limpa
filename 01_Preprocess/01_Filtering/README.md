# 01_Preprocess / 01_Filtering

Reads the DIA-NN report and removes rows. It never drops a sample, and it never filters on
missing values.

| | |
|---|---|
| **Script** | `a_script/01_filter.qmd` |
| **Reads** | `00_Input/report.parquet`, `00_Input/metadata.csv` |
| **Writes** | `c_data/precursors_filtered.rds` plus two CSVs recording what it removed |
| **Next** | `01_Preprocess/02_Quantification` |

## What it does, in order

**Read the report.** limpa's `readDIANN()` applies the q-value cutoffs and returns abundances on
the log2 scale.

**Name the samples.** Columns arrive named by MS run. One match against `metadata.csv` renames
them to our sample IDs and drops the one injection that was discarded.

**Check the samples that made it.** A short block reports how many peptides and how much signal
each sample carried, and fails if any of them looks like a failed injection. It removes nothing.
`metadata.csv` decides exclusions, and this checks that decision was right.

**Repair the search database fault.** The FASTA used for the search listed every protein twice,
once normally and once tagged as a contaminant. That makes creatine kinase and myoglobin look
ambiguous, and limpa's next two filters would delete them along with most of the muscle. So we
drop only the entries that are contaminant-only, then strip the tag from the rest. This has to
happen first: once the tag is gone, a contaminant-only entry is indistinguishable from a real
one.

**Remove peptides that map to more than one protein**, then **remove protein groups holding more
than one protein**. Both are limpa functions, run in the order its guide gives.

**Remove contaminants.** A needle biopsy bleeds and passes through skin, so a large share of the
signal is blood, plasma, skin and antibody protein. These go before quantification, because the
detection curve is fitted across the whole matrix and contaminants left in place would shape it.

## How contaminants are decided

A list of gene symbols in the notebook, plus four patterns for the antibody, keratin,
small-proline-rich and haemoglobin families. Symbols were checked against single-cell muscle
expression rather than whole-tissue values, because tissue measurements are themselves
blood-contaminated.

A protein group counts as contamination only when **every** gene name on it does. One protein can
carry several names, so testing just the first one would get it wrong in both directions.

Twelve proteins are blood-abundant but genuinely expressed in muscle, so they are kept on
purpose.

## The guard

One assertion at the end: a list of 28 proteins that must survive, mixing muscle markers with the
twelve blood-abundant keepers. It exists in order to fail. One line catches both a broken repair
and a contaminant rule that reaches too far.
