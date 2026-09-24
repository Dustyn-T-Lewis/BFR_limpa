# 01 · Preprocess

DIA-NN output to a protein table with uncertainty attached.

| Step | Runs | Writes |
|---|---|---|
| [`01_Filtering`](01_Filtering/README.md) | annotation repair, limpa's precursor filters, contaminant removal | `precursors_filtered.rds`, `01_filter.xlsx` |
| [`02_Quantification`](02_Quantification/README.md) | detection curve, protein roll-up, detection filter, cyclic loess | `proteins.rds`, `02_quantify.xlsx`, `quant_runs/` |

Run them in order: the second reads the first's `.rds`.
