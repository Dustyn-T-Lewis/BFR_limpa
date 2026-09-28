# 01 · Preprocess

This stage turns the DIA-NN output into a normalised protein table with a standard error on every
estimate. It has two steps, and they run in order because the second reads the first's `.rds`.

[`01_Filtering`](01_Filtering/README.md) repairs the search annotation, applies limpa's precursor
filters and removes contaminants. It writes `precursors_filtered.rds` and `01_filter.xlsx`.
[`02_Quantification`](02_Quantification/README.md) fits the detection curve, rolls precursors up
to proteins, applies the detection filter and normalises with cyclic loess. It writes
`proteins.rds`, `02_quantify.xlsx` and the checkpoints in `quant_runs/`.
