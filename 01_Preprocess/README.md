# 01 · Preprocess

This stage turns the DIA-NN output into a normalised protein table with a standard error on every
estimate, in two steps run in order.

[`01_Filtering`](01_Filtering/README.md) repairs the search annotation, applies limpa's precursor
filters and removes contaminants, writing `precursors_filtered.rds` and `01_filter.xlsx`.
[`02_Quantification`](02_Quantification/README.md) reads that `.rds`, fits the detection curve as a
diagnostic, rolls precursors up to proteins with the slope preset at 0.7, applies the detection
filter and normalises with cyclic loess, writing `proteins.rds`, `02_quantify.xlsx` and the
checkpoints in `quant_runs/`.
