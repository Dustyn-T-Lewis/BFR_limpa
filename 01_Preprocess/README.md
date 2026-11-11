# 01 · Preprocess

This stage turns the DIA-NN output into a normalised protein table with a standard error on every
estimate, in two steps run in order.

[`01_Filtering`](01_Filtering/README.md) repairs the search annotation, applies limpa's precursor
filters and removes contaminants, writing `precursors_filtered.rds`, the haemoglobin precursors it
sets aside (`haemoglobin_precursors.rds`) and `01_filter.xlsx`.
[`02_Quantification`](02_Quantification/README.md) reads those `.rds` files, fits the detection
curve as a diagnostic, rolls precursors up to proteins with the slope preset at 0.7, applies the
detection filter and normalises with cyclic loess. It adds a blood index per sample, the mean of
the HBB and HBA1 estimates, and writes `proteins.rds`, `proteins_adjusted.rds` (the blood index's
effect removed, for steps that take no design), `02_quantify.xlsx` and the checkpoints in
`quant_runs/`.
