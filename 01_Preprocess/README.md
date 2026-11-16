# 01 · Preprocess

Turns the DIA-NN report into a normalised protein matrix with a standard error on every value.

1. [`01_Filtering`](01_Filtering/README.md): repairs the search annotation, applies limpa's
   precursor filters and removes contaminants.
2. [`02_Quantification`](02_Quantification/README.md): rolls precursors up to proteins with
   `dpcQuant()` (slope 0.7), filters on detection, normalises with cyclic loess and adds a blood
   index per sample.
