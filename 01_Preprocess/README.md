# 01 · Preprocess

Turns the DIA-NN report into a normalised protein matrix with a standard error on every value.

1. [`01_Filtering`](01_Filtering/README.md): which precursors are trustworthy evidence for one
   human muscle protein? It repairs the cRAP annotation, applies limpa's precursor filters and
   removes blood and skin contaminants.
2. [`02_Quantification`](02_Quantification/README.md): how abundant is each protein, and how sure
   is that estimate? It rolls precursors up with `dpcQuant()`, filters on detection, normalises
   with cyclic loess and adds a blood index per sample.
