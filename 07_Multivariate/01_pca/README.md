# 01 · Multilevel PCA

PCA of all 131 samples after removing each participant's mean profile (mixOmics
`pca(multilevel = participant)`), so it shows how a person's four samples differ rather than how
people differ. Proteins are unscaled, as in the `02_Quantification` MDS plot.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds`.

```sh
quarto render 07_Multivariate/01_pca/a_script/01_pca.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_pca.html`: the walkthrough.
- `b_reports/01_pca.pdf`: samples on PC1 and PC2 by arm and time, each leg's pair joined, and the
  20 largest loadings on each component.
- `c_data/01_pca.xlsx`:
  - `variance`: proportion of variance per component, first 10;
  - `structure`: R² of timepoint, treatment and group on PC1 to PC4;
  - `scores`: per sample, metadata and PC1 to PC5 scores;
  - `loadings`: per protein, gene and PC1 to PC5 loadings.

PC1 (13.5%) carries sarcomere and Z-disc proteins (ACTN2, ACTN3, MYOZ1, MYOZ2, TTN, NEB, DES) and
PC2 (9.0%) extracellular matrix and fat (COL14A1, MYOC, FMOD, THBS4, PLIN1, FASN). Timepoint
explains at most 5% of each of the first three components and 11% of PC4, so within a person the
samples differ more by fibre and matrix content than by training. The PCA is descriptive and has no
test.
