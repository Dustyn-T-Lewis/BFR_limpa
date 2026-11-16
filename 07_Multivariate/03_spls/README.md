# 03 · sPLS regression on muscle growth

Asks whether a weighted panel of protein changes (T2 minus T1) predicts a leg's change in vastus
lateralis and rectus femoris cross-sectional area (CSA) and echo intensity (post minus pre). One
sPLS1 model per outcome, fitted on BFR legs, on HLRT legs and on both arms pooled.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds` and `00_Input/phenotype.csv`.

```sh
quarto render 07_Multivariate/03_spls/a_script/03_spls.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/03_spls.html`: the walkthrough.
- `b_reports/03_spls.pdf`: one page per outcome, held-out predicted against observed change for
  each fit (r from one cross-validation pass, Q² from `perf()`), and the 20 largest component 1
  loadings per fit shaded by stability.
- `c_data/03_spls.xlsx`:
  - `models`: per outcome, fit and component, legs, panel size, tuning MAE, and `perf()` MSEP,
    RMSEP, R², Q² and total Q²;
  - `panels`: per selected protein, outcome, fit and component, gene, loading and stability;
  - `predictions`: per leg, outcome and fit, arm, observed and held-out predicted change.

No fit predicts growth: component 1 Q² is below the package's 0.0975 threshold in all 12 fits and
below 0 in 10. The highest are 0.09 (HLRT, VL echo) and 0.00 (HLRT, VL CSA). In the BFR fits
held-out predictions correlate negatively with observed change for both echo outcomes (r −0.46 RF,
−0.28 VL), the pattern of a model fitting noise. The pooled panel size is tuned on folds that can
split a participant's two legs, so the pooled errors may lean optimistic. Growth is answered
protein by protein in `02_Differential_Expression/03_Protein_Association`, by set in
`03_Pathway_Enrichment/05_set_association` and by module in `04_Network/02_module_traits`.
