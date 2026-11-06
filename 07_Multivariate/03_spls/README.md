# 03 · sPLS regression on muscle growth

Asks whether a weighted panel of protein changes predicts how much a leg changed in vastus
lateralis and rectus femoris cross-sectional area (CSA) and echo intensity. One sparse PLS model
per outcome (sPLS1), protein change T2 minus T1 against outcome change post minus pre, fitted
three ways: BFR legs alone and HLRT legs alone (one leg per participant, so cross-validation folds
are folds of people), and both arms pooled as a sensitivity run.

`tune.spls()` picks a panel of 5 to 1,000 proteins per component by mean absolute error (mixOmics
divides each miss by the outcome's variance, so the value compares sizes within an outcome only),
`spls()` fits it and `perf()` cross-validates it, 5-fold repeated 50 times and seeded. A
participant's two legs grow together (r = 0.57 for VL CSA), so the pooled `perf()` holds out each
participant's legs together; `tune.spls()` cannot, so the pooled panel size may lean optimistic.
The panel size is chosen on the same data in every fit, so errors are somewhat optimistic; there
is no permutation reference. `perf()` does not return its predictions, so the figure's held-out
predictions come from one cross-validation pass of `spls()` and `predict()` on the same folds.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds` and `00_Input/phenotype.csv`.

```sh
quarto render 07_Multivariate/03_spls/a_script/03_spls.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/03_spls.html`: the walkthrough.
- `b_reports/03_spls.pdf`, one A4 page per outcome: held-out predicted against observed change for
  the BFR, HLRT and pooled fits with the correlation and component 1's Q²; below, the 20 largest
  component 1 loadings of each fit, shaded by stability.
- `c_data/03_spls.xlsx`:
  - `models`: one row per outcome, fit and component: legs, panel size, tuning MAE, and from
    `perf()` the MSEP and RMSEP (relative to the outcome's variance; 1 is no better than the
    mean), R², Q² and total Q²;
  - `panels`: one row per selected protein, outcome, fit and component: gene, loading and
    stability across the cross-validation refits;
  - `predictions`: one row per leg, outcome and fit: arm, observed change and held-out predicted
    change from component 1.

On current data no fit predicts growth: component 1 Q² is below the package's 0.0975 threshold
in all 12 fits and below 0 in 10 of them. The highest are 0.09 (HLRT, VL echo) and 0.00 (HLRT,
VL CSA). In the BFR fits held-out predictions correlate negatively with observed change for both
echo outcomes (r −0.46 RF, −0.28 VL), the pattern of a model fitting noise. This matches the
absence of change associations in stages 02 to 04, where the only CSA association is one
post-training level in `04_Network`.
