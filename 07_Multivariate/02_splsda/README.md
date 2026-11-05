# 02 · Multilevel sPLS-DA

Asks whether a weighted panel of proteins separates two groups in participants the model has not
seen. Four comparisons, each on two samples per participant with `multilevel = participant`
(each participant's mean removed first): T1 against T2 in the BFR leg and in the HLRT leg, and BFR
against HLRT on each leg's training change and on the post-training samples. Training is the
positive control; the arm comparisons ask whether any protein combination separates the arms,
which the per-protein tests of stage 02 cannot rule out.

For each comparison `tune.splsda()` picks a panel of 5 to 1,000 proteins per component by balanced
error rate (BER; 0.5 is chance), `splsda()` fits it on everyone and `perf()` cross-validates it.
Cross-validation is 5-fold, repeated 50 times, seeded, with folds of whole participants. The panel
size is chosen on the same data, so the reported errors are somewhat optimistic; there is no
permutation reference. A two-factor multilevel design on all four samples is not used: with one
sample per participant, leg and time, mixOmics computes its decomposition from the class labels.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 07_Multivariate/02_splsda/a_script/02_splsda.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_splsda.html`: the walkthrough.
- `b_reports/02_splsda.pdf`, one A4 page per comparison: samples on the two component scores with
  a line joining each participant's pair; the ROC curve of each of the 50 cross-validation repeats
  from component 1's held-out predictions, with the mean AUC; the 25 largest loadings on each
  component, shaded by stability.
- `c_data/02_splsda.xlsx`:
  - `models`: one row per comparison and component: participants, samples, panel size, tuning BER,
    `perf()` BER (mean and SD under the max distance, mean under the centroid and Mahalanobis
    distances) and cross-validated AUC with SD;
  - `panels`: one row per selected protein, comparison and component: gene, loading (its weight in
    the component score; the sign gives the direction) and stability (the share of the 250
    cross-validation refits that selected it);
  - `scores`: one row per sample and comparison: participant, group and the two component scores.

On current data training separates clearly: component 1 BER 0.014 and AUC 0.999 in the HLRT leg
(20 proteins), 0.18 and 0.87 in the BFR leg (1,000 proteins). The arm comparisons sit at chance:
BER 0.52 and AUC 0.48 on change, 0.46 and 0.58 after training. Choosing the panel size on the same
data pulls a tuned BER below 0.5 even without signal, so neither result shows a protein combination
that separates BFR from HLRT, consistent with stage 02.
