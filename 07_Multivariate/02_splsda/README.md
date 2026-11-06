# 02 · Multilevel sPLS-DA

Asks whether a weighted panel of proteins separates two groups in participants the model has not
seen. Five comparisons, each on two samples per participant with `multilevel = participant`
(each participant's mean removed first): T1 against T2 in the BFR leg and in the HLRT leg, and BFR
against HLRT before training, after training and on each leg's training change. Training is the
positive control; the arm comparisons ask whether any protein combination separates the arms,
which the per-protein tests of stage 02 cannot rule out. Each arm's training panel then scores the
other arm's legs (transfer), and the two training panels are annotated by WGCNA module and
over-represented gene sets.

For each comparison `tune.splsda()` picks a panel of 5 to 1,000 proteins per component by balanced
error rate (BER; 0.5 is chance), `splsda()` fits it on everyone and `perf()` cross-validates it.
Cross-validation is 5-fold, repeated 50 times, seeded, with folds of whole participants. The panel
size is chosen on the same data, so the reported errors are somewhat optimistic; there is no
permutation reference. A two-factor multilevel design on all four samples is not used: with one
sample per participant, leg and time, mixOmics computes its decomposition from the class labels.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds`, the `membership` and
`modules` sheets of `04_Network/01_modules/c_data/01_modules.xlsx`, and
`03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds`.

```sh
quarto render 07_Multivariate/02_splsda/a_script/02_splsda.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_splsda.html`: the walkthrough.
- `b_reports/02_splsda.pdf`, A4 portrait: one page per comparison (samples on the two component
  scores with a line joining each participant's pair; the ROC curve of each of the 50
  cross-validation repeats from component 1's held-out predictions, with the mean AUC; the 25
  largest loadings on each component, shaded by stability); a last page counts each training
  panel's component 1 proteins by main-network WGCNA module.
- `c_data/02_splsda.xlsx`:
  - `models`: one row per comparison and component: participants, samples, panel size, tuning BER,
    `perf()` BER (mean and SD under the max distance, mean under the centroid and Mahalanobis
    distances) and cross-validated AUC with SD;
  - `panels`: one row per selected protein, comparison and component: gene, loading (its weight in
    the component score; the sign gives the direction), stability (the share of the 250
    cross-validation refits that selected it), main-network WGCNA module and its BH top set;
  - `transfer`: one row per training model: the arm it scores, the AUC of its component 1 score
    for T1 against T2 in that arm, and that arm's own cross-validated AUC;
  - `panel_pathways`: one row per training panel and gene set with nominal p < 0.05: overlap,
    hypergeometric p and BH q (`fgsea::fora()` against every measured gene);
  - `scores`: one row per sample and comparison: participant, group and the two component scores.

On current data training separates clearly: component 1 BER 0.006 and AUC 1.000 in the HLRT leg
(20 proteins), 0.18 and 0.88 in the BFR leg (1,000 proteins). The arm comparisons sit near chance:
BER 0.39 and AUC 0.67 before training (5 proteins), 0.44 and 0.60 after training, 0.46 and 0.57 on
change. Choosing the panel size on the same data pulls a tuned BER below 0.5 even without signal,
so none shows a protein combination that separates BFR from HLRT; the pre-training value, the
lowest, sits on the baseline contrast that has no stage 02 hit. The training signature transfers:
the HLRT panel scores BFR legs at AUC 0.91 (BFR's own 0.88), the BFR panel HLRT legs at 0.98
(HLRT's own 1.00). The BFR panel draws most from the aerobic-respiration module M3 (182 proteins)
and the muscle module M1 (124) and is enriched for oxidative phosphorylation, mitochondrial protein
degradation and the TCA cycle (BH q < 0.001); the 20-protein HLRT panel draws 7 proteins from M3,
and its top set, arginine and proline metabolism, misses BH (q 0.058).
