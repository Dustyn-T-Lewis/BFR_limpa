# 07 · Multivariate

Stages 02 to 04 test proteins one at a time, curated sets and data-driven modules. This stage asks
the same three questions with mixOmics (Rohart et al. 2017, doi:10.1371/journal.pcbi.1005752),
where a panel of proteins is built to fit the outcome and judged by cross-validated prediction in
participants left out. It answers one thing the earlier stages cannot: whether any combination of
proteins separates the groups or predicts growth.

[`01_pca`](01_pca/README.md) draws a multilevel PCA of all samples with each participant's mean
removed.

[`02_splsda`](02_splsda/README.md) classifies T1 against T2 in each arm (the positive control) and
BFR against HLRT on change and after training, with multilevel sparse PLS-DA.

[`03_spls`](03_spls/README.md) predicts each leg's CSA and echo change from its protein change with
sparse PLS regression, per arm and pooled.

Every model uses 5-fold cross-validation repeated 50 times, with folds of whole participants. Panel
sizes are tuned on the same data and there is no permutation reference, so error rates lean
optimistic. Two-factor multilevel designs are not used: on this crossover they leak the labels.
