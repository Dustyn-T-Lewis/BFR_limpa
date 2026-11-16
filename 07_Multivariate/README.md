# 07 · Multivariate

Asks with mixOmics (Rohart et al. 2017, doi:10.1371/journal.pcbi.1005752) whether any combination
of proteins separates the groups or predicts growth, judged by cross-validated prediction in
participants left out.

- [`01_pca`](01_pca/README.md): multilevel PCA of all samples, each participant's mean removed.
- [`02_splsda`](02_splsda/README.md): multilevel sPLS-DA of T1 against T2 in each arm and of BFR
  against HLRT before, after and on change; transfer of each arm's training panel to the other.
- [`03_spls`](03_spls/README.md): sPLS regression of each leg's CSA and echo change on its protein
  change, per arm and pooled.

Cross-validation is 5-fold repeated 50 times with folds of whole participants; the pooled
regressions hold out one participant per fold. Panel sizes are tuned on the same data and there is
no permutation reference, so error rates lean optimistic. Two-factor multilevel designs are not
used: on this crossover they leak the labels.
