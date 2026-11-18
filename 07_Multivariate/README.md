# 07 · Multivariate

Asks with mixOmics (Rohart et al. 2017, doi:10.1371/journal.pcbi.1005752) whether any combination
of proteins separates the groups or predicts growth, judged by cross-validated prediction in
participants left out.

- [`01_pca`](01_pca/README.md): what drives the variation within a person, by multilevel PCA of
  all samples with each participant's mean removed.
- [`02_splsda`](02_splsda/README.md): whether a protein panel separates T1 from T2 in each arm, or
  BFR from HLRT before, after and on change, and whether each arm's training panel transfers to
  the other.
- [`03_spls`](03_spls/README.md): whether a panel of protein changes predicts each leg's CSA and
  echo change, per arm and pooled.

Cross-validation is 5-fold repeated 50 times with folds of whole participants; the pooled
regressions hold out one participant per fold. Multilevel designs are one-factor (participant):
two-factor designs leak the labels on this crossover. Each step README gives its inputs, outputs
and results.
