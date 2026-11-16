# 03 · Pathway Enrichment

Tests whether proteins that work together moved together, and scores every sample on every set.
The eight steps take about two minutes; each step README gives its inputs, outputs and results.

- [`00_build_gene_sets`](00_build_gene_sets/README.md): the one set list every later step reads,
  from a frozen MSigDB release and local GO Slim sets.
- [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md): fgsea and fry on every set and
  contrast; `collapsePathways` marks the non-redundant fgsea hits.
- [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md): protein volcanoes with the
  fgsea hits ringed.
- [`03_run_singscore`](03_run_singscore/README.md): a singscore per sample and set,
  `set_scores.rds`.
- [`04_set_classification`](04_set_classification/README.md): how well set scores separate the
  study groups.
- [`05_set_association`](05_set_association/README.md): whether set scores and protein slopes
  track the phenotype.
- [`06_leading_edge`](06_leading_edge/README.md): Sam's PR #26 analysis, whether both legs move
  the same pathways through the same proteins.
- [`07_pathway_maps`](07_pathway_maps/README.md): dotplots, enrichment maps, trees and protein
  networks of both training responses.

Figures run S1 to S22: S1 and S2 from `01_`, S3 and S4 from `02_`, S5 from `03_`, S6 to S10
from `04_`, S11 to S13 from `05_`, S14 to S18 from `06_` and S19 to S22 from `07_`. `01_`, `02_`,
`03_` and `06_` write one A4 PDF each; `07_` writes one per collection, and `04_` and `05_` one
per collection plus a chance summary.

## Results

- `01_`: training raises mitochondrial respiration, oxidative phosphorylation and mitochondrial
  gene expression in both legs: fgsea 433 BFR and 279 HLRT sets (102 and 62 after collapse), fry
  635 and 410. Between the arms fgsea puts cell motility, adhesion and cytoskeleton sets lower in
  BFR; fry confirms none.
- `04_`: scores separate pre from post in 30 sets with BFR and 112 with HLRT after BH,
  mitochondrial protein degradation among the strongest. Between the arms one set passes,
  HALLMARK_SPERMATOGENESIS on change (q 0.018).
- `05_`: no set tracks a CSA or echo change after BH or fry; the fgsea lead for RF CSA growth is
  mitochondrial translation, unconfirmed by fry. On volume load 59 sets pass fry (FDR 0.048) on
  the average slope.
- `06_`: across all 1,990 sets the two legs' NES correlate (Spearman rho 0.76), and in the 71
  pathways both legs move the leading edges overlap with median Jaccard index 0.55 (0.32 to 0.94).
- `07_`: figures only.

fgsea is competitive and assumes exchangeable proteins, which co-regulated sets violate; fry is
self-contained and does not, so between-leg fgsea calls without fry support are leads. BH runs
over all 1,990 sets per contrast in `01_`, within collection and task in `04_`, and within
collection, analysis and outcome in `05_`, so their counts answer different families of tests.
