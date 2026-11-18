# 03 · Pathway Enrichment

Asks whether proteins that work together moved together, and scores every sample on every set.
The eight steps run in order in about two minutes; each step README gives its inputs, outputs and
results.

- [`00_build_gene_sets`](00_build_gene_sets/README.md): which sets are tested? One set list for
  every later step, from a frozen MSigDB release and local GO Slim sets.
- [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md): which sets move on each contrast?
  fgsea and fry on every set; `collapsePathways` marks the non-redundant fgsea hits.
- [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md): how do the set hits sit on the
  protein results? Protein volcanoes with the fgsea hits ringed.
- [`03_run_singscore`](03_run_singscore/README.md): how high does each set score in each sample?
  A singscore per sample and set, `set_scores.rds`.
- [`04_set_classification`](04_set_classification/README.md): do set scores separate the study
  groups?
- [`05_set_association`](05_set_association/README.md): do set scores and protein slopes track
  the change in muscle size and echo intensity?
- [`06_leading_edge`](06_leading_edge/README.md): do both legs move the same pathways through the
  same proteins? Sam's PR #26 analysis on the fgsea leading edges.
- [`07_pathway_maps`](07_pathway_maps/README.md): how do the training responses cluster? Dotplots,
  enrichment maps, trees and protein networks of both legs.

Every step after `00_` writes one A4 figure PDF named like the step. Figures run S1 to S22: S1
and S2 from `01_`, S3 and S4 from `02_`, S5 from `03_`, S6 to S10 from `04_`, S11 to S13 from
`05_`, S14 to S18 from `06_` and S19 to S22 from `07_`. `04_` and `05_` open with their chance
summary; `04_`, `05_` and `07_` then give the collections in the order Hallmark, GO_Slim,
KEGG_Legacy, Reactome, GOBP, each page titled with its collection.

fgsea is competitive and assumes exchangeable proteins, which co-regulated sets violate; fry is
self-contained and does not, so fgsea calls without fry support are leads. BH runs over all sets
per contrast in `01_`, within collection and task in `04_`, and within collection, analysis and
outcome in `05_`, so their counts answer different families of tests.
