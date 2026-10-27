# 03 · Pathway Enrichment

This stage tests whether proteins that work together moved together, then scores every sample on
every set. The eight steps take about two minutes.

[`00_build_gene_sets`](00_build_gene_sets/README.md) builds the one set list every later step uses,
from a frozen MSigDB snapshot and GO Slim sets built locally, with the protein-to-gene map and the
size filter. [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) tests every set on every
contrast with fgsea and fry, and marks the non-redundant fgsea hits with `collapsePathways`.
[`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) draws protein volcanoes with the
pathway hits ringed. [`03_run_singscore`](03_run_singscore/README.md) scores every
sample on every set and saves the matrix as `set_scores.rds`.
[`04_set_classification`](04_set_classification/README.md) asks how well those scores separate the
study groups, and [`05_set_association`](05_set_association/README.md) whether they track the
phenotype. [`06_leading_edge`](06_leading_edge/README.md) carries Sam's core-enrichment analysis:
whether the pathways both legs move are driven by the same proteins, with his NES scatter of both
legs across every set.  [`07_pathway_maps`](07_pathway_maps/README.md) draws both
training responses as dotplots, enrichment maps, trees and protein networks.

Figures run S1 to S27 through the stage: S1–S7 from `01_`, S8–S9 from `02_`, S10 from `03_`,
S11–S15 from `04_`, S16–S18 from `05_`, S19–S23 from `06_` and S24–S27 from `07_`. Steps `01_` to
`03_` and `06_` write one PDF each; `07_` writes one per collection; `04_` and `05_` write one per
collection plus a chance summary. Each A4 page carries its S-number and a short title; encodings
and findings live in the step READMEs.

## What each step answers

| Step | Question | Answer |
|---|---|---|
| `01_run_fgsea_and_fry` | Which sets move with training? | fgsea 376 BFR and 281 HLRT sets (88 and 68 after collapse); fry 656 and 589. Between the arms fry finds none, so fgsea's between-arm sets are leads only. |
| `04_set_classification` | Does a set's score separate the groups? | 9 sets on BFR training and 49 on HLRT after BH; BFR against HLRT sits at chance. |
| `05_set_association` | Does a set's score change track the phenotype? | No set survives BH. |
| `06_leading_edge` | Do both legs move the same sets, through the same proteins? | Across all 1,990 sets the two legs' NES correlate (Spearman rho 0.83). Leading-edge overlap, untested: across the 71 pathways both legs move, a median 58% of leading-edge proteins are shared (41% to 100%). |
| `07_pathway_maps` | How do the pathways each leg moves cluster, and which proteins link them? | Figures only; the main sets of `01_` drawn as maps per collection. |

## Methods

fgsea is a competitive test. It asks whether a set sits toward one end of the protein ranking and
assumes exchangeable proteins, which co-regulated sets violate. fry is self-contained. It asks
whether the set moved at all under the fitted design and makes no such assumption. Both run on the
same 1,990 sets and both ship.

`01_` adjusts BH over all 1,990 sets in a contrast, for fgsea and fry alike. `04_` adjusts within
each collection and task, and `05_` within each collection, analysis and outcome, so a count from
`01_` and one from `04_` or `05_` answer different families of tests.
