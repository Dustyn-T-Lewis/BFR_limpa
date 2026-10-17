# 03 · Pathway Enrichment

This stage tests whether proteins that work together moved together, then scores every sample on
every set. The nine steps take about two minutes.

[`00_build_gene_sets`](00_build_gene_sets/README.md) builds the one set list every later step uses,
from a frozen MSigDB snapshot and GO Slim sets built locally, with the protein-to-gene map and the
size filter. [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) tests every set on every
contrast with fgsea and fry, and marks the non-redundant fgsea hits with `collapsePathways`.
[`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) draws protein volcanoes with the
pathway hits ringed, and [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md) plots each
set's BFR NES against its HLRT NES. [`04_run_singscore`](04_run_singscore/README.md) scores every
sample on every set and saves the matrix as `set_scores.rds`.
[`05_set_classification`](05_set_classification/README.md) asks how well those scores separate the
study groups, and [`06_set_association`](06_set_association/README.md) whether they track the
phenotype. [`07_leading_edge`](07_leading_edge/README.md) asks whether the pathways both legs move
are driven by the same proteins, and [`08_pathway_maps`](08_pathway_maps/README.md) draws both
training responses as dotplots, enrichment maps, trees and protein networks.

Figures run S1 to S26 through the stage: S1–S7 from `01_`, S8–S9 from `02_`, S10 from `03_`, S11
from `04_`, S12–S16 from `05_`, S17–S19 from `06_`, S20–S22 from `07_` and S23–S26 from `08_`.
Steps `01_` to `04_` and `07_` write one PDF each; `08_` writes one per collection;
`05_` and `06_` write one per collection plus a chance summary. Each A4 page has a supplement-style
caption with the title, one entry per lettered panel, the encodings and the source sheet. Findings
live in these READMEs, not on the figures.

## What each step answers

| Step | Question | Answer |
|---|---|---|
| `01_run_fgsea_and_fry` | Which sets move with training? | fgsea 376 BFR and 281 HLRT sets (88 and 68 after collapse); fry 656 and 589. Between the arms fry finds none, so fgsea's between-arm sets are leads only. |
| `03_enrich_scatter_fgsea` | Do both arms move the same sets? | Yes. Across Hallmark and GO Slim rho is 0.815, and all 25 significant sets move the same way in both. |
| `05_set_classification` | Does a set's score separate the groups? | 9 sets on BFR training and 49 on HLRT after BH; BFR against HLRT sits at chance. |
| `06_set_association` | Does a set's score change track the phenotype? | No set survives BH. |
| `07_leading_edge` | Do both legs move a shared pathway through the same proteins? | Overlap only, untested: across the 71 pathways both legs move, a median 58% of leading-edge proteins are shared (41% to 100%). |
| `08_pathway_maps` | How do the pathways each leg moves cluster, and which proteins link them? | Figures only; the main sets of `01_` drawn as maps per collection. |

## Methods

fgsea is a competitive test. It asks whether a set sits toward one end of the protein ranking and
assumes exchangeable proteins, which co-regulated sets violate. fry is self-contained. It asks
whether the set moved at all under the fitted design and makes no such assumption. Both run on the
same 1,990 sets and both ship.

`01_` adjusts BH over all 1,990 sets in a contrast, for fgsea and fry alike. `05_` adjusts within
each collection and task, and `06_` within each collection, analysis and outcome, so a count from
`01_` and one from `05_` or `06_` answer different families of tests.
