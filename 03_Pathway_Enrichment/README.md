# 03 · Pathway Enrichment

This stage tests whether proteins that work together moved together, then scores every sample on
every set. The seven steps take about 100 seconds.

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
phenotype.

Figures run S1 to S18 through the stage: S1–S7 from `01_`, S8–S9 from `02_`, S10 from `03_`, S11
from `04_`, S12–S16 from `05_` and S17–S18 from `06_`, with one PDF per collection in each. Each A4
page has a supplement-style caption with the title, one entry per lettered panel, the encodings and
the source sheet. Findings live in these READMEs, not on the figures.

## Methods

fgsea is competitive: it asks whether a set sits toward one end of the protein ranking, and assumes
exchangeable proteins, which co-regulated sets violate. fry is self-contained: it asks whether the
set moved at all under the fitted design, with no such assumption. Both run on the same 1,990 sets
and both ship.
