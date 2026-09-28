# 03 · Pathway Enrichment

This stage tests whether proteins that work together moved together, then scores every sample on
every set. The six steps ran in 104 seconds, 49 of them in `05_`.

[`00_build_gene_sets`](00_build_gene_sets/README.md) builds the one set list every later step
uses, from a frozen MSigDB snapshot and GO Slim sets built locally, with the protein-to-gene map
and the size filter. [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) tests every set on
every contrast with fgsea and fry, and marks the non-redundant fgsea hits with
`collapsePathways`. [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) draws protein
volcanoes with the pathway hits ringed, and [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md)
plots each set's BFR NES against its HLRT NES. [`04_run_singscore`](04_run_singscore/README.md)
scores every sample on every set and saves the matrix as `set_scores.rds`.
[`05_classify_and_associate_sets`](05_classify_and_associate_sets/README.md) asks how well those
scores separate the study groups and whether they track the phenotype.

Figures run S1 to S19 through the stage: S1–S7 from `01_`, S8–S9 from `02_`, S10–S11 from `03_`, S12
from `04_`, and S13–S19 from `05_` over 128 pages. Each A4 page has a supplement-style caption with
the title, one entry per lettered panel, the encodings and the source sheet. Findings live in these
READMEs, not on the figures.

## Methods

fgsea is competitive: it asks whether a set sits toward one end of the protein ranking, and assumes
exchangeable proteins, which co-regulated sets violate. fry is self-contained: it asks whether the
set moved at all under the fitted design, with no such assumption. Both run on the same 1,990 sets
and both ship.