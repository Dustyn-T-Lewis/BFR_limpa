# 03 · Pathway Enrichment

Tests whether proteins that work together moved together, then scores every sample on every set.

| Step | Runs | Writes |
|---|---|---|
| [`00_build_gene_sets`](00_build_gene_sets/README.md) | frozen MSigDB snapshot, protein-to-gene map, size filter, GO Slim sets | `gene_sets.rds`, workbook |
| [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) | fgsea and fry per contrast, `collapsePathways` | workbook, S1–S7 |
| [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) | protein volcanoes with pathway rings | S8–S9 |
| [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md) | BFR NES against HLRT NES | workbook, S10–S11 |
| [`04_run_singscore`](04_run_singscore/README.md) | per-sample set scores | `set_scores.rds`, workbook, S12 |
| [`05_classify_and_associate_sets`](05_classify_and_associate_sets/README.md) | set classification and phenotype association | workbook, S13–S19 over 128 pages |

The six steps ran in 104 seconds, 49 of them in `05_`. Figures are numbered S1 to S19 through
the stage. Each A4 page carries a caption beneath it in the supplement style: title, one entry per
lettered panel, encodings, and the workbook sheet the data come from. Findings live in these
READMEs, not on the figures.

## fgsea and fry answer different questions

fgsea is competitive: does a set sit at one end of the protein ranking? fry is self-contained:
did the set move at all under the fitted design? Both run on the same 1,990 sets and both ship.
fgsea assumes proteins are exchangeable, and co-regulated sets break that; fry does not assume it.
