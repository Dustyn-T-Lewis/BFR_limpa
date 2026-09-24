# 01_run_fgsea_and_fry

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| Reads | `gene_sets.rds`, `00_build_gene_sets.xlsx`, `fit.rds`, `design.rds`, `proteins.rds` |
| Writes | `01_run_fgsea_and_fry.xlsx`, `01_run_fgsea_and_fry_figures.pdf` (S1–S7) |

`topTable()` rebuilds the five contrasts from the saved fit, keeping `02_Differential`'s BH.
fgsea ranks proteins by moderated t, seeded. `collapsePathways` re-tests each significant set
against a stronger set's leading edge; survivors carry `main = TRUE`. fry reads `proteins$E` with
limpa's weights and no `block`, since participant is fixed in the design (residual within-leg
correlation 0.027).

## Outputs

The workbook's `set_tests` sheet has one row per set, contrast and method. `NES`, `leadingEdge`
and `main` are fgsea-only. `leadingEdge` holds `;`-joined gene symbols, which
`02_enrich_volcano_fgsea` matches to point labels. `protein_results` holds every protein and
contrast from the saved fit.

S1 is the ten strongest collapse survivors per contrast across all collections, S2–S6 the same
within each collection (KEGG_Legacy and Reactome have no interaction survivor, so three panels),
and S7 the significant sets per collection before and after collapse.
