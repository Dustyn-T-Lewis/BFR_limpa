# 01_run_fgsea_and_fry

Tests every set on every contrast with fgsea and fry, and marks non-redundant fgsea hits.

| | |
|---|---|
| **Reads** | `gene_sets.rds`, `fit.rds`, `design.rds`, `proteins.rds` |
| **Writes** | `set_tests.rds`, `set_tests.csv`, `01_run_fgsea_and_fry.xlsx`, 23 dot plots |

`topTable()` rebuilds the five contrasts from the saved fit, keeping `02_Differential`'s BH.
fgsea ranks proteins by moderated t, seeded. `collapsePathways` re-tests each significant set
against a stronger set's leading edge; survivors carry `main = TRUE`. fry reads `proteins$E` with
limpa's weights and no `block`, since participant is fixed in the design (residual within-leg
correlation 0.027).

The script prints the pre-training control first.

## Outputs

`set_tests` has one row per set, contrast and method. `NES`, `leadingEdge` and `main` are
fgsea-only. `leadingEdge` holds gene symbols, which `02_enrich_volcano_fgsea` matches to point
labels; it is a list column in the RDS and `;`-joined in the CSV.

`b_reports/` has one folder per collection plus `all_db/`, one dot plot per contrast: the ten
strongest collapse survivors. `KEGG_Legacy/` and `Reactome/` have three, with no interaction
survivor. `all_db/` adds `02_collapse_before_after`.
