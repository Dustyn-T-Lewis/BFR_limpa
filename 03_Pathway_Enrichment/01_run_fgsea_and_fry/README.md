# 03_Pathway_Enrichment / 01_run_fgsea_and_fry

Tests every set on every contrast with fgsea and fry, and marks the non-redundant fgsea hits.

- Reads: `00_build_gene_sets/c_data/gene_sets.rds`, `00_build_gene_sets.xlsx` (`set_catalog`,
  `protein_gene_map`), `02_Differential_Expression/02_Differential/c_data/fit.rds`,
  `02_Differential_Expression/01_Design/c_data/design.rds`,
  `01_Preprocess/02_Quantification/c_data/proteins.rds`
- Writes: `c_data/01_run_fgsea_and_fry.xlsx` (`overview`, `set_summary`, `protein_summary`,
  `set_tests`, `protein_results`; `02_` and `03_` read the last two),
  `b_reports/01_run_fgsea_and_fry_figures.pdf` (S1–S7)
- Run: `Rscript 03_Pathway_Enrichment/01_run_fgsea_and_fry/a_script/01_run_fgsea_and_fry.R`,
  34 seconds.

`topTable()` rebuilds the five contrasts from the saved fit, keeping `02_Differential`'s BH.
fgsea ranks proteins by moderated t, seeded. fry reads `proteins$E` with limpa's weights and no
`block`, since participant is fixed in the design (residual within-leg correlation 0.027).

`set_tests` has one row per set, contrast and method. `NES`, `leadingEdge` and `main` are
fgsea-only. `leadingEdge` holds `;`-joined gene symbols, which `02_enrich_volcano_fgsea` matches
to point labels.

## Collapse runs after testing, not before

All sets are tested and BH corrects within each contrast. `collapsePathways` then re-tests each
significant set against a stronger set's leading edge and marks survivors `main = TRUE`; it deletes
nothing. Deduplicating before testing lost two thirds of the discoveries (363 to 133 on BFR
training): a redundancy filter removes significant sets, not hopeless ones (Bourgon et al., PNAS
2010).

## fry finds nothing on any between-leg contrast

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 376 | 88 | 656 |
| HLRT_Post-Pre | 281 | 68 | 589 |
| BFR_Post-HLRT_Post | 159 | 52 | 0 |
| Modality_x_Time_Interaction | 26 | 10 | 0 |
| BFR_Pre-HLRT_Pre *(control)* | 1 | 1 | 0 |

fry's zeros on the between-leg contrasts make the interaction's 10 fgsea pathways leads, not
findings. The control's one fgsea hit,
`REACTOME_STRIATED_MUSCLE_CONTRACTION`, is the floor every count should be read against.

S1 shows the ten strongest collapse survivors per contrast across all collections, S2–S6 the same
within each collection (KEGG_Legacy and Reactome have no interaction survivor, so three panels),
and S7 the significant sets per collection before and after collapse.
