# 01_run_fgsea_and_fry

This step tests every set on every contrast with fgsea and fry and marks significant fgsea sets not
redundant with a stronger one. It reads `00_build_gene_sets/c_data/gene_sets.rds` and
`00_build_gene_sets.xlsx`, `02_Differential_Expression/02_Differential/c_data/fit.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds` and
`01_Preprocess/02_Quantification/c_data/proteins.rds`.

`c_data/01_run_fgsea_and_fry.xlsx` holds `set_tests`, one row per set, contrast and method, and
`protein_results`, every protein and contrast from the fit; `02_` and `03_` read both. Figures S1–S7
are in `b_reports/01_run_fgsea_and_fry_figures.pdf`.

```sh
Rscript 03_Pathway_Enrichment/01_run_fgsea_and_fry/a_script/01_run_fgsea_and_fry.R
```

`topTable()` rebuilds the five contrasts from the saved fit with `02_Differential`'s BH. fgsea ranks
proteins by moderated t, with a fixed seed. fry reads `proteins$E` with limpa's weights and no
`block`, since participant is fixed in the design (residual within-leg correlation 0.027). `NES`,
`leadingEdge` and `main` are fgsea-only; `leadingEdge` holds `;`-joined gene symbols that
`02_enrich_volcano_fgsea` matches to point labels.

Every set is tested, with BH within each contrast. `collapsePathways` then re-tests each significant
set against a stronger set's leading edge and marks survivors `main = TRUE`, deleting nothing.
Deduplicating before testing lost two thirds of the discoveries (363 to 133 on BFR training),
because a redundancy filter removes significant sets rather than hopeless ones (Bourgon et al., PNAS
2010).

Significant sets per contrast at an adjusted p below 0.05:

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 376 | 88 | 656 |
| HLRT_Post-Pre | 281 | 68 | 589 |
| BFR_Post-HLRT_Post | 159 | 52 | 0 |
| Modality_x_Time_Interaction | 26 | 10 | 0 |
| BFR_Pre-HLRT_Pre *(control)* | 1 | 1 | 0 |

fry finds nothing on the between-leg contrasts, so the interaction's 10 fgsea pathways are leads
rather than findings. The control's one fgsea hit, `REACTOME_STRIATED_MUSCLE_CONTRACTION`, is the
floor for every count.

S1 shows the ten strongest collapse survivors per contrast across all collections, S2–S6 the same
within each collection (KEGG_Legacy and Reactome have no interaction survivor, so three panels), and
S7 significant sets per collection before and after collapse.