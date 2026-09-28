# 01_run_fgsea_and_fry

This step tests every set on every contrast with fgsea and fry, and marks which significant fgsea
sets are not redundant with a stronger one. It reads the sets from
`00_build_gene_sets/c_data/gene_sets.rds` and the catalogue and protein map from
`00_build_gene_sets.xlsx`, the fit from `02_Differential_Expression/02_Differential/c_data/fit.rds`,
the design from `02_Differential_Expression/01_Design/c_data/design.rds`, and the protein table
from `01_Preprocess/02_Quantification/c_data/proteins.rds`.

Results go to `c_data/01_run_fgsea_and_fry.xlsx`. Its `set_tests` sheet has one row per set,
contrast and method, and `protein_results` has every protein and contrast from the fit. `02_` and
`03_` read those two sheets. Figures S1–S7 are in `b_reports/01_run_fgsea_and_fry_figures.pdf`.

```sh
Rscript 03_Pathway_Enrichment/01_run_fgsea_and_fry/a_script/01_run_fgsea_and_fry.R
```

`topTable()` rebuilds the five contrasts from the saved fit, keeping `02_Differential`'s BH. fgsea
ranks proteins by moderated t, with a fixed seed. fry reads `proteins$E` with limpa's weights and
no `block`, because participant is already a fixed term in the design (the residual within-leg
correlation is 0.027). In `set_tests`, `NES`, `leadingEdge` and `main` are filled for fgsea only.
`leadingEdge` holds the gene symbols joined by `;`, which `02_enrich_volcano_fgsea` matches to
point labels.

Every set is tested, and BH corrects within each contrast. `collapsePathways` then re-tests each
significant set against a stronger set's leading edge and marks the survivors `main = TRUE`; it
deletes nothing. Deduplicating the sets before testing lost two thirds of the discoveries (363 to
133 on BFR training), because a redundancy filter removes significant sets rather than hopeless
ones (Bourgon et al., PNAS 2010).

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
floor every count should be read against.

S1 shows the ten strongest collapse survivors in each contrast across all collections, and S2–S6
do the same within each collection. KEGG_Legacy and Reactome have no interaction survivor, so their
pages have three panels. S7 compares the significant sets per collection before and after
collapse.
