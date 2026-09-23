# 01_run_fgsea_and_fry

Tests every gene set two ways and reports both, then lets `collapsePathways` prune the redundant
fgsea hits.

| | |
|---|---|
| **Script** | `a_script/01_run_fgsea_and_fry.R` |
| **Reads** | `gene_sets.rds`, `fit.rds`, `design.rds`, `proteins.rds` |
| **Writes** | `c_data/set_tests.rds`, `01_run_fgsea_and_fry.xlsx`, `set_tests.csv`, `b_reports/` |

## What runs

`topTable()` rebuilds all five contrasts from the saved fit. BH stays as `02_Differential` applied
it, one adjustment per contrast; re-adjusting on the mapped subset would change what every FDR
means. `pi_score` is the Xiao et al. (2014, PMID 22321699) score — it controls no error rate and
selects nothing.

**fgsea is competitive**: it asks whether a set piles up at one end of the ranking relative to
every other protein. It runs on all 1,990 sets per contrast, seeded so a rerun reproduces.
`collapsePathways` then re-runs each significant set conditioned on a more significant one's
leading edge and keeps it only if it still stands alone. The `main` column marks survivors.

**fry is self-contained**: it asks whether a set moved at all, rotating the residuals of the fitted
model rather than shuffling gene labels, so correlation cannot inflate the null. It indexes rows of
`proteins$E` through the representative protein chosen in `00_build_gene_sets`, carries limpa's
per-observation weights, and takes no `block` argument because `01_Design` fixes participant in the
design and the residual within-leg correlation is 0.027. Both methods arrive in the same long
table, one row per set per contrast per method.

**Neither is simply better.** fry flags 656 of 1,990 sets on BFR training, because with a global
effect most sets moved a little. fgsea flags 376. Where they disagree is information.

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 376 | 88 | 656 |
| HLRT_Post-Pre | 281 | 68 | 589 |
| BFR_Post-HLRT_Post | 159 | 52 | 0 |
| Modality_x_Time_Interaction | 26 | 10 | 0 |
| BFR_Pre-HLRT_Pre *(control)* | 1 | 1 | 0 |

On the pre-training control fry calls 0 at a best FDR of exactly 1.00 while fgsea still calls 2.
The gap is what a competitive test costs when genes within a set are correlated, and it is not
specific to this implementation: across three independent implementations of fry and four of fgsea
on 996 shared sets, fry returned 0, 0, 0 at a p < 0.05 rate of 0.016 to 0.019, and fgsea returned
3, 3, 1, 14 at a rate of 0.057 to 0.068. Under a true null that rate should be 0.05.

## Read the control first

The script prints the negative control first. `BFR_Pre-HLRT_Pre` compares two legs of one person
before either was trained, so it cannot contain signal. fgsea calls one set significant there,
down from 14 before GO:CC and GO:MF were dropped. fry calls none.

## Handoff

```r
fg <- readRDS("03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds")
fg$protein_results # all proteins, all contrasts, BH FDR and pi_score, plus the plot label
fg$set_tests       # one row per set per contrast per method
fg$set_summary     # one row per contrast: sets, fgsea, collapsed, fry
```

`set_tests` is long: `method` is `fgsea` or `fry`, and both carry `p` and `padj` from the same
BH family, one per method per contrast. `NES`, `leadingEdge` and `main` are fgsea-only and are
`NA` on every fry row, so read them behind `filter(method == "fgsea")`.

`leadingEdge` is keyed on gene symbols, which is what lets `02_enrich_volcano_fgsea` match
`leadingEdge` against its point labels. Rank on anything else and the tick lines vanish without a
warning. `leadingEdge` is a list column in the RDS and a `;`-joined string in the workbook and CSV.

## Figures

`b_reports/` holds one directory per collection plus `all_db/`, one file per contrast inside each.
Every panel shows the ten strongest collapse survivors by adjusted p, the rule the volcano rings
use. Colour is `-log10(FDR)` inside a collection, the collection itself in `all_db/`, which also
carries `02_collapse_before_after.png`.

Labels come from `enrichVolcano::ev_clean_label`, so a dot plot, a volcano and a ROC strip name a
pathway the same way. `KEGG_Legacy/` and `Reactome/` hold three files rather than four, because
neither has a set surviving collapse on the interaction. 23 figures, each as PNG and PDF.
