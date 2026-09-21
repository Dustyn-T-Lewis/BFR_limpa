# 01_run_fgsea

Ranks proteins by moderated t and asks fgsea whether a set piles up at one end, then lets
`collapsePathways` prune the redundant hits.

| | |
|---|---|
| **Script** | `a_script/01_run_fgsea.R` |
| **Reads** | `gene_sets.rds`, `fit.rds`, `design.rds`, `config.yml` |
| **Writes** | `c_data/fgsea.rds`, `01_run_fgsea.xlsx`, `fgsea_results.csv` |

## What runs

`topTable()` rebuilds all five contrasts from the saved fit. BH stays as `02_Differential` applied
it, one adjustment per contrast; re-adjusting on the mapped subset would change what every FDR
means. `pi_score` is the Xiao et al. (2014, PMID 22321699) score — it controls no error rate and
selects nothing.

fgsea then runs on all 2,732 sets per contrast, seeded from `config.yml` so a rerun reproduces.
`collapsePathways` re-runs each significant set conditioned on a more significant one's leading
edge and keeps it only if it still stands alone. The `main` column marks survivors.

## Read the control first

The script prints the negative control before any other result. `BFR_vs_HLRT_at_T1` compares two
legs of one person before either was trained, so it cannot contain signal — and fgsea calls 14
sets significant there, 6 after collapse. `../README.md` carries the full table and what it does
and does not contaminate.

## Handoff

```r
fg <- readRDS("03_Pathway_Enrichment/01_run_fgsea/c_data/fgsea.rds")
fg$protein_results # all proteins, all contrasts, BH FDR and pi_score, plus the plot label
fg$fgsea_results   # all contrasts; `main` marks what survived collapse
```

`fgsea_results` is keyed on gene symbols, which is what lets `03_enrich_volcano_fgsea` match
`leadingEdge` against its point labels. Rank on anything else and the tick lines vanish without a
warning. `leadingEdge` is a list column in the RDS and a `;`-joined string in the workbook and CSV.
