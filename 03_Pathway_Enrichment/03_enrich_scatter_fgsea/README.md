# 03_enrich_scatter_fgsea

This step plots each set's BFR_Post-Pre NES against its HLRT_Post-Pre NES, from `set_tests` in
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`. `c_data/03_enrich_scatter_fgsea.xlsx` holds
the concordance table below, the discordant sets and both NES with adjusted p for every set; S10–S11
are in `b_reports/03_enrich_scatter_fgsea_figures.pdf`.

```sh
Rscript 03_Pathway_Enrichment/03_enrich_scatter_fgsea/a_script/03_enrich_scatter_fgsea.R
```

| Population | Sets | rho | Significant | Discordant |
|---|---:|---:|---:|---:|
| All collections | 1,990 | 0.83 | 418 | 8 |
| Collapse survivors | 126 | 0.86 | 126 | 3 |
| Hallmark and GO Slim | 100 | 0.82 | 25 | 0 |

No set significant in both arms changes sign: both modalities moved the same pathways by similar
amounts, which is why the interaction is empty. A set is discordant when its two NES have opposite
signs; all eight are significant in one arm only, with the other at p 0.43 to 0.99.

S10 shows all sets, collapse survivors and discordant sets. S11 shows Hallmark and GO Slim, whose
sets do not nest, then each concordant quadrant rescaled so every significant set is named. GO:BP's
overlapping branches would weight rho toward the largest branch.