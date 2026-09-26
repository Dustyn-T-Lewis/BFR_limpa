# 03_enrich_scatter_fgsea

Plots each set's BFR_Post-Pre NES against its HLRT_Post-Pre NES.

- Reads: `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx` (`set_tests`)
- Writes: `c_data/03_enrich_scatter_fgsea.xlsx` (`overview`, `concordance`, `discordant`,
  `nes_scatter`), `b_reports/03_enrich_scatter_fgsea_figures.pdf` (S10–S11)
- Run: `Rscript 03_Pathway_Enrichment/03_enrich_scatter_fgsea/a_script/03_enrich_scatter_fgsea.R`,
  3 seconds.

## Results

| Population | Sets | rho | Significant | Discordant |
|---|---:|---:|---:|---:|
| All collections | 1,990 | 0.83 | 418 | 8 |
| Collapse survivors | 126 | 0.86 | 126 | 3 |
| Hallmark and GO Slim | 100 | 0.82 | 25 | 0 |

No set significant in both arms differs in sign. Both modalities moved the same pathways by
similar amounts, which is why the interaction is empty.

S10: all sets, collapse survivors, discordant sets. S11: Hallmark and GO Slim, then each
concordant quadrant rescaled so every significant set is named. These two collections are used
because their sets do not nest; GO:BP's overlapping branches would weight rho toward whichever
branch is largest.

A set is discordant when its two NES differ in sign. All eight are significant in one arm only,
with the other arm at p 0.43 to 0.99.
