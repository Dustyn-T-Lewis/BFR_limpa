# 03_enrich_scatter_fgsea

Each set's BFR_Post-Pre NES against its HLRT_Post-Pre NES.

| | |
|---|---|
| Reads | `set_tests.rds` |
| Writes | 2 composites, `nes_scatter.csv`, `03_enrich_scatter_fgsea.xlsx` |

| Population | Sets | rho | Significant | Discordant |
|---|---:|---:|---:|---:|
| All collections | 1,990 | 0.83 | 418 | 8 |
| Collapse survivors | 126 | 0.86 | 126 | 3 |
| Hallmark and GO Slim | 100 | 0.82 | 25 | 0 |

`01_nes_concordance_all`: all sets, collapse survivors, discordant sets.
`02_nes_concordance_curated`: Hallmark and GO Slim, then each concordant quadrant rescaled so every
significant set is named. These two collections are used because their sets do not nest; GO:BP's
overlapping branches would weight rho toward whichever branch is largest.

A set is discordant when its two NES differ in sign. All eight are significant in one arm only,
with the other arm at p 0.43 to 0.99.
