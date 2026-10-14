# 03 · NES scatter, fgsea

Compares each set's BFR_Post-Pre NES with its HLRT_Post-Pre NES, to ask whether both modalities
move the same biology. No new test.

Inputs: `set_tests` from `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`.

```sh
quarto render 03_Pathway_Enrichment/03_enrich_scatter_fgsea/a_script/03_enrich_scatter_fgsea.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/03_enrich_scatter_fgsea.html`: the walkthrough.
- `b_reports/03_enrich_scatter_fgsea_figures.pdf`: S10, the NES scatter for Hallmark and GO Slim.
- `c_data/03_enrich_scatter_fgsea.xlsx`: the concordance table below, the discordant sets, and
  both NES with adjusted p for every set.

| Population | Sets | rho | Significant | Discordant |
|---|---:|---:|---:|---:|
| All collections | 1,990 | 0.83 | 418 | 8 |
| Collapse survivors | 126 | 0.86 | 126 | 3 |
| Hallmark and GO Slim | 100 | 0.82 | 25 | 0 |

A set is discordant when its two NES have opposite signs. None significant in both arms changes
sign: both modalities moved the same pathways by similar amounts.
