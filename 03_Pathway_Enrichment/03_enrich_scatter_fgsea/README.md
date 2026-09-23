# 03_enrich_scatter_fgsea

Puts the two training contrasts on one pair of axes. Whether blood flow restriction does something
high-load training does not is then a matter of looking, rather than an inference from an empty
interaction term.

| | |
|---|---|
| **Script** | `a_script/03_enrich_scatter_fgsea.R` |
| **Reads** | `set_tests.rds` |
| **Writes** | `b_reports/01_nes_scatter.png` and `.pdf`, `c_data/03_enrich_scatter_fgsea.xlsx`, `nes_scatter.csv` |

## What it shows

fgsea scores every set once per contrast, so each set has a `BFR_Post-Pre` NES and an
`HLRT_Post-Pre` NES. Plotting one against the other puts every set on the identity line if the
two modalities agree and off it if they do not.

**They agree.** Across 100 sets, Spearman rho is **0.82**, and **every one of the 25 significant
sets falls on the same side of zero in both contrasts**. 8 are significant in both, 14 in BFR
only, 3 in HLRT only. The interaction contrast returns almost nothing because there is no
interaction to find: both arms moved the same biology by the same amount.

Oxidative phosphorylation, mitochondrial gene expression and generation of precursor metabolites
rise in both arms. Muscle system process, myogenesis and mitotic cell cycle fall in both.

## Hallmark and GO Slim only

GO:BP's 1,366 sets are largely restatements of one another. They would bury the plot and pull the
rank correlation toward whichever branch of the ontology is densest. Hallmark and GO Slim are the
two collections whose members do not nest, which is what makes a correlation over them mean
something. Every set is still tested; this figure selects what to draw.

## Handoff

`nes_scatter.csv` carries one row per set with both NES values, both adjusted p-values and which
contrast it reached significance in. The workbook adds the concordance row and the quadrant counts.
