# 03_enrich_scatter_fgsea

Puts the two training contrasts on one pair of axes. Whether blood flow restriction does something
high-load training does not is then a matter of looking, rather than an inference from an empty
interaction term.

| | |
|---|---|
| **Script** | `a_script/03_enrich_scatter_fgsea.R` |
| **Reads** | `set_tests.rds` |
| **Writes** | two composites in `b_reports/`, `c_data/03_enrich_scatter_fgsea.xlsx`, `nes_scatter.csv` |

## What it shows

fgsea scores every set once per contrast, so each set has a `BFR_Post-Pre` NES and an
`HLRT_Post-Pre` NES. Plotting one against the other puts every set on the identity line if the
two modalities agree and off it if they do not.

**They agree, and the agreement strengthens as redundancy is removed.**

| Population | Sets | rho | Significant | Discordant |
|---|---:|---:|---:|---:|
| All collections | 1,990 | 0.83 | 418 | 8 |
| Survived collapsePathways | 126 | **0.86** | 126 | 3 |
| Hallmark and GO Slim | 100 | 0.82 | 25 | **0** |

The interaction contrast returns almost nothing because there is no interaction to find. Oxidative
phosphorylation, mitochondrial gene expression and generation of precursor metabolites rise in
both arms; muscle system process, myogenesis and mitotic cell cycle fall in both.

## Discordant sets

A set is discordant when the two contrasts put it on opposite sides of zero. Eight qualify, and
**none is significant in both arms**: each reaches significance in one arm while the other returns
p between 0.43 and 0.99, so the opposite sign rests on a null. Read them as one-armed results, not
as the two modalities disagreeing. Zero of the 239 sets significant in both arms disagree in sign.

## Two composites

`01_nes_concordance_all.png` carries three panels: every tested set, the 126 that survived
collapse, and the eight discordant sets on their own axes with names. `02_nes_concordance_curated.png`
carries the Hallmark and GO Slim cloud beside the two concordant quadrants, each scaled to its own
sets so all 25 can be named. Both also write PDF.

## Hallmark and GO Slim only

GO:BP's 1,366 sets are largely restatements of one another, so they dominate the cloud and pull
the rank correlation toward whichever branch of the ontology is densest. Hallmark and GO Slim are
the two collections whose members do not nest, which is what makes a correlation over them mean
something. The first composite reports every set anyway, so both views are on record.

## Handoff

`nes_scatter.csv` carries one row per set with both NES values, both adjusted p-values and which
contrast it reached significance in. The workbook adds the concordance row and the quadrant counts.
