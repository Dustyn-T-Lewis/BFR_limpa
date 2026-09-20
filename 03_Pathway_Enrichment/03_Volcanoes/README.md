# 03_Pathway_Enrichment / 03_Volcanoes

Draws the protein volcanoes with the enriched pathways ringed. Computes nothing.

| | |
|---|---|
| **Script** | `a_script/03_volcanoes.qmd` |
| **Reads** | `01_Gene_Sets/c_data/pathway_inputs.rds`, `design.rds` |
| **Figures** | `b_reports/figures/`: five FDR volcanoes and two pi-ranked views, PNG and PDF |

```sh
quarto render 03_Pathway_Enrichment/03_Volcanoes/a_script/03_volcanoes.qmd --output-dir ../b_reports
```

## Why this is its own stage

These were drawn inside `01_Gene_Sets` when the repo ran no pathway test. With no enrichment
table to pass, the notebook handed `enrichVolcano` an empty one and rebound the package's
validator in a private environment to stop it refusing — about 25 lines of environment surgery,
a version guard, and two provenance fields to document the workaround.

`01_` now runs fgsea, so there is a real table. The package is called as written and all of that
is gone. Splitting the stage was the fix.

## Two claims on one figure

Point colours and the count badges read protein-level BH FDR. The ring reads set-level fgsea FDR.
They are separate claims, and the ring is the weaker one — fgsea assumes the ranked proteins are
exchangeable, and it rings one pathway on the pre-training control, a contrast that cannot contain
signal. `02_Set_Tests` is where the set p-value you would quote comes from, and it finds nothing
on that contrast.

Read the control panel first for exactly that reason.

## The label trap

`volcano_ring()` matches leading-edge genes against the point labels, and a label carries its
accession when a symbol sits on more than one protein. The notebook translates the leading edges
into label space before drawing. Skip that and the tick lines draw nothing, with no warning and no
error — the plot simply comes out bare. All 3,396 leading-edge genes on `BFR_post_vs_pre`
currently translate and match.

## What the panels cannot do

The package rescales each contrast to fill its own plotting area, so one cloud looking taller or
wider than its neighbour says nothing at all. Read exact values from `01_Gene_Sets`'s workbook.
