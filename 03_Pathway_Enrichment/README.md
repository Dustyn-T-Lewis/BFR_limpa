# 03 · Pathway Enrichment

`01_Gene_Sets` is implemented. Whole-set tests, sample scores and hit-list
enrichment remain planned.

Testing protein by protein asks three thousand separate questions. Pathway analysis asks whether
a group of proteins that work together moved together. A small shift shared across forty
mitochondrial proteins is invisible one protein at a time and obvious as a set, which is why this
is the right stage when the protein-level result is thin.

| Sub-stage | Does / will do |
|---|---|
| [`01_Gene_Sets`](01_Gene_Sets/README.md) | freeze human MSigDB; filter on size and measured overlap; export the protein matrix, model and FDR/pi scores; draw enrichVolcano protein plots |
| `02_Set_Tests` | test each set as a whole, and score every set in every sample |
| `03_Enrichment` | for comparisons with hits, ask which processes those hits belong to |

```sh
quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd --output-dir ../b_reports
```

The first stage writes `b_reports/01_gene_sets.html`, PNG/PDF protein volcanoes,
and `c_data/` workbooks/RDS/CSV outputs. Its source cache pins human MSigDB
2026.1.Hs (Hallmark, Reactome, GO Biological Process) and supports offline
rerenders. The handoff includes the full limpa EList, weighted fit, and
participant design alongside the plain matrix.

Protein significance uses BH FDR. Pi-scores rank training-contrast labels;
they do not select hits. enrichVolcano draws the protein plots with no pathway
arcs until a later stage supplies valid enrichment results.

## What each would tell us

**The frozen collection** is not a result. Set membership changes between database releases, so
it gets cached and stamped with the version used. Without an overlap filter, one finding reappears
under six names.

**Testing whole sets** needs no single protein to be significant, so it is the analysis most
likely to say something about the interaction.

**Scoring sets per sample** gives one number per set per sample, which can be plotted against the
muscle measurements in `00_Input/phenotype.csv`. A set tracking the change in cross-sectional area
is worth more than a set with a small p-value.

**Enrichment on hit lists** only works where there are hits, and lowering a threshold to
manufacture one is not an option.

## Two tests we will not use

One assumes proteins vary independently, which mitochondrial and ribosomal sets do not, and it
reported significant sets on the negative control. Another handles correlation but silently
ignores the argument saying participants contributed four samples each. Whatever runs here must
accept the repeated-measures structure and be checked against the control first.
