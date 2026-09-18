# 03 · Pathway Enrichment

`01_Gene_Sets` is implemented. Whole-set tests, sample scores and hit-list
enrichment remain planned.

Testing protein by protein asks three thousand separate questions. Pathway analysis asks whether
a group of proteins that work together moved together. A small shift shared across forty
mitochondrial proteins is invisible one protein at a time and obvious as a set, which is why this
is the right stage when the protein-level result is thin.

| Sub-stage | Does / will do |
|---|---|
| [`01_Gene_Sets`](01_Gene_Sets/README.md) | freeze human MSigDB, filter on size and measured overlap, export the protein matrix with its fitted model and scores, draw the protein volcanoes |
| `02_Set_Tests` | test each set as a whole, and score every set in every sample |
| `03_Enrichment` | for comparisons with hits, ask which processes those hits belong to |

```sh
quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd --output-dir ../b_reports
```

That stage runs no pathway test. It freezes the gene sets, hands the next stage the protein
matrix together with the fitted model and the participant design, and draws the protein
volcanoes. The frozen snapshot is human MSigDB 2026.1.Hs, covering Hallmark, Reactome and GO
Biological Process, and every later render reads it from disk rather than the network.

Proteins are called on BH FDR. The pi score ranks labels on the two training contrasts and
selects nothing. The volcanoes carry no pathway ring until a later stage has a tested result to
put in one.

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
