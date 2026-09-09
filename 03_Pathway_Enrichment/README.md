# 03 · Pathway Enrichment

**Planned.** This folder holds this README.

Testing protein by protein asks three thousand separate questions. Pathway analysis asks whether
a group of proteins that work together moved together. A small shift shared across forty
mitochondrial proteins is invisible one protein at a time and obvious as a set, which is why this
is the right stage when the protein-level result is thin.

| Sub-stage | Will do |
|---|---|
| `01_Gene_Sets` | download and freeze a gene-set collection, filtered on size and overlap |
| `02_Set_Tests` | test each set as a whole, and score every set in every sample |
| `03_Enrichment` | for comparisons with hits, ask which processes those hits belong to |

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
