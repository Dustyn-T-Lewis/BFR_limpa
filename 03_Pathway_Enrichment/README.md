# 03 · Pathway Enrichment

**Planned.** This folder holds this README and nothing else.

Protein-by-protein testing asks three thousand separate questions and answers each on its own.
Pathway analysis asks a different one: did a group of proteins that work together move together?
A small shift shared across forty mitochondrial proteins is invisible one protein at a time and
obvious as a set. That makes this the right stage to run when the protein-level result is thin.

```
02_Differential_Expression results + the protein table
  01_Gene_Sets  freeze a set collection    -> a cached, version-stamped set list
  02_Set_Tests  test whole sets            -> per-set p-values and per-sample scores
  03_Enrichment ask what the hits share    -> enrichment tables
```

| Sub-stage | Will do |
|---|---|
| `01_Gene_Sets` | download and freeze the gene-set collection, filtered on size and on overlap between sets |
| `02_Set_Tests` | test each set as a whole, and score every set in every sample |
| `03_Enrichment` | for comparisons that produced a hit list, ask which processes those hits belong to |

## What each part would tell us

**The frozen set collection** is not a result. It exists so every number after it is
reproducible. Set membership changes between database releases, so rebuilding it silently would
move every downstream answer. It gets cached and stamped with the release used, and it is the one
data file worth committing.

Sets are filtered on size and on how much they overlap each other. Without the overlap filter,
one finding reappears under six names and looks like six findings.

**Testing whole sets** asks whether the proteins in a set moved as a group. It needs no single
protein to be significant on its own, so it is the analysis most likely to say something about the
interaction. A result here would read: restriction and heavy load changed this process
differently, even though no individual protein survived correction.

**Scoring sets per sample** turns the protein table into a few thousand set scores, one per
sample. Because they are per-sample, they can be plotted against the muscle measurements in
`00_Input/phenotype.csv`. A set whose score tracks the change in cross-sectional area is worth
more than a set with a small p-value, because it connects the proteomics to something we measured
on the leg.

**Enrichment on the hit lists** only works where there are hits. Where a comparison found none,
this cannot run, and lowering the threshold to manufacture a list is not an option.

## Two tests we will not use

One popular method assumes proteins vary independently, which mitochondrial and ribosomal sets
plainly do not, and it reports significant sets on the negative control. Another handles
correlation but silently ignores the argument that tells it participants contributed four samples
each. Both were tried in pilot work and both failed the control. Whatever runs here has to accept
the repeated-measures structure and be checked against the control first.

## What would make this convincing

A set flagged by both set-level methods, whose per-sample score also tracks a muscle measurement.
Three independent lines pointing the same way survives a reviewer who does not trust any single
p-value.

A null result is also an answer. It would mean the two training modes produced the same response
at this dose, which is worth reporting plainly.

## Packages

A gene-set source such as msigdbr, limma for the set-level tests, and singscore for the
per-sample scores. Package versions are not pinned, which is exactly why the set collection is
cached.

## Cost

Minutes once the sets are cached.
