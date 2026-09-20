# 03_Pathway_Enrichment / 00_Functions

Four functions the three pathway stages share. Runs nothing, reads nothing, writes nothing.

| | |
|---|---|
| **File** | `00_functions.R` |
| **Sourced by** | `01_Gene_Sets`, `02_Set_Tests`, `03_Volcanoes` |

```r
source(here::here("03_Pathway_Enrichment", "00_Functions", "00_functions.R"))
```

## Why this exists at all

Every other stage in this repo is self-contained: `a_script/` holds only the notebook, nothing is
sourced, and files on disk are the only handoff. This file breaks that rule once, deliberately.
The gene-set loader and the overlap reduction are needed by all three stages, and three copies
would drift apart. Nothing here has a side effect, so sourcing it only defines functions.

A function belongs here when two stages call it. One caller means it stays in the notebook.

## Loading sets

```r
load_sets(gene_sets, dbs, dedup = TRUE, scope = c("across", "within"), cutoff = 0.5)
```

Pass one database or five; they come back as a single named list. `scope` is the part worth
understanding.

`"within"` reduces each database separately, so a GO term cannot displace a Hallmark set covering
the same biology. That is what you want when reporting per database, and it reproduces the
1,040 sets in `gene_sets.rds` exactly.

`"across"` pools the requested databases first, so a Hallmark set and its GO twin compete and the
duplicate story survives once. That is what you want when several databases feed one figure.
Pooling all three collapses 12 more sets — Reactome's mitochondrial translation into GO's
mitochondrial gene expression at Jaccard 0.64, branched-chain amino acid catabolism into its GO
equivalent at 0.86. `dedup = FALSE` returns the raw union.

## Excluding sets by name

`drop_disease_terms()` returns both halves, `kept` and `excluded`, and the excluded rows carry the
pattern that caught them. It returns both because a name filter is a judgement rather than a
statistic: `REACTOME_INFLUENZA_INFECTION` is largely ribosome and translation machinery under a
misleading label, so removing it takes real biology with it. A reader has to be able to see what
went and decide whether that was fair. `patterns = NULL` keeps everything.

## The reduction itself

`reduce_sets()` walks the catalogue from largest measured membership down and drops any set whose
Jaccard overlap with something already kept reaches the cutoff. Ties go to the larger source set,
then the lower set ID, so the answer does not depend on row order. It moved here unchanged from
`01_gene_sets.qmd` and still returns the same 1,040 sets with the same `representative` and
`jaccard` columns — that equality is the test that it moved correctly.
