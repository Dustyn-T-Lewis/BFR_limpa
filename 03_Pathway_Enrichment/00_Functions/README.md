# 03_Pathway_Enrichment / 00_Functions

Nine functions the three pathway stages share. Defining them has no side effect; calling
`stage_inputs()` or `write_stage_outputs()` does.

| | |
|---|---|
| **File** | `00_functions.R` |
| **Sourced by** | `01_Gene_Sets`, `02_Set_Tests`, `03_Volcanoes` |

```r
source(here::here("03_Pathway_Enrichment", "00_Functions", "00_functions.R"))
```

## Why this exists at all

Every other stage in this repo is self-contained: `a_script/` holds only the notebook, nothing is
sourced, and files on disk are the only handoff. This file breaks that rule once, deliberately,
so the three pathway stages open the same way, close the same way, and cannot each invent their
own provenance record.

| Function | What it is for |
|---|---|
| `stage_inputs()` | resolve repo-relative paths, stop if one is missing, read them, return a manifest |
| `write_stage_outputs()` | one workbook, the RDS files, the CSVs, and one provenance record stamped into all of them |
| `load_sets()` | one or more databases as a single named list, deduped or not |
| `reduce_sets()` | the greedy Jaccard walk |
| `drop_disease_terms()` | split a results table into kept and name-excluded |
| `set_row_indices()` | gene symbols to row numbers of the protein matrix |
| `run_fgsea_all()` | fgsea for every collection against every contrast, seeded |
| `run_fry_all()` | fry for every contrast, one contrast vector at a time |
| `score_samples()` | singscore's per-sample score matrix |

Wrapping cost about as many lines as it saved: the notebooks lost 75 and this file gained 109.
What it bought is uniformity. Every workbook ends with `input_manifest` and `package_versions`
because `write_stage_outputs()` appends them, not because three notebooks each remembered to.

## Why singscore rather than GSVA

`score_samples()` uses singscore because a score must not depend on which other samples are in
the matrix. The design is paired within participant, so a participant's T1 score shifting when
the cohort changes would couple samples the contrast treats as fixed.

Measured on this matrix, 200 sets, dropping 51 of the 131 samples: a singscore value moves by 0,
a GSVA value by up to 0.34 on a range of 1.5. GSVA estimates a kernel density across the samples
present, which is what makes it cohort-relative. Its advantage is that reviewers recognise it.

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
