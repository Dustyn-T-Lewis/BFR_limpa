# 03_Pathway_Enrichment / 01_Gene_Sets

Freezes human gene sets, maps proteins onto genes, and draws the protein volcanoes. It runs no
pathway test.

| | |
|---|---|
| **Script** | `a_script/01_gene_sets.qmd` |
| **Reads** | `proteins.rds`, `fit.rds`, `design.rds` |
| **Writes** | `c_data/gene_sets.rds`, `pathway_inputs.rds`, `01_gene_sets.xlsx`, three CSVs |
| **Figures** | `b_reports/figures/`: five FDR volcanoes and two pi-ranked views, PNG and PDF |
| **Frozen source** | `c_data/cache/msigdb_2026.1.Hs_Hallmark-Reactome-GOBP.rds` and its `.md5` |

```sh
quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd --output-dir ../b_reports
```

## The freeze

Set membership changes between MSigDB releases, so one release is pinned and kept. The first
render fetches human MSigDB 2026.1.Hs through `msigdbr` and writes a snapshot with an md5 beside
it; every render after that reads the snapshot and checks the md5, so it needs neither the network
nor `msigdbr` installed. The filename carries the release and the collections, which are what make
one snapshot different from another.

Do not rebuild it casually. Restore the RDS and its checksum together when a check fails, and move
both aside first if a rebuild is genuinely wanted.

## Which sets are kept

A set needs 15 to 500 source genes and at least 15 of them measured here. The second bar counts
proteins this experiment actually detected, which is what governs power.

Many surviving sets then say almost the same thing. Within each database the stage walks the
qualifying sets from largest measured membership down and drops any whose Jaccard overlap with a
set already kept reaches 0.5. Databases reduce separately, so a GO term cannot erase a Hallmark set
covering the same biology. Every dropped set records which set displaced it and by how much.

1,837 sets qualify and 1,040 survive reduction. That cut is structural and fixed in advance: it
reads no p-value, and it does not make what survives independent.

## One protein per gene

Set testing needs one row per gene and the protein matrix does not supply that cleanly. A protein
with no gene symbol, or with several, is left out rather than split, since splitting would invent
measurements nobody made. Where several proteins share a symbol, the one with the most observed
precursors per sample represents it. That choice is made once, for all five contrasts, and reads no
fold change and no p-value.

Every protein stays in the matrix and in the volcanoes regardless. Only the gene-level view drops
them, and the `protein_gene_map` sheet records every decision.

## The volcanoes

[enrichVolcano](https://github.com/Dustyn-T-Lewis/enrichVolcano) draws them from an empty
enrichment table, so the protein volcano appears and no untested pathway is presented as enriched.

All five contrasts colour on BH FDR < 0.05, with no fold-change floor added. The two training
contrasts also get a view whose labels are ranked by pi, keeping the same FDR colours and counts —
so a pi label can land on a protein that never passed FDR. Pi ranks and selects nothing.

The package rescales each panel to its own data, so one cloud looking taller than its neighbour
means nothing. Read the numbers from the tables.

## Next-stage contract

`pathway_inputs.rds` holds what this stage derives, and nothing it merely read:

```r
x <- readRDS("03_Pathway_Enrichment/01_Gene_Sets/c_data/pathway_inputs.rds")
x$protein_results # All tested proteins, all contrasts, FDR and pi scores
x$protein_map     # Gene representative chosen per protein
x$gene_sets       # Retained sets of gene symbols
x$set_indices     # Row indices into proteins$E, not into anything built here
x$provenance      # Input checksums, settings, package versions, cache checksum
```

The matrix and the fitted model stay where the upstream stages wrote them. A set test opens
`01_Preprocess/02_Quantification/c_data/proteins.rds` and
`02_Differential_Expression/02_Differential/c_data/fit.rds` itself, the way every stage in this
repo reads its upstream. Copying them here would only let the copies go stale.

Whatever runs next has to keep the participant design and limpa's uncertainty, and it should report
the baseline control before anyone reads the interaction. The upstream approximation carries
forward: cyclic loess moved the abundances and left the standard errors behind.

Dependencies are the repo's usual set plus `enrichVolcano`, which needs `volcano_ring()` and its
`volc_sig_col` argument (verified with 0.3.0.9000). `msigdbr` is needed only to create a missing
snapshot.
