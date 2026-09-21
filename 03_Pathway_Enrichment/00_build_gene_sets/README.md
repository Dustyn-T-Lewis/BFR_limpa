# 00_build_gene_sets

Freezes one MSigDB release, maps the protein matrix onto gene symbols, and keeps the sets large
enough to test. `01_run_fgsea` and `02_run_singscore` both read the result, so membership is
decided once.

| | |
|---|---|
| **Script** | `a_script/00_build_gene_sets.R` |
| **Reads** | `proteins.rds`, `config.yml` |
| **Writes** | `c_data/gene_sets.rds`, `00_build_gene_sets.xlsx`, `c_data/cache/` |

## The freeze

The first run fetches ten collections through `msigdbr` and writes a snapshot with an md5 beside
it. Every run after that reads the snapshot and verifies the checksum, needing neither the network
nor `msigdbr`. The filename carries the release and the collection list.

Do not rebuild it casually. Restore the RDS and its checksum together when a check fails, and move
both aside first if a rebuild is genuinely wanted.

## Which sets qualify

14,636 frozen. Size is the only filter: 15 to 500 source genes, at least 15 measured here. The
measured bar does most of the work, because the universe is 3,039 genes.

| Collection | In MSigDB | Qualifying |
|---|---:|---:|
| Hallmark | 50 | 41 |
| KEGG_Legacy | 186 | 94 |
| KEGG_Medicus | 658 | 25 |
| Reactome | 1,839 | 430 |
| WikiPathways | 925 | 222 |
| PID | 196 | 67 |
| BioCarta | 292 | 11 |
| GO:BP | 7,538 | 1,366 |
| GO:CC | 1,080 | 231 |
| GO:MF | 1,872 | 245 |
| **Total** | **14,636** | **2,732** |

## One protein per gene

Set tests need one row per gene. A protein with no symbol, or several, is dropped rather than
split, because splitting invents measurements nobody made. Where several share a symbol, the one
with the most observed precursors represents it. Decided once for all five contrasts, reading no
fold change and no p-value.

Every protein stays in the matrix and in the volcano data; only the gene-level view drops them,
and `protein_gene_map` records every decision.

## Themes

Each GO set carries its nearest `GO.db` ancestor holding at least 300 genes. Striated muscle
contraction rolls to *muscle contraction*, oxidative phosphorylation to *cellular respiration*.
Grouping for reading only: every set is still tested alone, and non-GO collections get no theme.

## Handoff

```r
gs <- readRDS("03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds")
gs$sets          # 2,732 tested sets, as measured gene symbols
gs$set_catalog   # all 14,636, with size, qualification and GO theme
gs$protein_map   # gene representative per protein, plus the plot label
gs$gene_universe # the 3,039 measured symbols
```

`msigdbr` is needed only to create a missing snapshot.
