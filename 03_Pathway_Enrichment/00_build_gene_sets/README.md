# 00_build_gene_sets

Freezes one MSigDB release, maps the protein matrix onto gene symbols, and keeps the sets large
enough to test. `01_run_fgsea_and_fry` and `04_run_singscore` both read the result, so membership is
decided once.

| | |
|---|---|
| **Script** | `a_script/00_build_gene_sets.R` |
| **Reads** | `proteins.rds` |
| **Writes** | `c_data/gene_sets.rds`, `00_build_gene_sets.xlsx`, `c_data/cache/` |

## The freeze

The first run fetches four collections through `msigdbr` and writes a snapshot with an md5 beside
it. Every run after that reads the snapshot and verifies the checksum, needing neither the network
nor `msigdbr`. The filename carries the release and the collection list.

Do not rebuild it casually. Restore the RDS and its checksum together when a check fails, and move
both aside first if a rebuild is genuinely wanted.

## Which sets qualify

9,613 frozen from MSigDB, plus 71 GO Slim sets built here. Size is the only filter: 15 to 500 source genes, at least 15 measured here. The
measured bar does most of the work, because the universe is 3,039 genes.

| Collection | In source | Qualifying | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 41 | 44 |
| KEGG_Legacy | 186 | 94 | 26 |
| Reactome | 1,839 | 430 | 33 |
| GO:BP | 7,538 | 1,366 | 29 |
| GO Slim *(built here)* | 71 | 59 | 138 |
| **Total** | **9,684** | **1,990** | |

Six collections were dropped after measuring what each contributed. GO:CC and GO:MF cost the most
to keep: GO:CC is the ontology of physical complexes, whose members are co-regulated by
stoichiometry, which is exactly the between-gene correlation a competitive test assumes away.
Three of the six sets fgsea called significant on the pre-training control came from GO:CC.
Dropping the two extra ontologies takes the control from 14 significant to 2. PID has been retired
since 2012 and returned nothing; BioCarta qualified 11 sets; KEGG_Medicus is disease-centric and
qualified 25 of 658. WikiPathways duplicates Reactome without adding a curation standard of its
own.

## One protein per gene

Set tests need one row per gene. A protein with no symbol, or several, is dropped rather than
split, because splitting invents measurements nobody made. Where several share a symbol, the one
with the most observed precursors represents it. Decided once for all five contrasts, reading no
fold change and no p-value.

Every protein stays in the matrix and in the volcano data; only the gene-level view drops them,
and `protein_gene_map` records every decision.

## Themes

Each GO:BP set carries a term from the GO Consortium's generic slim, a curated 140-term list
frozen here as `c_data/cache/goslim_generic.obo` with an md5 checked on every run. The rule is the
one `GSEABase::goSlim` uses: a term belongs to every slim term it descends from, and the most
specific of those wins. 836 of the 1,366 GO:BP sets carry a slim label; the rest descend from no slim term and carry
`NA`.

| Theme | Sets |
|---|---:|
| anatomical structure development | 90 |
| signaling | 85 |
| immune system process | 66 |
| cell differentiation | 57 |
| lipid metabolic process | 39 |

Each slim term is also **a set**: every measured gene annotated to it or to any term beneath it in
GO:BP. 59 of the 71 terms with a measured gene pass the 15-to-500 rule.

Two details decide whether these sets are right. They are built from the frozen membership, not
from the qualifying subset, so a gene is not lost when the GO:BP set carrying it falls outside the
size filter. And the rule is read on measured size, because a slim term is broad by design and the
question is whether it is testable here. Built the other way, *protein folding* held 19 genes
instead of 129.

The `theme` column stays on every GO:BP row, so a result among the 1,366 traces to its slim
parent. A set gets one label, so the most specific ancestor takes it; the slim sets themselves use
every ancestor, which is what map2slim means.

## Handoff

```r
gs <- readRDS("03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds")
gs$sets          # 1,990 tested sets, as measured gene symbols
gs$set_catalog   # all 9,684, with size, qualification and GO Slim theme
gs$protein_map   # gene representative per protein, plus the plot label
gs$gene_universe # the 3,039 measured symbols
```

`msigdbr` is needed only to create a missing snapshot.
