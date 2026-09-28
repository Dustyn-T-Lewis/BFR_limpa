# 00_build_gene_sets

This step freezes the gene sets, maps proteins to genes and applies the size filter, so every later
step reads the same list. It reads the protein table from
`01_Preprocess/02_Quantification/c_data/proteins.rds`, and the frozen MSigDB snapshot and GO Slim
OBO from `c_data/cache/`, each stored beside its `.md5`.

The tested sets are saved as `c_data/gene_sets.rds`, which `01_run_fgsea_and_fry` and
`04_run_singscore` read. `c_data/00_build_gene_sets.xlsx` holds the set catalogue (`set_catalog`)
and the protein-to-gene map (`protein_gene_map`) that `01_`, `04_` and `05_` read, with counts per
collection and per mapping outcome. The step draws nothing, so it has no `b_reports/`.

```sh
Rscript 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.R
```

## Snapshot

The first run fetches four MSigDB collections (release 2026.1.Hs) through `msigdbr` and saves them
as an RDS with an md5 checksum. Later runs check the checksum and need neither a network
connection nor `msigdbr`. Restore the RDS and its `.md5` together, or move both aside to rebuild.

## Collections

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 41 | 44 |
| KEGG_Legacy | 186 | 94 | 26 |
| Reactome | 1,839 | 430 | 33 |
| GO:BP | 7,538 | 1,366 | 29 |
| GO Slim | 71 | 59 | 138 |
| Total | 9,684 | 1,990 | |

A set is tested when it has 15 to 500 source genes and at least 15 of them measured. No set is
excluded by name. Six other collections were measured and dropped, GO:CC chief among them. Its
sets are physical complexes, co-regulated by stoichiometry, and three of the control's original six
false positives came from it. Dropping GO:CC and GO:MF took the control from 14 hits to 1.

The GO Slim sets are built here from the GO Consortium's generic slim (`goslim_generic.obo`,
frozen with an md5). Each holds every measured gene annotated to the slim term or to any term below
it in GO:BP. Genes come from the full membership rather than the tested subset, and the size rule
reads measured size. Built the other way, *protein folding* kept 19 of its 129 genes.

## One protein per gene

Proteins with no symbol or with several symbols are left out of the set tests rather than split.
Where several proteins share a symbol, the one with the most observed precursors represents it.
`protein_gene_map` records each decision, and the measured universe is 3,039 symbols.
