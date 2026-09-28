# 00_build_gene_sets

This step freezes the gene sets, maps proteins to genes and applies the size filter, so every later
step reads one list. It reads `01_Preprocess/02_Quantification/c_data/proteins.rds` and the frozen
MSigDB snapshot and GO Slim OBO in `c_data/cache/`, each beside its `.md5`.

`01_run_fgsea_and_fry` and `04_run_singscore` read the tested sets, `c_data/gene_sets.rds`.
`c_data/00_build_gene_sets.xlsx` holds the set catalogue (`set_catalog`) and protein-to-gene map
(`protein_gene_map`) that `01_`, `04_`, `05_` and `06_` read, plus counts per collection and mapping
outcome. The step draws nothing, so it has no `b_reports/`.

```sh
Rscript 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.R
```

## Snapshot

The first run fetches four MSigDB collections (release 2026.1.Hs) through `msigdbr` and saves an RDS
with an md5; later runs check it and need neither network nor `msigdbr`. Restore the RDS and `.md5`
together, or move both aside to rebuild.

## Collections

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 41 | 44 |
| KEGG_Legacy | 186 | 94 | 26 |
| Reactome | 1,839 | 430 | 33 |
| GO:BP | 7,538 | 1,366 | 29 |
| GO Slim | 71 | 59 | 138 |
| Total | 9,684 | 1,990 | |

A set is tested with 15 to 500 source genes and at least 15 measured. No set is excluded by name.
Six other collections were measured and dropped, chiefly GO:CC: its sets are physical complexes,
co-regulated by stoichiometry, and three of the control's original six false positives came from it.
Dropping GO:CC and GO:MF took the control from 14 hits to 1.

GO Slim sets come from the GO Consortium's generic slim (`goslim_generic.obo`, frozen with an md5).
Each holds every measured gene annotated to the slim term or any GO:BP term below it, drawn from the
full membership rather than the tested subset, and the size rule reads measured size. Built the
other way, *protein folding* kept 19 of its 129 genes.

## One protein per gene

Proteins with no symbol or several symbols are left out of set tests rather than split. Where
proteins share a symbol, the one with the most observed precursors represents it. `protein_gene_map`
records each decision; the measured universe is 3,039 symbols.
