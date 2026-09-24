# 03_Pathway_Enrichment / 00_build_gene_sets

Freezes the gene sets, maps proteins to genes and applies the size filter. Every later step reads
this one list.

- Reads: `01_Preprocess/02_Quantification/c_data/proteins.rds`, `c_data/cache/` (the frozen
  MSigDB snapshot and the GO Slim OBO, each beside its `.md5`)
- Writes: `c_data/gene_sets.rds` (the tested sets, read by `01_run_fgsea_and_fry` and
  `04_run_singscore`), `c_data/00_build_gene_sets.xlsx` (`overview`, `collection_summary`,
  `set_catalog`, `protein_gene_map`, `mapping_summary`; `01_`, `04_` and `05_` read the catalogue
  and the map). It draws nothing, so it has no `b_reports/`.
- Run: `Rscript 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.R`, about
  five seconds.

## MSigDB 2026.1.Hs is frozen with a checksum

The first run fetches four MSigDB collections (release 2026.1.Hs) through `msigdbr` and writes an
RDS with an md5. Later runs verify the checksum and need neither network nor `msigdbr`. Restore
the RDS and its `.md5` together; move both aside to rebuild.

## 1,990 sets across five collections are tested

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 41 | 44 |
| KEGG_Legacy | 186 | 94 | 26 |
| Reactome | 1,839 | 430 | 33 |
| GO:BP | 7,538 | 1,366 | 29 |
| GO Slim | 71 | 59 | 138 |
| Total | 9,684 | 1,990 | |

A set is tested with 15 to 500 source genes and at least 15 measured. No set is excluded by name.
Six other collections were measured and dropped. GO:CC is the main one: its sets are physical
complexes, co-regulated by stoichiometry, and three of the control's original six false positives
came from it. Dropping GO:CC and GO:MF took the control from 14 hits to 1.

GO Slim sets are built here from the GO Consortium's generic slim (`goslim_generic.obo`, frozen
with an md5). Each holds every measured gene annotated to the slim term or any term below it in
GO:BP. Genes come from the full membership, not the tested subset, and the size rule reads
measured size. Built otherwise, *protein folding* kept 19 of its 129 genes.

## One protein represents each gene

Proteins with no symbol or several symbols are left out of set tests rather than split. Where
several proteins share a symbol, the one with the most observed precursors represents it.
`protein_gene_map` records each decision; the measured universe is 3,039 symbols.
