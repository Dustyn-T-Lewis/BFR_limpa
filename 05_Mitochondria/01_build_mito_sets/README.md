# 01_build_mito_sets

This step turns MitoCarta 3.0 into gene lists. It reads `05_Mitochondria/00_Input/Human.MitoCarta3.0.xls`
and the protein-to-gene map in
`03_Pathway_Enrichment/00_build_gene_sets/c_data/00_build_gene_sets.xlsx`, and writes
`c_data/01_build_mito_sets.xlsx` with one sheet, `membership`: one row per set and gene, flagged
when a representative protein carries that symbol.

```sh
Rscript 05_Mitochondria/01_build_mito_sets/a_script/01_build_mito_sets.R
```

The sets are the full MitoCarta list (1,136 genes, 660 measured) and the 149 MitoPathways. The
complex readouts in `02_test_mito` are the pathways "CI subunits" to "CV subunits": structural
subunits, mtDNA-encoded ones included. The parent pathways "Complex I" and so on also hold assembly
factors, which build a complex without tracking how much of it exists.
