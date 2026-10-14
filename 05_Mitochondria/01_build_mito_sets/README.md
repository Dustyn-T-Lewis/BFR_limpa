# 01 · Build mito sets

Turns MitoCarta 3.0 into gene sets, the full list and the 149 MitoPathways, and flags which genes
we measure.

Inputs: `05_Mitochondria/00_Input/Human.MitoCarta3.0.xls`, the `protein_gene_map` sheet of
`03_Pathway_Enrichment/00_build_gene_sets/c_data/00_build_gene_sets.xlsx`.

```sh
quarto render 05_Mitochondria/01_build_mito_sets/a_script/01_build_mito_sets.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_build_mito_sets.html`: the walkthrough.
- `c_data/01_build_mito_sets.xlsx`: the `membership` sheet, one row per set and gene, flagged when
  a representative protein carries the symbol. MitoCarta lists 1,136 genes, 660 of them measured.

The complex readouts downstream are the pathways "CI subunits" to "CV subunits", structural
subunits only. The parent pathways "Complex I" and so on also hold assembly factors.
