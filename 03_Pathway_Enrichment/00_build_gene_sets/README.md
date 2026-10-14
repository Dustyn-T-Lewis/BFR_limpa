# 00 · Build gene sets

Freezes one MSigDB release, maps each gene symbol to one protein and keeps the sets large enough
to test, so every later step reads one set list.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`, and the frozen MSigDB snapshot
(release 2026.1.Hs) and GO Slim OBO in `c_data/cache/`, each beside its `.md5`. With the cache in
place no network or `msigdbr` is needed; restore a file and its `.md5` together.

```sh
quarto render 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/00_build_gene_sets.html`: the walkthrough. No figures.
- `c_data/gene_sets.rds`: the 1,990 tested sets, read by `01_run_fgsea_and_fry`,
  `04_run_singscore`, `07_leading_edge` and `08_pathway_maps`.
- `c_data/00_build_gene_sets.xlsx`: the set catalogue (`set_catalog`) and protein-to-gene map
  (`protein_gene_map`), read by `01_`, `04_` to `07_`, plus counts per collection and mapping
  outcome.

A set is tested with 15 to 500 source genes and at least 15 measured; no set is excluded by name.
The measured universe is 3,039 symbols.

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 41 | 44 |
| Reactome | 1,839 | 430 | 33 |
| KEGG_Legacy | 186 | 94 | 26 |
| GO:BP | 7,538 | 1,366 | 29 |
| GO Slim | 71 | 59 | 138 |
