# 00 · Build gene sets

Freezes one MSigDB release, maps each gene symbol to one protein and keeps the sets large enough
to test, so every later step reads one set list.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`, and the frozen MSigDB snapshot
(release 2026.1.Hs) and GO Slim OBO in `c_data/cache/`, each beside its `.md5`. With the cache in
place no network or `msigdbr` is needed; restore a file and its `.md5` together. The GO Slim sets
also read the GO hierarchy from `GO.db`, which is not frozen: these sets were built with GO.db
3.23.1 (GO release 2026-01-23), and `06_leading_edge` reads the same package for child terms.

```sh
quarto render 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/00_build_gene_sets.html`: the walkthrough. No figures.
- `c_data/gene_sets.rds`: the 1,990 tested sets, read by `01_run_fgsea_and_fry`,
  `03_run_singscore`, `06_leading_edge` and `07_pathway_maps`.
- `c_data/go_bp_measured.rds`: the measured genes of every GO:BP set, tested or not, read by
  `06_leading_edge` for child terms.
- `c_data/00_build_gene_sets.xlsx`, read by `01_`, `03_` to `06_` and `05_Mitochondria`:
  - `collection_summary`: per collection, sets in the source, sets tested, median measured size.
  - `set_catalog`: one row per source set; `qualifies` marks the tested ones.
  - `protein_gene_map`: one row per protein; `selected` marks each symbol's representative.
  - `mapping_summary`: proteins per mapping outcome.

A set is tested with 15 to 500 source genes and at least 15 measured; no set is excluded by name.
The measured universe is 3,039 symbols. GO Slim keeps the BP slim terms that have GO:BP offspring:
the CC and MF terms drop out, and so do BP leaf terms.

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 41 | 44 |
| Reactome | 1,839 | 430 | 33 |
| KEGG_Legacy | 186 | 94 | 26 |
| GO:BP | 7,538 | 1,366 | 29 |
| GO Slim | 71 | 59 | 138 |
