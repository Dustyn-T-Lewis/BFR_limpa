# 00 · Build gene sets

Freezes one MSigDB release, maps each gene symbol to one protein and keeps the sets large enough
to test, so every later step reads one set list.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`; the frozen MSigDB snapshot
(2026.1.Hs) and GO Slim OBO in `c_data/cache/`, each beside its `.md5` (restore a file and its
`.md5` together). With the cache in place no network or `msigdbr` is needed. GO Slim also reads
`GO.db`, which is not frozen: these sets were built with GO.db 3.23.1 (GO release 2026-01-23).

```sh
quarto render 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/00_build_gene_sets.html`: the walkthrough.
- `c_data/gene_sets.rds`: the 1,990 tested sets, read by `01_`, `03_`, `05_`, `07_` and
  `04_Network/01_modules`.
- `c_data/go_bp_measured.rds`: measured genes of every GO:BP set, tested or not, read by `06_`.
- `c_data/00_build_gene_sets.xlsx`, read by `01_`, `03_` to `06_`, `04_Network` and
  `05_Mitochondria`:
  - `set_catalog`: one row per source set; `qualifies` marks the tested ones.
  - `protein_gene_map`: one row per protein; `selected` marks each symbol's representative.

A set is tested with 15 to 500 source genes and at least 15 measured; no set is excluded by name.
The measured universe is 3,040 symbols. GO Slim keeps the BP slim terms with GO:BP offspring.

| Collection | In source | Tested | Median measured |
|---|---:|---:|---:|
| Hallmark | 50 | 41 | 44 |
| Reactome | 1,839 | 430 | 33 |
| KEGG_Legacy | 186 | 94 | 26 |
| GO:BP | 7,538 | 1,366 | 29 |
| GO Slim | 71 | 59 | 138 |

Protein groups with no symbol or several symbols are dropped, not split.
