# 08_pathway_maps

This step draws the training response of both legs as pathway maps. It reads `protein_results` and
`set_tests` from `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx` and `gene_sets.rds` from
`00_build_gene_sets`, and writes one PDF per collection to `b_reports/` (`08_maps_Hallmark.pdf` and
so on). It computes no table, so it has no `c_data/`.

```sh
Rscript 03_Pathway_Enrichment/08_pathway_maps/a_script/08_pathway_maps.R
```

enrichplot's maps need a clusterProfiler result, so `compareCluster` runs GSEA on both training
rankings at once. It uses fgsea underneath, with the same sets and the same moderated t ranking as
`01_run_fgsea_and_fry`. Across all 3,980 set and contrast pairs its NES agree with 01's at
r ≥ 0.9998 in every collection, largest difference 0.094. The figures draw from it; every test and
table stays in 01.

Each map shows the main sets that 01 calls at FDR < 0.05 in either leg, one per cluster of
overlapping sets: 6 Hallmark, 5 GO Slim, 11 KEGG, 31 Reactome and the 30 strongest of GO:BP's 73.
Each PDF has four pages:

- S23, dotplot: BFR and HLRT side by side, colour NES, dot size the leading-edge share.
- S24, enrichment map: pathways as nodes, split into BFR and HLRT pies, joined when their proteins
  overlap.
- S25, tree: pathways clustered by shared proteins, each cluster named by its most frequent words.
- S26, network: the five strongest pathways and their leading-edge proteins. The proteins shared
  between pathways are named in S21 of `07_leading_edge`.

The step runs in about 30 seconds. `compareCluster` does not pass fgsea's `eps` through, so p
below 1e-10 is floored; that only shades the dots. The idea of richer pathway figures comes from
Sam's PR #26.
