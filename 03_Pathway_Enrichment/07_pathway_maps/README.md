# 07 · Pathway maps

Draws the training response of both legs as pathway maps, from a `compareCluster` GSEA run on the
same sets and rankings as `01_run_fgsea_and_fry`. clusterProfiler 4.20 runs that GSEA through
enrichit, whose multilevel code is derived from fgsea's but separate, so the pathways drawn are
chosen from 01's `set_tests` while their NES and leading edges come from this run. Figures only;
every test and table stays in 01. Adapted from Sam's PR #26.

Inputs: the `protein_results` and `set_tests` sheets of
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`, and
`00_build_gene_sets/c_data/gene_sets.rds`.

```sh
quarto render 03_Pathway_Enrichment/07_pathway_maps/a_script/07_pathway_maps.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/07_pathway_maps.html`: the walkthrough, with the NES check against 01.
- `b_reports/07_maps_<collection>.pdf`, one per collection, four pages: S19 dotplot (colour NES,
  dot size the leading-edge share), S20 enrichment map (pathways joined by shared proteins, pies
  split by leg), S21 pathway tree (clusters named by their most frequent words), S22 the five
  strongest pathways and their leading-edge proteins, pies showing which leg each protein leads
  in. These leading edges can differ slightly from the shared proteins in S15 of 06_leading_edge.

Each map shows the main sets 01 calls at FDR < 0.05 in either leg, at most 30 per collection: 6
Hallmark, 5 GO Slim, 11 KEGG, 30 of Reactome's 31 and 30 of GO:BP's 73. Across all 3,980 set
and contrast pairs, this run's NES agree with 01's at r ≥ 0.9998 in every collection; the largest
difference is 0.094, in Reactome. No `c_data/`.
