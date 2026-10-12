# 07_leading_edge

This step asks whether the pathways both legs move with training are driven by the same proteins.
It reads the fgsea results in `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx` and writes
`c_data/07_leading_edge.xlsx` and `b_reports/07_leading_edge.pdf` (S20).

```sh
Rscript 03_Pathway_Enrichment/07_leading_edge/a_script/07_leading_edge.R
```

A pathway's leading edge is the set of its proteins that fgsea ranks before the running-sum peak,
the ones carrying the enrichment. The step keeps every pathway called at FDR < 0.05 in both
training contrasts, in the same direction, and kept as a main set by `collapsePathways` in at
least one leg: 71 pathways (41 GO:BP, 20 Reactome, 7 KEGG, 2 GO Slim, 1 Hallmark). For each it
counts the leading-edge proteins shared by both legs and those in one leg only, and gives the
Jaccard index, shared over all.

The legs mostly share their drivers. The median Jaccard index is 0.58 and the lowest is 0.41, so
no pathway is driven by two separate groups of proteins. The lowest overlaps are in protein
ubiquitination and tRNA aminoacylation. Leading edges move with small changes in ranking, so these
counts describe and do not test.

The idea comes from Sam's PR #26, which ran a second GSEA through clusterProfiler, with GO cellular
component and molecular function, to get its leading edges. This step reads the fgsea run the
stage already has and keeps the stage's collections.
