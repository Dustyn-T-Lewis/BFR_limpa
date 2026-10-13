# 07_leading_edge

This step asks whether the pathways both legs move with training are driven by the same proteins.
It reads the fgsea results in `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx` and the set
catalogue and `gene_sets.rds` from `00_build_gene_sets`, and writes `c_data/07_leading_edge.xlsx`
and `b_reports/07_leading_edge.pdf` (S20 to S22).

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
ubiquitination and tRNA aminoacylation.

S20 draws the shared and one-leg counts for all 71 pathways on one page. S21 gives each pathway
its own Venn, areas to scale, with the shared proteins listed under the title, six to a page. S22
breaks the 27 GO:BP pathways that have tested child terms into those children (61 in all) and
counts the parent's leading-edge proteins of each leg inside each child, BFR to the left and HLRT
to the right. In 2 of the 66 parent and child pairs one leg holds at least five proteins and more
than twice the other's count, both times BFR.

Leading edges move with small changes in ranking, and child terms overlap, so none of these counts
is tested; they describe.

The questions and the Venn and GO decomposition views come from Sam's PR #26. His version ran a
second GSEA through clusterProfiler, with GO cellular component and molecular function, and tested
the child terms with Fisher's test; this step reads the fgsea run the stage already has, keeps the
stage's collections and leaves the overlapping child terms untested.
