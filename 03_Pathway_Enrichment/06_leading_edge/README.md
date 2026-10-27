# 06 · Leading edge

Asks whether the pathways both legs move with training are driven by the same proteins, by
comparing the fgsea leading edges of the two training contrasts. Sam's analysis from PR #26 (his
overlap measures, GO hierarchy breakdown and unique-protein split) on the pipeline's own fgsea
results and gene sets; his separate GSEA run, GO:CC and GO:MF, and his choice of duplicate genes by
largest |t| are not carried over.

Inputs: the `set_tests` and `protein_results` sheets of
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`, and the `set_catalog` sheet and
`go_bp_measured.rds` from `00_build_gene_sets/c_data/`.

```sh
quarto render 03_Pathway_Enrichment/06_leading_edge/a_script/06_leading_edge.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/06_leading_edge.html`: the walkthrough.
- `b_reports/06_leading_edge.pdf`, Sam's figures, each page titled only:
  - S14: overlap heatmap, one tile per pathway coloured by Jaccard index, labelled with the shared
    count, one panel per collection;
  - S15: one Venn per pathway (ggvenn), each leg's NES and adjusted p and the shared count against
    chance above, the shared proteins below, six to a page;
  - S16: NES in BFR against HLRT for every tested pathway, colour by collection, shape by where
    fgsea calls it at FDR < 0.05, quadrant counts and Spearman rho;
  - S17: GO:BP decomposition, one parent per panel (top 25 by the larger adjusted p), child terms
    [relation] by leg; fill is the share of the leading edge, size the count, a black outline a
    hypergeometric BH-adjusted p < 0.05; the subtitle gives the share in significant descendants;
  - S18: proteins unique to each leg by child term, one parent per panel (top 25 by divergence p),
    BFR-only left and HLRT-only right, shared count and Fisher asterisks per child.
  S17 and S18 put two panels on a page.
- `c_data/06_leading_edge.xlsx`, read by no later step:
  - `leading_edge`: one row per pathway: shared and one-leg counts, NES and q per leg, Jaccard
    index, overlap coefficient and the shared count expected by chance;
  - `proteins`: one row per pathway and leading-edge protein, with the legs it appears in;
  - `go_children`: one row per GO:BP parent, child term and leg: proteins in the child, share of
    the leading edge, fold, hypergeometric p and q;
  - `go_coverage`: one row per parent and leg: share of the leading edge in significant
    descendant terms, and the proteins outside them;
  - `unique_children`: one row per parent and child: BFR-only and HLRT-only proteins in the child,
    Fisher p and q, the parent's divergence and its shuffle p and q, and `edge_pct` per leg.

The overlap measures, Venn statistics and child tests describe: leading edges move with small
ranking changes, children overlap, and the Fisher and shuffle tests treat proteins as
exchangeable.
