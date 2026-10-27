# 07 · Leading edge

Asks whether the pathways both legs move with training are driven by the same proteins, by
comparing the fgsea leading edges of the two training contrasts. Sam's analysis from PR #26 (his
overlap measures, GO hierarchy breakdown and unique-protein split) on the pipeline's own fgsea
results and gene sets; his separate GSEA run, GO:CC and GO:MF, and his choice of duplicate genes by
largest |t| are not carried over.

Inputs: the `set_tests` and `protein_results` sheets of
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`, and the `set_catalog` sheet and
`go_bp_measured.rds` from `00_build_gene_sets/c_data/`.

```sh
quarto render 03_Pathway_Enrichment/07_leading_edge/a_script/07_leading_edge.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/07_leading_edge.html`: the walkthrough.
- `b_reports/07_leading_edge.pdf`: S20, shared and one-leg leading-edge counts for every pathway;
  S21, one Venn per pathway with NES, q, shared against chance and the shared proteins; S22, GO:BP
  parents broken into child terms (share of the leading edge, count, hypergeometric q); S23, the
  proteins unique to each leg by child term.
- `c_data/07_leading_edge.xlsx`, read by no later step:
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
