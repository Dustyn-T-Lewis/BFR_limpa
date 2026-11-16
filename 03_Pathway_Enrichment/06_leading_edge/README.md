# 06 · Leading edge

Asks whether the pathways both legs move with training are driven by the same proteins, by
comparing the fgsea leading edges of the two training contrasts. It runs Sam's PR #26 analysis
(overlap measures, GO hierarchy breakdown, unique-protein split) on the pipeline's own fgsea
results and sets; his separate GSEA run, GO:CC and GO:MF, and duplicate choice by |t| are not used.

Inputs: the `set_tests` and `protein_results` sheets of
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`; the `set_catalog` sheet and
`go_bp_measured.rds` from `00_build_gene_sets/c_data/`.

```sh
quarto render 03_Pathway_Enrichment/06_leading_edge/a_script/06_leading_edge.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/06_leading_edge.html`: the walkthrough.
- `b_reports/06_leading_edge.pdf`:
  - S14: overlap heatmap, one tile per pathway coloured by Jaccard index, labelled with the
    shared count, one panel per collection;
  - S15: one Venn per pathway, each leg's NES and adjusted p and the shared count against chance
    above, the shared proteins below, six to a page;
  - S16: NES in BFR against HLRT for every tested pathway, colour by collection, shape by where
    fgsea calls it at FDR < 0.05, with quadrant counts and Spearman rho;
  - S17: GO:BP decomposition for the 25 parents with the smallest adjusted p in the weaker leg,
    up to 15 child terms [relation] by leg; fill is the share of the leading edge, size the
    count, a black outline hypergeometric BH p < 0.05; the subtitle gives the share in
    significant descendants;
  - S18: proteins unique to each leg by child term for the top 25 parents by divergence p,
    BFR-only left and HLRT-only right, an asterisk at Fisher BH p < 0.05.
  S17 and S18 put two panels on a page.
- `c_data/06_leading_edge.xlsx`, read by no later step:
  - `leading_edge`: one row per pathway: shared and one-leg counts, NES and q per leg, Jaccard
    index, overlap coefficient, shared count expected by chance.
  - `proteins`: one row per pathway and leading-edge protein, with the legs it appears in.
  - `go_children`: one row per GO:BP parent, child and leg: proteins in the child, share of the
    leading edge, fold, hypergeometric p and q.
  - `go_coverage`: one row per parent and leg: share of the leading edge in significant
    descendants, and the proteins outside them.
  - `unique_children`: one row per parent and child: BFR-only and HLRT-only proteins in the
    child, Fisher p and q, the parent's divergence with its shuffle p and q, `edge_pct` per leg.

Both legs move the same pathways the same way, mostly through the same proteins. Across all 1,990
sets the two legs' NES correlate (Spearman rho 0.76), and all 71 pathways both legs move go the
same direction in each: 30 rise (mitochondrial respiration, translation, protein import, TCA
cycle) and 41 fall (muscle contraction, myogenesis, proteasome). Their leading edges overlap with
median Jaccard index 0.55 (0.32 to 0.94). The closest overlap is in falling muscle sets
(dystrophin glycoprotein complex 0.94, striated muscle contraction 0.86); in aerobic respiration
and electron transport the BFR leading edge holds more proteins of its own (79 shared, 47 BFR
only, 26 HLRT only).

The overlaps and tests are descriptive: leading edges move with small ranking changes, children
overlap, and the Fisher and shuffle tests treat proteins as exchangeable.
