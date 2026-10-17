# 07 · Leading edge

Asks whether the pathways both legs move with training are driven by the same proteins, by
comparing the fgsea leading edges of the two training contrasts. Adapted from Sam's PR #26.

Inputs: the `set_tests` sheet of `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`, and the
`set_catalog` sheet and `gene_sets.rds` from `00_build_gene_sets/c_data/`.

```sh
quarto render 03_Pathway_Enrichment/07_leading_edge/a_script/07_leading_edge.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/07_leading_edge.html`: the walkthrough.
- `b_reports/07_leading_edge.pdf`: S20, shared and one-leg leading-edge counts for every pathway;
  S21, one Venn per pathway with the shared proteins listed; S22, GO:BP parents broken into their
  tested child terms.
- `c_data/07_leading_edge.xlsx`, read by no later step:
  - `leading_edge`: one row per pathway: shared and one-leg counts, Jaccard index, NES per leg.
  - `proteins`: one row per pathway and leading-edge protein, with the legs it appears in.
  - `go_children`: one row per GO:BP parent, tested child term and leg, with that leg's count.

71 pathways qualify (41 GO:BP, 20 Reactome, 7 KEGG, 2 GO Slim, 1 Hallmark). The two legs' leading
edges overlap with a median Jaccard index of 0.58, the lowest 0.41 (protein polyubiquitination).
27 GO:BP parents have 61 tested children in 66 pairs; in 2 pairs one leg holds at least five
proteins and more than twice the other's count, both times BFR. None of these counts is tested,
and the Jaccard index has no chance reference here.
