# 02 · Enrichment volcanoes, fgsea

Draws protein volcanoes with the collapse-surviving fgsea sets ringed; computes nothing. The
baseline contrast is not drawn.

Inputs: `protein_results` and `set_tests` from `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`.

```sh
quarto render 03_Pathway_Enrichment/02_enrich_volcano_fgsea/a_script/02_enrich_volcano_fgsea.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_enrich_volcano_fgsea.html`: the walkthrough.
- `b_reports/02_enrich_volcano_fgsea_figures.pdf`: S8 (the interaction, both training responses
  and post-training BFR against HLRT, five lowest-FDR proteins named) and S9 (the training
  responses with labels ranked by Π, which selects nothing, each followed by the protein's mean
  observed precursors).

No table is written, so there is no `c_data/`. fry finds no set on S8's panels A and D, the
interaction and BFR_Post-HLRT_Post, so their rings are fgsea leads only.
