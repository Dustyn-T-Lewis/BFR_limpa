# 02 · Enrichment volcanoes, fgsea

Draws protein volcanoes with the collapse-surviving fgsea sets ringed; computes nothing. The
baseline contrast is not drawn.

Inputs: `protein_results` and `set_tests` from
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`.

```sh
quarto render 03_Pathway_Enrichment/02_enrich_volcano_fgsea/a_script/02_enrich_volcano_fgsea.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_enrich_volcano_fgsea.html`: the walkthrough.
- `b_reports/02_enrich_volcano_fgsea_figures.pdf`: S3 (the interaction, both training responses
  and post-training BFR against HLRT, five lowest-FDR proteins named) and S4 (the training
  responses with labels ranked by pi-score, which selects nothing, each followed by the protein's
  mean observed precursors). Points are coloured when protein BH FDR < 0.05; the ring shows up to
  twelve collapse survivors, lowest adjusted p first, coloured by NES, with ticks on their
  leading-edge proteins.

No table is written, so there is no `c_data/`.

The training volcanoes are ringed by rising mitochondrial sets in both legs, led by aerobic
respiration and oxidative phosphorylation. NAMPT tops both legs (log2 +0.69 with BFR, +1.01 with
HLRT) and MYLK2 falls in both; the BFR panel also names SUCLA2, HK1 and IDH3B, the HLRT panel
ALDH1B1, ART3 and GSTK1. The interaction and post-training panels name no protein, since none
reaches BH 0.05, and their rings hold the between-leg fgsea leads. The interaction ring holds cell
motility, adhesion and cytoskeleton sets, lower in BFR. The post-training ring holds muscle
contraction and cytoplasmic translation sets lower in BFR and amino acid metabolism sets higher.
fry confirms none of the sets on those two panels (see `01_run_fgsea_and_fry`).
