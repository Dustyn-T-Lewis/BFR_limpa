# 02_enrich_volcano_fgsea

This step draws protein volcanoes with the collapse-surviving fgsea sets ringed around them. It
computes nothing. It reads `protein_results` and `set_tests` from
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx` and writes figures S8–S9 to
`b_reports/02_enrich_volcano_fgsea_figures.pdf`. Having no table to write, it has no `c_data/`.

```sh
Rscript 03_Pathway_Enrichment/02_enrich_volcano_fgsea/a_script/02_enrich_volcano_fgsea.R
```

Point colour shows each protein's BH FDR, and the ring shows each set's fgsea FDR, with red for up
and blue for down. The ring holds the twelve collapse survivors with the lowest adjusted p, in
either direction. S8 has four panels: the interaction, both training responses, and post-training
BFR against HLRT. S9 repeats the two training responses with the labels ranked by Π. The control
is not drawn.

Leading-edge genes are matched to the label of the protein fgsea ranked for that symbol, because
the label carries an accession when two proteins share a symbol. When two ringed sets shorten to
the same display name, the collection is appended.
