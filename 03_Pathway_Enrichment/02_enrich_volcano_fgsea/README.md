# 02_enrich_volcano_fgsea

This step draws protein volcanoes with the collapse-surviving fgsea sets ringed, computing nothing.
It reads `protein_results` and `set_tests` from
`01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx` and writes S8–S9 to
`b_reports/02_enrich_volcano_fgsea_figures.pdf`; with no table to write, it has no `c_data/`.

```sh
Rscript 03_Pathway_Enrichment/02_enrich_volcano_fgsea/a_script/02_enrich_volcano_fgsea.R
```

Point colour shows protein BH FDR and the ring set fgsea FDR, red for up and blue for down. The ring
holds the twelve collapse survivors with the lowest adjusted p in either direction. S8 has four
panels (the interaction, both training responses, post-training BFR against HLRT); S9 repeats the
training responses labelled by Π. The control is not drawn.

Leading-edge genes match the label of the protein fgsea ranked for that symbol, since a label
carries an accession when two proteins share a symbol. When two ringed sets shorten to the same
name, the collection is appended.