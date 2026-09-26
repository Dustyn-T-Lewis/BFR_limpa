# 02_enrich_volcano_fgsea

Draws protein volcanoes with the collapse-surviving fgsea sets ringed. Computes nothing.

- Reads: `01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx` (`protein_results`, `set_tests`)
- Writes: `b_reports/02_enrich_volcano_fgsea_figures.pdf` (S8–S9). It writes no table, so it has no
  `c_data/`.
- Run: `Rscript 03_Pathway_Enrichment/02_enrich_volcano_fgsea/a_script/02_enrich_volcano_fgsea.R`,
  3 seconds.

Point colour reads protein BH FDR; the ring reads set fgsea FDR. Red is up, blue down. The ring
holds the twelve collapse survivors with the lowest adjusted p, in either direction. S8 has four
panels (interaction, both training responses, post-training BFR against HLRT); S9 repeats the two
training responses with labels ranked by Π. The control is not drawn.

Leading-edge genes are translated to the label of the protein fgsea ranked for that symbol, since
a label carries its accession when two proteins share a symbol. When two ringed sets reduce to one
display name, the collection is appended.
