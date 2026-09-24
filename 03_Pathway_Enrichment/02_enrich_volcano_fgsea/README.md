# 02_enrich_volcano_fgsea

Protein volcanoes with collapse-surviving pathways ringed. Computes nothing.

| | |
|---|---|
| Reads | `set_tests.rds` |
| Writes | 6 volcanoes |

Point colour reads protein BH FDR; the ring reads set fgsea FDR. Red is up, blue down. Four FDR
panels (interaction, both training responses, post-training BFR against HLRT) and two Π-ranked
repeats. The control is not drawn.

Leading-edge genes are translated to point labels before drawing, since a label carries its
accession when two proteins share a symbol. When two ringed sets reduce to one display name, the
collection is appended.
