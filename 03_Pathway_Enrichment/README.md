# 03 · Pathway Enrichment

Tests whether proteins that work together moved together, then scores every sample on every set.

| Step | Runs | Writes |
|---|---|---|
| [`00_build_gene_sets`](00_build_gene_sets/README.md) | frozen MSigDB snapshot, protein-to-gene map, size filter, GO Slim sets | `gene_sets.rds` |
| [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) | fgsea and fry per contrast, `collapsePathways` | `set_tests.rds`, 23 dot plots |
| [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) | protein volcanoes with pathway rings | 6 volcanoes |
| [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md) | BFR NES against HLRT NES | 2 composites |
| [`04_run_singscore`](04_run_singscore/README.md) | per-sample set scores | `singscore.rds`, 2 figures |
| [`05_classify_and_associate_sets`](05_classify_and_associate_sets/README.md) | set classification and phenotype association | `set_results.rds`, 7 figures, 96 pages |

```sh
for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea \
         03_enrich_scatter_fgsea 04_run_singscore 05_classify_and_associate_sets; do
  Rscript 03_Pathway_Enrichment/$s/a_script/$s.R
done
```

About three minutes in total, most of it drawing the `05_` pages. Each substage writes every
figure as PDF, and a PNG when the figure is one page, and bundles the PDFs into one
`<substage>_figures.pdf`.
Figure titles name the figure, subtitles give method and counts, captions state the encodings and
the source table. Findings live here, not on the figures.

## Methods

**fgsea** is competitive: does a set sit at one end of the protein ranking? **fry** is
self-contained: did the set move at all under the fitted design? Both run on the same 1,990 sets
and both ship. fgsea assumes proteins are exchangeable, and co-regulated sets break that; fry does
not assume it.

**singscore** gives one rank-based score per set per sample. It carries no p-value and never sees
a contrast, and a sample's score does not change with the cohort, which the paired design needs.

**Order: pool, test, then collapse.** All sets are tested and BH corrects within each contrast.
`collapsePathways` then marks non-redundant hits in the `main` column and deletes nothing.
Deduplicating before testing was measured and lost two thirds of the discoveries (363 to 133 on
BFR training): a redundancy filter removes significant sets, not hopeless ones (Bourgon et al.,
PNAS 2010).

No set is excluded by name.

## Results

| Contrast | fgsea | after collapse | fry |
|---|---:|---:|---:|
| BFR_Post-Pre | 376 | 88 | 656 |
| HLRT_Post-Pre | 281 | 68 | 589 |
| BFR_Post-HLRT_Post | 159 | 52 | 0 |
| Modality_x_Time_Interaction | 26 | 10 | 0 |
| BFR_Pre-HLRT_Pre *(control)* | 1 | 1 | 0 |

fry finds nothing on any between-leg contrast, so the interaction's 10 fgsea pathways are leads,
not findings. The control's one fgsea hit, `REACTOME_STRIATED_MUSCLE_CONTRACTION`, is the floor
every count should be read against.

**Classification by set score**, nominal hits over chance, per collection:

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 1.4 | 3.8 | 3.4 | 4.0 | 7.3 |
| pre vs post, HLRT | 4.8 | 4.0 | 4.2 | 3.7 | 4.7 |
| post, BFR vs HLRT | 0.95 | 0.64 | 0.37 | 1.16 | 0.67 |
| change, BFR vs HLRT | 0.95 | 0.00 | 0.14 | 0.53 | 0.67 |
| baseline *(control)* | 0.95 | 0.43 | 0.47 | 0.48 | 0.00 |

Training clears chance in every collection; no between-leg task does. Under BH, 9 sets survive on
BFR training and 49 on HLRT, and none elsewhere. No set, and no protein
(`02_Differential_Expression/03_Phenotype`), tracks the phenotype after correction. With 32 pairs
the study resolves a paired correlation near 0.5, so this is a null at that resolution.

**Concordance.** BFR and HLRT NES correlate at rho 0.83 over all sets, 0.86 over collapse
survivors and 0.82 over Hallmark and GO Slim. No set significant in both arms differs in sign.
Both modalities moved the same pathways by similar amounts, which is why the interaction is empty.
