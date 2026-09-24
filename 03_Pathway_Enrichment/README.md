# 03 · Pathway Enrichment

Tests whether proteins that work together moved together, then scores every sample on every set.

| Step | Runs | Writes |
|---|---|---|
| [`00_build_gene_sets`](00_build_gene_sets/README.md) | frozen MSigDB snapshot, protein-to-gene map, size filter, GO Slim sets | `gene_sets.rds`, workbook |
| [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) | fgsea and fry per contrast, `collapsePathways` | workbook, S1–S7 |
| [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) | protein volcanoes with pathway rings | S8–S9 |
| [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md) | BFR NES against HLRT NES | workbook, S10–S11 |
| [`04_run_singscore`](04_run_singscore/README.md) | per-sample set scores | `set_scores.rds`, workbook, S12 |
| [`05_classify_and_associate_sets`](05_classify_and_associate_sets/README.md) | set classification and phenotype association | workbook, S13–S19 over 128 pages |

```sh
for s in 00_build_gene_sets 01_run_fgsea_and_fry 02_enrich_volcano_fgsea \
         03_enrich_scatter_fgsea 04_run_singscore 05_classify_and_associate_sets; do
  Rscript 03_Pathway_Enrichment/$s/a_script/$s.R
done
```

About two minutes in total, half of it in `05_`. Each substage writes at most two things:

- `c_data/<substage>.xlsx`, an `overview` sheet listing the others, then one sheet per table. Later
  substages read their tables from these workbooks.
- `b_reports/<substage>_figures.pdf`, A4 pages numbered S1 to S19 through the stage, each with a
  supplement-style caption beneath: title, one entry per lettered panel, encodings, and the sheet
  the data come from. Findings live here, not on the figures.

An `.rds` is written only where the next substage needs an R object: the set list, and the
singscore matrix, whose exact rank ties a trip through Excel would break.

## Methods

fgsea is competitive: does a set sit at one end of the protein ranking? fry is
self-contained: did the set move at all under the fitted design? Both run on the same 1,990 sets
and both ship. fgsea assumes proteins are exchangeable, and co-regulated sets break that; fry does
not assume it.

singscore gives one rank-based score per set per sample. It carries no p-value and never sees
a contrast, and a sample's score does not change with the cohort, which the paired design needs.

The order is pool, test, then collapse. All sets are tested and BH corrects within each contrast.
`collapsePathways` then marks non-redundant hits in the `main` column and deletes nothing.
Deduplicating before testing lost two thirds of the discoveries (363 to 133 on BFR training): a
redundancy filter removes significant sets, not hopeless ones (Bourgon et al., PNAS 2010).

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

Classification by set score, nominal hits over chance, per collection:

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

BFR and HLRT NES correlate at rho 0.83 over all sets, 0.86 over collapse survivors and 0.82 over
Hallmark and GO Slim. No set significant in both arms differs in sign.
Both modalities moved the same pathways by similar amounts, which is why the interaction is empty.
