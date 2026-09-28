# 05_classify_and_associate_sets

This step measures how well each set's score separates study groups and whether it tracks the
phenotype. It reads `04_run_singscore/c_data/set_scores.rds`, the catalogue in
`00_build_gene_sets/c_data/00_build_gene_sets.xlsx`,
`01_Preprocess/02_Quantification/c_data/proteins.rds` and `00_Input/phenotype.csv`. It writes
`c_data/05_classify_and_associate_sets.xlsx`, led by the chance comparison below, and S12–S18, 128
pages, to `b_reports/05_classify_and_associate_sets_figures.pdf`.

```sh
Rscript 03_Pathway_Enrichment/05_classify_and_associate_sets/a_script/05_classify_and_associate_sets.R
```

Five paired tasks: pre against post in each arm, and BFR against HLRT at baseline, after training
and in change. AUC comes from `pROC` with `direction = "<"` fixed, since left free it flips sets
below chance, and p from the paired Wilcoxon test. Phenotype associations use Spearman `cor.test`,
pooled over 65 legs and differential over 32 participants.

Nominal hits as a multiple of chance, each collection read against its own chance count (5% of its
sets), with BH within collection and task:

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 1.4 | 3.8 | 3.4 | 4.0 | 7.3 |
| pre vs post, HLRT | 4.8 | 4.0 | 4.2 | 3.7 | 4.7 |
| post, BFR vs HLRT | 0.95 | 0.64 | 0.37 | 1.16 | 0.67 |
| change, BFR vs HLRT | 0.95 | 0.00 | 0.14 | 0.53 | 0.67 |
| baseline *(control)* | 0.95 | 0.43 | 0.47 | 0.48 | 0.00 |

After BH, 9 sets survive on BFR training, 49 on HLRT and none elsewhere. No set or protein tracks
the phenotype after correction (`02_Differential_Expression/03_Phenotype`).

S12–S15 draw an ROC panel for every set reaching nominal p on each task except the control; S16–S17
a scatter for every nominal set-outcome pair, pooled then differential; S18 nominal hits against
chance by collection. Panels run twelve to an A4 page, smallest p first, with BH q in each header.
Read every panel against S18: at a ratio of 1, nominal hits are what chance returns.