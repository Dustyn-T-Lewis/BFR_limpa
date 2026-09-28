# 05_set_classification

This step measures how well each set's singscore separates the study groups. It reads
`04_run_singscore/c_data/set_scores.rds`, the catalogue in
`00_build_gene_sets/c_data/00_build_gene_sets.xlsx` and the sample targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds`. It writes `c_data/05_set_classification.xlsx`,
led by the chance comparison below, one ROC PDF per collection (`05_classification_Hallmark.pdf` and
so on) and `05_chance_by_collection.pdf`.

```sh
Rscript 03_Pathway_Enrichment/05_set_classification/a_script/05_set_classification.R
```

There are five paired tasks: pre against post in each arm, and BFR against HLRT at baseline, after
training and in change. AUC comes from `pROC` with `direction = "<"` fixed, since left free it flips
sets below chance, and p comes from the paired Wilcoxon test. Each collection is read against its own
chance count, 5% of its sets, and BH runs within collection and task.

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 1.4 | 3.8 | 3.4 | 4.0 | 7.3 |
| pre vs post, HLRT | 4.8 | 4.0 | 4.2 | 3.7 | 4.7 |
| post, BFR vs HLRT | 0.95 | 0.64 | 0.37 | 1.16 | 0.67 |
| change, BFR vs HLRT | 0.95 | 0.00 | 0.14 | 0.53 | 0.67 |
| baseline *(control)* | 0.95 | 0.43 | 0.47 | 0.48 | 0.00 |

The table gives nominal hits as a multiple of chance. After BH, 9 sets survive on BFR training, 49 on
HLRT and none elsewhere. Protein-level classification is `02_Differential_Expression/02_Differential`.

S12–S15 draw an ROC panel for every set reaching nominal p on each task except the control, twelve
to an A4 page, smallest p first, with BH q in each header. S16 plots nominal hits against chance by
collection; a ratio of 1 means the nominal hits are what chance returns.
