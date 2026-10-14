# 05 · Set classification

Measures how well each set's singscore separates the study groups in five paired tasks: pre against
post in each arm, and BFR against HLRT at baseline, after training and in change.

Inputs: `04_run_singscore/c_data/set_scores.rds`, the `set_catalog` sheet of
`00_build_gene_sets/c_data/00_build_gene_sets.xlsx`, and the targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 03_Pathway_Enrichment/05_set_classification/a_script/05_set_classification.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/05_set_classification.html`: the walkthrough.
- `b_reports/05_classification_<collection>.pdf`, one per collection: S12 to S15, an ROC panel for
  every set reaching nominal p on each task except the baseline, smallest p first, BH q in each
  header.
- `b_reports/05_chance_by_collection.pdf`: S16, nominal hits against chance by collection and task.
- `c_data/05_set_classification.xlsx`: `chance_expectation`, `set_auc` (AUC, paired Wilcoxon p and
  BH q per set and task) and `set_catalog`. Read by no later step.

Nominal hits as a multiple of chance (5% of the sets tested):

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 1.43 | 3.83 | 3.44 | 4.01 | 7.33 |
| pre vs post, HLRT | 4.76 | 4.04 | 4.23 | 3.73 | 4.67 |
| post, BFR vs HLRT | 0.95 | 0.64 | 0.37 | 1.16 | 0.67 |
| change, BFR vs HLRT | 0.95 | 0.00 | 0.14 | 0.53 | 0.67 |
| baseline | 0.95 | 0.43 | 0.47 | 0.48 | 0.00 |

After BH, 9 sets survive on BFR training, 49 on HLRT and none elsewhere. Protein-level
classification is `02_Differential_Expression/04_Protein_Classification`.
