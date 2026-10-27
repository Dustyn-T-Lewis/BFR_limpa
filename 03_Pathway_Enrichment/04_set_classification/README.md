# 04 · Set classification

Measures how well each set's singscore separates the study groups in five paired tasks: pre against
post in each arm, and BFR against HLRT at baseline, after training and in change.

Inputs: `03_run_singscore/c_data/set_scores.rds`, the `set_catalog` sheet of
`00_build_gene_sets/c_data/00_build_gene_sets.xlsx`, and the targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds`.

```sh
quarto render 03_Pathway_Enrichment/04_set_classification/a_script/04_set_classification.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/04_set_classification.html`: the walkthrough.
- `b_reports/04_classification_<collection>.pdf`, one per collection: S11 to S14, an ROC panel for
  every set reaching nominal p on each task except the baseline, smallest p first, BH q in each
  header. Red curves separate higher in the second group (post, or BFR), blue lower; shading is
  the area under the curve and the dashed line chance. p is the paired Wilcoxon signed-rank test,
  q is BH within collection and task.
- `b_reports/04_chance_by_collection.pdf`: S15, nominal hits against chance by collection and task.
  Bar length is observed nominal hits over the null expectation, red above 1 and grey below; labels
  give observed of tested.
- `c_data/04_set_classification.xlsx`, one sheet per table, read by no later step:
  - `chance_expectation`: one row per task and collection, with nominal hits against chance and BH
    hits.
  - `set_auc`: one row per set and task, with AUC, paired Wilcoxon p and BH q, smallest p first.
  - `set_catalog`: one row per tested set, with its collection and measured size.

Nominal hits as a multiple of chance (5% of the sets tested):

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 1.46 | 3.83 | 3.44 | 4.01 | 7.46 |
| pre vs post, HLRT | 4.88 | 4.04 | 4.23 | 3.73 | 4.75 |
| post, BFR vs HLRT | 0.98 | 0.64 | 0.37 | 1.16 | 0.68 |
| change, BFR vs HLRT | 0.98 | 0.00 | 0.14 | 0.53 | 0.68 |
| baseline | 0.98 | 0.43 | 0.47 | 0.48 | 0.00 |

After BH, 9 sets survive on BFR training, 49 on HLRT and none elsewhere. Protein-level
classification is `02_Differential_Expression/04_Protein_Classification`.
