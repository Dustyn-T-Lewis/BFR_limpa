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
- `b_reports/04_classification_<collection>.pdf`, one per collection: S6 to S9, one per task
  except the baseline, with an ROC panel for every set reaching nominal p, smallest p first, BH q
  in each header. Red curves separate higher in the second group (post, or BFR), blue lower;
  shading is the area under the curve and the dashed line chance. p is the paired Wilcoxon
  signed-rank test, q is BH within collection and task.
- `b_reports/04_chance_by_collection.pdf`: S10, nominal hits against chance by collection and task.
  Bar length is observed nominal hits over the null expectation, red above 1 and grey below; labels
  give observed of tested.
  Features co-vary, so counts spread widely around chance; the shuffle range is in
  `02_Differential_Expression/03_Protein_Association`.
- `c_data/04_set_classification.xlsx`, one sheet per table, read by no later step:
  - `chance_expectation`: one row per task and collection, with nominal hits against chance and BH
    hits.
  - `set_auc`: one row per set and task, with AUC, paired Wilcoxon p and BH q, smallest p first.
  - `set_catalog`: one row per tested set, with its collection and measured size.

Nominal hits as a multiple of chance (5% of the sets tested):

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 2.44 | 3.40 | 3.72 | 4.07 | 7.46 |
| pre vs post, HLRT | 4.39 | 4.47 | 5.21 | 4.10 | 4.75 |
| post, BFR vs HLRT | 0.98 | 0.85 | 0.33 | 1.04 | 0.68 |
| change, BFR vs HLRT | 1.95 | 0.64 | 0.28 | 0.82 | 1.69 |
| baseline | 1.46 | 1.06 | 0.51 | 1.26 | 0.68 |

Set scores separate pre from post training in both arms: 30 sets pass BH with BFR and 112 with
HLRT. With HLRT the strongest sets are rises in tRNA aminoacylation, mitochondrial protein
degradation and mitochondrial RNA metabolism (AUC 0.69 to 0.71). With BFR mitochondrial protein
degradation also rises (AUC 0.74), while growth factor response and NF-κB signalling sets fall
(AUC 0.22 to 0.29). Between the arms one set passes BH, HALLMARK_SPERMATOGENESIS on change (AUC
0.68, q 0.018). The nominal leads on change are NF-κB signalling and programmed cell death, both
lower in BFR's change (AUC 0.27), and after training tRNA processing sits higher in BFR (AUC 0.67,
q 0.14). The between-arm rows of the table sit near chance, so the single BH hit among the many
sets tested and these leads need confirmation. Protein-level classification is
`02_Differential_Expression/04_Protein_Classification`.
