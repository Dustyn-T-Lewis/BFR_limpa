# 05_classify_and_associate_sets

This step measures how well each set's score separates the study groups, and whether it tracks the
phenotype. It reads the scores from `04_run_singscore/c_data/set_scores.rds`, the set catalogue
from `00_build_gene_sets/c_data/00_build_gene_sets.xlsx`, the sample targets from
`01_Preprocess/02_Quantification/c_data/proteins.rds` and the outcomes from `00_Input/phenotype.csv`.
Tables go to `c_data/05_classify_and_associate_sets.xlsx`, whose first table is the chance
comparison below, and figures S13–S19, 128 pages of them, to
`b_reports/05_classify_and_associate_sets_figures.pdf`.

```sh
Rscript 03_Pathway_Enrichment/05_classify_and_associate_sets/a_script/05_classify_and_associate_sets.R
```

There are five tasks: pre against post in each arm, and BFR against HLRT at baseline, after
training and in change. All are paired. AUC comes from `pROC` with `direction = "<"` fixed, since
left free it flips sets that fall below chance, and p comes from the paired Wilcoxon test.
Associations with the phenotype use Spearman `cor.test`, pooled over 65 legs and differential over
32 participants.

The table gives nominal hits as a multiple of chance. Each collection is read against its own
chance count, 5% of its sets, and BH runs within each collection and task.

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 1.4 | 3.8 | 3.4 | 4.0 | 7.3 |
| pre vs post, HLRT | 4.8 | 4.0 | 4.2 | 3.7 | 4.7 |
| post, BFR vs HLRT | 0.95 | 0.64 | 0.37 | 1.16 | 0.67 |
| change, BFR vs HLRT | 0.95 | 0.00 | 0.14 | 0.53 | 0.67 |
| baseline *(control)* | 0.95 | 0.43 | 0.47 | 0.48 | 0.00 |

After BH, 9 sets survive on BFR training and 49 on HLRT, and none on any other task. No set tracks
the phenotype after correction, and no protein does either
(`02_Differential_Expression/03_Phenotype`).

S13–S16 draw an ROC panel for every set reaching nominal p on each task except the control.
S17–S18 draw a scatter for every set and outcome pair reaching nominal p, pooled and then
differential, and S19 shows nominal hits against chance by collection. Panels run twelve to an A4
page, smallest p first, with the BH q in each header. Read every panel against
S19: where the ratio is 1, the nominal hits are what chance returns.
