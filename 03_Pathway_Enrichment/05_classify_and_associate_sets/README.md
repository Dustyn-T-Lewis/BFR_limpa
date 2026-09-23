# 05_classify_and_associate_sets

How well each set's score separates study groups, and whether it tracks the phenotype.

| | |
|---|---|
| **Reads** | `gene_sets.rds`, `set_tests.rds`, `singscore.rds`, `proteins.rds`, `phenotype.csv` |
| **Writes** | `set_results.rds`, `05_classify_and_associate_sets.xlsx`, 7 figures over 96 pages |

Five tasks: pre against post in each arm, BFR against HLRT at baseline, after training, and in
change. All are paired, so AUC comes from `pROC` with `direction = "<"` fixed (unfixed, it flips
below-chance sets) and p from the paired Wilcoxon test. Associations use Spearman `cor.test`,
pooled over 65 legs and differential over 32 participants.

Each collection is read against its own chance count (5% of its sets), with BH within collection
and task. The table is in `../README.md`.

Figures: an ROC panel for every set reaching nominal p on each task (control excluded), a scatter
for every set-outcome pair reaching nominal p, and nominal hits over chance by collection. Panels
run 16 to a page, by collection then p, with BH q in each header; paged figures are PDF only.
Most nominal hits are expected by chance; read them against the chance figure.

```r
sr <- readRDS("03_Pathway_Enrichment/05_classify_and_associate_sets/c_data/set_results.rds")
sr$chance_expectation # nominal against chance, per collection
sr$set_auc            # set x task: AUC, paired p, BH
sr$set_association    # set x outcome, pooled and differential
sr$by_arm             # correlations within each arm
```
