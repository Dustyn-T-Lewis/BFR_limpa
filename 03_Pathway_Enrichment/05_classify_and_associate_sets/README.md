# 05_classify_and_associate_sets

How well each set's score separates study groups, and whether it tracks the phenotype.

| | |
|---|---|
| Reads | `00_build_gene_sets.xlsx`, `set_scores.rds`, `proteins.rds`, `phenotype.csv` |
| Writes | `05_classify_and_associate_sets.xlsx`, `05_classify_and_associate_sets_figures.pdf` (S13–S19, 128 pages) |

Five tasks: pre against post in each arm, BFR against HLRT at baseline, after training, and in
change. All are paired, so AUC comes from `pROC` with `direction = "<"` fixed (unfixed, it flips
below-chance sets) and p from the paired Wilcoxon test. Associations use Spearman `cor.test`,
pooled over 65 legs and differential over 32 participants.

Each collection is read against its own chance count (5% of its sets), with BH within collection
and task. The table is in `../README.md`.

Figures: S13–S16 draw an ROC panel for every set reaching nominal p on each task (control
excluded), S17–S18 a scatter for every set-outcome pair reaching nominal p, pooled then
differential, and S19 nominal hits over chance by collection. Panels run twelve to an A4 page, by
collection then p, with BH q in each header. Most nominal hits are expected by chance; read them
against S19.

The workbook holds `chance_expectation`, `set_auc`, `set_association`, `set_by_arm` and
`set_catalog`, listed with their row counts on its `overview` sheet.
