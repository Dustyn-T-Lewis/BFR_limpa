# 05_classify_and_associate_sets

Measures how well each set's score separates study groups, and whether it tracks the phenotype.

- Reads: `00_build_gene_sets/c_data/00_build_gene_sets.xlsx` (`set_catalog`),
  `04_run_singscore/c_data/set_scores.rds`, `01_Preprocess/02_Quantification/c_data/proteins.rds`,
  `00_Input/phenotype.csv`
- Writes: `c_data/05_classify_and_associate_sets.xlsx` (`overview`, `chance_expectation`,
  `set_auc`, `set_association`, `set_by_arm`, `set_catalog`),
  `b_reports/05_classify_and_associate_sets_figures.pdf` (S13–S19, 128 pages)
- Run: `Rscript 03_Pathway_Enrichment/05_classify_and_associate_sets/a_script/05_classify_and_associate_sets.R`,
  49 seconds.

Five tasks: pre against post in each arm, BFR against HLRT at baseline, after training, and in
change. All are paired, so AUC comes from `pROC` with `direction = "<"` fixed (unfixed, it flips
below-chance sets) and p from the paired Wilcoxon test. Associations use Spearman `cor.test`,
pooled over 65 legs and differential over 32 participants.

## Only training clears chance

Nominal hits over chance, per collection. Each collection is read against its own chance count
(5% of its sets), with BH within collection and task.

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| pre vs post, BFR | 1.4 | 3.8 | 3.4 | 4.0 | 7.3 |
| pre vs post, HLRT | 4.8 | 4.0 | 4.2 | 3.7 | 4.7 |
| post, BFR vs HLRT | 0.95 | 0.64 | 0.37 | 1.16 | 0.67 |
| change, BFR vs HLRT | 0.95 | 0.00 | 0.14 | 0.53 | 0.67 |
| baseline *(control)* | 0.95 | 0.43 | 0.47 | 0.48 | 0.00 |

Under BH, 9 sets survive on BFR training and 49 on HLRT, and none elsewhere. No set tracks the
phenotype after correction, and no protein does either (`02_Differential_Expression/03_Phenotype`).

S13–S16 draw an ROC panel for every set reaching nominal p on each task (control excluded),
S17–S18 a scatter for every set-outcome pair reaching nominal p, pooled then differential, and S19
nominal hits over chance by collection. Panels run twelve to an A4 page, by collection then p,
with BH q in each header. Read every panel against S19: at a ratio of 1, the nominal hits are what
chance returns.
