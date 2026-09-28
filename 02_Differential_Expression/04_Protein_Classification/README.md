# 02_Differential_Expression / 04_Protein_Classification

This step measures how well each protein separates the study groups, using the same five paired
tasks as `03_Pathway_Enrichment/05_set_classification`. It reads
`01_Preprocess/02_Quantification/c_data/proteins.rds`, the DE results in
`02_Differential/c_data/02_differential.xlsx` and the correlations in
`03_Protein_Association/c_data/03_protein_association.xlsx`. It writes
`c_data/04_protein_classification.xlsx` one ROC PDF per task
(`04_classification_pre_vs_post_BFR.pdf` and so on) and `04_chance.pdf`.

```sh
Rscript 02_Differential_Expression/04_Protein_Classification/a_script/04_protein_classification.R
```

AUC comes from `pROC` with `direction = "<"` fixed, and p from the paired Wilcoxon test, with BH
within task. The baseline control is reported but not drawn.

| Task | Nominal | Expected | Ratio | BH < 0.05 |
|---|---:|---:|---:|---:|
| Pre to post, BFR | 477 | 152 | 3.14 | 53 |
| Pre to post, HLRT | 664 | 152 | 4.37 | 220 |
| Baseline BFR vs HLRT *(control)* | 140 | 152 | 0.92 | 0 |
| Post BFR vs HLRT | 172 | 152 | 1.13 | 0 |
| Change, BFR vs HLRT | 120 | 152 | 0.79 | 0 |

No protein separates BFR from HLRT beyond chance. The best protein on the untrained baseline legs
reaches AUC 0.82, higher than the best after training (0.78).

The `candidates` sheet lists the 45 proteins that pass DE at FDR 0.05 in any contrast and also reach
nominal p for association and classification. Nine of them also separate the arms after training or
in change: NAMPT, GSTK1, ATP1A1, MSN, NPEPPS, SUCLG1, ACTR2, ART3 and ALDH1B1. Each protein gets
eight association tests, so about 44 of the 130 DE hits would clear that bar by chance; treat the
list as leads.
