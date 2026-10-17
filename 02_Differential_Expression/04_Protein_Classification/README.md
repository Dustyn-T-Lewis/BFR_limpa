# 04 · Protein classification

Measures how well each protein separates the groups in five paired tasks (AUC from `pROC` with
`direction = "<"`, p from the paired Wilcoxon signed-rank test, BH within task), the same tasks as
`03_Pathway_Enrichment/05_set_classification`.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential/c_data/02_differential.xlsx`, `03_Protein_Association/c_data/03_protein_association.xlsx`.

```sh
quarto render 02_Differential_Expression/04_Protein_Classification/a_script/04_protein_classification.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/04_protein_classification.html`: the walkthrough, with previews and the chance figure.
- `b_reports/04_classification_<task>.pdf`: ROC curves for every protein at nominal p, one file per
  task except baseline, which is reported but not drawn.
- `b_reports/04_chance.pdf`: nominal hits against chance per task.
- `c_data/04_protein_classification.xlsx`, one sheet per table:
  - `chance_expectation`: one row per task, with nominal hits against chance and BH hits.
  - `protein_auc`: one row per protein and task, with AUC, paired Wilcoxon p and BH q, smallest p
    first.
  - `candidates`: one row per protein that is DE at FDR 0.05 and reaches nominal p for association
    and classification.

| Task | Nominal | Expected | Ratio | BH < 0.05 |
|---|---:|---:|---:|---:|
| Pre to post, BFR | 477 | 152 | 3.14 | 53 |
| Pre to post, HLRT | 664 | 152 | 4.37 | 220 |
| Baseline BFR vs HLRT | 140 | 152 | 0.92 | 0 |
| Post BFR vs HLRT | 172 | 152 | 1.13 | 0 |
| Change, BFR vs HLRT | 120 | 152 | 0.79 | 0 |

No protein separates BFR from HLRT beyond chance. The `candidates` sheet lists 45 proteins that pass
DE at FDR 0.05 in any contrast and reach nominal p for association and classification; nine also
separate the arms after training or in change (NAMPT, GSTK1, ATP1A1, MSN, NPEPPS, SUCLG1, ACTR2,
ART3, ALDH1B1). Each protein gets eight association tests, so chance alone would put about 44 of
the 130 DE hits past the association bar, close to the 45 listed. Treat the list, the nine
included, as leads.
