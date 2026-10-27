# 04 · Protein classification

Measures how well each protein separates a participant's BFR leg from their HLRT leg at baseline,
after training and in change (AUC from `pROC` with `direction = "<"`, p from the paired Wilcoxon
signed-rank test, BH within task). Pre to post within an arm is left to `02_Differential`;
`03_Pathway_Enrichment/04_set_classification` keeps those two tasks for sets.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential/c_data/02_differential.xlsx`, `03_Protein_Association/c_data/03_protein_association.xlsx`.

```sh
quarto render 02_Differential_Expression/04_Protein_Classification/a_script/04_protein_classification.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/04_protein_classification.html`: the walkthrough, with previews and the chance figure.
- `b_reports/04_classification_<task>.pdf`: ROC curves for every protein at nominal p, one file per
  task except baseline, which is reported but not drawn. Red curves separate higher in BFR, blue
  lower; shading is the area under the curve and the dashed line is chance. Panel headers give AUC,
  paired Wilcoxon p and BH q.
- `b_reports/04_chance.pdf`: nominal hits against chance per task. Bar length is nominal hits
  divided by the chance count; red bars exceed 1, and labels give observed of tested.
  Features co-vary, so under no effect the count swings widely around 5% (shuffled vl_csa gave
  63 to 308 protein hits against 152); only a large departure from 1 means anything.
- `c_data/04_protein_classification.xlsx`, one sheet per table:
  - `chance_expectation`: one row per task, with nominal hits against chance and BH hits.
  - `protein_auc`: one row per protein and task, with AUC, paired Wilcoxon p and BH q, smallest p
    first.
  - `candidates`: one row per protein that is DE at FDR 0.05 and reaches nominal p for association
    and classification.

| Task | Nominal | Expected | Ratio | BH < 0.05 |
|---|---:|---:|---:|---:|
| Baseline BFR vs HLRT | 140 | 152 | 0.92 | 0 |
| Post BFR vs HLRT | 172 | 152 | 1.13 | 0 |
| Change, BFR vs HLRT | 120 | 152 | 0.79 | 0 |

No protein separates BFR from HLRT beyond chance. The `candidates` sheet lists the proteins that
pass DE at FDR 0.05 in any contrast and reach nominal p for association and for separating the arms
after training or in change (NAMPT, GSTK1, ATP1A1, MSN, NPEPPS, SUCLG1, ACTR2, ART3, ALDH1B1). Each
protein gets eight association tests, so about a third of DE hits clear the association bar by
chance. Treat the list as leads.
