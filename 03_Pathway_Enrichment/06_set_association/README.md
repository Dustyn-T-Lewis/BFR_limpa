# 06 · Set association

Correlates each set's singscore change with the same leg's change in muscle size and echo
intensity, pooled over 65 legs and within participant (BFR minus HLRT) over 32 participants.

Inputs: `04_run_singscore/c_data/set_scores.rds`, the `set_catalog` sheet of
`00_build_gene_sets/c_data/00_build_gene_sets.xlsx`, the targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds`, and `00_Input/phenotype.csv`.

```sh
quarto render 03_Pathway_Enrichment/06_set_association/a_script/06_set_association.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/06_set_association.html`: the walkthrough.
- `b_reports/06_association_<collection>.pdf`, one per collection: S17 (pooled) and S18
  (differential), every set-outcome pair reaching nominal p, smallest p first, with the
  within-arm correlations inset.
- `b_reports/06_chance_by_collection.pdf`: S19, nominal hits against chance by collection.
- `c_data/06_set_association.xlsx`: `chance_expectation`, `set_association` (Spearman r, p and
  BH q), `set_by_arm` (the same correlation inside each arm, descriptive) and `set_catalog`. Read
  by no later step.

No set tracks the phenotype after correction, and no protein does either
(`02_Differential_Expression/03_Protein_Association`).
