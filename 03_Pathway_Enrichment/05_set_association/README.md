# 06 · Set association

Correlates each set's singscore change with the same leg's change in muscle size and echo
intensity, pooled over 65 legs and within participant (BFR minus HLRT) over 32 participants.

Inputs: `03_run_singscore/c_data/set_scores.rds`, the `set_catalog` sheet of
`00_build_gene_sets/c_data/00_build_gene_sets.xlsx`, the targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds`, and `00_Input/phenotype.csv`.

```sh
quarto render 03_Pathway_Enrichment/05_set_association/a_script/05_set_association.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/05_set_association.html`: the walkthrough.
- `b_reports/05_association_<collection>.pdf`, one per collection: every set-outcome pair
  reaching nominal p, smallest p first. S17 (pooled) plots each leg, with the within-arm
  correlations inset; S18 (differential) plots the BFR minus HLRT differences that were tested,
  one point per participant.
- `b_reports/05_chance_by_collection.pdf`: S19, nominal hits against chance by collection.
- `c_data/05_set_association.xlsx`, read by no later step:
  - `chance_expectation`: one row per analysis, outcome and collection; nominal hits against chance.
  - `set_association`: one row per set, analysis and outcome; Spearman r, p and BH q, smallest p
    first.
  - `set_by_arm`: one row per set, arm and outcome; the same correlation inside each arm.
  - `set_catalog`: one row per tested set, with its collection and measured size.

No set tracks the phenotype after correction, and no protein does either
(`02_Differential_Expression/03_Protein_Association`).
