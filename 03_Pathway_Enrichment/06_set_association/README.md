# 06_set_association

This step correlates each set's singscore change with the same leg's change in muscle size and
quality. It reads `04_run_singscore/c_data/set_scores.rds`, the catalogue in
`00_build_gene_sets/c_data/00_build_gene_sets.xlsx`, the sample targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds` and `00_Input/phenotype.csv`. It writes
`c_data/06_set_association.xlsx` one scatter PDF per collection (`06_association_Hallmark.pdf` and
so on) and `06_chance_by_collection.pdf`.

```sh
Rscript 03_Pathway_Enrichment/06_set_association/a_script/06_set_association.R
```

Correlations are Spearman `cor.test`, pooled over 65 legs and differential over 32 participants,
where the differential analysis correlates the BFR minus HLRT change in score with the same
difference in outcome. The `set_by_arm` sheet repeats each correlation inside BFR and inside HLRT,
for description only. BH runs within collection and outcome.

No set tracks the phenotype after correction, and no protein does either
(`02_Differential_Expression/03_Protein_Association`).

S17 draws every set-outcome pair reaching nominal p in the pooled analysis and S18 in the
differential one, twelve to an A4 page, smallest p first, with the within-arm correlations inset.
S19 plots nominal hits against chance by collection.
