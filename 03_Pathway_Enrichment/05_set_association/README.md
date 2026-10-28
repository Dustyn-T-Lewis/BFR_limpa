# 05 · Set association

Correlates each set's singscore change with the same leg's change in muscle size and echo
intensity, pooled over 65 legs and within participant (BFR minus HLRT) over 32 participants. It
also tests sets on the protein association slope from
`02_Differential_Expression/03_Protein_Association`: fgsea on its t ranking for leads, fry on the
same designs for confirmation.

Inputs:

- `03_run_singscore/c_data/set_scores.rds`;
- the `set_catalog` and `protein_gene_map` sheets and `gene_sets.rds` from
  `00_build_gene_sets/c_data/`;
- the targets in `01_Preprocess/02_Quantification/c_data/proteins.rds`;
- `00_Input/phenotype.csv`;
- `association_t.rds` and `association_designs.rds` from
  `02_Differential_Expression/03_Protein_Association/c_data/`.

```sh
quarto render 03_Pathway_Enrichment/05_set_association/a_script/05_set_association.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/05_set_association.html`: the walkthrough.
- `b_reports/05_association_<collection>.pdf`, one per collection: every set-outcome pair
  reaching nominal p, smallest p first. S11 (pooled) plots each leg, with a line fitted in each
  arm and the within-arm correlations inset; S12 (differential) plots the BFR minus HLRT
  differences that were tested, one point per participant, with a least-squares line. Each panel
  header gives the tested Spearman r, its p and BH q within collection and outcome.
- `b_reports/05_chance_by_collection.pdf`: S13, nominal hits against chance by collection. Bar
  length is observed nominal hits over the null expectation, red above 1 and grey below; labels
  give observed of tested.
  Features co-vary, so counts spread widely around chance; the shuffle range is in
  `02_Differential_Expression/03_Protein_Association`.
- `c_data/05_set_association.xlsx`, read by no later step:
  - `chance_expectation`: one row per analysis, outcome and collection; nominal hits against chance.
  - `set_association`: one row per set, analysis and outcome; Spearman r, p and BH q, smallest p
    first.
  - `set_by_arm`: one row per set, arm and outcome; the same correlation inside each arm.
  - `association_sets`: one row per set, analysis and outcome; fgsea NES, p, BH q and leading edge,
    and fry direction, p and FDR. fgsea treats proteins as independent, so a set counts when fry
    confirms it.
  - `set_catalog`: one row per tested set, with its collection and measured size.

No set's singscore tracks the phenotype after correction, and no protein does either
(`02_Differential_Expression/03_Protein_Association`). On the protein association ranking, fgsea
calls up to 149 sets per outcome (mitochondrial translation and respiration top pooled rf_csa),
and fry confirms none (smallest FDR 0.748).
