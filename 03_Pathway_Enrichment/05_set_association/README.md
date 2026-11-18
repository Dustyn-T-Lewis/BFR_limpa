# 05 · Set association

Correlates each set's singscore change with the same leg's change in muscle size and echo
intensity, over all 65 legs (`all_legs`) and within participant (`within_participant`, BFR minus
HLRT, 32 participants). It also tests sets on the protein association rankings of
`02_Differential_Expression/03_Protein_Association`: `average` (protein-on-outcome slope averaged
over BFR and HLRT) and `difference` (BFR slope minus HLRT slope), for the four change outcomes and
volume load. fgsea on the t ranking gives leads; fry on the same design and contrast confirms.

Inputs: `03_run_singscore/c_data/set_scores.rds`; `gene_sets.rds` and the `set_catalog` and
`protein_gene_map` sheets from `00_build_gene_sets/c_data/`; the targets in
`01_Preprocess/02_Quantification/c_data/proteins.rds`; `00_Input/phenotype.csv`;
`association_t.rds` and `association_designs.rds` from
`02_Differential_Expression/03_Protein_Association/c_data/`.

```sh
quarto render 03_Pathway_Enrichment/05_set_association/a_script/05_set_association.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/05_set_association.html`: the walkthrough.
- `b_reports/05_set_association.pdf`: S13 first, nominal hits over the chance expectation by
  collection, red above 1, labels giving observed of tested; then per collection (Hallmark,
  GO_Slim, KEGG_Legacy, Reactome, GOBP) every set-outcome pair at nominal p, smallest p first.
  S11 (all legs) plots each leg with a line per arm and the within-arm correlations inset; S12
  (within participant) plots the tested BFR minus HLRT differences, one point per participant.
  Panel headers give Spearman r, p and BH q within collection and outcome.
- `c_data/05_set_association.xlsx`, read by no later step:
  - `chance_expectation`: one row per analysis, outcome and collection: nominal hits against chance.
  - `set_association`: one row per set, analysis and outcome: Spearman r, p and BH q.
  - `set_by_arm`: one row per set, arm and outcome: the same correlation inside each arm.
  - `association_sets`: one row per set, ranking (`coefficient`) and outcome: fgsea NES, p, q and
    leading edge, fry direction, p and FDR.

Volume load is the one ranking fry confirms. On the average slope 59 sets reach fry FDR 0.048
(all at the same BH value): protein import into mitochondria and glucose and glycogen synthesis
fall, and cell junction assembly and immune regulation rise, in legs that lifted more. fgsea's
strongest volume-load leads, aerobic respiration and oxidative phosphorylation (NES −2.8), stop
short in fry (FDR 0.089 and 0.097). The difference ranking on volume load gives 215 fgsea sets
and no fry set (smallest FDR 0.093).

On the change outcomes fry confirms no set. On the average slope, legs whose RF CSA grew more show
larger rises in mitochondrial gene expression and translation (fgsea NES 3.0, 136 sets; fry FDR
0.82). The difference rankings give 4 to 231 fgsea sets per outcome; the closest to fry is VL echo
(148 sets, smallest fry FDR 0.23). No singscore change tracks an outcome after BH; the closest is
prostanoid metabolism against the BFR minus HLRT difference in RF echo (r −0.66, q 0.081).

fgsea counts proteins as independent and these sets overlap heavily, so its calls without fry
support are leads.
