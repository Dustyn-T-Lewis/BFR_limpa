# 01 · Plasticity

Compares the BFR and HLRT training responses as wholes: how far each leg moved (magnitude), whether
the legs moved the same proteins the same way (direction), whether rankings, sets and modules agree
(coherence), which BFR against HLRT leads remain, and whether any per-leg response summary tracks
growth. Magnitude uses a cross-fit over 200 participant splits (as in the YvO Figure 4C) with a
200-swap arm-label null, each leg's directional change over the proteins a time-only limpa model
calls at BH, and fry on the interaction.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds` and `proteins_adjusted.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`, `fit.rds` and the `DEP_matrix` sheet of
`02_Differential_Expression/02_Differential/c_data/02_differential.xlsx`,
`03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds`, the `set_tests` sheet of
`03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/01_run_fgsea_and_fry.xlsx`, the
`eigengene_tests` and `modules` sheets of `04_Network/01_modules/c_data/01_modules.xlsx`, the `mpe`
sheet of `05_Mitochondria/02_test_mito/c_data/02_test_mito.xlsx`, the `slow_fraction` sheet of
`06_Composition/02_fibre_type/c_data/02_fibre_type.xlsx`, the `models` and `panels` sheets of
`07_Multivariate/02_splsda/c_data/02_splsda.xlsx`, and `00_Input/phenotype.csv`.

```sh
quarto render 08_Plasticity/01_plasticity/a_script/01_plasticity.qmd --output-dir ../b_reports
```

The run takes about ten minutes, most of it the swap null.

Outputs:

- `b_reports/01_plasticity.html`: the walkthrough.
- `b_reports/01_plasticity.pdf`, five A4 pages: magnitude (hit counts, directional change, cross-fit
  bars); direction (BFR against HLRT scatter, top sets per quadrant); coherence (RRHO2 map, NES and
  module scatters, HLRT hits in the BFR ranking); classification and leads (AUCs, interaction sets,
  volcano, interaction of the training-responsive sets); response against growth (3 x 4 scatter
  grid).
- `c_data/01_plasticity.xlsx`:
  - `counts`: per arm and threshold (FDR or nominal 0.05), proteins up and down;
  - `magnitude`: per choosing arm and set, proteins, median |log2 change| per arm, ratio,
    `split_low` and `split_high` (middle 95% of splits, not a confidence interval), swap p;
  - `legs`: per leg, directional change, MPE change, slow-fraction change, outcome changes and
    whether the participant has both legs;
  - `direction`: per protein at FDR 0.05 in either arm, logFC and FDR per arm, hit and quadrant;
  - `quadrant_sets`: the ten sets with the smallest p per quadrant, size, overlap, p, BH q;
  - `coherence`: NES rho, module r and the two fry tests of HLRT hits in BFR training;
  - `classification`: per sPLS-DA comparison, component 1 participants, panel size, BER, AUC, SD;
  - `leads`: per protein at nominal interaction p < 0.05, interaction and per-arm statistics and
    whether it sits in the sPLS-DA change panel;
  - `tests`: the directional-change test (proteins, mean per arm, ratio, paired p) and the two fry
    interaction tests;
  - `growth`: per response, outcome and analysis (BFR, HLRT, within participant), rho, n, p.

The two arms run the same programme, and BFR may run it slightly less far. Of 262 proteins at FDR
0.05 in either arm (44 BFR only, 51 both, 167 HLRT only), 94% change in the same direction in both
arms: 125 rise in both and 121 fall in both, against 8 in each discordant quadrant. The shared rise
is mitochondrial (mitochondrial protein degradation, TCA cycle and oxidative phosphorylation, BH q
3e-8 to 1e-7); no other quadrant has a set at BH. Set NES agree across arms (Spearman rho 0.85 over
530 sets), as do main-network module changes (r 0.92 over 7 modules), and HLRT's training hits move
the same way in BFR (fry, up set p 2e-11, down set p 7e-14). BFR moves the training-responsive
proteome about 15% less than HLRT: over the 435 proteins the time-only model calls at BH, the mean
directional change is 0.258 log2 with BFR and 0.304 with HLRT (ratio 0.85, paired p 0.13).
Cross-fitted, BFR's median change is 0.75 to 0.94 of HLRT's (swap p 0.21 to 0.35), and fry on the
interaction leans the same way without passing (rising set median -0.032, p 0.45; falling set
+0.049, p 0.15). None of these passes, so the smaller BFR response is a lean, not a finding.

Between the arms, sPLS-DA sits near chance (AUC 0.67 before training, 0.60 after, 0.58 on change)
against 0.88 and 1.00 for training within each arm. The interaction has 126 leads at p < 0.05 (54
larger with BFR, 72 with HLRT), led by MSN, KLHL30, ARL8A, RPL36 and FAHD2A, and 45 of the 50
change-panel proteins are among them. fgsea calls 10 interaction sets at BH, all rising less or
falling more with BFR (vesicle-mediated transport, cell motility, cell adhesion,
epithelial-mesenchymal transition, apical junction, q 0.0008 to 0.005), but fry passes none (min q
0.96), so these are leads.

No response summary tracks growth within an arm (all |rho| <= 0.34, p >= 0.05). Pooled within
participant, the leg with the larger directional change had the smaller VL echo intensity change
(rho -0.51, p 0.0014); the other 11 within-participant tests give |rho| <= 0.22 and p >= 0.21. That
echo correlation is one of 36 tests and unconfirmed by the per-arm tests, so it is a lead.
