# 01 · Plasticity

Compares the BFR and HLRT training responses as wholes, in five parts.

1. Magnitude. Hit counts per arm and direction; a cross-fitted magnitude (as in the YvO Figure 4C):
   in each of 200 random splits half the participants choose the proteins (all, the top 500 or 250,
   or as many as the choosing arm's BH hits) by one arm's response, and the other half measure both
   arms' median |log2 change| on them, once choosing by HLRT and once by BFR. The swap p compares
   the BFR to HLRT ratio with 200 data sets in which each participant's arm labels are swapped at
   random. A signature score projects each leg's change onto the training direction learned from
   the other participants, over the proteins whose average training change across both arms passes
   BH 0.05; BFR against HLRT is a paired t and a 5,000-flip sign test. fry tests the interaction
   over the same proteins, split into those that rise and those that fall.
2. Direction. BFR against HLRT training change for proteins at FDR 0.05 in either arm, with sign
   concordance, quadrant counts and quadrant over-representation (`fgsea::fora()`, Hallmark, GO
   Slim and Reactome, universe every measured gene). The Deming slope is drawn as a descriptive
   line only.
3. Coherence. RRHO2 of the two training rankings (signed -log10 p); fgsea NES per set and
   main-network module eigengene change, BFR against HLRT; fry of HLRT's up and down training hits
   in the BFR training contrast.
4. Classification and leads. sPLS-DA AUCs from `07_Multivariate/02_splsda`; the interaction
   volcano with the sPLS-DA change panel circled; the interaction sets with the smallest fgsea p.
5. Response against growth. Per leg, the signature score, MPE change and slow-fraction change
   (T2 minus T1) against the change in VL and RF CSA and echo intensity (post minus pre). Spearman
   rho within each arm, and pooled within participant after centring each participant's two legs,
   with p from 5,000 random swaps of the response between a participant's legs.

Changes for the cross-fit, signature score and growth come from `proteins_adjusted.rds` (blood
effect removed); fry uses `proteins.rds` with the stage 02 design. The cross-fit and signature test
use the 32 participants with all four biopsies; growth uses all 65 legs, scoring the one leg of
the participant missing a biopsy on the direction of all 32.

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
- `b_reports/01_plasticity.pdf`, five A4 pages: magnitude (hit counts, signature scores, cross-fit
  bars labelled with ratio, split interval and swap p); direction (scatter with Deming and identity
  lines, top three sets per quadrant); coherence (RRHO2 map, NES scatter, module scatter, rank
  histogram of HLRT hits in the BFR ranking); classification and leads (AUCs, interaction sets,
  volcano, interaction of the average-response sets); response against growth (3 x 4 scatter grid
  with rho and p per arm and within participant).
- `c_data/01_plasticity.xlsx`:
  - `counts`: one row per arm and threshold (FDR or nominal 0.05): proteins up and down;
  - `magnitude`: one row per choosing arm and set (all proteins once, `selected_on` "none"):
    proteins, BFR and HLRT median |log2 change| over splits, ratio, `split_low` and `split_high`
    (the middle 95% of splits, not a confidence interval) and swap p;
  - `legs`: one row per leg: arm, signature score, MPE change, slow-fraction change, the four
    outcome changes and whether the participant has both legs;
  - `direction`: one row per protein at FDR 0.05 in either arm: BFR and HLRT logFC and FDR, which
    arm passes, and quadrant;
  - `quadrant_sets`: the ten sets with the smallest p per quadrant: size, overlap, p and BH q;
  - `coherence`: one row per summary statistic: NES rho, module r and the two fry tests of HLRT
    hits in BFR training;
  - `classification`: one row per sPLS-DA comparison, component 1: participants, panel size, BER,
    AUC and its SD;
  - `leads`: one row per protein at nominal interaction p < 0.05: interaction logFC, p and FDR, BFR
    and HLRT logFC and p, and whether it sits in the sPLS-DA change panel;
  - `tests`: the signature test (proteins, mean scores, ratio, paired p, sign-flip p) and the two
    fry interaction tests (proteins, median average training change, median interaction,
    direction, p);
  - `growth`: one row per response, outcome and analysis (BFR, HLRT, within participant): rho, n
    and p.

The two arms run the same programme, and BFR runs it less far. Of 262 proteins at FDR 0.05 in either
arm (44 BFR only, 51 both, 167 HLRT only), 94% change in the same direction in both arms: 125 rise
in both and 121 fall in both, against 8 in each discordant quadrant. The shared rise is
mitochondrial (mitochondrial protein degradation, TCA cycle and oxidative phosphorylation, BH q 3e-8
to 1e-7); no other quadrant has a set at BH. Set NES agree across arms (Spearman rho 0.85 over 530
sets), as do main-network module changes (r 0.83 over 42 modules), and HLRT's training hits move the
same way in BFR (fry, up set p 2e-11, down set p 7e-14). By signature score BFR moves 0.77 as far as
HLRT (paired p 0.031, sign-flip p 0.030). Cross-fitted, BFR's median change is 0.75 to 0.94 of
HLRT's across sets and both choosing arms (swap p 0.21 to 0.35), and fry on the interaction leans
the same way without passing (rising set median -0.037, p 0.36; falling set +0.050, p 0.13).

Between the arms, sPLS-DA sits near chance (AUC 0.67 before training, 0.60 after, 0.58 on change)
against 0.88 and 1.00 for training within each arm. The interaction has 126 leads at p < 0.05 (54
larger with BFR, 72 with HLRT), led by MSN, KLHL30, ARL8A, RPL36 and FAHD2A, and 45 of the 50
change-panel proteins are among them. fgsea calls 10 interaction sets at BH, all rising less or
falling more with BFR (vesicle-mediated transport, cell motility, cell adhesion,
epithelial-mesenchymal transition, apical junction, q 0.0008 to 0.005), but fry passes none (min q
0.96), so these are leads.

No response summary tracks growth within an arm (all |rho| <= 0.32, p >= 0.07). Pooled within
participant, the leg with the larger signature score had the smaller VL echo intensity change (rho
-0.49, p 0.003); the other 11 within-participant tests give |rho| <= 0.28 and p >= 0.12. That single
echo correlation is one of 36 tests and is unconfirmed by the per-arm tests, so treat it as a lead.
