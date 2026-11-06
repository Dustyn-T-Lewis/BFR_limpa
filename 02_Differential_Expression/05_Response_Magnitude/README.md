# 05 · Response magnitude

Asks whether BFR training moved the proteome less far than HLRT training, every protein moving the
same way but less, which no single-protein interaction test can see. Hit counts (95 against 218)
cannot answer it: a difference in significance is not a significant difference. Three tests do,
all within the 32 participants with all four biopsies:

1. Cross-fitted magnitude, as in the YvO Figure 4C. In each of 200 random splits, half the
   participants choose the proteins (all, the top 500 or 250, or as many as the choosing arm's BH
   hits) by one arm's response, and the other half measure both arms' median |log2 change| on
   them, so selection does not inflate the choosing arm. Run once choosing by HLRT and once by
   BFR. The swap p compares the BFR to HLRT ratio with 200 data sets in which each participant's
   arm labels are swapped at random.
2. Signature score. For each participant the training direction is the mean change across both
   arms in the other participants, over the proteins whose average training change across both
   arms passes BH 0.05 (a set fixed before any comparison). Each leg's change is projected onto it;
   BFR against HLRT is a paired t and a 5,000-flip sign test on one number per participant.
3. fry on the interaction contrast (BFR change minus HLRT change) over the same proteins, split
   into those that rise and those that fall with training, in the limpa model with the blood
   covariate. A smaller BFR response pushes the rising set negative and the falling set positive.

Changes come from `proteins_adjusted.rds` (blood effect removed); fry uses `proteins.rds` with the
stage 02 design.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds` and `proteins_adjusted.rds`,
`01_Design/c_data/design.rds`, the `contrast_summary` sheet of
`02_Differential/c_data/02_differential.xlsx`.

```sh
quarto render 02_Differential_Expression/05_Response_Magnitude/a_script/05_response_magnitude.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/05_response_magnitude.html`: the walkthrough.
- `b_reports/05_response_magnitude.pdf`, one A4 page: cross-fitted median |log2 change| per arm
  and set, chosen by HLRT and by BFR, each pair labelled with its BFR to HLRT ratio and swap p;
  below, each participant's BFR and HLRT signature scores.
- `c_data/05_response_magnitude.xlsx`:
  - `magnitude`: one row per choosing arm and set: proteins, BFR and HLRT median |log2 change|
    over splits, ratio with the middle 95% of splits, and swap p;
  - `signature_scores`: one row per participant: BFR and HLRT score;
  - `tests`: the signature test (proteins, mean scores, ratio, paired p, swap p) and the two fry
    tests (proteins, median average training change, median interaction, direction, p).

BFR training moves the proteome in the same direction as HLRT training and about a quarter less
far. The signature score puts BFR at 0.77 of HLRT (paired p 0.031, swap p 0.035). Cross-fitted,
BFR's median change is 0.75 to 0.94 of HLRT's across sets and both choosing arms, and HLRT changes
more even when BFR chooses the proteins. fry leans the same way: the rising set's interaction is
−0.04 against an average rise of 0.27 (p 0.36), the falling set's +0.05 against an average fall of
0.23 (p 0.13). This is an aggregate difference that single-protein tests cannot resolve, and only
the signature score passes 0.05; the cross-fitted ratios (swap p 0.23 to 0.34) and fry do not.
