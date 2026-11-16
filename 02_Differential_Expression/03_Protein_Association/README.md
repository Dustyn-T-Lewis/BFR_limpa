# 03 · Protein association

Tests whether each protein's training change tracks the same leg's outcome: change in VL and RF
cross-sectional area (CSA) and echo intensity, and total volume load lifted. Each slope is a
coefficient in a limpa model (`dpcDE()`, `eBayes(robust = TRUE)`), with BH within each outcome and
coefficient. The average and difference t values rank proteins for gene set tests in
`03_Pathway_Enrichment/05_set_association`.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds` (models) and `proteins_adjusted.rds`
(figures), `00_Input/phenotype.csv`.

```sh
quarto render 02_Differential_Expression/03_Protein_Association/a_script/03_protein_association.qmd --output-dir ../b_reports
```

## Models

```
arm:    ~ 0 + leg_id + treatment:post + treatment:post:x + blood
within: ~ 0 + leg_id + participant:post + bfr_post + post:x + blood
```

`post` is 1 at T2; `x` is the leg's outcome, centred on the mean leg.

| Term | What it does |
|---|---|
| `leg_id` | own baseline per leg, so only change within a leg counts |
| `treatment:post` | mean training change per arm, so an arm difference cannot pose as a slope |
| `treatment:post:x` | one slope on x per arm |
| `participant:post` | mean training change per participant (within model) |
| `bfr_post` | extra training change in the BFR leg (within model) |
| `post:x` | one slope on x, read only from a participant's two legs (within model) |
| `blood` | blood index, as in `02_Differential` |

| Coefficient | Model | Question |
|---|---|---|
| BFR | arm | Across BFR legs, does protein change track the outcome? |
| HLRT | arm | Across HLRT legs, does it? |
| difference | arm | Is the BFR slope different from the HLRT slope? |
| average | arm | Averaged over both arms, does protein change track the outcome? |
| within | within | Within a participant, does the leg that changed more change more in protein? |

The arm model uses all 65 legs with both samples (BFR 32, HLRT 33); the within slope comes from
the 64 legs of the 32 participants with both legs. Volume load has no within slope: dose differs by
arm by design, so a participant's leg difference is the arm difference.

## Outputs

- `b_reports/03_protein_association.html`: the walkthrough.
- `b_reports/03_association_<outcome>.pdf`: one page per outcome, the three smallest p for each
  coefficient. Arm coefficients plot each leg's protein change against its outcome with a line per
  arm; the within row plots BFR minus HLRT for both, one point per participant. Panel headers give
  the tested coefficient's t, p and BH q.
- `c_data/association_t.rds`: moderated t, one row per protein, one column per ranking
  (`average_<outcome>`, `difference_<outcome>`), for fgsea in `05_set_association`.
- `c_data/association_designs.rds`: per ranking, the arm design and its contrast vector, for fry
  in `05_set_association`.
- `c_data/03_protein_association.xlsx`:
  - `outcomes`: one row per outcome and coefficient; legs used, formula, covariates, nominal and
    BH counts, best FDR, top five genes with direction;
  - `top_leads`: the 25 smallest p per outcome and coefficient; slope, t, p and FDR;
  - `associations`: every protein, outcome and coefficient; slope, t, p and FDR.

## Results

Volume load is the only outcome with BH hits. Across BFR legs, CA14 (t −4.9, q 0.019) and TOMM40
(t −4.6, q 0.022) fall more and HSPA12A (t 4.7, q 0.019) rises more in legs that lifted more. The
average slope keeps CA14 (q 0.004) and TOMM40 (q 0.018). HSPA12A is the one BFR versus HLRT
difference (t 4.8, q 0.034): it rises with dose in BFR legs and not in HLRT legs. Nothing tracks
volume load across HLRT legs.

No protein tracks a change in CSA or echo after BH. Leads by smallest p (+ rises, − falls with
the outcome; for the difference, + means the BFR slope is the higher):

| Outcome | BFR | HLRT | Difference (BFR - HLRT) |
|---|---|---|---|
| vl_csa | TJP1−, PPP2R5E−, MRPS16+ | PTCD3+, DDT+, JPH1− | PPP2R5E−, MRPS16+, GNA13− |
| vl_echo | PHKG1+, COA4+, AP2A2− | MEMO1−, UBE2L3+, COL4A1− | SCARB2−, GLG1−, FASN− |
| rf_csa | TJP1−, MRPS16+, DNAJA2− | NMNAT1+, MTHFS+, AGMAT+ | PICALM−, MPC1−, MTHFS− |
| rf_echo | ITGA1−, PLCB1−, FKBP5− | NUBPL−, PSMC3+, ATG7+ | MRPL24+, NUBPL+, SPCS2+ |

The closest change lead is AGMAT on RF CSA (average, t 4.4, q 0.11); SMTNL1 falls in the leg that
grew more in RF CSA within participant (t −4.3, q 0.40). Nominal hits sit near the 152 expected by
chance from 3,042 proteins, except VL echo in HLRT legs (305) and its BFR versus HLRT difference
(361), where the best q is 0.23.

Nominal p < 0.05 and BH < 0.05 per outcome and coefficient:

| Outcome | BFR | HLRT | Difference | Average | Within |
|---|---:|---:|---:|---:|---:|
| vl_csa | 81 / 0 | 155 / 0 | 84 / 0 | 138 / 0 | 67 / 0 |
| vl_echo | 186 / 0 | 305 / 0 | 361 / 0 | 164 / 0 | 101 / 0 |
| rf_csa | 119 / 0 | 85 / 0 | 88 / 0 | 125 / 0 | 129 / 0 |
| rf_echo | 43 / 0 | 179 / 0 | 100 / 0 | 106 / 0 | 116 / 0 |
| volume_load_t | 415 / 3 | 96 / 0 | 295 / 1 | 411 / 2 | not fitted |

Proteins move together, so nominal counts swing widely under no association (63 to 308 in 40
within-arm shuffles of VL CSA, before the blood covariate), and the volume-load slopes compare
people, so strength, size and sex ride along with dose.
