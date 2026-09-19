# 00 · Input

Study data. Nothing here runs.

| File | What one row is | Read by |
|---|---|---|
| `report.parquet` | one precursor in one MS run | `01_Preprocess/01_Filtering` |
| `metadata.csv` | one MS sample, 131 rows | `01_Preprocess/01_Filtering` |
| `phenotype.csv` | one leg, 70 rows | nothing yet |

## Get report.parquet

Too large for git, so it is attached to a release. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet
```

## Two things to get right

**Treatment is not the leg letter.** Run names contain `D` and `E` (right and left). Which leg
got restriction was counterbalanced, so the letter gives the side only. `metadata.csv` carries
`treatment` as its own column, done by hand.

**The same sample has three names**: the MS run name, our sample ID, and participant plus leg.
`metadata.csv` holds all three and `01_Filtering` reconciles them once.

`phenotype.csv` stays at one row per leg because training volume is one number for the whole
leg. Stored once per MS sample it would be repeated, and summing the column would report twice
the real volume.

## Correction (2026-09-19): vasto lateral pre-CSA, 8 of 70 rows

The coefficient of variation between the two pre-training vasto lateral CSA measurements was
too high. A few of the underlying ultrasound images were re-analyzed, and `vl_csa_pre_cm2`
changed for 8 legs (participants S07, S10, S11, S18 ×2, S21, S24, S29). Nothing else in the file
changed — not post, not echo, not rectus femoris, not training volume.

The baseline difference between the BFR and HLRT legs (paired, same person, opposite legs) was
already there before this fix and is still there after it: BFR legs average about 1.15 cm²
smaller than HLRT legs at pre (paired t-test, p≈0.02 either way). The fix did not create or
remove this imbalance, and the trained groups still show equivalent percent gains (p≈0.72), so
it does not change the study's conclusions.

This has not been discussed with the group yet. Anyone who prefers the previous numbers should
say so — `git log -p -- 00_Input/phenotype.csv` shows exactly which 8 values moved and by how
much.
