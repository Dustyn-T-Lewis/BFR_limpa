# 00 · Input

Study data. 

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

## Notes

**Treatment is not the leg letter.** Run names contain `D` and `E` (right and left). Which leg
got restriction was counterbalanced, so the letter gives the side only. `metadata.csv` carries
`treatment` as its own column, done by hand.

**The same sample has three names**: the MS run name, our sample ID, and participant plus leg.
`metadata.csv` holds all three and `01_Filtering` reconciles them once.

`phenotype.csv` stays at one row per leg because training volume is one number for the whole
leg. Stored once per MS sample it would be repeated, and summing the column would report twice
the real volume.

### Correction (2026-09-19): vasto lateral pre-CSA, 8 of 70 rows

The CV between the two pre-training VL CSA measurements was high. Ultrasound images were re-analyzed, 
`vl_csa_pre_cm2` changed for 8 legs (participants S07, S10, S11, S18 ×2, S21, S24, S29). Nothing else 
in the file changed — pre VL CSA only

