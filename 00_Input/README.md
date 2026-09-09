# 00 · Input

Study data. Nothing here runs.

| File | What one row is | Used by |
|---|---|---|
| `report.parquet` | one peptide in one MS run | `01_Preprocess` |
| `metadata.csv` | one MS sample, 131 rows | `01_Preprocess` |
| `phenotype.csv` | one leg, 70 rows | later stages |

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
`metadata.csv` holds all three and `01_Preprocess` reconciles them once.

`phenotype.csv` stays at one row per leg because training volume is one number for the whole
leg. If we stored it once per MS sample it would be repeated, and summing the column would
report twice the real volume.
