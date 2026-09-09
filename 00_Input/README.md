# 00 · Input

The study data. Nothing here runs, and nothing computed belongs here. Anything a script
produces goes in that stage's `c_data/` folder.

| File | What one row is | Used by |
|---|---|---|
| `report.parquet` | one peptide measured in one MS run | `01_Preprocess` |
| `metadata.csv` | one MS sample, 131 rows | `01_Preprocess` |
| `phenotype.csv` | one leg, 70 rows | later stages |

## Get report.parquet first

It is too large for git, so it is attached to a release. From the repo root:

```sh
curl -L -o 00_Input/report.parquet \
  https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/download/data-v1/report.parquet
```

Or click `report.parquet` on the
[data-v1 release page](https://github.com/Dustyn-T-Lewis/BFR_limpa/releases/tag/data-v1).

`metadata.csv` says which sample is which: participant, leg, treatment, timepoint, sex, and the
original MS run name.

`phenotype.csv` holds the muscle measurements taken alongside the biopsies: cross-sectional
area and echo intensity before and after, sessions completed, and total training volume. No
stage reads it yet. It is the only route from the proteomics to a physical outcome, so it will
matter when we ask whether a protein change tracks a muscle change.

## Two things that are easy to get wrong

**Treatment is not the leg letter.** Run names contain `D` and `E` (right and
left). Which leg got restriction was counterbalanced across participants, so the letter tells you
the side only. `metadata.csv` carries `treatment` as its own column, simply done by hand.

**The same sample has three names.** The MS run name, our analysis ID, and the participant-plus-leg
pair. `metadata.csv` carries all three, and `01_Preprocess` reconciles them in one step so no
later stage has to.

## Why phenotype.csv has one row per leg

Training volume is one number for the whole leg. If we store it once per MS sample, it gets
repeated, so adding up the column would report twice the real volume.
