# 00 · Input

The study data. Nothing here runs, and nothing computed belongs here. Anything a script
produces goes in that stage's `c_data/` folder.

| File | What one row is | Used by |
|---|---|---|
| `report.parquet` | one peptide measured in one MS run | `01_Preprocess` |
| `metadata.csv` | one MS sample, 131 rows | `01_Preprocess` |
| `phenotype.csv` | one leg, 70 rows | later stages |

`report.parquet` is the DIA-NN search output, 355 MB, too large for git. It is not in the repo.
Ask for a copy and put it here before running anything.

`metadata.csv` says which sample is which: participant, leg, treatment, timepoint, sex, and the
original MS run name.

`phenotype.csv` holds the muscle measurements taken alongside the biopsies: cross-sectional
area and echo intensity before and after, sessions completed, and total training volume. No
stage reads it yet. It is the only route from the proteomics to a physical outcome, so it will
matter when we ask whether a protein change tracks a muscle change.

## Two things that are easy to get wrong

**Treatment is not the leg letter.** Run names contain `D` and `E`, Portuguese for right and
left. Which leg got restriction was counterbalanced across participants, so the letter tells you
the side and nothing else. `metadata.csv` carries `treatment` as its own column, worked out once
by hand.

**The same sample has three names.** The MS run name, our analysis ID, and the participant-plus-leg
pair. `metadata.csv` carries all three, and `01_Preprocess` reconciles them in one step so no
later stage has to.

## Why phenotype.csv has one row per leg

Training volume is one number for the whole leg. Store it once per MS sample instead and it gets
repeated, so adding up the column would report twice the real volume.
