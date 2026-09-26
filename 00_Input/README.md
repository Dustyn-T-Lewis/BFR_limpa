# 00 · Input

Study data.

| File | What one row is | Read by |
|---|---|---|
| `report.parquet` | one precursor in one MS run | `01_Preprocess/01_Filtering` |
| `metadata.csv` | one MS sample, 131 rows | `01_Preprocess/01_Filtering` |
| `phenotype.csv` | one leg, 70 rows | `02_Differential_Expression/03_Phenotype`, `03_Pathway_Enrichment/05_classify_and_associate_sets` |

`report.parquet` is too large for git and sits on the `data-v1` release. The download command is
in the root README.

## Notes

Treatment is not the leg letter. Run names contain `D` and `E` (right and left). Which leg
got restriction was counterbalanced, so the letter gives the side only. `metadata.csv` carries
`treatment` as its own column, done by hand.

The same sample has three names: the MS run name, our sample ID, and participant plus leg.
`metadata.csv` holds all three and `01_Filtering` reconciles them once.

`phenotype.csv` stays at one row per leg because training volume is one number for the whole
leg. Stored once per MS sample it would be repeated, and summing the column would report twice
the real volume.

### Correction (2026-09-19): vastus lateralis pre-CSA, 8 of 70 rows

The CV between the two pre-training VL CSA measurements was high, so the ultrasound images were
re-analysed. `vl_csa_pre_cm2` changed for 8 legs (S07, S10, S11, S18 ×2, S21, S24, S29). Nothing
else in the file changed.

