# 00 · Input

Study data. No script runs here. Anything a script produces belongs in that stage's
`c_data` folder instead.

| File | What one row is | Read by |
|---|---|---|
| `report.parquet` | one precursor measured in one MS run | `01_Preprocess/01_Filtering` |
| `metadata.csv` | one MS sample, 131 rows | `01_Preprocess/01_Filtering` |
| `phenotype.csv` | one leg, 70 rows | nothing yet |

`report.parquet` is 355 MB. That is too large for git, so `.gitignore` lists it. Clone this
repo and you will not get the file. Copy it into this folder by hand before you run
`01_Filtering`. `readDIANN()` reads the DIA-NN v2 parquet format directly, and needs the
`nanoparquet` package installed to do it.

## The same sample has three names

`01_Filtering` reconciles them in one `match()` against `metadata$run`. No other stage does
this.

| Name | Example | Where it appears |
|---|---|---|
| DIA-NN run | `CL_MDR_2E-T1` | column names in `report.parquet`, and `metadata$run` |
| Sample ID | `S02_L_T1` | `metadata$sample_id`, and every matrix column downstream |
| Participant and leg | `S02` and `L` | `metadata`, and the join into `phenotype.csv` |

The match ignores case, because the instrument wrote `Cl_MDR_7E-T2` where the sheet has
`CL_MDR_7E-T2`.

`D` and `E` in a run name are Portuguese for right and left, *direita* and *esquerda*. They
say which leg, not which treatment. Leg assignment was counterbalanced across participants,
so you cannot read treatment off the letter. `metadata.csv` carries `treatment`, `timepoint`
and `group` as their own columns, written and checked by hand.

## Why phenotype.csv has one row per leg

Two reasons. It holds 140 ultrasound measurements, and a file with one row per MS sample
holds only 131 of them. And `volume_load_total_kg` counts the whole leg, so one row per
sample would repeat that number, and `sum()` would report twice the training volume.

## Deleted on 2026-09-07

`git log` still holds all three.

- `precursors.csv`. 53,852 rows by 141 columns, 45 MB. The same numbers `report.parquet`
  holds, written wide instead of long. No script read it.
- `BFR_meta.csv`. 130 rows in the old `2N_T1_D` naming, with leg measurements repeated onto
  sample rows. `metadata.csv` replaced it. The leg mapping survives there: the run name
  carries D or E, and the `leg` column carries R or L.
- `BFR_phenotype_data.xlsx`. The two-sheet workbook behind `phenotype.csv`. Its
  `mCSA_ECHO_PA_FL` sheet also held pennation angle and fascicle length, which
  `phenotype.csv` does not carry. Those two measurements now exist only in `git log` and at
  the collaborating site. Ask them for the workbook if anyone wants a fascicle analysis.

## Unanswered

Nobody recorded which DIA-NN version wrote `report.parquet`, which FASTA it searched, or
what FDR setting it used. `Unblinding.xlsx` was opened in this folder and is now gone.
`metadata.csv` carries the leg-to-treatment map, so the analysis runs without any of that.
Ask the collaborating site for the search settings before this goes to a journal. Their
answer will take longer to arrive than the remaining code work.
