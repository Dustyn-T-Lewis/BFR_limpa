# 00 · Input

Proteomics measurements, sample information and muscle measurements from the training study.

| File | Contents | One row | Read by |
|---|---|---|---|
| `report.parquet` | DIA-NN precursor measurements and protein annotations | one precursor in one MS run | `01_Filtering`, `06_Composition/01_fibre_type_export` |
| `metadata.csv` | sample identifiers, participant, leg, treatment, timepoint, group and sex | one MS sample | `01_Filtering`, `06_Composition/01_fibre_type_export` |
| `phenotype.csv` | before-and-after muscle size and echo intensity, session count and training volume | one participant’s leg | `03_Protein_Association`, `03_Pathway_Enrichment/05_set_association` |

Treatment cannot be inferred from leg side: blood-flow restriction was counterbalanced between left
and right legs. Use the `treatment` column in `metadata.csv`.

`metadata.csv` links MS run names to sample IDs and identifies each sample’s participant, leg and
timepoint. In run names, `D` means right and `E` means left. Seven run names end in `r`; for one of
them, S08’s pre-training sample, the report also holds an earlier injection (`CL_MDR_8E-T1`) that
the metadata does not list, so `01_Filtering` drops it.

DIA-NN wrote the report already filtered at a precursor `Q.Value` and `Lib.Q.Value` of 0.01; only
the protein-group cutoff in `readDIANN()` removes further rows.

`phenotype.csv` stores both timepoints on the same row, with total training volume recorded per
leg. It covers 35 participants: S22 and S34 have no MS runs, and the association steps join
phenotype onto the proteomics legs, so they drop out.

`report.parquet` is stored in the `data-v1` release because it is too large for Git. See the
[project README](../README.md) for the download command.

## Correction (2026-09-19): vastus lateralis pre-CSA, 8 of 70 rows

The CV between the two pre-training VL CSA measurements was high, so the ultrasound images were
re-analysed. `vl_csa_pre_cm2` changed for 8 legs (S07, S10, S11, S18 ×2, S21, S24, S29). Nothing
else in the file changed.
