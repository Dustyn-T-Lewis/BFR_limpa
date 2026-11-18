# 02 · Fibre type

Estimates each biopsy's slow (type I) fibre fraction with the FibeRtypeR method (Van der Stede et
al. 2026, J Physiol 604:2624, doi:10.1113/JP290082), rebuilt in R, and tests whether training or
arm shifted it. Markers are the 77 slow and 68 fast proteome markers of the paper's Table S1; the
reference is sex-specific, the median abundance of the pre-training slow and fast fibre pools of
each sex in Moesgaard et al. 2026 (Nat Commun, doi:10.1038/s41467-026-77275-9, PRIDE PXD056664).
The two-fraction least-squares fit has a closed form. Each sample is fitted against its own sex's
reference on the 99 markers both references and the matrix hold.

Two choices differ from the paper: the input is the limpa matrix rather than DIA-NN LFQ, because it
measures more markers in every sample, and each sample is shifted so its marker median equals the
reference's, because the scales differ. The deposit holds 48 pre-training pools where the paper
reports 46.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`, `references/fibertyper_reference.csv`.

The reference CSV was built once, outside the pipeline, from the deposit's DIA-NN 1.8 `pg_matrix`
(zeros as missing, median-normalised per run): per marker, the median log2 abundance over the
pre-training pools of each sex and fibre type. Columns: `gene`, `marker_type`, `table_s1_column`,
`slow_F`, `fast_F`, `slow_M`, `fast_M`; empty where no pool measured the marker.

```sh
quarto render 06_Composition/02_fibre_type/a_script/02_fibre_type.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_fibre_type.html`: the report.
- `b_reports/02_fibre_type.pdf`: each leg's slow fraction before and after training by arm.
- `c_data/02_fibre_type.xlsx`:
  - `slow_fraction`: per sample, participant, leg, treatment, timepoint, sex and slow fraction.
  - `contrasts`: the five stage 02 contrasts in percentage points, with SE, df, p and 95% interval.
  - `group_summary`: mean, SD and n of the slow fraction per arm and timepoint.

Results: training did not shift the slow fraction: BFR +1.0 points (95% interval -1.7 to 3.6,
p 0.47), HLRT +1.3 (-1.4 to 3.9, p 0.34), interaction -0.3 (-4.0 to 3.4, p 0.87). The median
biopsy is 42% slow (range 22% to 70%), and the estimate behaves like fibre type: legs agree at T1
(r 0.69), one leg's two biopsies agree (r 0.77), participants hold 77% of the variance and women
are 4.2 points slower. No histology exists for these biopsies, so the level is unvalidated and
only the ranking is used.
