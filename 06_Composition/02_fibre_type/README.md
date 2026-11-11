# 02 · Fibre type

Estimates each biopsy's slow (type I) fibre fraction with the FibeRtypeR method (Van der Stede et
al. 2026, J Physiol 604:2624, doi:10.1113/JP290082) and tests whether training or arm shifted it.
The method is rebuilt in R from the paper: the 77 slow and 68 fast proteome markers of its Table S1
(Shared and Proteome columns), sex-specific reference profiles (the median abundance of the
pre-training slow and fast fibre pools of each sex in Moesgaard et al. 2026, Nat Commun,
doi:10.1038/s41467-026-77275-9, PRIDE PXD056664), and a non-negative least-squares fit with the two
fractions summing to one, which for two fractions has a closed form. Each sample is fitted against
its own sex's reference, on the 99 markers that both references and the matrix hold.

Two choices differ from the paper. The input is the limpa matrix rather than DIA-NN LFQ, because
limpa's complete matrix measures more markers in every sample. Each sample is shifted so its median
over the markers equals the reference's, because limpa's scale differs from the deposit's; the paper
says only "log2-normalized". The deposit also holds 48 pre-training pools (12 per sex and fibre
type) where the paper reports 46.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`, and `references/fibertyper_reference.csv`.

The reference CSV was built once, outside the pipeline, from the deposit's sample sheet and log2
protein matrix (DIA-NN 1.8 `pg_matrix`, zeros as missing, median-normalised per run): for each of
the 145 markers, the median log2 abundance over the pre-training pools of each sex and fibre type,
missing values dropped. Columns: `gene`, `marker_type`, `table_s1_column`, `slow_F`, `fast_F`,
`slow_M`, `fast_M`. A marker no pool of one sex and type measured is left empty.

```sh
quarto render 06_Composition/02_fibre_type/a_script/02_fibre_type.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_fibre_type.html`: the walkthrough.
- `b_reports/02_fibre_type.pdf`: one page, each leg's slow fraction before and after training by
  arm, with the arm mean.
- `c_data/02_fibre_type.xlsx`:
  - `slow_fraction`: one row per sample: participant, leg, treatment, timepoint, sex and slow
    fraction;
  - `contrasts`: the five stage 02 contrasts on the slow fraction in percentage points: estimate,
    standard error, df, p and 95% interval;
  - `group_summary`: mean, SD and n of the slow fraction per arm and timepoint;
  - `checks`: agreement between legs and between timepoints, the share of variance between
    participants and the sex difference.

The median biopsy is 42% slow (range 22% to 70%), and the estimate behaves like fibre type: left
and right legs agree at T1 (r 0.69), one leg's two biopsies agree (r 0.77), participants hold 77%
of the variance and women are 4.2 points slower. Training did not shift it: BFR +1.0 points (95%
interval −1.7 to 3.6, p 0.47), HLRT +1.3 (−1.4 to 3.9, p 0.34), interaction −0.3 (−4.0 to 3.4,
p 0.87). No histology exists for these biopsies, so the level of the fraction is unvalidated and
only its ranking is used.
