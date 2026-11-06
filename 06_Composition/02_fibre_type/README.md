# 02 · Fibre type

Estimates each biopsy's slow (type I) fibre fraction with the FibeRtypeR method (Van der Stede et
al. 2026, J Physiol 604:2624, doi:10.1113/JP290082) and tests whether training or arm shifted it.
The method is rebuilt in R from the paper: the 77 slow and 68 fast proteome markers of its Table S1
(Shared and Proteome columns), reference profiles from the pre-training slow and fast fibre pools of
Moesgaard et al. 2026 (Nat Commun, doi:10.1038/s41467-026-77275-9, PRIDE PXD056664), and a
non-negative least-squares fit with the two fractions summing to one, which for two fractions has a
closed form. The input is the limpa matrix, which measures 100 of the 145 markers in every sample;
the paper used DIA-NN MaxLFQ.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`, and `references/fibertyper_reference.csv`
(per marker: type, Table S1 column, and mean log2 abundance in the 24 slow and 24 fast pools; built
from the PRIDE deposit, median-normalised per pool).

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

The median biopsy is 40% slow (range 18% to 67%), and the estimate behaves like fibre type: left
and right legs agree at T1 (r 0.71), one leg's two biopsies agree (r 0.78), participants hold 78%
of the variance and women are 5.2 points slower. Training did not shift it: BFR +1.0 points (95%
interval −1.7 to 3.8, p 0.45), HLRT +1.5 (−1.2 to 4.2, p 0.28), interaction −0.4 (−4.3 to 3.4,
p 0.83). It is not a covariate: participant already absorbs most of it, and in a scratch test
adding it to the stage 02 design changed hit counts by under 10% and fold changes not at all
(r ≥ 0.996). No histology exists for these biopsies, so the level of the fraction is unvalidated
and only its ranking is used.
