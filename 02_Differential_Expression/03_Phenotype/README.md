# 03_Phenotype

Correlates every protein's pre-to-post change with how much the muscle actually changed.

| | |
|---|---|
| **Notebook** | `a_script/03_phenotype.qmd` |
| **Reads** | `proteins.rds`, `phenotype.csv` |
| **Writes** | `b_reports/03_phenotype.html`, `c_data/phenotype.rds`, `03_phenotype.xlsx` |

## Why it is not part of 02_Differential

A contrast compares group means. It cannot ask whether a protein tracks a per-leg outcome, because
the outcome never enters the design. This is the question the ultrasound data exists to answer and
nothing else in the pipeline asks it.

Two framings. **Pooled** uses all 65 legs and correlates each protein's change with that leg's
phenotype change. **Differential** uses the 32 participants with both legs and asks whether the
BFR-minus-high-load difference in a protein tracks the same difference in phenotype, which removes
every between-participant confound.

## Nothing survives, and the null is well behaved

3,042 proteins return 152 nominal hits at p < 0.05 when nothing is there, so nominal counts are
read against that rather than reported.

| Analysis | Outcome | n | Nominal | Ratio to chance | BH < 0.05 |
|---|---|---:|---:|---:|---:|
| pooled | vl_csa | 65 | 183 | 1.20 | 0 |
| pooled | vl_echo | 65 | 127 | 0.83 | 0 |
| pooled | rf_csa | 65 | 134 | 0.88 | 0 |
| pooled | rf_echo | 65 | 113 | 0.74 | 0 |
| differential | vl_csa | 32 | 146 | 0.96 | 0 |
| differential | vl_echo | 32 | 130 | 0.85 | 0 |
| differential | rf_csa | 32 | 168 | 1.10 | 0 |
| differential | rf_echo | 32 | 116 | 0.76 | 0 |

Every ratio sits between 0.74 and 1.20 and the best BH value across all 24,336 correlations is
0.281. This is a null, cleanly. With 32 pairs the study resolves a paired correlation of about
0.5, so it is a null at this resolution and not evidence that no protein tracks hypertrophy.

Spearman throughout, since neither a protein change nor an echo intensity change is expected to be
normal. Classification is deliberately absent: `02_Differential` already asks which proteins
changed, with the fixed participant term and 95 residual df, which beats a paired Wilcoxon on 32
pairs.
