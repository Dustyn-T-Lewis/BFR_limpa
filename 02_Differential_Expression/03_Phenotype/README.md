# 02_Differential_Expression / 03_Phenotype

Correlates each protein's pre-to-post change with the leg's change in muscle size and quality.

| | |
|---|---|
| Script | `a_script/03_phenotype.qmd` |
| Reads | `proteins.rds`, `00_Input/phenotype.csv` |
| Writes | `c_data/phenotype.rds`, `c_data/03_phenotype.xlsx` |

A contrast compares group means and cannot use a per-leg outcome, so this is a separate question.
Pooled: 65 legs, protein change against phenotype change. Differential: 32 participants,
the BFR-minus-HLRT difference in a protein against the same difference in phenotype, which removes
between-participant confounds. Spearman throughout; BH within each analysis and outcome.

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

Chance is 152 nominal hits from 3,042 proteins. Nothing survives BH (best 0.281). With 32 pairs
the resolvable paired correlation is about 0.5, so this is a null at that resolution.

Classification is left to `02_Differential`, whose fitted model beats a paired Wilcoxon on 32
pairs.
