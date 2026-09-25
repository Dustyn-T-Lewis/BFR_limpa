# 03_Phenotype

Correlates each protein's pre-to-post change with the leg's change in muscle size and quality.

- Reads: `01_Preprocess/02_Quantification/c_data/proteins.rds`, `00_Input/phenotype.csv`
- Writes: `c_data/03_phenotype.xlsx` (`overview`, `summary`, `correlations`),
  `b_reports/03_phenotype.html`
- Run: `quarto render 02_Differential_Expression/03_Phenotype/a_script/03_phenotype.qmd --output-dir ../b_reports`,
  5 seconds.

## No protein tracks the phenotype after correction

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
