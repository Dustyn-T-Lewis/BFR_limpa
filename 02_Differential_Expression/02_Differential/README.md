# 02 · Differential expression

Fits the protein estimates with `dpcDE()` and tests the five contrasts with
`eBayes(robust = TRUE)` and BH adjustment within each contrast. `dpcDE()` uses the standard errors
as precision weights and, as in limpa's default, adds no sample quality weights. The design
carries the haemoglobin blood index as a covariate; 951 of 3,042 proteins have a blood slope at
BH 0.05. With it, BFR and HLRT training give 95 and 218 hits; without it 86 and 161, all 86 and
153 of the 161 shared. When the fitted-slope quantification run exists, the notebook also refits
on it: 100 and 221 training hits (84 and 200 shared with the primary). No contrast between legs
has a hit under any version.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`, `01_Design/c_data/design.rds`.

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_differential.html`: the walkthrough, with diagnostic plots.
- `c_data/fit.rds`: the fitted contrasts, read by `03_Pathway_Enrichment/01_run_fgsea_and_fry`,
  `04_Network/02_module_traits` and `05_Mitochondria/02_test_mito`.
- `c_data/02_differential.xlsx`:
  - `DEP_matrix`: one row per protein, with log2 fold change, p and BH-adjusted p per contrast,
    and the blood slope (log2 protein per log2 haemoglobin) with its BH-adjusted p;
  - `contrast_summary`: discoveries per contrast at FDR 0.05 and 0.10, and at FDR 0.05 without
    the blood covariate;
  - `null_summary`: per contrast, estimated share of true nulls, median |t| and share of p below
    0.2.
