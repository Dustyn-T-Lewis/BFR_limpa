# 02 · Differential expression

Fits the protein estimates with `dpcDE()`, which uses the standard errors as precision weights, and
tests the five contrasts with `eBayes(robust = TRUE)` and BH within each contrast. The design
carries the blood index. When the fitted-slope quantification run exists, the notebook refits on it.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`, `01_Design/c_data/design.rds`.

```sh
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/02_differential.html`: the report, with diagnostic plots.
- `c_data/fit.rds`: the fitted contrasts, for `03_Pathway_Enrichment/01_run_fgsea_and_fry`,
  `04_Network/02_module_traits`, `05_Mitochondria/02_test_mito` and `08_Plasticity`.
- `c_data/02_differential.xlsx`:
  - `protein_results`: per protein, log2 fold change, p and BH p per contrast, and the blood slope
    (log2 protein per log2 haemoglobin) with its BH p.
  - `contrast_summary`: discoveries per contrast at FDR 0.05 and 0.10, and at 0.05 without blood.

Results: training changes 95 proteins with BFR (58 up, 37 down) and 218 with HLRT (106 up, 112
down) at FDR 0.05. The 51 shared hits all move the same way: RBP1, CA14, NAMPT and the IDH3
subunits rise, MYLK2 and SMTNL1 fall. 167 pass only with HLRT, among them COL14A1 (log2 +2.01),
BGN and PC (+1.25, +1.24) and a fall in MYH1 (-1.32); 44 pass only with BFR, 27 of them in
MitoCarta (CHCHD2, TACO1, TIMM44, SLC25A3). Most one-arm hits lean the same way in the other arm
(153 of 167, 42 of 44). No between-leg contrast or interaction has a BH hit (smallest interaction q
0.68). 951 of 3,042 proteins have a blood slope at BH 0.05. Without the blood covariate the
training hits are 86 and 161 (86 and 153 shared); on the fitted-slope run, 100 and 221 (84 and
200 shared).

The interaction's 126 nominal hits (p < 0.05) are below the 152 expected by chance, so they are
leads only: BGN, PC, RPL36, MSN and MUSTN1 rise more with HLRT; ARL8A and FAHD2A more with BFR.
