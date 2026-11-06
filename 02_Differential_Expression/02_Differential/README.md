# 02 · Differential expression

Fits the protein estimates with `dpcDE()` and tests the five contrasts with
`eBayes(robust = TRUE)` and BH adjustment within each contrast. `dpcDE()` uses the standard errors
as precision weights and, as in limpa's default, adds no sample quality weights. The design
carries the haemoglobin blood index as a covariate; 951 of 3,042 proteins have a blood slope at
BH 0.05. When the fitted-slope quantification run exists, the notebook also refits on it.

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

Training changes 95 proteins with BFR (58 up, 37 down) and 218 with HLRT (106 up, 112 down) at
FDR 0.05. The 51 shared hits all move the same way in both legs: RBP1, CA14, NAMPT and the IDH3
subunits rise, MYLK2 and SMTNL1 fall. 167 proteins pass only with HLRT, among them COL14A1 (log2
+2.01), BGN and PC (+1.25 and +1.24) and a fall in MYH1 (−1.32); 44 pass only with BFR, 27 of
them in MitoCarta (CHCHD2, TACO1, TIMM44, SLC25A3). Most one-arm hits lean the same way in the
other arm (153 of 167 and 42 of 44). Without the blood covariate the training hits are 86 and
161, all 86 and 153 of the 161 shared with the primary; on the fitted-slope run they are 100 and
221 (84 and 200 shared).

The interaction gives 126 proteins at nominal p < 0.05, the exploratory leads for a difference
between the arms. In 72 the HLRT change is the higher, led by BGN (BFR −0.34, HLRT +1.25), PC
(−0.28, +1.24), RPL36 (−0.31, +0.67), MSN (−0.05, +0.46) and MUSTN1 (−0.32, +0.34). In the other
54 the BFR change is the higher, such as ARL8A (+0.34, −0.39) and FAHD2A (+0.43, −0.10). None
passes BH (smallest q 0.68), no between-leg contrast has a hit under any version, and 126 is below
the 152 expected by chance, so these are leads only.
