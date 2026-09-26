# 01_Preprocess / 02_Quantification

Fits the detection curve, rolls precursors up to proteins, filters on detection and normalises.

- Reads: `01_Preprocess/01_Filtering/c_data/precursors_filtered.rds`,
  `c_data/quant_runs/proteins_slope_0.7.rds`, and `c_data/quant_runs/proteins_slope_fitted.rds`
  when present
- Writes: `c_data/proteins.rds`, `c_data/02_quantify.xlsx` (`overview`, `dpc_parameters`,
  `protein_quality`), `b_reports/02_quantify.html`
- Run: `quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports`,
  9 seconds.

`dpcQuant()` costs about 100 minutes per call, so `a_script/02_quantify_run.R` computes it and
`02_quantify.qmd` only loads the result. A render cannot start a long run. Run the script when the
precursor matrix changes:

```sh
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity
```

The first writes `c_data/quant_runs/proteins_slope_0.7.rds`, the primary run at the preset slope.
The second writes only `proteins_slope_fitted.rds`, the fitted-slope sensitivity run. They can run
at the same time. Both are committed, so a fresh clone renders without running either. Nothing ties
a checkpoint to the precursors it came from, so rerunning `01_Filtering` means rerunning this too.
