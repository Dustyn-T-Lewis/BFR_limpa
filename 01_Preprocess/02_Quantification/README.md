# 01_Preprocess / 02_Quantification

This step fits the detection curve, rolls precursors up to proteins, filters on detection and
normalises. It reads `01_Preprocess/01_Filtering/c_data/precursors_filtered.rds` and writes
`c_data/proteins.rds`, which every later stage reads, plus `c_data/02_quantify.xlsx` (curve
parameters and per-protein coverage, abundance and standard error) and `b_reports/02_quantify.html`.

```sh
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

`dpcQuant()` takes about 100 minutes per call, so `a_script/02_quantify_run.R` runs it outside the
notebook, which only loads the result. Rerun the script when the precursor matrix changes:

```sh
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity
```

The first command writes the primary run at the preset slope,
`c_data/quant_runs/proteins_slope_0.7.rds`; the second writes only the fitted-slope sensitivity run,
`c_data/quant_runs/proteins_slope_fitted.rds`, which the notebook reads when present. They can run
in parallel. Both are committed, so a fresh clone renders without either. Nothing ties a checkpoint
to its precursors, so rerunning `01_Filtering` means rerunning this too.

Protein groups are not capped. Titin's holds 3,245 precursors, where the median group holds 7.
