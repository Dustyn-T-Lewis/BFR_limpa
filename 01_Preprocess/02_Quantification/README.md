# 01_Preprocess / 02_Quantification

This step fits the detection curve, rolls precursors up to proteins, filters proteins on detection
and normalises them. It starts from `01_Preprocess/01_Filtering/c_data/precursors_filtered.rds`
and writes `c_data/proteins.rds`, the protein table every later stage reads.
`c_data/02_quantify.xlsx` holds the fitted curve's intercept and slope and each protein's
coverage, abundance and standard error. The rendered notebook is `b_reports/02_quantify.html`.

```sh
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

The roll-up, `dpcQuant()`, takes about 100 minutes per call. It runs outside the notebook, in
`a_script/02_quantify_run.R`, and the notebook only loads what it saved, so a render never starts
a long run. Run the script again whenever the precursor matrix changes:

```sh
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity
```

The first command writes `c_data/quant_runs/proteins_slope_0.7.rds`, the primary run at the preset
slope. The second writes only `c_data/quant_runs/proteins_slope_fitted.rds`, the sensitivity run
at the fitted slope, which the notebook reads when it is present. The two can run at the same
time. Both files are committed, so a fresh clone renders without running either. Nothing ties a
checkpoint to the precursors it came from, so rerunning `01_Filtering` means rerunning this too.
