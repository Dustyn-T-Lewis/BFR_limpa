# 02 · Quantification

Estimates protein abundance from the filtered precursors with limpa, filters proteins on detection
and normalises them with cyclic loess.

Inputs: `01_Filtering/c_data/precursors_filtered.rds`.

`dpcQuant()` takes about 100 minutes, so a script runs it and saves checkpoints in
`c_data/quant_runs/`; rerun it whenever the filtered precursors change (the checkpoints do not
record their input). `--sensitivity` writes the fitted-slope run instead of the primary 0.7 run.

```sh
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_quantify.html`: the walkthrough, with diagnostic plots.
- `c_data/proteins.rds`: abundance, standard errors and annotations, read by every later stage.
- `c_data/02_quantify.xlsx`: samples, curve parameters, detection counts, normalisation summary and
  per-protein quality.
