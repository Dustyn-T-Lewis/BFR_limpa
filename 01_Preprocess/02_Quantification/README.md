# 02 · Quantification

Estimates protein abundance from the filtered precursors with limpa, filters proteins on detection
and normalises them with cyclic loess.

Inputs: `01_Filtering/c_data/precursors_filtered.rds`.

`dpcQuant()` takes about 100 minutes, so a script runs it and saves checkpoints in
`c_data/quant_runs/`; rerun it whenever the filtered precursors change. The primary run presets the
detection-curve slope at 0.7; `--sensitivity` writes the run with the slope fitted to these data.

```sh
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_quantify.html`: the walkthrough, with diagnostic plots.
- `c_data/proteins.rds`: abundance, standard errors and annotations, read by stages 02 to 06.
- `c_data/02_quantify.xlsx`:
  - `samples`: sample identifiers and study metadata;
  - `dpc_parameters`: the diagnostic detection-curve fit (the primary run uses slope 0.7);
  - `detection_summary`: proteins before and after the detection filter, and its threshold;
  - `protein_quality`: per protein, precursor count, proportion observed, mean abundance, mean
    standard error and mean observed precursors per sample.
