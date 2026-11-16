# 02 · Quantification

Estimates protein abundance from the filtered precursors with `dpcQuant()`, filters proteins on
detection, normalises with cyclic loess and adds a blood index per sample (mean log2 HBB and HBA1,
centred).

Inputs: `01_Filtering/c_data/precursors_filtered.rds`, `haemoglobin_precursors.rds` and the
`sample_contamination` sheet of `01_filter.xlsx`.

`dpcQuant()` takes about 100 minutes, so a script writes checkpoints to `c_data/quant_runs/`; rerun
it when the filtered precursors change. `--sensitivity` writes the run with the fitted slope.

```sh
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_quantify.html`: the report, with diagnostic plots.
- `c_data/proteins.rds`: abundance, standard errors, annotations and sample table (with blood index
  and contamination load), for stages 02 to 08.
- `c_data/proteins_adjusted.rds`: the same with the blood effect removed (`removeBatchEffect()`,
  group and participant protected), for steps that take no design.
- `c_data/02_quantify.xlsx`:
  - `samples`: sample identifiers, study metadata, blood index and contamination load.
  - `dpc_parameters`: the diagnostic detection-curve fit.
  - `detection_summary`: proteins before and after the detection filter, and its threshold.
  - `protein_quality`: per protein, precursors, proportion observed, mean abundance, mean standard
    error and mean observed precursors per sample.

Results: the detection filter (32 samples) keeps 3,042 of 3,451 protein groups. The fitted
detection-curve slope is 0.478, below the 0.7 to 0.9 the limpa FAQ calls typical, so the primary
run presets 0.7 and the fitted slope is a sensitivity run. The two haemoglobins correlate at
r 0.90, and the blood index tracks contamination load (Spearman 0.89).
