# 01 · Filtering

Repairs the search annotation, removes shared and compound assignments, and removes blood and skin
contaminants before quantification. Nothing is filtered on missing values.

Inputs: `00_Input/report.parquet`, `00_Input/metadata.csv`.

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_filter.html`: the report.
- `c_data/precursors_filtered.rds`: filtered precursors, for `02_Quantification`.
- `c_data/haemoglobin_precursors.rds`: HBB and HBA1 precursors, for the blood index.
- `c_data/01_filter.xlsx`:
  - `filtering_summary`: precursors, protein groups and % of starting signal after each filter.
  - `annotation_exclusions`: cRAP and decoy categories (overlapping), with their signal.
  - `removed_contaminants`: each removed gene-symbol string, its precursors and signal.
  - `sample_contamination`: per sample, % of signal the contaminant filter removed.

Results: 47,201 precursors in 3,451 protein groups remain, holding 59.7% of the starting signal.
The contaminant filter removes 27.9 points of it, mostly HBA1 (9.4%), HBB (9.3%) and ALB (8.0%);
per sample it removes a median 28% (range 7% to 79%). The contaminant list was picked by
judgement, not by a rule; the notebook says why a scripted HPA rule failed to replace it.
