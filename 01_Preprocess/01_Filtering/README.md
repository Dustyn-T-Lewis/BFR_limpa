# 01 · Filtering

Turns the DIA-NN report into a filtered precursor matrix: repairs the search annotation, removes
shared and compound assignments, and removes blood and skin contaminants before quantification.
Nothing is filtered on missing values.

Inputs: `00_Input/report.parquet`, `00_Input/metadata.csv`.

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_filter.html`: the walkthrough, with previews at each step.
- `c_data/precursors_filtered.rds`: filtered precursors, read by `02_Quantification`.
- `c_data/01_filter.xlsx`:
  - `samples`: sample identifiers and study metadata;
  - `filtering_summary`: precursors, protein groups and percentage of starting signal after each
    filter;
  - `annotation_exclusions`: the overlapping cRAP and decoy categories, with their signal;
  - `removed_contaminants`: each removed gene-symbol string, its precursors and signal;
  - `sample_contamination`: per sample, the percentage of signal the contaminant filter removed.

The contaminant list was picked by judgement, not by a rule; the notebook says why a scripted HPA
rule failed to replace it.
