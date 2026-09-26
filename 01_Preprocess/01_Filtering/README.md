# 01_Preprocess / 01_Filtering

Reads the DIA-NN report and removes precursor rows. It drops no sample and filters nothing on
missing values.

- Reads: `00_Input/report.parquet`, `00_Input/metadata.csv`
- Writes: `c_data/precursors_filtered.rds` (read by `02_Quantification`), `c_data/01_filter.xlsx`
  (`overview`, `filter_log`, `contaminants_removed`), `b_reports/01_filter.html`
- Run: `quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports`,
  12 seconds.

The notebook gives the reason for each filter beside its code.
