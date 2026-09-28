# 01_Preprocess / 01_Filtering

This step reads `00_Input/report.parquet` and `00_Input/metadata.csv` and removes precursor rows,
dropping no sample and filtering nothing on missing values.

`02_Quantification` reads the result, `c_data/precursors_filtered.rds`. `c_data/01_filter.xlsx`
records the rows, groups and signal left after each filter and every contaminant removed;
`b_reports/01_filter.html` gives the reason for each filter beside its code.

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports
```
