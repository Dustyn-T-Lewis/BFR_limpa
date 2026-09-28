# 01_Preprocess / 01_Filtering

This step reads the DIA-NN report and the sample metadata from `00_Input/report.parquet` and
`00_Input/metadata.csv`, and removes precursor rows. It drops no sample and filters nothing on
missing values.

The filtered precursors are saved as `c_data/precursors_filtered.rds`, which `02_Quantification`
reads. `c_data/01_filter.xlsx` records the rows, groups and signal left after each filter, and
every contaminant removed with its signal. The rendered notebook, `b_reports/01_filter.html`,
gives the reason for each filter beside its code.

```sh
quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports
```
