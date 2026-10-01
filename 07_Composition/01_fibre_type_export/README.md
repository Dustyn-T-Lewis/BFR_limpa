# 01_fibre_type_export

This step writes the file that FibeRtypeR (Van der Stede et al. 2026, J Physiol 604:2624) reads to
estimate the slow and fast fibre fractions of each biopsy. It reads `Run`, `Genes` and
`Genes.MaxLFQ.Unique` from `00_Input/report.parquet` and the run names in `00_Input/metadata.csv`,
and writes `c_data/fibertyper_input.csv`: one row per gene, first column `identifier`, one column
per `sample_id`. It is the only CSV in the repo, because the app accepts nothing else.

```sh
Rscript 07_Composition/01_fibre_type_export/a_script/01_fibre_type_export.R
```

The app wants raw LFQ intensities and normalises them itself. `Genes.MaxLFQ.Unique` is DIA-NN's
gene-level LFQ built from unique peptides, the quantification the authors found most accurate.
DIA-NN writes 0 on precursor rows that do not count toward it, so a 0 is read as not quantified
and left blank. The file holds 4,317 genes; 34.5% of cells are blank and 1,176 genes are complete
in all 131 biopsies. Of the 145 proteins in the published signature (77 slow, 68 fast), 118 are
present.

## Upload

1. Open https://ccgg.ugent.be/shiny/fibertyper/ and choose proteomics.
2. Load `c_data/fibertyper_input.csv` and say the identifiers are gene names.
3. If the app asks for sex, the proteomics reference is sex-specific; `sex` is a column of
   `00_Input/metadata.csv`.
4. Run it and download the table as .xlsx to
   `07_Composition/02_fibre_type_import/c_data/fibertyper_output.xlsx`.
