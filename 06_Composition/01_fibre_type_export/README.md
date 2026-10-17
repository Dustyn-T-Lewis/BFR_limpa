# 01 · Fibre type export

Writes the file FibeRtypeR (Van der Stede et al. 2026, J Physiol 604:2624) reads to estimate the
slow and fast fibre fractions of each biopsy: DIA-NN's MaxLFQ intensities from gene-specific
peptides, as the app asks; DIA-NN normalised them across runs and the app normalises again.

Inputs: `Run`, `Genes` and `Genes.MaxLFQ.Unique` from `00_Input/report.parquet`,
`00_Input/metadata.csv`.

```sh
quarto render 06_Composition/01_fibre_type_export/a_script/01_fibre_type_export.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_fibre_type_export.html`: the walkthrough.
- `c_data/fibertyper_input.csv`: one row per gene, first column `identifier`, one column per
  `sample_id`; unquantified cells blank. 4,317 genes, 34.5% of cells blank, 1,176 genes complete in
  all 131 biopsies. The only CSV output, because the app accepts nothing else.

The search FASTA lists 43 human proteins twice, because the appended cRAP database repeats them
(MB, CKM, FABP3 and CYCS among them). DIA-NN writes zero MaxLFQ for them in every run, so they are
absent from the upload. `01_Preprocess/01_Filtering` repairs them at precursor level for the main
analysis; this file comes straight from DIA-NN's gene values.

## Upload

1. Open https://ccgg.ugent.be/shiny/fibertyper/ and choose proteomics.
2. Load `c_data/fibertyper_input.csv` and say the identifiers are gene names.
3. If the app asks for sex, the proteomics reference is sex-specific; `sex` is a column of
   `00_Input/metadata.csv`.
4. Run it and download the table as .xlsx to
   `06_Composition/02_fibre_type_import/c_data/fibertyper_output.xlsx`.
