# The FibeRtypeR upload: DIA-NN's gene-level LFQ from unique peptides, one row per gene and one
# column per biopsy. The app wants raw LFQ and normalises it itself, and its authors found
# unique-peptide quantification the most accurate (Van der Stede et al. 2026, J Physiol 604:2624).

pacman::p_load(here, dplyr, tidyr, readr)

out <- here("07_Composition", "01_fibre_type_export", "c_data")
metadata <- read_csv(here("00_Input", "metadata.csv"), show_col_types = FALSE)
report <- nanoparquet::read_parquet(
  here("00_Input", "report.parquet"),
  col_select = c("Run", "Genes", "Genes.MaxLFQ.Unique")
)

# Runs match metadata$run case-insensitively, as in 01_filter.qmd, which absorbs Cl_MDR_7E-T2
# and drops the discarded S08 injection. The gene value repeats on every precursor row, and
# DIA-NN writes 0 on rows whose precursor is not unique to the gene, so 0 means not quantified.
lfq <- report |>
  filter(nzchar(Genes), Genes.MaxLFQ.Unique > 0) |>
  distinct() |>
  mutate(sample_id = metadata$sample_id[match(tolower(Run), tolower(metadata$run))]) |>
  filter(!is.na(sample_id))
stopifnot(!anyDuplicated(select(lfq, sample_id, Genes)))

upload <- lfq |>
  select(identifier = Genes, sample_id, Genes.MaxLFQ.Unique) |>
  pivot_wider(names_from = sample_id, values_from = Genes.MaxLFQ.Unique) |>
  select(identifier, all_of(metadata$sample_id))
stopifnot(!anyDuplicated(upload$identifier), ncol(upload) == nrow(metadata) + 1)

# The one CSV in the repo: the app reads nothing else. Missing values stay blank.
write_csv(upload, file.path(out, "fibertyper_input.csv"), na = "")
message(nrow(upload), " genes x ", ncol(upload) - 1, " samples")
sessionInfo()
