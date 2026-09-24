# Run by hand when the precursor matrix changes; about 100 minutes per dpcQuant() call.
# 02_quantify.qmd only loads what this writes, so a render never starts one.
#   Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
#   Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity

pacman::p_load(here, limpa)

y <- readRDS(here("01_Preprocess", "01_Filtering", "c_data", "precursors_filtered.rds"))
out <- here("01_Preprocess", "02_Quantification", "c_data", "quant_runs")

if ("--sensitivity" %in% commandArgs(trailingOnly = TRUE)) {
  proteins <- dpcQuant(y, "Protein.Group", dpc = dpc(y$E))
  file <- file.path(out, "proteins_slope_fitted.rds")
} else {
  proteins <- dpcQuant(y, "Protein.Group", dpc.slope = 0.7)
  file <- file.path(out, "proteins_slope_0.7.rds")
}
saveRDS(proteins, file)
message("wrote ", file)
