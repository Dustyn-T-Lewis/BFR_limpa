# Run by hand when the precursor matrix changes; about 100 minutes per dpcQuant() call.
# 02_quantify.qmd only loads what this writes, so a render can never start one.
#   Rscript 01_Preprocess/02_Quantification/a_script/dpc_quant.R
#   Rscript 01_Preprocess/02_Quantification/a_script/dpc_quant.R --sensitivity

library(here)
library(limpa)

precursors_file <- here("01_Preprocess", "01_Filtering", "c_data",
                        "precursors_filtered.rds")
out <- here("01_Preprocess", "02_Quantification", "c_data", "dpcQuant")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

y <- readRDS(precursors_file)
# Recorded, not enforced. The notebook loads whatever is here; this says which precursor
# matrix it came from, for anyone reconstructing a run after the fact.
stamp <- unname(tools::md5sum(precursors_file))

write_run <- function(object, file) {
  saveRDS(list(object = object, input_md5 = stamp), file)
  message("wrote ", file)
}

write_run(
  dpcQuant(y, "Protein.Group", dpc.slope = 0.7),
  file.path(out, "proteins_slope_0.7.rds")
)

if ("--sensitivity" %in% commandArgs(trailingOnly = TRUE)) {
  write_run(
    dpcQuant(y, "Protein.Group", dpc = dpc(y$E)),
    file.path(out, "proteins_slope_fitted.rds")
  )
}
