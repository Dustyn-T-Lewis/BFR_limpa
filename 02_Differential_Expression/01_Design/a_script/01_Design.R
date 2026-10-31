# ---02_Differential_Expression / 01_Design ----------

# Builds the design matrix and the five contrasts, checks both, writes design.rds.
# Fits no model. Reads no protein abundance.

library(here)
library(limpa)
library(limma)
library(dplyr)
library(tibble)
library(readr)

design_dir <- here("02_Differential_Expression", "01_Design", "c_data")
dir.create(design_dir, recursive = TRUE, showWarnings = FALSE)

# --- Input ------------------------
proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))
targets <- proteins$targets 
dim(proteins) # 3042 proteins x 131 samples

# 1 What experiment did we actually run?
with(targets, table(treatment, timepoint))
with(targets, table(group))
# sample_id must be unique: it is the key used to align design rows to matrix columns.
stopifnot(!anyDuplicated(targets$sample_id))
# group must be exactly treatment_timepoint. Catches a hand-edited metadata cell.
stopifnot(identical(as.character(targets$group),
                    paste(targets$treatment, targets$timepoint, sep = "_")))
# Which participants are missing a cell, and is that the gap we expect?
cells <- count(targets, participant, name = "n_samples")
table(cells$n_samples) # S06: 3 samples, no BFR_T2 biopsy
filter(cells, n_samples != 4)     

# The negative control is only a control if leg assignment was randomized. If BFR were
# always the right leg, BFR_vs_HLRT_at_T1 would be measuring limb dominance, not error.
leg_balance <- targets |>
  distinct(participant, leg, treatment) |>
  filter(treatment == "BFR") |>
  count(leg, name = "n_participants")
leg_balance   # expect roughly half L, half R              

# --- 02 Design Matrix --------------------------
targets$group <- factor(targets$group, 
                        levels = c("BFR_T1", "BFR_T2", "HLRT_T1", "HLRT_T2"))
targets$participant <- factor(targets$participant)

design <- model.matrix(~0 + group + participant, data = targets)

colnames(design) <- sub("^group", "", colnames(design))   # groupBFR_T1 -> BFR_T1
rownames(design) <- targets$sample_id                     # required for the check in stage 02

c(samples      = nrow(design),
  coefficients = ncol(design),
  residual_df  = nrow(design) - ncol(design))             # 131, 36, 95

design[1:6, 1:6]

# 3 Three checks, each guarding a silent failure

# (a) Column names must be legal R names. makeContrasts() PARSES its arguments as code,
#     so a column called "2E-T1" would be read as a subtraction.
stopifnot(identical(make.names(colnames(design)), colnames(design)))

# (b) Full rank: every column carries information no combination of the others carries.
#     If rank < ncol, some coefficient is not estimable and limma may return NA rather
#     than an error. This is the check that would fire if you ever added `+ sex`.
stopifnot(qr(design)$rank == ncol(design))

# (c) Alignment: design rows are the same samples, in the same order, as matrix columns.
#     limma compares only their COUNTS. Reorder one and every result is silently wrong.
stopifnot(identical(rownames(design), colnames(proteins$E)))

# Looking the the intra block correlation between samples (PRE and POST) taken from the same leg
leg_id   <- paste(targets$participant, targets$leg, sep = "_")   # one level per leg
leg_corr <- duplicateCorrelation(proteins$E, design, block = leg_id)$consensus.correlation

round(leg_corr, 3) # 0.021 is close to 0 which suggests that samples from coming from the same leg are as different as 
# samples coming from two different legs within the same person.

# 5: The five contrasts
contrasts <- makeContrasts(
  BFR_post_vs_pre   = BFR_T2 - BFR_T1,                          # descriptive
  HLRT_post_vs_pre  = HLRT_T2 - HLRT_T1,                        # descriptive
  BFR_vs_HLRT_at_T1 = BFR_T1 - HLRT_T1,                         # NEGATIVE CONTROL
  BFR_vs_HLRT_at_T2 = BFR_T2 - HLRT_T2,                         # descriptive
  interaction       = (BFR_T2 - BFR_T1) - (HLRT_T2 - HLRT_T1),  # PRIMARY
  levels = design
)

contrasts[rowSums(abs(contrasts)) > 0, ]     # show only rows that carry a weight

# The five are not five independent questions. Assert the identity rather than trusting
# the comment: interaction IS the difference of the two between-leg contrasts.
stopifnot(all.equal(unname(contrasts[, "interaction"]),
                    unname(contrasts[, "BFR_vs_HLRT_at_T2"] -
                           contrasts[, "BFR_vs_HLRT_at_T1"])))

# 5: Writing the design object
design_obj <- list(
  design           = design,
  contrasts        = contrasts,
  primary          = "interaction",
  negative_control = "BFR_vs_HLRT_at_T1",
  leg_correlation  = leg_corr,
  built_on         = Sys.Date(),
  limma_version    = as.character(packageVersion("limma")),
  limpa_version    = as.character(packageVersion("limpa"))
)
saveRDS(design_obj, file.path(design_dir, "design.rds"))

write_csv(
  tibble(level = "within leg, participant fixed (unweighted matrix correlation)",
         correlation = round(leg_corr, 3),
         samples = nrow(design), coefficients = ncol(design),
         residual_df = nrow(design) - ncol(design)),
  file.path(design_dir, "design_summary.csv")
)

cat("design written:", nrow(design), "samples,", ncol(design), "coefficients,",
    nrow(design) - ncol(design), "residual df\n")
