# --- 02_Differential_Expression / 02_Differential ----------
# limpa's sequence for a contrast-based analysis: dpcDE -> contrasts.fit -> eBayes -> topTable.
# dpcDE() is the ONLY entry point that reads the standard errors dpcQuant() produced, so it
# replaces lmFit()

library(here)
library(limpa)
library(limma)
library(dplyr)
library(purrr)
library(readr)
library(tibble)

de_dir <- here("02_Differential_Expression", "02_Differential", "c_data")
dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)

# 1: Input, and the assertion
proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))
d        <- readRDS(here("02_Differential_Expression", "01_Design", "c_data", "design.rds"))

# limma compares only the COUNT of design rows to matrix columns, never their
# identity or order. Are the row names of the design equal to the column names of the expression matrix?
stopifnot(identical(rownames(d$design), colnames(proteins$E)))
cat("leg correlation from design stage:", round(d$leg_correlation, 3), "\n")
dim(proteins)

# --- Fit ---------------------------
# sample.weights = TRUE estimates one quality weight per biopsy. Samples differed in how much
# blood they carried; removing the contaminant ROWS afterwards does not undo what heavy
# contamination did to the sample underneath, because that scan time is gone.
fit <- dpcDE(proteins, d$design, sample.weights = TRUE, plot = TRUE)
# Extract the sample weights BEFORE contrasts.fit(), which does not carry $targets through.
stopifnot(!is.null(fit$targets$sample.weight))
sample_weights <- tibble(
  sample_id = colnames(proteins),
  group     = proteins$targets$group,
  weight    = fit$targets$sample.weight
)

fit <- contrasts.fit(fit, d$contrasts)
fit <- eBayes(fit, robust = TRUE)

summary(fit$df.prior)      # a VECTOR under robust = TRUE, one value per protein
fit$df.total[1:5]

# Did the weighting stay even across the four cells?
# Weights downweight noisy samples; none is dropped. What matters is whether the
# downweighting falls EVENLY across the four groups. 
weight_balance <- sample_weights |>
  summarise(n = n(),
            median_weight = round(median(weight), 2),
            lowest_weight = round(min(weight), 2),
            .by = group)
weight_balance
arrange(sample_weights, weight) |> head(5)            # the five worst biopsies


# What effect could the primary contrast have detected?
primary <- d$primary  # "interaction"

se <- fit$stdev.unscaled[, primary] * sqrt(fit$s2.post)

# (a) Nominal: two-sided alpha = 0.05, 80% power. A LOWER bound - ignores multiplicity.
mde_nominal <- se * (qt(0.975, fit$df.total) + qt(0.80, fit$df.total))
print(mde_nominal)
# (b) Multiplicity-adjusted: the rank-1 BH bar, alpha/m.
# Accounting for multiplicity roughly doubles the effect you need.
alpha_bh <- 0.05 / nrow(fit)
mde_bh   <- se * (qt(1 - alpha_bh / 2, fit$df.total) + qt(0.80, fit$df.total))

detectable <- tibble(
  contrast   = primary,
  quantile   = c("median", "90th pct"),
  mde_log2_nominal = round(quantile(mde_nominal, c(0.5, 0.9)), 3),
  mde_log2_bh      = round(quantile(mde_bh,      c(0.5, 0.9)), 3)
) |>
  mutate(fold_nominal = round(2^mde_log2_nominal, 2),
         fold_bh      = round(2^mde_log2_bh, 2))
detectable

# --- Results ----------------
results <- map(colnames(d$contrasts), \(cn) {
  topTable(fit, coef = cn, number = Inf, adjust.method = "BH", sort.by = "p") |>
    rownames_to_column("protein") |>
    mutate(contrast = cn, .before = 1)
}) |>
  list_rbind() |>
  mutate(pi_score = P.Value^abs(logFC))

contrast_summary <- results |>
  summarise(n_tested  = n(),
            n_fdr_05  = sum(adj.P.Val < 0.05),
            n_fdr_10  = sum(adj.P.Val < 0.10),
            min_adj_p = signif(min(adj.P.Val), 3),
            median_abs_logFC = round(median(abs(logFC)), 4),
            .by = contrast)
contrast_summary

# Is a contrast empty, or is it being tested conservatively?
# Three numbers together, because each misleads alone.
#   prop_true_null  estimates the fraction of tests that are truly null
#   median_abs_t    is 0.674 under a calibrated null; below that means conservative
#   p_below_0.2     is 0.20 under a uniform null; below that means the distribution
#                   rises toward 1 rather than being flat
#
# Read every contrast against BFR_vs_HLRT_at_T1, NOT against theory. That contrast is a true
# null by construction, so whatever it returns is what a null looks like in THIS design.
null_calibration <- results |>
  summarise(prop_true_null = round(propTrueNull(P.Value), 3),
            median_abs_t   = round(median(abs(t)), 3),
            p_below_0.2    = round(mean(P.Value < 0.2), 3),
            .by = contrast)
null_calibration

treat_summary <- map(c(1.1, 1.15, 1.2), \(fc) {
  ft <- treat(fit, lfc = log2(fc), robust = TRUE)
  tibble(fold_change = fc,
         lfc_floor   = round(log2(fc), 3),
         contrast    = colnames(fit$coefficients),
         n_q_05      = map_int(colnames(fit$coefficients),
                               \(cn) sum(topTreat(ft, coef = cn, number = Inf)$adj.P.Val < 0.05)))
}) |>
  list_rbind()
treat_summary


# =========================================================================================
# 8 · The pi score - ranking only, never selection
# =========================================================================================
# pi_score = P.Value^|logFC|, which equals 10^(-pi) for the score of Xiao et al. (2014),
# so it ranks exactly as the published -log10(p) * |log2FC| does. Smaller = higher ranked.
# It controls NO error rate, so it orders within a contrast and selects nothing.
#
# It rewards large |logFC|, and here the largest fold changes come from proteins with the
# fewest detected precursors, whose values lean on DPC recovery. The observation count is
# printed alongside so that is visible rather than hidden.
#
# Applied only to the two contrasts that carry signal. Ranking the interaction would produce
# a table of the twenty largest fold changes in a contrast whose prop_true_null sits at its
# ceiling, which invites reading structure into noise whatever the caption says.
mean_obs <- rowMeans(proteins$other$n.observations)

pi_ranking <- results |>
  filter(contrast %in% c("BFR_post_vs_pre", "HLRT_post_vs_pre")) |>
  slice_min(pi_score, n = 10, by = contrast) |>
  mutate(mean_obs = round(mean_obs[protein], 2)) |>
  select(contrast, protein, Genes, logFC, adj.P.Val, pi_score, NPrec, mean_obs)
pi_ranking


# Write - BEFORE the control is checked, so a failure leaves the evidence
write_csv(results,          file.path(de_dir, "protein_contrasts_long.csv"))
write_csv(contrast_summary, file.path(de_dir, "contrast_summary.csv"))
write_csv(null_calibration, file.path(de_dir, "null_calibration.csv"))
write_csv(treat_summary,    file.path(de_dir, "treat_summary.csv"))
write_csv(detectable,       file.path(de_dir, "detectable_effect.csv"))
write_csv(sample_weights,   file.path(de_dir, "sample_weights.csv"))
write_csv(weight_balance,   file.path(de_dir, "weight_balance.csv"))
saveRDS(fit,                file.path(de_dir, "fit.rds"))

# Stale output from the previous contrast naming would be silently read by downstream
# scripts. Name them rather than delete them.
stale <- list.files(de_dir, pattern = "^DE_.*\\.csv$|^de_summary\\.csv$", full.names = FALSE)
if (length(stale)) {
  warning("Stale files from the old contrast names are still present: ",
          paste(stale, collapse = ", "), ". Delete them before running 05_Figures.")
}


# =========================================================================================
# 10 · Negative control
# =========================================================================================
# Two randomly assigned, untrained legs of the same person on the same day should not differ.
# If they do, the cause is upstream - annotation, contamination, a sample swap - and that
# failure is the most informative result this pipeline can produce.
nc <- filter(results, contrast == d$negative_control)

# Asserted because sum(logical(0)) is 0: a contrast name matching nothing would look like a
# clean control while testing nothing at all. This is exactly how a rename breaks a check.
stopifnot(nrow(nc) == nrow(proteins))

n_control <- sum(nc$adj.P.Val < 0.05)
cat("\n=== NEGATIVE CONTROL (", d$negative_control, ") ===\n", sep = "")
cat("proteins at FDR < 0.05:", n_control,
    if (n_control == 0) " -> PASS\n" else " -> FAIL, investigate upstream before using any result\n")
if (n_control > 0) warning("Negative control is not empty. Do not relax the threshold.")

slice_head(nc, n = 5) |> select(protein, Genes, logFC, P.Value, adj.P.Val)

proteins <- readRDS(here("01_Preprocess","02_Quantification","c_data","proteins.rds"))
sd(apply(proteins$E, 2, median))

