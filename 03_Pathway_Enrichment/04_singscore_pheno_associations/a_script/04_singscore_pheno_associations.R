# Ask whether a pathway's movement tracks the phenotype, and whether it tracks it differently
# under BFR than under high load.
#
# Every participant contributed one BFR leg and one HLRT leg, so the differential question is
# asked within participant: does the leg that gained more pathway score also gain more muscle?
# One number per participant removes every between-participant confound and, at n = 32, has
# roughly double the power of comparing two independent groups of legs.
#
# 2,732 sets x 4 outcomes is 10,928 tests per analysis. Nominal p is the screening statistic,
# BH is reported beside it, and signal_check gives the count expected by chance at each
# threshold. A row is interesting only when observed clears expected by a visible margin.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})
here::i_am("config.yml")

cfg <- yaml::read_yaml(here("config.yml"))
stage <- here("03_Pathway_Enrichment", "04_singscore_pheno_associations")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports", "figures")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  singscore = "03_Pathway_Enrichment/02_run_singscore/c_data/singscore.rds",
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  proteins = "01_Preprocess/02_Quantification/c_data/proteins.rds",
  phenotype = "00_Input/phenotype.csv",
  metadata = "00_Input/metadata.csv"
)
paths <- map_chr(inputs, here::here)
if (!all(file.exists(paths))) {
  stop(
    "Run 02_run_singscore first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
scores <- readRDS(paths[["singscore"]])$scores
set_catalog <- readRDS(paths[["gene_sets"]])$set_catalog
targets <- readRDS(paths[["proteins"]])$targets
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
metadata <- readr::read_csv(paths[["metadata"]], show_col_types = FALSE)
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

# The design this analysis leans on, asserted rather than assumed. The one asymmetry is real:
# S06's left leg has T1 only, never acquired, and that leg is BFR, which is why one treatment
# contributes 32 legs and the other 33.
targets$leg_id <- paste(targets$participant, targets$leg, sep = "_")
legs_per_participant <- targets |>
  distinct(participant, leg, treatment) |>
  summarise(treatments = paste(sort(treatment), collapse = "+"), .by = participant)
treatment_agreement <- inner_join(
  distinct(targets, participant, leg, treatment),
  distinct(phenotype, participant, leg, treatment),
  by = c("participant", "leg"), suffix = c("_ms", "_pheno")
)
stopifnot(
  identical(colnames(scores), targets$sample_id),
  all(legs_per_participant$treatments == "BFR+HLRT"),
  nrow(anti_join(distinct(targets, participant, leg), phenotype,
    by = c("participant", "leg")
  )) == 0,
  all(treatment_agreement$treatment_ms == treatment_agreement$treatment_pheno),
  !"reinjection" %in% names(metadata)
)
message(
  "design: ", ncol(scores), " samples, ", n_distinct(targets$participant), " participants, ",
  n_distinct(targets$leg_id), " legs, treatment agrees ",
  sum(treatment_agreement$treatment_ms == treatment_agreement$treatment_pheno), "/",
  nrow(treatment_agreement)
)

# One change score per set per leg. A leg missing either timepoint drops out here.
sample_by_timepoint <- targets |>
  select(leg_id, timepoint, sample_id) |>
  pivot_wider(names_from = timepoint, values_from = sample_id)
legs <- targets |>
  distinct(leg_id, participant, leg, treatment) |>
  left_join(sample_by_timepoint, by = "leg_id") |>
  filter(!is.na(T1), !is.na(T2)) |>
  arrange(participant, treatment)
dropped <- setdiff(unique(targets$leg_id), legs$leg_id)
if (length(dropped)) message("legs without both timepoints: ", paste(dropped, collapse = ", "))

delta_score <- scores[, legs$T2] - scores[, legs$T1]
colnames(delta_score) <- legs$leg_id

# Vastus lateralis is the biopsied muscle and the primary outcome. Rectus femoris was not
# biopsied, so a set tracking VL but not RF is more convincing than one tracking both. Echo
# intensity indexes tissue quality rather than size.
outcomes <- c(
  vl_csa = "vl_csa_post_cm2 - vl_csa_pre_cm2",
  vl_echo = "vl_echo_post_au - vl_echo_pre_au",
  rf_csa = "rf_csa_post_cm2 - rf_csa_pre_cm2",
  rf_echo = "rf_echo_post_au - rf_echo_pre_au"
)
leg_phenotype <- legs |>
  left_join(select(phenotype, -treatment), by = c("participant", "leg")) |>
  mutate(!!!set_names(map(outcomes, \(e) rlang::parse_expr(e)), names(outcomes)))

# Spearman is Pearson on ranks, so one cor() call covers every set at once. The t
# approximation is what cor.test(exact = FALSE) uses.
spearman_by_set <- function(matrix_by_leg, outcome) {
  usable <- !is.na(outcome)
  matrix_by_leg <- matrix_by_leg[, usable, drop = FALSE]
  outcome <- outcome[usable]
  n <- length(outcome)
  r <- as.vector(stats::cor(t(matrix_by_leg), outcome, method = "spearman"))
  statistic <- r * sqrt((n - 2) / (1 - r^2))
  tibble(
    set_id = rownames(matrix_by_leg), n = n, r = r,
    p = 2 * stats::pt(-abs(statistic), n - 2)
  ) |>
    mutate(fdr = p.adjust(p, "BH"))
}
label_sets <- function(x) {
  left_join(x, select(set_catalog, set_id, database, pathway, theme), by = "set_id") |>
    relocate(database, pathway, theme, .after = set_id)
}

# 1. Pooled: does this pathway track adaptation at all, ignoring treatment.
pooled <- imap(outcomes, \(expression, name) {
  spearman_by_set(delta_score, leg_phenotype[[name]]) |> mutate(outcome = name, .before = 1)
}) |>
  list_rbind() |>
  label_sets() |>
  arrange(p)

# 2. By treatment: the same correlation computed inside each arm. Descriptive only; two
# separate confidence intervals are not a test of their difference.
by_treatment <- map(set_names(c("BFR", "HLRT")), function(arm) {
  keep <- leg_phenotype$treatment == arm
  imap(outcomes, \(expression, name) {
    spearman_by_set(delta_score[, keep, drop = FALSE], leg_phenotype[[name]][keep]) |>
      mutate(outcome = name, .before = 1)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "treatment") |>
  pivot_wider(
    id_cols = c(outcome, set_id), names_from = treatment,
    values_from = c(n, r, p), names_glue = "{treatment}_{.value}"
  ) |>
  label_sets()

# 3. Differential, the study question. Within each participant holding both legs complete,
# the BFR-minus-HLRT difference in score change against the same difference in outcome.
paired <- legs |>
  summarise(n_legs = n(), .by = participant) |>
  filter(n_legs == 2) |>
  pull(participant)
bfr <- filter(legs, participant %in% paired, treatment == "BFR") |> arrange(participant)
hlrt <- filter(legs, participant %in% paired, treatment == "HLRT") |> arrange(participant)
stopifnot(identical(bfr$participant, hlrt$participant))

paired_score <- delta_score[, bfr$leg_id] - delta_score[, hlrt$leg_id]
colnames(paired_score) <- bfr$participant
paired_outcome <- map(set_names(names(outcomes)), function(name) {
  value <- set_names(leg_phenotype[[name]], leg_phenotype$leg_id)
  value[bfr$leg_id] - value[hlrt$leg_id]
})
message("paired participants: ", length(paired))

differential <- imap(paired_outcome, \(value, name) {
  spearman_by_set(paired_score, value) |> mutate(outcome = name, .before = 1)
}) |>
  list_rbind() |>
  label_sets() |>
  arrange(p)

# What a table of this size returns under the null. Read this before reading any row above.
signal_check <- map(
  list(pooled = pooled, differential = differential),
  function(table) {
    map(c(0.05, 0.01, 0.001), function(threshold) {
      summarise(table,
        threshold = threshold, tests = n(),
        observed = sum(p < threshold, na.rm = TRUE),
        expected = round(n() * threshold), .by = outcome
      )
    }) |>
      list_rbind()
  }
) |>
  list_rbind(names_to = "analysis") |>
  mutate(ratio = round(observed / pmax(expected, 1), 2)) |>
  arrange(analysis, outcome, threshold)
print(as.data.frame(filter(signal_check, threshold == 0.01)))

message("strongest differential associations:")
print(as.data.frame(
  differential |>
    slice_head(n = 8) |>
    transmute(outcome, database,
      pathway = substr(pathway, 1, 44),
      n, r = round(r, 3), p = signif(p, 2), fdr = signif(fdr, 2)
    )
))

# One figure per outcome: the three sets with the largest paired correlation, one point per
# participant. The origin is where a participant's two legs moved together.
walk(names(outcomes), function(name) {
  top <- differential |>
    filter(outcome == name, is.finite(r)) |>
    slice_max(abs(r), n = 3) |>
    mutate(panel = paste0(
      sub("^[A-Z0-9]+_", "", pathway), "\nr = ", round(r, 2),
      ", p = ", signif(p, 2)
    ))
  plot_data <- map(seq_len(nrow(top)), function(i) {
    tibble(
      panel = top$panel[i], participant = colnames(paired_score),
      d_score = paired_score[top$set_id[i], ], d_outcome = paired_outcome[[name]]
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = top$panel))
  figure <- ggplot(plot_data, aes(d_score, d_outcome)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
    geom_smooth(
      method = "lm", formula = y ~ x, se = TRUE, colour = "#1f5e8c",
      fill = "#1f5e8c", alpha = 0.12, linewidth = 0.6
    ) +
    geom_point(size = 2, alpha = 0.85, colour = "#16202b") +
    facet_wrap(~panel, scales = "free_x") +
    labs(
      x = "BFR minus high load, change in set score",
      y = paste0("BFR minus high load, change in ", name),
      title = paste0("Differential association with ", name),
      caption = paste0(
        "one point per participant, n = ", length(paired),
        "; nominal p, uncorrected across ", nrow(filter(differential, outcome == name)), " sets"
      )
    ) +
    theme_minimal(base_size = 11) +
    theme(strip.text = element_text(size = 8.5), plot.title = element_text(face = "bold"))
  ggsave(file.path(figure_dir, paste0("differential_", name, ".png")),
    figure,
    width = 9, height = 3.8, dpi = 200, bg = "white"
  )
  message("drew differential_", name)
})

packages <- c("here", "limpa", "singscore", "dplyr", "tidyr", "purrr", "ggplot2")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
design_check <- tibble(
  samples = ncol(scores), participants = n_distinct(targets$participant),
  legs = n_distinct(targets$leg_id), legs_with_both_timepoints = nrow(legs),
  paired_participants = length(paired),
  treatment_agreement = paste0(
    sum(treatment_agreement$treatment_ms == treatment_agreement$treatment_pheno),
    "/", nrow(treatment_agreement)
  ),
  legs_dropped = if (length(dropped)) paste(dropped, collapse = ", ") else "none"
)

saveRDS(
  list(
    differential = differential, pooled = pooled, by_treatment = by_treatment,
    signal_check = signal_check, design_check = design_check,
    paired_score = paired_score, paired_outcome = paired_outcome,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config = cfg,
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "pheno_associations.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    signal_check = signal_check,
    design_check = design_check,
    differential = differential,
    by_treatment = by_treatment,
    pooled = pooled,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "04_singscore_pheno_associations.xlsx")
)
message("wrote pheno_associations.rds and 04_singscore_pheno_associations.xlsx")
