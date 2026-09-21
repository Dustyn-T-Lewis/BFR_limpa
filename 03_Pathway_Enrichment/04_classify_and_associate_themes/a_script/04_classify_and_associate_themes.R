# Group the tested sets into GO themes, then answer two questions about them:
# how well a theme separates the study's groups, and whether it tracks the phenotype.
#
# A theme is this pipeline's module: a GO ancestor holding at least theme_min_sets qualifying
# members, scored as the mean singscore of those members. One readable row replaces a dozen
# nested terms saying the same thing.
#
# Every comparison here is paired. Pre versus post is the same leg twice; BFR versus high load
# is two legs of one person. So AUC is the effect-size descriptor, taken from the rank sums,
# and the p-value comes from the paired Wilcoxon signed-rank test. An unpaired p would answer
# for a design this study did not run.
#
# Nominal p only. Nothing is corrected. chance_expectation carries the count a table this size
# returns under the null, and baseline_BFR_vs_HLRT is the empirical floor because no signal
# can exist there.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})
here::i_am("config.yml")

cfg <- yaml::read_yaml(here("config.yml"))
stage <- here("03_Pathway_Enrichment", "04_classify_and_associate_themes")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports", "figures")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  fgsea = "03_Pathway_Enrichment/01_run_fgsea/c_data/fgsea.rds",
  singscore = "03_Pathway_Enrichment/02_run_singscore/c_data/singscore.rds",
  proteins = "01_Preprocess/02_Quantification/c_data/proteins.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here::here)
if (!all(file.exists(paths))) {
  stop(
    "Run 00 through 02 first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
set_catalog <- readRDS(paths[["gene_sets"]])$set_catalog
fgsea_results <- readRDS(paths[["fgsea"]])$fgsea_results
set_score <- readRDS(paths[["singscore"]])$scores
targets <- readRDS(paths[["proteins"]])$targets
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

# ---- themes ------------------------------------------------------------------------------

membership <- set_catalog |>
  filter(qualifies, !is.na(theme)) |>
  select(theme, set_id, database, pathway, measured_size) |>
  add_count(theme, name = "sets_in_theme") |>
  filter(sets_in_theme >= cfg$theme_min_sets) |>
  arrange(desc(sets_in_theme), theme, set_id)
themes <- distinct(membership, theme, sets_in_theme)
stopifnot(nrow(themes) > 0, all(membership$set_id %in% rownames(set_score)))

theme_score <- membership |>
  split(~theme) |>
  map(\(part) colMeans(set_score[part$set_id, , drop = FALSE])) |>
  do.call(what = rbind)
theme_score <- theme_score[themes$theme, ]
message(nrow(themes), " themes covering ", nrow(membership), " sets")

# fgsea_results already carries a per-set theme column, so drop it before joining the
# membership table or the join silently returns theme.x and theme.y.
theme_nes <- fgsea_results |>
  select(-theme) |>
  inner_join(select(membership, theme, set_id), by = "set_id") |>
  summarise(
    nes = mean(NES), nes_sd = sd(NES),
    members_nominal = sum(pval < cfg$theme_nominal_p), .by = c(theme, contrast)
  )

# ---- the five classification tasks -------------------------------------------------------

targets$leg_id <- paste(targets$participant, targets$leg, sep = "_")
by_timepoint <- targets |>
  select(leg_id, participant, leg, treatment, timepoint, sample_id) |>
  pivot_wider(names_from = timepoint, values_from = sample_id)
legs <- filter(by_timepoint, !is.na(T1), !is.na(T2)) |> arrange(participant, treatment)

# Legs of the same participant, one per treatment, at whichever timepoint is asked for.
leg_pairs <- function(column) {
  by_timepoint |>
    filter(!is.na(.data[[column]])) |>
    select(participant, treatment, sample = all_of(column)) |>
    pivot_wider(names_from = treatment, values_from = sample) |>
    filter(!is.na(BFR), !is.na(HLRT))
}
baseline <- leg_pairs("T1")
post <- leg_pairs("T2")
delta_pairs <- legs |>
  select(participant, treatment, leg_id) |>
  pivot_wider(names_from = treatment, values_from = leg_id) |>
  filter(!is.na(BFR), !is.na(HLRT))

delta_score <- theme_score[, legs$T2] - theme_score[, legs$T1]
colnames(delta_score) <- legs$leg_id
delta_set <- set_score[, legs$T2] - set_score[, legs$T1]
colnames(delta_set) <- legs$leg_id

# Each task names the matrix it reads and the two paired column sets. positive is the group an
# AUC above 0.5 favours, which is what makes the direction readable.
tasks <- list(
  pre_vs_post_BFR = list(
    matrix = "score", positive = filter(legs, treatment == "BFR")$T2,
    negative = filter(legs, treatment == "BFR")$T1,
    label = "Pre to post, BFR", favours = "post"
  ),
  pre_vs_post_HLRT = list(
    matrix = "score", positive = filter(legs, treatment == "HLRT")$T2,
    negative = filter(legs, treatment == "HLRT")$T1,
    label = "Pre to post, HL", favours = "post"
  ),
  baseline_BFR_vs_HLRT = list(
    matrix = "score", positive = baseline$BFR, negative = baseline$HLRT,
    label = "Baseline BFR vs HL (control)", favours = "BFR"
  ),
  post_BFR_vs_HLRT = list(
    matrix = "score", positive = post$BFR, negative = post$HLRT,
    label = "Post BFR vs HL", favours = "BFR"
  ),
  delta_BFR_vs_HLRT = list(
    matrix = "delta", positive = delta_pairs$BFR, negative = delta_pairs$HLRT,
    label = "Change, BFR vs HL", favours = "BFR"
  )
)

# AUC from the rank sums, vectorised over every row at once; p from the paired signed-rank.
classify <- function(values, task) {
  positive <- values[, task$positive, drop = FALSE]
  negative <- values[, task$negative, drop = FALSE]
  n_pos <- ncol(positive)
  n_neg <- ncol(negative)
  ranks <- t(apply(cbind(positive, negative), 1, rank))
  auc <- (rowSums(ranks[, seq_len(n_pos), drop = FALSE]) - n_pos * (n_pos + 1) / 2) /
    (n_pos * n_neg)
  tibble(
    feature = rownames(values), n_pairs = n_pos, auc = auc,
    p_paired = map_dbl(seq_len(nrow(values)), \(i) {
      suppressWarnings(
        stats::wilcox.test(positive[i, ], negative[i, ], paired = TRUE)$p.value
      )
    })
  )
}
run_tasks <- function(score_matrix, delta_matrix) {
  # The argument is `spec`, not `task`: mutate() creates a `task` column first, and a later
  # `task$label` in the same call would read that column instead of the argument.
  imap(tasks, function(spec, name) {
    values <- if (spec$matrix == "delta") delta_matrix else score_matrix
    classify(values, spec) |>
      mutate(task = name, task_label = spec$label, favours = spec$favours, .before = 1)
  }) |>
    list_rbind()
}
theme_auc <- run_tasks(theme_score, delta_score) |> rename(theme = feature)
set_auc <- run_tasks(set_score, delta_set) |> rename(set_id = feature)

auc_counts <- theme_auc |>
  summarise(
    features = n(), nominal = sum(p_paired < cfg$theme_nominal_p),
    expected = round(n() * cfg$theme_nominal_p, 1),
    max_auc = round(max(auc), 2), .by = c(task, task_label, n_pairs)
  ) |>
  mutate(analysis = "classification", .before = 1)
print(as.data.frame(select(auc_counts, task, n_pairs, features, nominal, expected, max_auc)))

# ---- phenotype association ---------------------------------------------------------------

outcomes <- c(
  vl_csa = "vl_csa_post_cm2 - vl_csa_pre_cm2",
  vl_echo = "vl_echo_post_au - vl_echo_pre_au",
  rf_csa = "rf_csa_post_cm2 - rf_csa_pre_cm2",
  rf_echo = "rf_echo_post_au - rf_echo_pre_au"
)
leg_phenotype <- legs |>
  left_join(select(phenotype, -treatment), by = c("participant", "leg")) |>
  mutate(!!!set_names(map(outcomes, rlang::parse_expr), names(outcomes)))

# Spearman is Pearson on ranks, so one cor() call covers every feature at once.
spearman_by_row <- function(values, outcome) {
  usable <- !is.na(outcome)
  values <- values[, usable, drop = FALSE]
  outcome <- outcome[usable]
  n <- length(outcome)
  r <- as.vector(stats::cor(t(values), outcome, method = "spearman"))
  tibble(
    feature = rownames(values), n = n, r = r,
    p = 2 * stats::pt(-abs(r * sqrt((n - 2) / (1 - r^2))), n - 2)
  )
}

# Pooled across all legs, then within each arm, then the differential: every participant
# trained one leg each way, so BFR minus high load is taken inside the person.
paired_score <- delta_score[, delta_pairs$BFR] - delta_score[, delta_pairs$HLRT]
colnames(paired_score) <- delta_pairs$participant
paired_set <- delta_set[, delta_pairs$BFR] - delta_set[, delta_pairs$HLRT]
paired_outcome <- map(set_names(names(outcomes)), function(name) {
  value <- set_names(leg_phenotype[[name]], leg_phenotype$leg_id)
  value[delta_pairs$BFR] - value[delta_pairs$HLRT]
})

associate <- function(values, paired_values) {
  pooled <- imap(outcomes, \(expression, name) {
    spearman_by_row(values, leg_phenotype[[name]]) |>
      mutate(analysis = "pooled", outcome = name)
  }) |>
    list_rbind()
  differential <- imap(paired_outcome, \(value, name) {
    spearman_by_row(paired_values, value) |>
      mutate(analysis = "differential", outcome = name)
  }) |>
    list_rbind()
  bind_rows(pooled, differential) |> relocate(analysis, outcome)
}
theme_association <- associate(delta_score, paired_score) |> rename(theme = feature)
set_association <- associate(delta_set, paired_set) |> rename(set_id = feature)

by_arm <- map(set_names(c("BFR", "HLRT")), function(arm) {
  keep <- leg_phenotype$treatment == arm
  imap(outcomes, \(expression, name) {
    spearman_by_row(delta_score[, keep, drop = FALSE], leg_phenotype[[name]][keep]) |>
      mutate(outcome = name)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "treatment") |>
  rename(theme = feature)

association_counts <- theme_association |>
  summarise(
    features = n(), nominal = sum(p < cfg$theme_nominal_p),
    expected = round(n() * cfg$theme_nominal_p, 1),
    .by = c(analysis, outcome, n)
  ) |>
  rename(n_pairs = n) |>
  mutate(analysis = paste0("association: ", analysis), .before = 1)

chance_expectation <- bind_rows(
  select(auc_counts, analysis, comparison = task_label, n_pairs, features, nominal, expected),
  select(association_counts, analysis,
    comparison = outcome, n_pairs, features,
    nominal, expected
  )
) |>
  mutate(ratio = round(nominal / pmax(expected, 0.1), 2))
print(as.data.frame(chance_expectation))

# ---- one figure per comparison, significant themes only -----------------------------------

# Every figure here draws only what reached nominal p. The full map of every theme against
# every task lives in the workbook, not on a figure.

roc_curve <- function(values, positive_columns, negative_columns) {
  positive <- values[positive_columns]
  negative <- values[negative_columns]
  cuts <- sort(unique(c(positive, negative)), decreasing = TRUE)
  tibble(
    fpr = c(0, map_dbl(cuts, \(t) mean(negative >= t)), 1),
    tpr = c(0, map_dbl(cuts, \(t) mean(positive >= t)), 1)
  )
}

draw_roc_figure <- function(task_name) {
  spec <- tasks[[task_name]]
  hits <- theme_auc |>
    filter(task == task_name, p_paired < cfg$theme_nominal_p) |>
    arrange(p_paired)
  if (!nrow(hits)) {
    message("no theme reaches nominal p for ", task_name)
    return(invisible(NULL))
  }
  values <- if (spec$matrix == "delta") delta_score else theme_score
  curves <- pmap(hits, function(theme, auc, p_paired, ...) {
    roc_curve(values[theme, ], spec$positive, spec$negative) |>
      mutate(
        # AUC below 0.5 means the theme separates the other way, which is a result and not
        # a failure, so the curve is drawn in the down colour rather than the up one.
        direction = if_else(auc >= 0.5, "higher", "lower"),
        panel = paste0(
          stringr::str_trunc(theme, 34),
          "\nAUC ", sprintf("%.2f", auc), "    p = ", signif(p_paired, 2)
        )
      )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel)))

  columns <- min(3, nrow(hits))
  rows <- ceiling(nrow(hits) / columns)
  figure <- ggplot(curves, aes(fpr, tpr, colour = direction, fill = direction)) +
    geom_abline(linetype = "22", linewidth = 0.35, colour = "grey60") +
    geom_ribbon(aes(ymin = 0, ymax = tpr), alpha = 0.16, colour = NA) +
    geom_step(linewidth = 0.8, direction = "hv") +
    scale_colour_manual(values = c(higher = "#B2182B", lower = "#2166AC"), guide = "none") +
    scale_fill_manual(values = c(higher = "#B2182B", lower = "#2166AC"), guide = "none") +
    facet_wrap(~panel, ncol = columns) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    scale_x_continuous(breaks = c(0, 0.5, 1)) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(
      x = "1 - specificity", y = "sensitivity",
      title = spec$label,
      subtitle = paste0(
        nrow(hits), " of ", nrow(themes), " themes reach nominal p    |    ",
        spec$n_pairs_label, "    |    red: higher in ", spec$favours,
        "    blue: lower in ", spec$favours
      ),
      caption = paste0(
        "AUC from rank sums; p from the paired Wilcoxon signed-rank test. Nominal p, ",
        "uncorrected \u2014 ", round(nrow(themes) * cfg$theme_nominal_p, 1),
        " of ", nrow(themes), " is what chance returns."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 8, lineheight = 1.3, margin = margin(3, 3, 5, 3)),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(6, "mm"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40")
    )
  name <- paste0("roc_", task_name)
  ggsave(file.path(figure_dir, paste0(name, ".png")), figure,
    width = 1.2 + columns * 2.6, height = 1.9 + rows * 2.6, dpi = 220, bg = "white"
  )
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), figure,
    width = 1.2 + columns * 2.6, height = 1.9 + rows * 2.6, bg = "white"
  )
  message("drew ", name, ": ", nrow(hits), " panels")
  nrow(hits)
}
tasks <- imap(tasks, function(spec, name) {
  n <- length(spec$positive)
  spec$n_pairs_label <- paste0(n, " paired ", if (grepl("^pre", name)) "legs" else "participants")
  spec
})
roc_drawn <- map_int(set_names(names(tasks)), \(name) as.integer(draw_roc_figure(name) %||% 0L))

# ---- association scatters, one figure per analysis -----------------------------------------

draw_association_figure <- function(which_analysis, title, subtitle) {
  hits <- theme_association |>
    filter(analysis == which_analysis, p < cfg$theme_nominal_p) |>
    arrange(p)
  if (!nrow(hits)) {
    message("no theme reaches nominal p for ", which_analysis)
    return(invisible(NULL))
  }
  points <- pmap(hits, function(analysis, outcome, theme, n, r, p) {
    arms <- by_arm |>
      filter(theme == !!theme, outcome == !!outcome) |>
      mutate(text = sprintf("%-5s n=%d  r=%+.2f  p=%.3f", treatment, n, r, p))
    tibble(
      panel = paste0(
        stringr::str_trunc(theme, 34), "   vs   ", outcome,
        "\nr = ", sprintf("%+.2f", r), "    p = ", signif(p, 2)
      ),
      d_score = delta_score[theme, ],
      d_outcome = leg_phenotype[[outcome]],
      treatment = leg_phenotype$treatment,
      caption = paste(arms$text, collapse = "\n")
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel)))

  columns <- min(3, nrow(hits))
  rows <- ceiling(nrow(hits) / columns)
  figure <- ggplot(points, aes(d_score, d_outcome, colour = treatment, fill = treatment)) +
    geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
    geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
    geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.12, linewidth = 0.6) +
    geom_point(size = 1.9, alpha = 0.9) +
    geom_text(
      data = distinct(points, panel, caption), inherit.aes = FALSE,
      aes(x = -Inf, y = Inf, label = caption), family = "mono",
      hjust = -0.05, vjust = 1.25, size = 2.3, lineheight = 1.25, colour = "grey25"
    ) +
    facet_wrap(~panel, scales = "free", ncol = columns) +
    scale_y_continuous(expand = expansion(mult = c(0.06, 0.3))) +
    scale_colour_manual(values = c(BFR = "#B2182B", HLRT = "#2166AC"), name = NULL) +
    scale_fill_manual(values = c(BFR = "#B2182B", HLRT = "#2166AC"), name = NULL) +
    labs(
      x = "change in theme score, T2 - T1, one point per leg",
      y = "change in phenotype", title = title,
      subtitle = paste0(
        nrow(hits), " of ", nrow(theme_association) / 2,
        " theme-outcome pairs reach nominal p    |    ", subtitle
      ),
      caption = paste0(
        "Lines are fitted within each arm; the header r is the ", which_analysis,
        " correlation that selected the panel. Nominal p, uncorrected \u2014 ",
        round(nrow(themes) * cfg$theme_nominal_p, 1), " per outcome is what chance returns."
      )
    ) +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 7.6, lineheight = 1.3, margin = margin(3, 3, 5, 3)),
      panel.spacing = unit(5, "mm"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40"),
      legend.position = "top"
    )
  name <- paste0("association_", which_analysis)
  ggsave(file.path(figure_dir, paste0(name, ".png")), figure,
    width = 1.2 + columns * 3.1, height = 2.2 + rows * 2.9, dpi = 220, bg = "white"
  )
  ggsave(file.path(figure_dir, paste0(name, ".pdf")), figure,
    width = 1.2 + columns * 3.1, height = 2.2 + rows * 2.9, bg = "white"
  )
  message("drew ", name, ": ", nrow(hits), " panels")
  nrow(hits)
}
association_drawn <- c(
  pooled = draw_association_figure(
    "pooled", "Training response against phenotype, all legs",
    "65 legs, both arms pooled"
  ),
  differential = draw_association_figure(
    "differential", "BFR minus high load, within participant",
    "32 paired participants"
  )
)

# Everything nominally significant is on a figure, asserted rather than assumed.
stopifnot(
  sum(roc_drawn) == sum(theme_auc$p_paired < cfg$theme_nominal_p),
  sum(association_drawn) == sum(theme_association$p < cfg$theme_nominal_p)
)

# ---- one workbook -------------------------------------------------------------------------

packages <- c("here", "fgsea", "singscore", "dplyr", "tidyr", "purrr", "ggplot2", "patchwork")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
member_detail <- membership |>
  left_join(
    fgsea_results |>
      filter(contrast == "BFR_post_vs_pre") |>
      select(set_id, bfr_nes = NES, bfr_p = pval),
    by = "set_id"
  ) |>
  arrange(desc(sets_in_theme), theme, bfr_p)

sheets <- list(
  chance_expectation = chance_expectation,
  theme_auc = arrange(theme_auc, p_paired),
  theme_association = arrange(theme_association, p),
  theme_by_arm = by_arm,
  theme_training_nes = theme_nes,
  theme_membership = member_detail,
  set_auc = arrange(set_auc, p_paired),
  set_association = arrange(set_association, p),
  theme_scores = rownames_to_column(as.data.frame(theme_score), "theme"),
  input_manifest = manifest,
  package_versions = versions
)
descriptions <- c(
  chance_expectation = "Nominal hits against the count chance returns. Read this first.",
  theme_auc = "How well each theme separates each task. AUC from ranks, p from paired test.",
  theme_association = "Theme against phenotype change: pooled, then within participant.",
  theme_by_arm = "The same correlation computed inside BFR and inside high load, descriptive.",
  theme_training_nes = "Mean member NES per contrast, and how many members reached nominal p.",
  theme_membership = "Every set behind every theme, with its own NES. Drill from a theme here.",
  set_auc = "The same classification at set level, all qualifying sets, all five tasks.",
  set_association = "The same phenotype correlations at set level.",
  theme_scores = "The theme by sample score matrix the analyses above were computed on.",
  input_manifest = "Which files were read and their checksums.",
  package_versions = "Package versions at the time of the run."
)
read_me <- tibble(
  sheet = names(sheets),
  rows = map_int(sheets, nrow),
  holds = unname(descriptions[names(sheets)])
)
writexl::write_xlsx(
  c(list(read_me = read_me), sheets),
  file.path(out, "04_classify_and_associate_themes.xlsx")
)
saveRDS(
  c(
    list(
      themes = themes, membership = membership, theme_score = theme_score,
      theme_auc = theme_auc, theme_association = theme_association, by_arm = by_arm,
      set_auc = set_auc, set_association = set_association,
      theme_nes = theme_nes, chance_expectation = chance_expectation,
      provenance = list(
        created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config = cfg,
        inputs = manifest, packages = versions
      )
    )
  ),
  file.path(out, "theme_results.rds"),
  compress = "xz"
)
message("wrote 04_classify_and_associate_themes.xlsx (", length(sheets) + 1, " sheets)")
