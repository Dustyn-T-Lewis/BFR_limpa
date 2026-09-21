# Group the tested sets into GO themes, score each theme per sample, and lay out how the
# themes move with training and how they couple to the phenotype.
#
# A theme is this pipeline's module: a GO ancestor holding at least theme_min_sets qualifying
# member sets. Its score is the mean singscore of its members, so one readable row replaces a
# dozen nested terms saying the same thing.
#
# Nominal p only at this stage. Nothing here is corrected and nothing is gated on q. A dashed
# cell border marks p < theme_nominal_p and bold text marks |r| at or above theme_bold_r,
# which are reading aids, not claims.

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
stage <- here("03_Pathway_Enrichment", "05_classify_pathway_themes")
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
scores <- readRDS(paths[["singscore"]])$scores
targets <- readRDS(paths[["proteins"]])$targets
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

# ---- themes and their scores ------------------------------------------------------------

membership <- set_catalog |>
  filter(qualifies, !is.na(theme)) |>
  select(theme, set_id, database, pathway, measured_size) |>
  add_count(theme, name = "sets_in_theme") |>
  filter(sets_in_theme >= cfg$theme_min_sets) |>
  arrange(desc(sets_in_theme), theme, set_id)
themes <- distinct(membership, theme, sets_in_theme)
stopifnot(nrow(themes) > 0, all(membership$set_id %in% rownames(scores)))

# Mean member score, so a theme's value is on the same scale as a set's.
theme_score <- membership |>
  split(~theme) |>
  map(\(part) colMeans(scores[part$set_id, , drop = FALSE])) |>
  do.call(what = rbind)
theme_score <- theme_score[themes$theme, ]
message(nrow(themes), " themes covering ", nrow(membership), " sets")

# ---- how each theme moved with training --------------------------------------------------

# Mean NES across member sets, which says which way the theme went and how consistently.
# fgsea_results already carries a per-set theme column, so drop it before joining the
# membership table or the join silently returns theme.x and theme.y.
theme_nes <- fgsea_results |>
  select(-theme) |>
  inner_join(select(membership, theme, set_id), by = "set_id") |>
  summarise(
    nes = mean(NES), nes_sd = sd(NES),
    members_nominal = sum(pval < cfg$theme_nominal_p), .by = c(theme, contrast)
  )

# ---- leg-level change, and the two phenotype questions -----------------------------------

targets$leg_id <- paste(targets$participant, targets$leg, sep = "_")
sample_by_timepoint <- targets |>
  select(leg_id, timepoint, sample_id) |>
  pivot_wider(names_from = timepoint, values_from = sample_id)
legs <- targets |>
  distinct(leg_id, participant, leg, treatment) |>
  left_join(sample_by_timepoint, by = "leg_id") |>
  filter(!is.na(T1), !is.na(T2)) |>
  arrange(participant, treatment)
delta_theme <- theme_score[, legs$T2] - theme_score[, legs$T1]
colnames(delta_theme) <- legs$leg_id

outcomes <- c(
  vl_csa = "vl_csa_post_cm2 - vl_csa_pre_cm2",
  vl_echo = "vl_echo_post_au - vl_echo_pre_au",
  rf_csa = "rf_csa_post_cm2 - rf_csa_pre_cm2",
  rf_echo = "rf_echo_post_au - rf_echo_pre_au"
)
leg_phenotype <- legs |>
  left_join(select(phenotype, -treatment), by = c("participant", "leg")) |>
  mutate(!!!set_names(map(outcomes, rlang::parse_expr), names(outcomes)))

# Spearman is Pearson on ranks, so one cor() call covers every theme at once.
spearman_by_theme <- function(matrix_by_column, outcome) {
  usable <- !is.na(outcome)
  matrix_by_column <- matrix_by_column[, usable, drop = FALSE]
  outcome <- outcome[usable]
  n <- length(outcome)
  r <- as.vector(stats::cor(t(matrix_by_column), outcome, method = "spearman"))
  tibble(
    theme = rownames(matrix_by_column), n = n, r = r,
    p = 2 * stats::pt(-abs(r * sqrt((n - 2) / (1 - r^2))), n - 2)
  )
}

pooled <- imap(outcomes, \(expression, name) {
  spearman_by_theme(delta_theme, leg_phenotype[[name]]) |> mutate(outcome = name)
}) |>
  list_rbind()

by_arm <- map(set_names(c("BFR", "HLRT")), function(arm) {
  keep <- leg_phenotype$treatment == arm
  imap(outcomes, \(expression, name) {
    spearman_by_theme(delta_theme[, keep, drop = FALSE], leg_phenotype[[name]][keep]) |>
      mutate(outcome = name)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "treatment")

# Every participant trained one leg under each condition, so the differential question is
# asked inside the person: BFR minus HLRT, in both the score change and the outcome change.
paired <- legs |>
  summarise(n_legs = n(), .by = participant) |>
  filter(n_legs == 2) |>
  pull(participant)
bfr <- filter(legs, participant %in% paired, treatment == "BFR") |> arrange(participant)
hlrt <- filter(legs, participant %in% paired, treatment == "HLRT") |> arrange(participant)
stopifnot(identical(bfr$participant, hlrt$participant))

paired_score <- delta_theme[, bfr$leg_id] - delta_theme[, hlrt$leg_id]
colnames(paired_score) <- bfr$participant
paired_outcome <- map(set_names(names(outcomes)), function(name) {
  value <- set_names(leg_phenotype[[name]], leg_phenotype$leg_id)
  value[bfr$leg_id] - value[hlrt$leg_id]
})
differential <- imap(paired_outcome, \(value, name) {
  spearman_by_theme(paired_score, value) |> mutate(outcome = name)
}) |>
  list_rbind()

# ---- panel A: the association map ---------------------------------------------------------

# Rows ordered by clustering the whole profile, so themes that behave alike sit together.
profile <- theme_nes |>
  select(theme, contrast, nes) |>
  pivot_wider(names_from = contrast, values_from = nes) |>
  column_to_rownames("theme")
theme_order <- rownames(profile)[hclust(dist(scale(profile)))$order]
theme_label <- set_names(
  paste0(themes$theme, "  (", themes$sets_in_theme, ")"), themes$theme
)

cell_style <- function(x, value_column) {
  mutate(x,
    nominal = !is.na(p) & p < cfg$theme_nominal_p,
    bold = abs(.data[[value_column]]) >= cfg$theme_bold_r,
    theme = factor(theme, levels = rev(theme_order)),
    label = sprintf("%+.2f", .data[[value_column]])
  )
}
draw_map <- function(data, value_column, scale_low, scale_high, scale_title,
                     block_title, show_rows) {
  ggplot(data, aes(column, theme, fill = .data[[value_column]])) +
    geom_tile(colour = "white", linewidth = 0.6) +
    geom_tile(
      data = filter(data, nominal), colour = "black", linewidth = 0.55,
      linetype = "22", fill = NA
    ) +
    geom_text(
      aes(label = label, fontface = if_else(bold, "bold", "plain")),
      size = 2.7, colour = "grey15"
    ) +
    facet_grid(~block, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(
      low = scale_low, mid = "white", high = scale_high, midpoint = 0,
      name = scale_title, limits = c(-1, 1) * max(abs(data[[value_column]]), na.rm = TRUE)
    ) +
    scale_y_discrete(labels = theme_label) +
    labs(x = NULL, y = NULL, title = block_title) +
    theme_minimal(base_size = 9) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 40, hjust = 1),
      axis.text.y = if (show_rows) element_text(size = 8) else element_blank(),
      strip.text = element_text(face = "bold", size = 8.5),
      plot.title = element_text(face = "bold", size = 10),
      legend.position = "bottom", legend.key.height = unit(3, "mm"),
      legend.key.width = unit(12, "mm"), legend.title = element_text(size = 8)
    )
}

nes_cells <- theme_nes |>
  filter(contrast %in% c("BFR_post_vs_pre", "HLRT_post_vs_pre", "interaction")) |>
  mutate(
    column = recode(contrast,
      BFR_post_vs_pre = "BFR", HLRT_post_vs_pre = "High load", interaction = "Interaction"
    ),
    column = factor(column, levels = c("BFR", "High load", "Interaction")),
    block = "Training response (mean NES)", p = NA_real_
  ) |>
  cell_style("nes") |>
  mutate(nominal = FALSE)

# Blocks and columns are factors, because facet_grid and the x axis otherwise sort
# alphabetically and put the differential before the pooled, and rf before vl.
pooled_block <- paste0("Coupled to phenotype, all legs (n=", max(pooled$n), ")")
differential_block <- paste0(
  "BFR minus high load, within participant (n=", max(differential$n), ")"
)
r_cells <- bind_rows(
  mutate(pooled, block = pooled_block),
  mutate(differential, block = differential_block)
) |>
  mutate(
    column = factor(outcome, levels = names(outcomes)),
    block = factor(block, levels = c(pooled_block, differential_block))
  ) |>
  cell_style("r")

panel_a <- draw_map(
  nes_cells, "nes", "#2166AC", "#B2182B", "mean NES",
  "A   Themes by training response", TRUE
) +
  draw_map(
    r_cells, "r", "#1b7837", "#762a83", "Spearman r",
    "      and by phenotype coupling", FALSE
  ) +
  plot_layout(widths = c(1, 2.6))

# ---- panel B: the coupling itself, for the themes that move most --------------------------

highlight <- differential |>
  filter(is.finite(r)) |>
  slice_max(abs(r), n = 6) |>
  arrange(desc(abs(r)))

coupling <- map(seq_len(nrow(highlight)), function(i) {
  row <- highlight[i, ]
  stats <- by_arm |>
    filter(theme == row$theme, outcome == row$outcome) |>
    mutate(text = sprintf("%s (n=%d): r=%+.2f, p=%.3f", treatment, n, r, p))
  tibble(
    facet = paste0(
      sub("^(.{46}).*", "\\1...", row$theme), "  vs  ", row$outcome,
      "\nwithin participant: r = ", sprintf("%+.2f", row$r), ", p = ", signif(row$p, 2)
    ),
    leg_id = colnames(delta_theme),
    d_score = delta_theme[row$theme, ],
    d_outcome = leg_phenotype[[row$outcome]],
    treatment = leg_phenotype$treatment,
    caption = paste(stats$text, collapse = "\n")
  )
}) |>
  list_rbind() |>
  mutate(facet = factor(facet, levels = unique(facet)))

panel_b <- ggplot(coupling, aes(d_score, d_outcome, colour = treatment, fill = treatment)) +
  geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.12, linewidth = 0.6) +
  geom_point(size = 1.6, alpha = 0.9) +
  geom_text(
    data = distinct(coupling, facet, caption), inherit.aes = FALSE,
    aes(x = -Inf, y = Inf, label = caption),
    hjust = -0.04, vjust = 1.15, size = 2.4, lineheight = 1.15, colour = "grey20"
  ) +
  facet_wrap(~facet, scales = "free", nrow = 2) +
  scale_colour_manual(values = c(BFR = "#B2182B", HLRT = "#2166AC"), name = NULL) +
  scale_fill_manual(values = c(BFR = "#B2182B", HLRT = "#2166AC"), name = NULL) +
  labs(
    x = "change in theme score (T2 - T1), one point per leg",
    y = "change in phenotype",
    title = "B   Theme-phenotype coupling, by training mode",
    caption = paste0(
      "Lines are per-arm fits on all legs; the header r is the within-participant ",
      "BFR-minus-high-load correlation across ", length(paired),
      " pairs. Nominal p, uncorrected."
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 7.2, lineheight = 1.2),
    plot.title = element_text(face = "bold", size = 10),
    plot.caption = element_text(hjust = 0, size = 7, colour = "grey35"),
    legend.position = "top"
  )

composite <- panel_a / panel_b + plot_layout(heights = c(1.5, 1))
ggsave(file.path(figure_dir, "theme_classification.png"), composite,
  width = 13.5, height = 13, dpi = 200, bg = "white"
)
ggsave(file.path(figure_dir, "theme_classification.pdf"), composite,
  width = 13.5, height = 13, bg = "white"
)
message("drew theme_classification")

# ---- supplement: every member set behind each drawn theme ---------------------------------

member_detail <- membership |>
  left_join(
    differential |>
      filter(outcome == "vl_csa") |>
      select(theme, theme_r = r, theme_p = p),
    by = "theme"
  ) |>
  left_join(
    fgsea_results |>
      filter(contrast == "BFR_post_vs_pre") |>
      select(set_id, bfr_nes = NES, bfr_p = pval),
    by = "set_id"
  ) |>
  arrange(desc(sets_in_theme), theme, bfr_p)

supplement <- member_detail |>
  filter(theme %in% highlight$theme) |>
  mutate(theme = factor(theme, levels = rev(unique(highlight$theme))))
supp_figure <- ggplot(supplement, aes(bfr_nes, theme)) +
  geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey70") +
  geom_point(aes(size = measured_size, colour = bfr_p < cfg$theme_nominal_p), alpha = 0.75) +
  scale_colour_manual(
    values = c(`TRUE` = "#B2182B", `FALSE` = "grey65"),
    labels = c(`TRUE` = "nominal p < 0.05", `FALSE` = "not nominal"), name = NULL
  ) +
  scale_size_continuous(range = c(1, 5), name = "measured genes") +
  scale_y_discrete(labels = theme_label) +
  labs(
    x = "member-set NES, BFR post vs pre", y = NULL,
    title = "Member sets behind each drawn theme",
    subtitle = "one point per qualifying set; a theme is only as coherent as its members"
  ) +
  theme_minimal(base_size = 9) +
  theme(plot.title = element_text(face = "bold"), legend.position = "bottom")
ggsave(file.path(figure_dir, "theme_members.png"), supp_figure,
  width = 9, height = 5.5, dpi = 200, bg = "white"
)
message("drew theme_members")

# ---- outputs ------------------------------------------------------------------------------

packages <- c("here", "fgsea", "singscore", "dplyr", "tidyr", "purrr", "ggplot2", "patchwork")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
nominal_counts <- bind_rows(
  mutate(pooled, analysis = "pooled"), mutate(differential, analysis = "differential")
) |>
  summarise(
    themes = n(), nominal = sum(p < cfg$theme_nominal_p),
    expected = round(n() * cfg$theme_nominal_p), .by = c(analysis, outcome)
  ) |>
  mutate(ratio = round(nominal / pmax(expected, 1), 2))
print(as.data.frame(nominal_counts))

saveRDS(
  list(
    themes = themes, membership = membership, theme_score = theme_score,
    theme_nes = theme_nes, pooled = pooled, by_arm = by_arm, differential = differential,
    nominal_counts = nominal_counts,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config = cfg,
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "pathway_themes.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    nominal_counts = nominal_counts,
    themes = themes,
    theme_training_nes = theme_nes,
    theme_differential = arrange(differential, p),
    theme_pooled = arrange(pooled, p),
    theme_by_arm = by_arm,
    theme_membership = member_detail,
    theme_scores = rownames_to_column(as.data.frame(theme_score), "theme"),
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "05_classify_pathway_themes.xlsx")
)
message("wrote pathway_themes.rds and 05_classify_pathway_themes.xlsx")
