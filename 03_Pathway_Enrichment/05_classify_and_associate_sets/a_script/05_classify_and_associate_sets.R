# Answer two questions about every tested set: how well it separates the study's groups, and
# whether it tracks the phenotype. The unit is the set, with no collapse and no grouping, and
# results are read per database so each collection carries its own chance expectation. Four
# independently curated collections agreeing is stronger than any single grouping of them.
#
# Every comparison here is paired. Pre versus post is the same leg twice; BFR versus high load
# is two legs of one person. AUC comes from pROC as the effect-size descriptor and the p-value
# from the paired Wilcoxon signed-rank test. An unpaired p would answer for a design this study
# did not run.
#
# Nominal p is read against chance_expectation, which carries the count a table that size
# returns under the null, per database. BH within database and task is reported beside it.
# baseline_BFR_vs_HLRT is the empirical floor, because no signal can exist there.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(purrr)
  library(ggplot2)
})

stage <- here("03_Pathway_Enrichment", "05_classify_and_associate_sets")
out <- file.path(stage, "c_data")
figure_dir <- file.path(stage, "b_reports")
for (path in c(out, figure_dir)) dir.create(path, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  set_tests = "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds",
  singscore = "03_Pathway_Enrichment/04_run_singscore/c_data/singscore.rds",
  proteins = "01_Preprocess/02_Quantification/c_data/proteins.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 00 through 02 first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
set_catalog <- readRDS(paths[["gene_sets"]])$set_catalog
fgsea_results <- readRDS(paths[["set_tests"]])$set_tests |> filter(method == "fgsea")
set_score <- readRDS(paths[["singscore"]])$scores
targets <- readRDS(paths[["proteins"]])$targets
phenotype <- readr::read_csv(paths[["phenotype"]], show_col_types = FALSE)
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

catalog <- set_catalog |>
  filter(qualifies) |>
  select(set_id, database, pathway, theme, measured_size)
stopifnot(setequal(catalog$set_id, rownames(set_score)))
set_score <- set_score[catalog$set_id, ]
collection_sizes <- count(catalog, database)
message(nrow(catalog), " sets across ", nrow(collection_sizes), " collections")

# The same labels the volcanoes and dot plots carry, with the line breaks flattened for a strip.
clean_name <- function(x) gsub("\n", " ", enrichVolcano::ev_clean_label(x))

# ---- the five classification tasks -------------------------------------------------------

targets$leg_id <- paste(targets$participant, targets$leg, sep = "_")
by_timepoint <- targets |>
  select(leg_id, participant, leg, treatment, timepoint, sample_id) |>
  pivot_wider(names_from = timepoint, values_from = sample_id)
legs <- filter(by_timepoint, !is.na(T1), !is.na(T2)) |> arrange(participant, treatment)

# One row per participant holding both of their legs, for whichever column identifies them.
pair_by_treatment <- function(data, column) {
  data |>
    filter(!is.na(.data[[column]])) |>
    select(participant, treatment, value = all_of(column)) |>
    pivot_wider(names_from = treatment, values_from = value) |>
    filter(!is.na(BFR), !is.na(HLRT))
}
baseline <- pair_by_treatment(by_timepoint, "T1")
post <- pair_by_treatment(by_timepoint, "T2")
delta_pairs <- pair_by_treatment(legs, "leg_id")

delta_set <- set_score[, legs$T2] - set_score[, legs$T1]
colnames(delta_set) <- legs$leg_id

# Each task names the matrix it reads and the two paired column sets. `favours` is the group an
# AUC above 0.5 points to, which is what makes the direction readable on the figures.
tasks <- list(
  pre_vs_post_BFR = list(
    matrix = "score", positive = filter(legs, treatment == "BFR")$T2,
    negative = filter(legs, treatment == "BFR")$T1,
    label = "Pre to post, BFR", favours = "post", unit = "legs"
  ),
  pre_vs_post_HLRT = list(
    matrix = "score", positive = filter(legs, treatment == "HLRT")$T2,
    negative = filter(legs, treatment == "HLRT")$T1,
    label = "Pre to post, HL", favours = "post", unit = "legs"
  ),
  baseline_BFR_vs_HLRT = list(
    matrix = "score", positive = baseline$BFR, negative = baseline$HLRT,
    label = "Baseline BFR vs HL (control)", favours = "BFR", unit = "participants"
  ),
  post_BFR_vs_HLRT = list(
    matrix = "score", positive = post$BFR, negative = post$HLRT,
    label = "Post BFR vs HL", favours = "BFR", unit = "participants"
  ),
  delta_BFR_vs_HLRT = list(
    matrix = "delta", positive = delta_pairs$BFR, negative = delta_pairs$HLRT,
    label = "Change, BFR vs HL", favours = "BFR", unit = "participants"
  )
)

# pROC::roc() auto-orients by default: on a case whose true directional AUC is 0.194 it returns
# 0.806. direction = "<" pins it, without which every below-chance theme flips and the red/blue
# encoding on the figures inverts silently.
fit_roc <- function(values, spec) {
  labels <- rep(c("neg", "pos"), c(length(spec$negative), length(spec$positive)))
  pROC::roc(labels, values[c(spec$negative, spec$positive)],
    levels = c("neg", "pos"), direction = "<", quiet = TRUE
  )
}
classify <- function(values, spec) {
  tibble(feature = rownames(values), n_pairs = length(spec$positive)) |>
    mutate(
      auc = apply(values, 1, \(row) as.numeric(pROC::auc(fit_roc(row, spec)))),
      p_paired = apply(values, 1, \(row) {
        suppressWarnings(
          stats::wilcox.test(row[spec$positive], row[spec$negative], paired = TRUE)$p.value
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
      mutate(task = name, task_label = spec$label, .before = 1)
  }) |>
    list_rbind()
}
set_auc <- run_tasks(set_score, delta_set) |>
  rename(set_id = feature) |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p_paired, "BH"), .by = c(task, database))

# Each collection gets its own denominator, so a 41-set collection is not read against a
# 1,366-set one. BH is applied inside the same two-way split for the same reason.
auc_counts <- set_auc |>
  summarise(
    features = n(), nominal = sum(p_paired < 0.05), expected = round(n() * 0.05, 1),
    fdr_sig = sum(fdr < 0.05), max_auc = round(max(auc), 2),
    .by = c(task, task_label, database, n_pairs)
  ) |>
  mutate(analysis = "classification", .before = 1)
print(as.data.frame(
  auc_counts |>
    mutate(ratio = round(nominal / expected, 2)) |>
    select(task, database, features, nominal, expected, ratio, fdr_sig) |>
    arrange(task, database)
))

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


paired_set <- delta_set[, delta_pairs$BFR] - delta_set[, delta_pairs$HLRT]
colnames(paired_set) <- delta_pairs$participant
paired_outcome <- map(set_names(names(outcomes)), function(name) {
  value <- set_names(leg_phenotype[[name]], leg_phenotype$leg_id)
  value[delta_pairs$BFR] - value[delta_pairs$HLRT]
})

# cor.test computes the exact Spearman p at these sample sizes. The t approximation this
# replaced was off by up to 9e-4, which is enough to move a result across the 0.05 line the
# figures report against.
spearman_by_row <- function(values, outcome) {
  usable <- !is.na(outcome)
  values <- values[, usable, drop = FALSE]
  outcome <- outcome[usable]
  # Ties make cor.test fall back from the exact p to its approximation and warn each time.
  # Reading the two fields off the htest rather than tidying it runs ten times faster over the
  # 32,000 tests this stage performs, and returns the same numbers to the bit.
  fits <- suppressWarnings(apply(values, 1, \(row) {
    test <- stats::cor.test(row, outcome, method = "spearman")
    c(r = unname(test$estimate), p = test$p.value)
  }))
  tibble(
    feature = rownames(values), n = length(outcome),
    r = unname(fits["r", ]), p = unname(fits["p", ])
  )
}
associate <- function(values, paired_values) {
  imap(outcomes, \(expression, name) {
    bind_rows(
      mutate(spearman_by_row(values, leg_phenotype[[name]]), analysis = "pooled"),
      mutate(spearman_by_row(paired_values, paired_outcome[[name]]), analysis = "differential")
    ) |>
      mutate(outcome = name)
  }) |>
    list_rbind() |>
    relocate(analysis, outcome)
}
set_association <- associate(delta_set, paired_set) |>
  rename(set_id = feature) |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(analysis, outcome, database))

by_arm <- map(set_names(c("BFR", "HLRT")), function(arm) {
  keep <- leg_phenotype$treatment == arm
  imap(outcomes, \(expression, name) {
    spearman_by_row(delta_set[, keep, drop = FALSE], leg_phenotype[[name]][keep]) |>
      mutate(outcome = name)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "treatment") |>
  rename(set_id = feature)

chance_expectation <- bind_rows(
  select(auc_counts, analysis,
    comparison = task_label, database, n_pairs,
    features, nominal, expected, fdr_sig
  ),
  set_association |>
    summarise(
      features = n(), nominal = sum(p < 0.05), expected = round(n() * 0.05, 1),
      fdr_sig = sum(fdr < 0.05), .by = c(analysis, outcome, database, n)
    ) |>
    mutate(analysis = paste0("association: ", analysis)) |>
    select(analysis,
      comparison = outcome, database, n_pairs = n,
      features, nominal, expected, fdr_sig
    )
) |>
  mutate(ratio = round(nominal / pmax(expected, 0.1), 2)) |>
  relocate(ratio, .after = expected)

# ---- figures: one per comparison, significant results only --------------------------------

# A task can reach nominal p in hundreds of sets, so a figure shows the strongest few from each
# collection rather than everything. Reading the collections side by side is the point: the
# full map of every set against every task and outcome lives in the workbook.
per_database <- 2

# The computed geometry assumes square panels, which the ROC and association grids are. A figure
# whose facets are not square passes its own width and height.
save_faceted <- function(plot, name, n_panels, columns, panel = 2.6, header = 2.1,
                         width = 1.2 + columns * panel,
                         height = header + ceiling(n_panels / columns) * panel) {
  figure <- plot +
    theme_minimal(base_size = 10) +
    theme(
      strip.text = element_text(size = 7.6, lineheight = 1.3, margin = margin(3, 3, 5, 3)),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(5, "mm"),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, colour = "grey30"),
      plot.caption = element_text(hjust = 0, size = 7.5, colour = "grey40"),
      legend.position = "top"
    )
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 220, bg = "white"
    )
  })
  message("drew ", name, ": ", n_panels, " panels")
  n_panels
}
chance_line <- paste0(
  "Nominal p, uncorrected. Each collection is read against its own expectation: ",
  paste(sprintf(
    "%s %.1f of %d", collection_sizes$database,
    collection_sizes$n * 0.05, collection_sizes$n
  ), collapse = ", "), "."
)

draw_roc_figure <- function(task_name) {
  spec <- tasks[[task_name]]
  hits <- set_auc |>
    filter(task == task_name, p_paired < 0.05) |>
    slice_min(p_paired, n = per_database, by = database, with_ties = FALSE) |>
    arrange(database, p_paired)
  if (!nrow(hits)) {
    message("no set reaches nominal p for ", task_name)
    return(0L)
  }
  values <- if (spec$matrix == "delta") delta_set else set_score
  curves <- pmap(hits, function(set_id, database, pathway, auc, p_paired, ...) {
    coordinates <- pROC::coords(fit_roc(values[set_id, ], spec), "all")
    tibble(
      fpr = 1 - coordinates$specificity, tpr = coordinates$sensitivity,
      # An AUC below 0.5 means the set separates the other way, which is a result about
      # direction rather than a failure, so it is coloured instead of hidden.
      direction = if_else(auc >= 0.5, "higher", "lower"),
      panel = paste0(
        database, ": ", stringr::str_trunc(clean_name(pathway), 32),
        "\nAUC ", sprintf("%.2f", auc), "    p = ", signif(p_paired, 2)
      )
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel))) |>
    arrange(panel, fpr, tpr)

  columns <- min(4, nrow(hits))
  save_faceted(
    ggplot(curves, aes(fpr, tpr, colour = direction, fill = direction)) +
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
        x = "1 - specificity", y = "sensitivity", title = spec$label,
        subtitle = sprintf(
          "singscore per set, %d paired %s, AUC with direction fixed", length(spec$positive),
          spec$unit
        ),
        caption = sprintf(
          paste(
            "Strongest %d sets per collection of %d reaching nominal p. Red separates higher in",
            "%s, blue lower. Shading is area under the curve, dashed line chance.",
            "p from the paired Wilcoxon signed-rank test. %s"
          ),
          per_database, sum(set_auc$task == task_name & set_auc$p_paired < 0.05),
          spec$favours, chance_line
        )
      ),
    paste0("roc_", task_name), nrow(hits), columns
  )
}
# The baseline control is tested and reported, in chance_expectation and in the README, but
# it gets no figure: it exists to say what the method returns when nothing is there, which is
# a number to read rather than a panel to present.
drawn_tasks <- setdiff(names(tasks), "baseline_BFR_vs_HLRT")
roc_drawn <- map_int(set_names(drawn_tasks), draw_roc_figure)

draw_association_figure <- function(which_analysis, title, subtitle) {
  hits <- set_association |>
    filter(analysis == which_analysis, p < 0.05) |>
    slice_min(p, n = per_database, by = database, with_ties = FALSE) |>
    arrange(database, p)
  if (!nrow(hits)) {
    message("no set reaches nominal p for ", which_analysis)
    return(0L)
  }
  points <- pmap(hits, function(set_id, database, pathway, outcome, n, r, p, ...) {
    arms <- by_arm |>
      filter(set_id == !!set_id, outcome == !!outcome) |>
      mutate(text = sprintf("%-5s n=%d  r=%+.2f  p=%.3f", treatment, n, r, p))
    tibble(
      panel = paste0(
        database, ": ", stringr::str_trunc(clean_name(pathway), 28), "   vs   ", outcome,
        "\nr = ", sprintf("%+.2f", r), "    p = ", signif(p, 2)
      ),
      d_score = delta_set[set_id, ], d_outcome = leg_phenotype[[outcome]],
      treatment = leg_phenotype$treatment, caption = paste(arms$text, collapse = "\n")
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel)))

  columns <- min(4, nrow(hits))
  save_faceted(
    ggplot(points, aes(d_score, d_outcome, colour = treatment, fill = treatment)) +
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
        x = "change in set score, T2 - T1, one point per leg",
        y = "change in phenotype", title = title,
        subtitle = sprintf("Spearman, %s analysis, %s", which_analysis, subtitle),
        caption = sprintf(
          paste(
            "Strongest %d set-outcome pairs per collection of %d reaching nominal p.",
            "Lines are fitted within each arm; the header r is the %s correlation",
            "that selected the panel. %s"
          ),
          per_database,
          sum(set_association$analysis == which_analysis & set_association$p < 0.05),
          which_analysis, chance_line
        )
      ),
    paste0("association_", which_analysis), nrow(hits), columns,
    panel = 3.1, header = 2.4
  )
}
invisible(draw_association_figure(
  "pooled", "Training response against phenotype, all legs", "65 legs, both arms pooled"
))
invisible(draw_association_figure(
  "differential", "BFR minus high load, within participant", "32 paired participants"
))

# Every collection with a nominal hit is represented on each figure that was drawn.
stopifnot(all(map_lgl(drawn_tasks, \(task_name) {
  available <- filter(set_auc, task == task_name, p_paired < 0.05)
  roc_drawn[[task_name]] == min(nrow(available), per_database * n_distinct(available$database))
})))
# Whether a collection clears chance is the headline, so it gets a panel of its own rather
# than living only in a sheet.
chance_figure <- chance_expectation |>
  filter(analysis == "classification") |>
  mutate(comparison = factor(comparison, levels = map_chr(tasks, "label")))
save_faceted(
  ggplot(chance_figure, aes(ratio, database, fill = ratio > 1)) +
    geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_text(aes(label = sprintf("%d of %d", nominal, features)),
      hjust = -0.12, size = 2.5, colour = "grey25"
    ) +
    facet_wrap(~comparison, ncol = 1) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey72"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(
      x = "observed nominal hits / chance expectation", y = NULL,
      title = "Nominal hits relative to chance, by collection",
      subtitle = sprintf(
        "paired Wilcoxon per set, %d collections, uncorrected p", nrow(collection_sizes)
      ),
      caption = paste(
        "Bar length is observed nominal hits divided by the count that collection returns under",
        "the null. Red clears 1, grey does not. Labels give observed of tested.",
        "Table: c_data/05_classify_and_associate_sets.xlsx, chance_expectation sheet."
      )
    ) +
    theme(panel.grid.major.y = element_blank()),
  "chance_by_database", n_distinct(chance_figure$comparison), 1,
  width = 6.5, height = 8
)

# ---- one workbook -------------------------------------------------------------------------

packages <- c("here", "fgsea", "pROC", "dplyr", "purrr", "ggplot2", "enrichVolcano")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
set_detail <- catalog |>
  left_join(
    fgsea_results |>
      filter(contrast == "BFR_Post-Pre") |>
      select(set_id, bfr_nes = NES, bfr_p = p),
    by = "set_id"
  ) |>
  arrange(database, bfr_p)

sheets <- list(
  chance_expectation = chance_expectation,
  set_auc = arrange(set_auc, p_paired),
  set_association = arrange(set_association, p),
  set_by_arm = by_arm,
  set_catalog = set_detail,
  set_scores = rownames_to_column(as.data.frame(set_score), "set_id"),
  input_manifest = manifest,
  package_versions = versions
)
descriptions <- c(
  chance_expectation = "Nominal hits against chance, per collection. Read this first.",
  set_auc = "How well each set separates each task. AUC from ranks, p from paired test.",
  set_association = "Set against phenotype change: pooled, then within participant.",
  set_by_arm = "The same correlation computed inside BFR and inside high load, descriptive.",
  set_catalog = "Every tested set with its collection, GO Slim theme and training NES.",
  set_scores = "The set by sample score matrix the analyses above were computed on.",
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
  file.path(out, "05_classify_and_associate_sets.xlsx")
)
saveRDS(
  list(
    catalog = catalog, set_auc = set_auc, set_association = set_association,
    by_arm = by_arm, chance_expectation = chance_expectation,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "set_results.rds"),
  compress = "xz"
)
combined <- file.path(figure_dir, "05_classify_and_associate_sets_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message(
  "wrote 05_classify_and_associate_sets.xlsx (", length(sheets) + 1, " sheets) and a ",
  length(pages), "-page figure PDF"
)
