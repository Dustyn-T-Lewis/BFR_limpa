# Per set, uncollapsed: how well it separates the study groups (pROC AUC as effect size, paired
# Wilcoxon signed-rank p) and whether it tracks phenotype. Every comparison is paired: pre/post is
# one leg twice, BFR/HLRT two legs of one person; an unpaired p fits a design not run.
# Nominal p is read per database against chance_expectation, with BH within database and task
# beside it; agreement across independent collections beats any single grouping of them.
# baseline_BFR_vs_HLRT is the empirical floor, because no signal can exist there.

pacman::p_load(here, dplyr, tibble, tidyr, purrr, stringr, ggplot2, readr, readxl, writexl)

stage <- here("03_Pathway_Enrichment", "05_classify_and_associate_sets")
set_catalog <- read_excel(
  here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data", "00_build_gene_sets.xlsx"),
  "set_catalog"
)
set_score <- readRDS(here("03_Pathway_Enrichment", "04_run_singscore", "c_data", "set_scores.rds"))
targets <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))$targets
phenotype <- read_csv(here("00_Input", "phenotype.csv"), show_col_types = FALSE)

catalog <- set_catalog |>
  filter(qualifies) |>
  select(set_id, database, pathway, measured_size)
set_score <- set_score[catalog$set_id, ]
collection_sizes <- count(catalog, database)
message(nrow(catalog), " sets across ", nrow(collection_sizes), " collections")


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

# Each task names its matrix and two paired column sets. `favours` is the group an AUC above 0.5
# points to, so the figures can state direction.
tasks <- list(
  pre_vs_post_BFR = list(
    matrix = "score", positive = filter(legs, treatment == "BFR")$T2,
    negative = filter(legs, treatment == "BFR")$T1,
    label = "Pre to post, BFR", favours = "post", unit = "legs"
  ),
  pre_vs_post_HLRT = list(
    matrix = "score", positive = filter(legs, treatment == "HLRT")$T2,
    negative = filter(legs, treatment == "HLRT")$T1,
    label = "Pre to post, HLRT", favours = "post", unit = "legs"
  ),
  baseline_BFR_vs_HLRT = list(
    matrix = "score", positive = baseline$BFR, negative = baseline$HLRT,
    label = "Baseline BFR vs HLRT (control)", favours = "BFR", unit = "participants"
  ),
  post_BFR_vs_HLRT = list(
    matrix = "score", positive = post$BFR, negative = post$HLRT,
    label = "Post BFR vs HLRT", favours = "BFR", unit = "participants"
  ),
  delta_BFR_vs_HLRT = list(
    matrix = "delta", positive = delta_pairs$BFR, negative = delta_pairs$HLRT,
    label = "Change, BFR vs HLRT", favours = "BFR", unit = "participants"
  )
)

# pROC::roc() auto-orients by default: a case with true directional AUC 0.194 returns 0.806.
# direction = "<" pins it; without it every below-chance set flips and the figures' red/blue
# encoding silently inverts.
fit_roc <- function(values, spec) {
  pROC::roc(
    controls = values[spec$negative], cases = values[spec$positive],
    direction = "<", quiet = TRUE
  )
}
# The argument is `spec`, not `task`: tibble() builds a `task` column first, and a later
# `task$label` in the same call would read that column instead of the argument.
set_auc <- imap(tasks, function(spec, name) {
  values <- if (spec$matrix == "delta") delta_set else set_score
  tibble(
    task = name, task_label = spec$label, set_id = rownames(values),
    n_pairs = length(spec$positive),
    auc = apply(values, 1, \(row) as.numeric(pROC::auc(fit_roc(row, spec)))),
    p_paired = apply(values, 1, \(row) {
      suppressWarnings(wilcox.test(row[spec$positive], row[spec$negative], paired = TRUE)$p.value)
    })
  )
}) |>
  list_rbind() |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p_paired, "BH"), .by = c(task, database))


# ---- phenotype association ---------------------------------------------------------------

outcomes <- c("vl_csa", "vl_echo", "rf_csa", "rf_echo")
# A duplicated phenotype row would silently misalign every leg below.
leg_phenotype <- legs |>
  left_join(
    select(phenotype, -treatment),
    by = c("participant", "leg"), relationship = "one-to-one"
  ) |>
  mutate(
    vl_csa = vl_csa_post_cm2 - vl_csa_pre_cm2,
    vl_echo = vl_echo_post_au - vl_echo_pre_au,
    rf_csa = rf_csa_post_cm2 - rf_csa_pre_cm2,
    rf_echo = rf_echo_post_au - rf_echo_pre_au
  )

paired_set <- delta_set[, delta_pairs$BFR] - delta_set[, delta_pairs$HLRT]
colnames(paired_set) <- delta_pairs$participant
paired_outcome <- map(set_names(outcomes), function(name) {
  value <- set_names(leg_phenotype[[name]], leg_phenotype$leg_id)
  value[delta_pairs$BFR] - value[delta_pairs$HLRT]
})

# cor.test gives the exact Spearman p at these sample sizes. The t approximation is off by up
# to 9e-4, enough to move a result across the 0.05 line the figures report against.
spearman_by_row <- function(values, outcome) {
  usable <- !is.na(outcome)
  values <- values[, usable, drop = FALSE]
  outcome <- outcome[usable]
  # Ties make cor.test fall back from the exact p to its approximation and warn each time.
  # Reading two fields off the htest, not tidying it, is ten times faster over the 32,000 tests
  # here and returns the same numbers to the bit.
  fits <- suppressWarnings(apply(values, 1, \(row) {
    test <- cor.test(row, outcome, method = "spearman")
    c(r = unname(test$estimate), p = test$p.value)
  }))
  tibble(
    set_id = rownames(values), n = length(outcome),
    r = unname(fits["r", ]), p = unname(fits["p", ])
  )
}
set_association <- map(set_names(outcomes), \(name) {
  bind_rows(
    pooled = spearman_by_row(delta_set, leg_phenotype[[name]]),
    differential = spearman_by_row(paired_set, paired_outcome[[name]]),
    .id = "analysis"
  ) |>
    mutate(outcome = name)
}) |>
  list_rbind() |>
  relocate(analysis, outcome) |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(analysis, outcome, database))

by_arm <- map(set_names(c("BFR", "HLRT")), function(arm) {
  keep <- leg_phenotype$treatment == arm
  map(set_names(outcomes), \(name) {
    spearman_by_row(delta_set[, keep, drop = FALSE], leg_phenotype[[name]][keep]) |>
      mutate(outcome = name)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "treatment")

# Each collection gets its own denominator, so a 41-set collection is not read against a
# 1,366-set one. BH is applied inside the same two-way split for the same reason.
chance_expectation <- bind_rows(
  transmute(set_auc,
    analysis = "classification", comparison = task_label, database, n_pairs,
    p = p_paired, fdr
  ),
  transmute(set_association,
    analysis = paste0("association: ", analysis), comparison = outcome, database,
    n_pairs = n, p, fdr
  )
) |>
  summarise(
    features = n(), nominal = sum(p < 0.05), expected = round(n() * 0.05, 1),
    ratio = round(nominal / expected, 2), fdr_sig = sum(fdr < 0.05),
    .by = c(analysis, comparison, database, n_pairs)
  )
print(as.data.frame(filter(chance_expectation, analysis == "classification")))


# ---- figures: S13 to S19, A4 pages ---------------------------------------------------------

figure_theme <- theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 6.8, lineheight = 1.2, margin = margin(3, 3, 4, 3)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(4, "mm"),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    plot.caption.position = "plot",
    legend.position = "top"
  )
supplement <- function(number, title, text, page = 1, n_pages = 1) {
  pages <- if (n_pages > 1) sprintf(" (page %d of %d)", page, n_pages) else ""
  str_wrap(sprintf("S%d Figure%s. %s. %s", number, pages, title, text), 115)
}
chance_line <- paste0(
  "Nominal p, uncorrected. Each collection is read against its own expectation: ",
  paste(sprintf(
    "%s %.1f of %d", collection_sizes$database,
    collection_sizes$n * 0.05, collection_sizes$n
  ), collapse = ", "), "."
)
set_label <- function(database, pathway) {
  str_wrap(paste0(database, ": ", enrichVolcano::clean_label(pathway, width = 1000)), 40)
}

# Every set reaching nominal p gets a panel, twelve to a page, by collection then p. `draw` builds
# the plot from one page's rows: a paginating facet would build every panel for every page.
paginate <- function(data, draw, number, title, text, scales = "fixed") {
  page_of <- (as.integer(data$panel) - 1) %/% 12 + 1
  n_pages <- max(page_of)
  map(seq_len(n_pages), \(page) {
    draw(droplevels(data[page_of == page, ])) +
      facet_wrap(~panel, ncol = 3, nrow = 4, scales = scales) +
      labs(caption = supplement(number, title, text, page, n_pages)) +
      figure_theme
  })
}

roc_figure <- function(task_name, number) {
  spec <- tasks[[task_name]]
  hits <- set_auc |>
    filter(task == task_name, p_paired < 0.05) |>
    arrange(database, p_paired)
  values <- if (spec$matrix == "delta") delta_set else set_score
  curves <- pmap(hits, function(set_id, database, pathway, auc, p_paired, fdr, ...) {
    coordinates <- pROC::coords(fit_roc(values[set_id, ], spec), "all")
    tibble(
      fpr = 1 - coordinates$specificity, tpr = coordinates$sensitivity,
      # An AUC below 0.5 separates the other way: a direction, not a failure, so it is
      # coloured, not hidden.
      direction = if_else(auc >= 0.5, "higher", "lower"),
      panel = paste0(
        set_label(database, pathway),
        "\nAUC ", sprintf("%.2f", auc), "   p ", signif(p_paired, 2), "   q ", signif(fdr, 2)
      )
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel))) |>
    arrange(panel, fpr, tpr)

  paginate(curves, \(page) {
    ggplot(page, aes(fpr, tpr, colour = direction, fill = direction)) +
      geom_abline(linetype = "22", linewidth = 0.35, colour = "grey60") +
      geom_ribbon(aes(ymin = 0, ymax = tpr), alpha = 0.16, colour = NA) +
      geom_step(linewidth = 0.7, direction = "hv") +
      scale_colour_manual(
        values = c(higher = "#B2182B", lower = "#2166AC"), aesthetics = c("colour", "fill"),
        guide = "none"
      ) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      labs(x = "1 - specificity", y = "sensitivity")
  }, number, paste("ROC curves for set scores,", spec$label), paste(
    sprintf(
      "ROC curves for all %d sets whose singscore separates the groups at nominal p, across %d",
      nrow(hits), length(spec$positive)
    ),
    sprintf("paired %s, by collection then p. AUC has its direction fixed: red", spec$unit),
    sprintf("separates higher in %s, blue lower.", spec$favours),
    "Shading is area under the curve, dashed line chance. p from the paired Wilcoxon",
    "signed-rank test, q from BH within collection and task.", chance_line,
    "Data: set_auc sheet of 05_classify_and_associate_sets.xlsx."
  ))
}

association_figure <- function(which_analysis, number, title, population) {
  hits <- set_association |>
    filter(analysis == which_analysis, p < 0.05) |>
    arrange(database, p)
  points <- pmap(hits, function(set_id, database, pathway, outcome, n, r, p, fdr, ...) {
    arms <- by_arm |>
      filter(set_id == !!set_id, outcome == !!outcome) |>
      mutate(text = sprintf("%-5s n=%d  r=%+.2f  p=%.3f", treatment, n, r, p))
    tibble(
      panel = paste0(
        set_label(database, pathway), "\nvs ", outcome,
        "   r = ", sprintf("%+.2f", r), "   p ", signif(p, 2), "   q ", signif(fdr, 2)
      ),
      d_score = delta_set[set_id, ], d_outcome = leg_phenotype[[outcome]],
      treatment = leg_phenotype$treatment, caption = paste(arms$text, collapse = "\n")
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel)))

  paginate(points, \(page) {
    ggplot(page, aes(d_score, d_outcome, colour = treatment, fill = treatment)) +
      geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey85") +
      geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey85") +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE, alpha = 0.12, linewidth = 0.5) +
      geom_point(size = 1.2, alpha = 0.9) +
      geom_text(
        data = distinct(page, panel, caption), inherit.aes = FALSE,
        aes(x = -Inf, y = Inf, label = caption), family = "mono",
        hjust = -0.05, vjust = 1.25, size = 1.9, lineheight = 1.2, colour = "grey25"
      ) +
      scale_y_continuous(expand = expansion(mult = c(0.06, 0.35))) +
      scale_colour_manual(
        values = c(BFR = "#B2182B", HLRT = "#2166AC"), aesthetics = c("colour", "fill"),
        name = NULL
      ) +
      labs(x = "change in set score, T2 - T1, one point per leg", y = "change in phenotype")
  }, number, title, paste(
    sprintf(
      "All %d set-outcome pairs whose Spearman correlation reaches nominal p in the %s",
      nrow(hits), which_analysis
    ),
    sprintf("analysis (%s), by collection then p.", population),
    "Lines are fitted within each arm; the header r is the", which_analysis, "correlation and",
    "the inset gives it within each arm. q from BH within collection and outcome.", chance_line,
    "Data: set_association and set_by_arm sheets of 05_classify_and_associate_sets.xlsx."
  ), scales = "free")
}

# The baseline control is reported in chance_expectation and the README but not drawn: it shows
# what the method returns when nothing is there, a number to read, not a panel to present.
drawn_tasks <- setdiff(names(tasks), "baseline_BFR_vs_HLRT")
figures <- c(
  imap(drawn_tasks, \(task_name, i) roc_figure(task_name, 12 + i)) |> list_flatten(),
  association_figure(
    "pooled", 17, "Training response against phenotype, all legs",
    sprintf("%d legs, both arms pooled", nrow(legs))
  ),
  association_figure(
    "differential", 18, "BFR minus HLRT against phenotype, within participant",
    sprintf("%d paired participants", nrow(delta_pairs))
  )
)

# Whether a collection clears chance is the headline, so it gets its own figure, not just a sheet.
chance_figure <- chance_expectation |>
  filter(analysis == "classification") |>
  mutate(comparison = factor(comparison, levels = map_chr(tasks, "label")))
figures <- c(figures, list(
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
      caption = supplement(19, "Nominal hits relative to chance, by collection", paste(
        "Paired Wilcoxon per set across", nrow(collection_sizes), "collections, uncorrected p.",
        "Bar length is observed nominal hits divided by the count that collection returns under",
        "the null; red clears 1, grey does not. Labels give observed of tested. Data:",
        "chance_expectation sheet of 05_classify_and_associate_sets.xlsx."
      ))
    ) +
    figure_theme +
    theme(panel.grid.major.y = element_blank())
))

pdf(
  file.path(stage, "b_reports", "05_classify_and_associate_sets_figures.pdf"),
  width = 8.27, height = 11.69
)
walk(figures, print)
invisible(dev.off())


# ---- one workbook -------------------------------------------------------------------------

sheets <- list(
  chance_expectation = chance_expectation,
  set_auc = arrange(set_auc, p_paired),
  set_association = arrange(set_association, p),
  set_by_arm = by_arm,
  set_catalog = catalog
)
overview <- tibble(
  sheet = names(sheets),
  rows = map_int(sheets, nrow),
  columns = map_int(sheets, ncol),
  description = c(
    "Nominal hits against chance, per collection. Read this first",
    "How well each set separates each task. AUC from ranks, p from the paired test",
    "Set against phenotype change: pooled, then within participant",
    "The same correlation computed inside BFR and inside HLRT, descriptive",
    "Every tested set with its collection and measured size"
  )
)
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "05_classify_and_associate_sets.xlsx")
)
message("wrote 05_classify_and_associate_sets.xlsx and a ", length(figures), "-page figure PDF")
sessionInfo()
