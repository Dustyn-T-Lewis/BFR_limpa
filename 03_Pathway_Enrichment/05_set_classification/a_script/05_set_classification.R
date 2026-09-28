# Per set and task: AUC and paired Wilcoxon p from the singscore matrix. Every comparison is paired.
# Nominal p is read against chance_expectation per collection, BH beside it. baseline_BFR_vs_HLRT
# is the empirical floor.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, readxl, writexl)

stage <- here("03_Pathway_Enrichment", "05_set_classification")
set_catalog <- read_excel(
  here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data", "00_build_gene_sets.xlsx"),
  "set_catalog"
)
set_score <- readRDS(here("03_Pathway_Enrichment", "04_run_singscore", "c_data", "set_scores.rds"))
targets <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))$targets

catalog <- set_catalog |>
  filter(qualifies) |>
  select(set_id, database, pathway, measured_size)
values <- set_score[catalog$set_id, ]
collection_sizes <- count(catalog, database)
message(nrow(catalog), " sets across ", nrow(collection_sizes), " collections")

by_timepoint <- targets |>
  mutate(leg_id = paste(participant, leg, sep = "_")) |>
  select(leg_id, participant, leg, treatment, timepoint, sample_id) |>
  pivot_wider(names_from = timepoint, values_from = sample_id)
legs <- drop_na(by_timepoint, T1, T2) |> arrange(participant, treatment)

# One row per participant holding both of their legs, for whichever column identifies them.
pair_by_treatment <- function(data, column) {
  data |>
    select(participant, treatment, value = all_of(column)) |>
    pivot_wider(names_from = treatment, values_from = value) |>
    drop_na(BFR, HLRT)
}
baseline <- pair_by_treatment(by_timepoint, "T1")
post <- pair_by_treatment(by_timepoint, "T2")
delta_pairs <- pair_by_treatment(legs, "leg_id")

delta <- values[, legs$T2] - values[, legs$T1]
colnames(delta) <- legs$leg_id

# Each task holds its matrix and two paired column sets. `favours` is the group an AUC above 0.5
# points to, so the figures can state direction.
arm <- split(legs, legs$treatment)
tasks <- list(
  pre_vs_post_BFR = list(
    values = values, positive = arm$BFR$T2, negative = arm$BFR$T1,
    label = "Pre to post, BFR", favours = "post", unit = "legs"
  ),
  pre_vs_post_HLRT = list(
    values = values, positive = arm$HLRT$T2, negative = arm$HLRT$T1,
    label = "Pre to post, HLRT", favours = "post", unit = "legs"
  ),
  baseline_BFR_vs_HLRT = list(
    values = values, positive = baseline$BFR, negative = baseline$HLRT,
    label = "Baseline BFR vs HLRT (control)", favours = "BFR", unit = "participants"
  ),
  post_BFR_vs_HLRT = list(
    values = values, positive = post$BFR, negative = post$HLRT,
    label = "Post BFR vs HLRT", favours = "BFR", unit = "participants"
  ),
  delta_BFR_vs_HLRT = list(
    values = delta, positive = delta_pairs$BFR, negative = delta_pairs$HLRT,
    label = "Change, BFR vs HLRT", favours = "BFR", unit = "participants"
  )
)
# The baseline control is reported in chance_expectation but not drawn: it shows what the method
# returns when nothing is there.
drawn_tasks <- setdiff(names(tasks), "baseline_BFR_vs_HLRT")

# pROC::roc() auto-orients: a true AUC of 0.194 comes back 0.806. direction = "<" pins it.
fit_roc <- function(values, spec) {
  pROC::roc(
    controls = values[spec$negative], cases = values[spec$positive],
    direction = "<", quiet = TRUE
  )
}
paired_test <- function(spec) {
  tibble(
    auc = apply(spec$values, 1, \(row) as.numeric(pROC::auc(fit_roc(row, spec)))),
    p_paired = apply(spec$values, 1, \(row) {
      suppressWarnings(wilcox.test(row[spec$positive], row[spec$negative], paired = TRUE)$p.value)
    })
  )
}
# The argument is `spec`, not `task`: tibble() builds a `task` column first, and a later
# `task$label` in the same call would read that column instead of the argument.
set_auc <- imap(tasks, function(spec, name) {
  tibble(
    task = name, task_label = spec$label, set_id = rownames(spec$values),
    n_pairs = length(spec$positive), paired_test(spec)
  )
}) |>
  list_rbind() |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p_paired, "BH"), .by = c(task, database))

# Each collection gets its own denominator, so a 41-set collection is not read against a
# 1,366-set one. BH is applied inside the same two-way split for the same reason.
chance_expectation <- set_auc |>
  transmute(
    analysis = "classification", comparison = task_label, database, n_pairs, p = p_paired, fdr
  ) |>
  summarise(
    features = n(), nominal = sum(p < 0.05), expected = round(n() * 0.05, 1),
    ratio = round(nominal / expected, 2), fdr_sig = sum(fdr < 0.05),
    .by = c(analysis, comparison, database, n_pairs)
  )

figure_theme <- theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 6.8, lineheight = 1.2, margin = margin(3, 3, 4, 3)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(4, "mm"),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    plot.caption.position = "plot",
    legend.position = "top"
  )
caption_text <- function(title, text, page = 1, n_pages = 1) {
  pages <- if (n_pages > 1) sprintf(" (page %d of %d)", page, n_pages) else ""
  str_wrap(sprintf("%s%s. %s", title, pages, text), 115)
}
# Twelve panels to a page. `draw` builds the plot from one page's rows: a paginating facet would
# build every panel for every page.
paginate <- function(data, draw, title, text, scales = "fixed") {
  pages <- split(data, (as.integer(data$panel) - 1) %/% 12)
  imap(unname(pages), \(page, i) {
    draw(droplevels(page)) +
      facet_wrap(~panel, ncol = 3, nrow = 4, scales = scales) +
      labs(caption = caption_text(title, text, i, length(pages))) +
      figure_theme
  })
}
chance_page <- function(data, y, facet, title, text) {
  ggplot(data, aes(ratio, .data[[y]], fill = ratio > 1)) +
    geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey40") +
    geom_col(width = 0.65) +
    geom_text(
      aes(label = sprintf("%d of %d", nominal, features)),
      hjust = -0.12, size = 2.5, colour = "grey25"
    ) +
    facet_wrap(vars(.data[[facet]]), ncol = 1) +
    scale_fill_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "grey72"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.22))) +
    labs(
      x = "observed nominal hits / chance expectation", y = NULL,
      caption = caption_text(title, text)
    ) +
    figure_theme +
    theme(panel.grid.major.y = element_blank())
}
write_pdf <- function(pages, name) {
  pdf(file.path(stage, "b_reports", name), width = 8.27, height = 11.69)
  walk(pages, print)
  invisible(dev.off())
}
chance_line <- function(db) {
  n <- collection_sizes$n[collection_sizes$database == db]
  sprintf("Nominal p, uncorrected: chance alone returns %.1f of the %d %s sets.", n * 0.05, n, db)
}
set_label <- function(database, pathway) {
  str_wrap(paste0(database, ": ", enrichVolcano::clean_label(pathway, width = 1000)), 40)
}

roc_figure <- function(task_name, number, db) {
  spec <- tasks[[task_name]]
  hits <- set_auc |>
    filter(task == task_name, database == db, p_paired < 0.05) |>
    arrange(p_paired)
  if (nrow(hits) == 0) {
    return(list())
  }
  curves <- pmap(hits, function(set_id, database, pathway, auc, p_paired, fdr, ...) {
    coordinates <- pROC::coords(fit_roc(spec$values[set_id, ], spec), "all")
    tibble(
      fpr = 1 - coordinates$specificity, tpr = coordinates$sensitivity,
      # An AUC below 0.5 separates the groups the other way, so it is coloured by direction.
      direction = if_else(auc >= 0.5, "higher", "lower"),
      panel = sprintf(
        "%s\nAUC %.2f   p %s   q %s",
        set_label(database, pathway), auc, signif(p_paired, 2), signif(fdr, 2)
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
  }, sprintf("S%d Figure", number), paste(
    paste0(spec$label, ", ", db, "."),
    sprintf(
      "ROC curves for all %d %s sets whose singscore separates the groups at nominal p, across %d",
      nrow(hits), db, length(spec$positive)
    ),
    sprintf("paired %s, smallest p first. AUC has its direction fixed: red", spec$unit),
    sprintf("separates higher in %s, blue lower.", spec$favours),
    "Shading is area under the curve, dashed line chance. p from the paired Wilcoxon",
    "signed-rank test, q from BH within collection and task.", chance_line(db),
    "Data: set_auc sheet of 05_set_classification.xlsx."
  ))
}

# One file per collection, so each can be read on its own.
walk(unique(catalog$database), \(db) {
  write_pdf(
    list_flatten(imap(drawn_tasks, \(task_name, i) roc_figure(task_name, 11 + i, db))),
    paste0("05_classification_", db, ".pdf")
  )
})
write_pdf(list(chance_page(
  mutate(chance_expectation, comparison = factor(comparison, map_chr(tasks, "label"))),
  "database", "comparison", "S16 Figure", paste(
    "Nominal hits relative to chance, by collection. Paired Wilcoxon per set across",
    nrow(collection_sizes), "collections, uncorrected p. Bar length is observed nominal hits",
    "divided by the count that collection returns under the null; red clears 1, grey does not.",
    "Labels give observed of tested. Data: chance_expectation sheet of 05_set_classification.xlsx."
  )
)), "05_chance_by_collection.pdf")

sheets <- list(
  chance_expectation = chance_expectation,
  set_auc = arrange(set_auc, p_paired),
  set_catalog = catalog
)
overview <- data.frame(
  sheet = names(sheets), rows = map_int(sheets, nrow), columns = map_int(sheets, ncol),
  description = c(
    "Nominal hits against chance, per collection. Read this first.",
    "How well each set separates each task. AUC from ranks, p from paired test.",
    "Every tested set with its collection and measured size."
  )
)
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "05_set_classification.xlsx")
)
sessionInfo()
