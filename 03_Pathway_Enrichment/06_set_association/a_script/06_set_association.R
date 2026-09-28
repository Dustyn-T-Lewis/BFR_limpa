# Per set: Spearman correlation of the singscore change with each phenotype change, pooled over
# legs and within participant. Nominal p is read against chance_expectation per collection, BH
# beside it.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, readr, readxl, writexl)

stage <- here("03_Pathway_Enrichment", "06_set_association")
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
delta_pairs <- pair_by_treatment(legs, "leg_id")

delta <- values[, legs$T2] - values[, legs$T1]
colnames(delta) <- legs$leg_id

outcomes <- c("vl_csa", "vl_echo", "rf_csa", "rf_echo")
# A duplicated phenotype row would silently misalign every leg below. Echo intensity falls as
# muscle composition improves.
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
paired <- delta[, delta_pairs$BFR] - delta[, delta_pairs$HLRT]
colnames(paired) <- delta_pairs$participant
paired_outcome <- map(set_names(outcomes), \(name) {
  value <- set_names(leg_phenotype[[name]], leg_phenotype$leg_id)
  value[delta_pairs$BFR] - value[delta_pairs$HLRT]
})

# cor.test gives the exact Spearman p at these sample sizes; the t approximation is off by up to
# 9e-4. Ties make it warn and fall back to the approximation, so the warnings are suppressed.
spearman_by_row <- function(values, outcome) {
  fits <- suppressWarnings(apply(values, 1, \(row) {
    test <- cor.test(row, outcome, method = "spearman")
    c(r = unname(test$estimate), p = test$p.value)
  }))
  tibble(
    set_id = rownames(values), n = sum(!is.na(outcome)),
    r = unname(fits["r", ]), p = unname(fits["p", ])
  )
}
set_association <- map(set_names(outcomes), \(name) {
  bind_rows(
    pooled = spearman_by_row(delta, leg_phenotype[[name]]),
    differential = spearman_by_row(paired, paired_outcome[[name]]),
    .id = "analysis"
  ) |>
    mutate(outcome = name)
}) |>
  list_rbind() |>
  relocate(analysis, outcome) |>
  left_join(select(catalog, set_id, database, pathway), by = "set_id") |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(analysis, outcome, database))
by_arm <- map(set_names(c("BFR", "HLRT")), \(treatment) {
  keep <- leg_phenotype$treatment == treatment
  map(set_names(outcomes), \(name) {
    spearman_by_row(delta[, keep], leg_phenotype[[name]][keep]) |>
      mutate(outcome = name)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "treatment")

# Each collection gets its own denominator, so a 41-set collection is not read against a
# 1,366-set one. BH is applied inside the same two-way split for the same reason.
chance_expectation <- set_association |>
  transmute(
    analysis = paste0("association: ", analysis), comparison = outcome, database,
    n_pairs = n, p, fdr
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

association_figure <- function(which_analysis, number, title, population, db) {
  hits <- set_association |>
    filter(analysis == which_analysis, database == db, p < 0.05) |>
    arrange(p)
  if (nrow(hits) == 0) {
    return(list())
  }
  points <- pmap(hits, function(set_id, database, pathway, outcome, r, p, fdr, ...) {
    arms <- by_arm |>
      filter(set_id == !!set_id, outcome == !!outcome) |>
      mutate(text = sprintf("%-5s n=%d  r=%+.2f  p=%.3f", treatment, n, r, p))
    tibble(
      panel = sprintf(
        "%s\nvs %s   r = %+.2f   p %s   q %s",
        set_label(database, pathway), outcome, r, signif(p, 2), signif(fdr, 2)
      ),
      change = delta[set_id, ], d_outcome = leg_phenotype[[outcome]],
      treatment = leg_phenotype$treatment, caption = paste(arms$text, collapse = "\n")
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel)))
  paginate(points, \(page) {
    ggplot(page, aes(change, d_outcome, colour = treatment, fill = treatment)) +
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
  }, sprintf("S%d Figure", number), paste(
    paste0(title, ", ", db, "."),
    sprintf(
      "All %d %s set-outcome pairs whose Spearman correlation reaches nominal p in the %s",
      nrow(hits), db, which_analysis
    ),
    sprintf("analysis (%s), smallest p first.", population),
    "Lines are fitted within each arm; the header r is the", which_analysis, "correlation and",
    "the inset gives it within each arm. q from BH within collection and outcome.", chance_line(db),
    "Data: set_association and set_by_arm sheets of 06_set_association.xlsx."
  ), scales = "free")
}

# One file per collection, so each can be read on its own.
walk(unique(catalog$database), \(db) {
  write_pdf(
    c(
      association_figure(
        "pooled", 17, "Training response against phenotype, all legs",
        sprintf("%d legs, both arms pooled", nrow(legs)), db
      ),
      association_figure(
        "differential", 18, "BFR minus HLRT, within participant",
        sprintf("%d paired participants", nrow(delta_pairs)), db
      )
    ),
    paste0("06_association_", db, ".pdf")
  )
})
write_pdf(list(chance_page(
  mutate(chance_expectation, panel = paste0(sub("association: ", "", analysis), ", ", comparison)),
  "database", "panel", "S19 Figure", paste(
    "Nominal hits relative to chance, by collection. Spearman per set, uncorrected p. Bar length",
    "is observed nominal hits divided by the count that collection returns under the null; red",
    "clears 1, grey does not. Labels give observed of tested. Data: chance_expectation sheet of",
    "06_set_association.xlsx."
  )
)), "06_chance_by_collection.pdf")

sheets <- list(
  chance_expectation = chance_expectation,
  set_association = arrange(set_association, p),
  set_by_arm = by_arm,
  set_catalog = catalog
)
overview <- data.frame(
  sheet = names(sheets), rows = map_int(sheets, nrow), columns = map_int(sheets, ncol),
  description = c(
    "Nominal hits against chance, per collection. Read this first.",
    "Set against phenotype change: pooled, then within participant.",
    "The same correlation computed inside BFR and inside HLRT, descriptive.",
    "Every tested set with its collection and measured size."
  )
)
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "06_set_association.xlsx")
)
sessionInfo()
