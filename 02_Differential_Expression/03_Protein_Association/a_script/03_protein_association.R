# Per protein: Spearman correlation of the protein change with each phenotype change, pooled over
# legs and within participant, mirroring 06_set_association on the protein matrix. BH within
# analysis and outcome.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, readr, readxl, writexl)

stage <- here("02_Differential_Expression", "03_Protein_Association")
proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))
targets <- proteins$targets
values <- proteins$E
genes <- set_names(coalesce(proteins$genes$Genes, rownames(proteins)), rownames(proteins))
phenotype <- read_csv(here("00_Input", "phenotype.csv"), show_col_types = FALSE)

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
    protein = rownames(values), n = sum(!is.na(outcome)),
    r = unname(fits["r", ]), p = unname(fits["p", ])
  )
}
results <- map(set_names(outcomes), \(name) {
  bind_rows(
    pooled = spearman_by_row(delta, leg_phenotype[[name]]),
    differential = spearman_by_row(paired, paired_outcome[[name]]),
    .id = "analysis"
  )
}) |>
  list_rbind(names_to = "outcome") |>
  mutate(fdr = p.adjust(p, "BH"), .by = c(analysis, outcome)) |>
  mutate(gene = proteins$genes$Genes[match(protein, rownames(proteins$genes))]) |>
  relocate(analysis, outcome, protein, gene)
by_arm <- map(set_names(c("BFR", "HLRT")), \(treatment) {
  keep <- leg_phenotype$treatment == treatment
  map(set_names(outcomes), \(name) {
    spearman_by_row(delta[, keep], leg_phenotype[[name]][keep]) |>
      mutate(outcome = name)
  }) |>
    list_rbind()
}) |>
  list_rbind(names_to = "treatment")

summary_table <- results |>
  summarise(
    proteins = n(), nominal = sum(p < 0.05), expected = round(n() * 0.05, 1),
    fdr_sig = sum(fdr < 0.05), best_fdr = signif(min(fdr), 3),
    .by = c(analysis, outcome, n)
  ) |>
  mutate(ratio = round(nominal / expected, 2)) |>
  arrange(analysis, outcome)

# Training hits are BH < 0.05 in either Post-Pre contrast. They vary more across legs, which gives
# them more power in a pooled correlation, so an excess here is a lead, not a finding.
training_hits <- read_excel(
  here("02_Differential_Expression", "02_Differential", "c_data", "02_differential.xlsx"),
  "DEP_matrix"
) |>
  filter(`BFR_Post-Pre_adj.P.Val` < 0.05 | `HLRT_Post-Pre_adj.P.Val` < 0.05) |>
  pull(protein)
training_overlap <- results |>
  mutate(training = protein %in% training_hits) |>
  summarise(
    training_hits = sum(training),
    nominal = sum(training & p < 0.05),
    expected = round(sum(training) * mean(p < 0.05), 1),
    fisher_p = signif(
      fisher.test(table(p < 0.05, training), alternative = "greater")$p.value, 2
    ),
    .by = c(analysis, outcome)
  ) |>
  arrange(analysis, outcome)

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
chance_text <- sprintf(
  "chance alone returns %.0f of %s proteins.", nrow(values) * 0.05,
  format(nrow(values), big.mark = ",")
)

association_figure <- function(which_analysis, outcome_name, title, population) {
  hits <- results |>
    filter(analysis == which_analysis, outcome == outcome_name, p < 0.05) |>
    arrange(p) |>
    mutate(gene = coalesce(gene, protein))
  points <- pmap(hits, function(protein, gene, outcome, r, p, fdr, ...) {
    arms <- by_arm |>
      filter(protein == !!protein, outcome == !!outcome) |>
      mutate(text = sprintf("%-5s n=%d  r=%+.2f  p=%.3f", treatment, n, r, p))
    tibble(
      panel = sprintf(
        "%s\nvs %s   r = %+.2f   p %s   q %s",
        paste(gene, protein, sep = "  "), outcome, r, signif(p, 2), signif(fdr, 2)
      ),
      change = delta[protein, ], d_outcome = leg_phenotype[[outcome]],
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
      labs(x = "change in protein, log2, T2 - T1, one point per leg", y = "change in phenotype")
  }, paste0(title, ", ", outcome_name), paste(
    sprintf(
      "All %d proteins whose Spearman correlation with %s reaches nominal p in the %s",
      nrow(hits), outcome_name, which_analysis
    ),
    sprintf("analysis (%s), smallest p first;", population), chance_text,
    "Lines are fitted within each arm; the header r is the", which_analysis, "correlation and",
    "the inset gives it within each arm. q from BH within analysis and outcome. Data:",
    "correlations and protein_by_arm sheets of 03_protein_association.xlsx."
  ), scales = "free")
}

# One file per outcome, so each can be read on its own.
walk(outcomes, \(outcome_name) {
  write_pdf(
    c(
      association_figure(
        "pooled", outcome_name, "Training response against phenotype, all legs",
        sprintf("%d legs, both arms pooled", nrow(legs))
      ),
      association_figure(
        "differential", outcome_name, "BFR minus HLRT, within participant",
        sprintf("%d paired participants", nrow(delta_pairs))
      )
    ),
    paste0("03_association_", outcome_name, ".pdf")
  )
})
write_pdf(list(chance_page(
  rename(summary_table, features = proteins), "outcome", "analysis",
  "Nominal hits relative to chance", paste(
    "Spearman per protein, uncorrected p;", chance_text, "Bar length is observed nominal hits",
    "divided by that count; red clears 1, grey does not. Labels give observed of tested. Data:",
    "summary sheet of 03_protein_association.xlsx."
  )
)), "03_chance.pdf")

sheets <- list(
  summary = summary_table,
  training_overlap = training_overlap,
  correlations = arrange(results, fdr),
  protein_by_arm = by_arm
)
overview <- data.frame(
  sheet = names(sheets), rows = map_int(sheets, nrow), columns = map_int(sheets, ncol),
  description = c(
    "Hits against chance per outcome",
    "Training hits at nominal p against chance",
    "Spearman r, p and FDR per protein",
    "The same correlation computed inside BFR and inside HLRT, descriptive"
  )
)
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "03_protein_association.xlsx")
)
sessionInfo()
