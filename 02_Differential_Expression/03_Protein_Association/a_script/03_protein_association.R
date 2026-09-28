# Each protein's pre-to-post change against the same leg's change in muscle size and quality.
# Pooled uses all legs; differential uses participants with both legs and correlates the BFR minus
# HLRT differences, so anything constant within a participant cancels. BH within analysis and
# outcome.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, readr, readxl, writexl)

stage <- here("02_Differential_Expression", "03_Protein_Association")
proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))
phenotype <- read_csv(here("00_Input", "phenotype.csv"), show_col_types = FALSE)

# A leg needs both timepoints for a change score; S06's left leg has T1 only.
legs <- proteins$targets |>
  mutate(leg_id = paste(participant, leg, sep = "_")) |>
  select(leg_id, participant, leg, treatment, timepoint, sample_id) |>
  pivot_wider(names_from = timepoint, values_from = sample_id) |>
  drop_na(T1, T2) |>
  arrange(participant, treatment)
leg_pairs <- legs |>
  select(participant, treatment, leg_id) |>
  pivot_wider(names_from = treatment, values_from = leg_id) |>
  drop_na(BFR, HLRT)

delta <- proteins$E[, legs$T2] - proteins$E[, legs$T1]
colnames(delta) <- legs$leg_id
paired <- delta[, leg_pairs$BFR] - delta[, leg_pairs$HLRT]
colnames(paired) <- leg_pairs$participant

# Echo intensity falls as muscle composition improves.
leg_phenotype <- legs |>
  left_join(select(phenotype, -treatment), by = c("participant", "leg")) |>
  mutate(
    vl_csa = vl_csa_post_cm2 - vl_csa_pre_cm2,
    vl_echo = vl_echo_post_au - vl_echo_pre_au,
    rf_csa = rf_csa_post_cm2 - rf_csa_pre_cm2,
    rf_echo = rf_echo_post_au - rf_echo_pre_au
  )
outcomes <- c("vl_csa", "vl_echo", "rf_csa", "rf_echo")
paired_outcome <- map(set_names(outcomes), \(name) {
  value <- set_names(leg_phenotype[[name]], leg_phenotype$leg_id)
  value[leg_pairs$BFR] - value[leg_pairs$HLRT]
})

# Ties make cor.test fall back to its approximation and warn on every one of 24,336 calls.
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

chance_page <- ggplot(summary_table, aes(ratio, outcome, fill = analysis)) +
  geom_vline(xintercept = 1, linewidth = 0.4, colour = "grey40") +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  scale_fill_manual(values = c(pooled = "#B2182B", differential = "#2166AC"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = "nominal hits / chance expectation", y = NULL,
    caption = str_wrap(paste(
      "Nominal hits relative to chance, per outcome.", format(nrow(proteins), big.mark = ","),
      "proteins give about",
      round(nrow(proteins) * 0.05), "nominal hits by chance, so a ratio near 1 is what the tests",
      "return when nothing is there. Data: summary sheet of 03_protein_association.xlsx."
    ), 110)
  ) +
  theme_minimal(base_size = 10) +
  theme(
    panel.grid.major.y = element_blank(), legend.position = "top",
    plot.caption = element_text(hjust = 0, size = 8.5), plot.caption.position = "plot"
  )

# Every pair at nominal p, twelve to a page, one block per analysis and outcome, smallest p first.
blocks <- results |>
  filter(p < 0.05) |>
  arrange(analysis, outcome, p) |>
  mutate(panel = sprintf(
    "%s  %s\n%s   r = %+.2f   p %s   q %s",
    coalesce(gene, protein), protein, outcome, r, signif(p, 2), signif(fdr, 2)
  )) |>
  split(~ analysis + outcome, drop = TRUE, sep = ", ")
blocks <- blocks[paste(rep(c("pooled", "differential"), each = 4), outcomes, sep = ", ")]
points_for <- function(hit) {
  pooled <- hit$analysis == "pooled"
  tibble(
    panel = hit$panel,
    change_protein = if (pooled) delta[hit$protein, ] else paired[hit$protein, ],
    change_outcome = if (pooled) leg_phenotype[[hit$outcome]] else paired_outcome[[hit$outcome]],
    arm = if (pooled) leg_phenotype$treatment else "BFR minus HLRT"
  )
}
pages <- imap(blocks, \(block, name) {
  points <- map(seq_len(nrow(block)), \(i) points_for(block[i, ])) |>
    list_rbind() |>
    mutate(panel = factor(panel, unique(panel)))
  chunks <- split(points, (as.integer(points$panel) - 1) %/% 12)
  imap(unname(chunks), \(page, i) {
    ggplot(droplevels(page), aes(change_protein, change_outcome, colour = arm)) +
      geom_point(size = 1.1, alpha = 0.85) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
      facet_wrap(~panel, ncol = 3, nrow = 4, scales = "free") +
      scale_colour_manual(
        values = c(BFR = "#B2182B", HLRT = "#2166AC", `BFR minus HLRT` = "grey30"), name = NULL
      ) +
      labs(
        x = "change in protein, log2", y = "change in outcome",
        caption = str_wrap(sprintf(
          paste(
            "%s (page %d of %d). %d proteins reach nominal p; chance alone gives about %d of %s.",
            "Spearman r and p; q is BH within analysis and outcome. Data: correlations sheet of",
            "03_protein_association.xlsx."
          ),
          name, i, length(chunks), nrow(block), round(nrow(proteins) * 0.05),
          format(nrow(proteins), big.mark = ",")
        ), 110)
      ) +
      theme_minimal(base_size = 9) +
      theme(
        strip.text = element_text(size = 7), legend.position = "top",
        plot.caption = element_text(hjust = 0, size = 8.5), plot.caption.position = "plot"
      )
  })
}) |>
  list_flatten()
pdf(
  file.path(stage, "b_reports", "03_protein_association_figures.pdf"),
  width = 8.27, height = 11.69
)
walk(c(list(chance_page), pages), print)
invisible(dev.off())

sheets <- list(
  summary = summary_table,
  training_overlap = training_overlap,
  correlations = arrange(results, fdr)
)
overview <- data.frame(
  sheet = names(sheets), rows = map_int(sheets, nrow), columns = map_int(sheets, ncol),
  description = c(
    "Hits against chance per outcome",
    "Training hits at nominal p against chance",
    "Spearman r, p and FDR per protein"
  )
)
write_xlsx(
  c(list(overview = overview), sheets), file.path(stage, "c_data", "03_protein_association.xlsx")
)
sessionInfo()
