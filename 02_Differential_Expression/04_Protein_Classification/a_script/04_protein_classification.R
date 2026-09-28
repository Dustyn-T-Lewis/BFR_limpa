# Per protein and task: AUC and paired Wilcoxon p, mirroring 05_set_classification on the protein
# matrix. Nominal p is read against chance_expectation, BH beside it. baseline_BFR_vs_HLRT is the
# empirical floor.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, readxl, writexl)

stage <- here("02_Differential_Expression", "04_Protein_Classification")
proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))
de <- read_excel(
  here("02_Differential_Expression", "02_Differential", "c_data", "02_differential.xlsx"),
  "DEP_matrix"
)
correlations <- read_excel(
  here(
    "02_Differential_Expression", "03_Protein_Association", "c_data",
    "03_protein_association.xlsx"
  ),
  "correlations"
)
genes <- set_names(coalesce(proteins$genes$Genes, rownames(proteins)), rownames(proteins))

by_timepoint <- proteins$targets |>
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

delta <- proteins$E[, legs$T2] - proteins$E[, legs$T1]
colnames(delta) <- legs$leg_id

# Each task holds its matrix and two paired column sets. `favours` is the group an AUC above 0.5
# points to, so the figures can state direction.
arm <- split(legs, legs$treatment)
tasks <- list(
  pre_vs_post_BFR = list(
    values = proteins$E, positive = arm$BFR$T2, negative = arm$BFR$T1,
    label = "Pre to post, BFR", favours = "post", unit = "legs"
  ),
  pre_vs_post_HLRT = list(
    values = proteins$E, positive = arm$HLRT$T2, negative = arm$HLRT$T1,
    label = "Pre to post, HLRT", favours = "post", unit = "legs"
  ),
  baseline_BFR_vs_HLRT = list(
    values = proteins$E, positive = baseline$BFR, negative = baseline$HLRT,
    label = "Baseline BFR vs HLRT (control)", favours = "BFR", unit = "participants"
  ),
  post_BFR_vs_HLRT = list(
    values = proteins$E, positive = post$BFR, negative = post$HLRT,
    label = "Post BFR vs HLRT", favours = "BFR", unit = "participants"
  ),
  delta_BFR_vs_HLRT = list(
    values = delta, positive = delta_pairs$BFR, negative = delta_pairs$HLRT,
    label = "Change, BFR vs HLRT", favours = "BFR", unit = "participants"
  )
)

# pROC::roc() auto-orients: a true AUC of 0.194 comes back 0.806. direction = "<" pins it.
fit_roc <- function(values, spec) {
  pROC::roc(
    controls = values[spec$negative], cases = values[spec$positive],
    direction = "<", quiet = TRUE
  )
}
# The argument is `spec`, not `task`: tibble() builds a `task` column first, and a later
# `task$label` in the same call would read that column instead of the argument.
protein_auc <- imap(tasks, function(spec, name) {
  values <- spec$values
  tibble(
    task = name, task_label = spec$label, protein = rownames(values),
    gene = genes[rownames(values)], n_pairs = length(spec$positive),
    auc = apply(values, 1, \(row) as.numeric(pROC::auc(fit_roc(row, spec)))),
    p_paired = apply(values, 1, \(row) {
      suppressWarnings(wilcox.test(row[spec$positive], row[spec$negative], paired = TRUE)$p.value)
    })
  )
}) |>
  list_rbind() |>
  mutate(fdr = p.adjust(p_paired, "BH"), .by = task)

chance_expectation <- protein_auc |>
  summarise(
    features = n(), nominal = sum(p_paired < 0.05), expected = round(n() * 0.05, 1),
    ratio = round(nominal / expected, 2), fdr_sig = sum(fdr < 0.05),
    .by = c(task_label, n_pairs)
  )

# DE at FDR 0.05 in any contrast, association and classification at nominal p. Each protein gets
# eight association tests, so about a third of any list clears the association bar by chance.
de_hits <- de |>
  select(protein, ends_with("adj.P.Val")) |>
  pivot_longer(-protein, names_to = "contrast", values_to = "de_fdr") |>
  filter(de_fdr < 0.05) |>
  summarise(
    de_contrasts = paste(sub("_adj.P.Val", "", contrast), collapse = "; "),
    best_de_fdr = min(de_fdr), .by = protein
  )
associated <- correlations |>
  filter(p < 0.05) |>
  summarise(
    association = paste0(analysis, " ", outcome, " r ", sprintf("%+.2f", r), collapse = "; "),
    best_association_p = min(p), .by = protein
  )
classified <- protein_auc |>
  filter(p_paired < 0.05, task != "baseline_BFR_vs_HLRT") |>
  summarise(
    classification = paste0(task_label, " AUC ", sprintf("%.2f", auc), collapse = "; "),
    separates_arms = any(task %in% c("post_BFR_vs_HLRT", "delta_BFR_vs_HLRT")),
    .by = protein
  )
candidates <- de_hits |>
  inner_join(associated, by = "protein") |>
  inner_join(classified, by = "protein") |>
  mutate(gene = genes[protein], .after = protein) |>
  arrange(desc(separates_arms), best_association_p)

figure_theme <- theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(size = 6.8, lineheight = 1.2, margin = margin(3, 3, 4, 3)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(4, "mm"),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    plot.caption.position = "plot",
    legend.position = "top"
  )

# Every protein reaching nominal p gets a panel, twelve to a page, smallest p first.
roc_figure <- function(task_name) {
  spec <- tasks[[task_name]]
  hits <- protein_auc |>
    filter(task == task_name, p_paired < 0.05) |>
    arrange(p_paired)
  curves <- pmap(hits, function(protein, gene, auc, p_paired, fdr, ...) {
    coordinates <- pROC::coords(fit_roc(spec$values[protein, ], spec), "all")
    tibble(
      fpr = 1 - coordinates$specificity, tpr = coordinates$sensitivity,
      # An AUC below 0.5 separates the groups the other way, so it is coloured by direction.
      direction = if_else(auc >= 0.5, "higher", "lower"),
      panel = sprintf(
        "%s  %s\nAUC %.2f   p %s   q %s", gene, protein, auc, signif(p_paired, 2), signif(fdr, 2)
      )
    )
  }) |>
    list_rbind() |>
    mutate(panel = factor(panel, levels = unique(panel))) |>
    arrange(panel, fpr, tpr)
  pages <- split(curves, (as.integer(curves$panel) - 1) %/% 12)
  text <- paste(
    sprintf(
      "ROC curves for all %d proteins that separate the groups at nominal p, across %d paired %s,",
      nrow(hits), length(spec$positive), spec$unit
    ),
    sprintf(
      "smallest p first; chance alone returns %.0f of %s. AUC has its direction fixed: red",
      nrow(proteins) * 0.05, format(nrow(proteins), big.mark = ",")
    ),
    sprintf("separates higher in %s, blue lower.", spec$favours),
    "Shading is area under the curve, dashed line chance. p from the paired Wilcoxon",
    "signed-rank test, q from BH within task. Data: protein_auc sheet of",
    "04_protein_classification.xlsx."
  )
  imap(unname(pages), \(page, i) {
    ggplot(droplevels(page), aes(fpr, tpr, colour = direction, fill = direction)) +
      geom_abline(linetype = "22", linewidth = 0.35, colour = "grey60") +
      geom_ribbon(aes(ymin = 0, ymax = tpr), alpha = 0.16, colour = NA) +
      geom_step(linewidth = 0.7, direction = "hv") +
      facet_wrap(~panel, ncol = 3, nrow = 4) +
      scale_colour_manual(
        values = c(higher = "#B2182B", lower = "#2166AC"), aesthetics = c("colour", "fill"),
        guide = "none"
      ) +
      coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
      scale_x_continuous(breaks = c(0, 0.5, 1)) +
      scale_y_continuous(breaks = c(0, 0.5, 1)) +
      labs(
        x = "1 - specificity", y = "sensitivity",
        caption = str_wrap(
          sprintf("%s (page %d of %d). %s", spec$label, i, length(pages), text), 115
        )
      ) +
      figure_theme
  })
}

# One file per task. The baseline control is reported in chance_expectation but not drawn: it
# shows what the method returns when nothing is there.
walk(setdiff(names(tasks), "baseline_BFR_vs_HLRT"), \(task_name) {
  pdf(
    file.path(stage, "b_reports", paste0("04_classification_", task_name, ".pdf")),
    width = 8.27, height = 11.69
  )
  walk(roc_figure(task_name), print)
  invisible(dev.off())
})

sheets <- list(
  chance_expectation = chance_expectation,
  protein_auc = arrange(protein_auc, p_paired),
  candidates = candidates
)
overview <- data.frame(
  sheet = names(sheets), rows = map_int(sheets, nrow), columns = map_int(sheets, ncol),
  description = c(
    "Nominal hits against chance, per task. Read this first.",
    "How well each protein separates each task. AUC from ranks, p from paired test.",
    "DE at FDR 0.05, associated and classified at nominal p."
  )
)
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "04_protein_classification.xlsx")
)
sessionInfo()
