# For pathways that both legs move the same way with training, do the same proteins drive them?
# Compares fgsea leading edges, the proteins ranked before the running-sum peak. Adapted from Sam's
# PR #26 (sc-norton), on the fgsea run 01_run_fgsea_and_fry already made. Leading edges shift with
# small changes in ranking, so the overlap describes and does not test.

pacman::p_load(here, dplyr, tidyr, stringr, ggplot2, readxl, writexl)

stage <- here("03_Pathway_Enrichment", "07_leading_edge")
set_tests <- read_excel(
  here("03_Pathway_Enrichment", "01_run_fgsea_and_fry", "c_data", "01_run_fgsea_and_fry.xlsx"),
  "set_tests"
)
arms <- c(BFR = "BFR_Post-Pre", HLRT = "HLRT_Post-Pre")

# Significant in both, same direction, and a main set after collapsePathways in at least one leg,
# so one representative stands for each cluster of overlapping sets.
fgsea <- set_tests |>
  filter(method == "fgsea", contrast %in% arms) |>
  mutate(arm = names(arms)[match(contrast, arms)]) |>
  select(set_id, database, pathway, arm, NES, padj, main, leadingEdge)
shared_sets <- fgsea |>
  summarise(
    keep = all(padj < 0.05) && n_distinct(sign(NES)) == 1 && any(main), .by = set_id
  ) |>
  filter(keep) |>
  pull(set_id)

proteins <- fgsea |>
  filter(set_id %in% shared_sets) |>
  separate_longer_delim(leadingEdge, ";") |>
  rename(gene = leadingEdge) |>
  summarise(
    in_arm = paste(sort(arm), collapse = " and "), .by = c(set_id, database, pathway, gene)
  ) |>
  mutate(in_arm = recode(in_arm, "BFR and HLRT" = "both"))
leading_edge <- proteins |>
  count(set_id, database, pathway, in_arm) |>
  pivot_wider(names_from = in_arm, values_from = n, values_fill = 0) |>
  mutate(jaccard = both / (both + BFR + HLRT)) |>
  left_join(
    fgsea |>
      select(set_id, arm, NES) |>
      pivot_wider(names_from = arm, values_from = NES, names_prefix = "NES_"),
    by = "set_id"
  ) |>
  arrange(jaccard)

overlap_page <- leading_edge |>
  mutate(label = str_trunc(enrichVolcano::clean_label(pathway, width = 1000), 60)) |>
  mutate(label = factor(label, unique(label))) |>
  pivot_longer(c(both, BFR, HLRT), names_to = "in_arm", values_to = "proteins") |>
  mutate(in_arm = factor(in_arm, c("BFR", "both", "HLRT"))) |>
  ggplot(aes(proteins, label, fill = in_arm)) +
  geom_col(width = 0.75) +
  facet_grid(database ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(
    values = c(BFR = "#B2182B", both = "grey55", HLRT = "#2166AC"), name = "leading edge in"
  ) +
  labs(
    x = "leading-edge proteins", y = NULL,
    caption = str_wrap(paste(
      "S20 Figure. Leading-edge overlap between the legs. Every pathway that fgsea calls at FDR",
      sprintf("< 0.05 in both training contrasts, in the same direction (%d,", nrow(leading_edge)),
      "each a main set after collapse in at least one leg). Bars split each pathway's",
      "leading-edge proteins into those shared by both legs and those in one leg only; lowest",
      "overlap at the bottom of each collection. Leading edges move with small ranking changes,",
      "so this describes and does not test. Data: leading_edge sheet of 07_leading_edge.xlsx."
    ), 115)
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 6), strip.text.y = element_text(angle = 0),
    legend.position = "top", panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
    plot.caption.position = "plot"
  )

pdf(file.path(stage, "b_reports", "07_leading_edge.pdf"), width = 8.27, height = 11.69)
print(overlap_page)
invisible(dev.off())

sheets <- list(leading_edge = leading_edge, proteins = proteins)
overview <- data.frame(
  sheet = names(sheets), rows = sapply(sheets, nrow),
  description = c(
    "One row per pathway: leading-edge proteins in both legs or one, Jaccard index, NES per leg.",
    "One row per pathway and leading-edge protein, with the legs it appears in."
  )
)
write_xlsx(c(list(overview = overview), sheets), file.path(stage, "c_data", "07_leading_edge.xlsx"))
sessionInfo()
