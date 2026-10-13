# For pathways that both legs move the same way with training, do the same proteins drive them?
# Compares fgsea leading edges, the proteins ranked before the running-sum peak. Adapted from Sam's
# PR #26 (sc-norton), on the fgsea run 01_run_fgsea_and_fry already made. Leading edges shift with
# small changes in ranking, so the overlap describes and does not test.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, patchwork, readxl, writexl)

stage <- here("03_Pathway_Enrichment", "07_leading_edge")
set_tests <- read_excel(
  here("03_Pathway_Enrichment", "01_run_fgsea_and_fry", "c_data", "01_run_fgsea_and_fry.xlsx"),
  "set_tests"
)
set_catalog <- read_excel(
  here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data", "00_build_gene_sets.xlsx"),
  "set_catalog"
)
sets <- readRDS(here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data", "gene_sets.rds"))
arms <- c(BFR = "BFR_Post-Pre", HLRT = "HLRT_Post-Pre")
arm_colours <- c(BFR = "#B2182B", HLRT = "#2166AC")

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
  arrange(jaccard) |>
  mutate(label = enrichVolcano::clean_label(pathway, width = 1000))

# For a GO:BP pathway both legs move, which of its tested child terms hold each leg's leading-edge
# proteins. Counts only; child terms overlap, so no test across them.
go_ids <- set_catalog |>
  filter(database == "GOBP", qualifies) |>
  select(set_id, source_id)
parents <- inner_join(filter(leading_edge, database == "GOBP"), go_ids, by = "set_id")
go_children_of <- AnnotationDbi::mget(parents$source_id, GO.db::GOBPCHILDREN, ifnotfound = NA)
by_arm <- proteins |>
  mutate(arm = if_else(in_arm == "both", "BFR;HLRT", in_arm)) |>
  separate_longer_delim(arm, ";")
go_children <- parents |>
  transmute(parent_id = set_id, parent = label, child_go = map(source_id, \(id) {
    unname(go_children_of[[id]])
  })) |>
  unnest(child_go) |>
  inner_join(rename(go_ids, child_id = set_id, child_go = source_id), by = "child_go") |>
  mutate(child = enrichVolcano::clean_label(
    set_catalog$pathway[match(child_id, set_catalog$set_id)], width = 1000
  )) |>
  cross_join(tibble(arm = names(arms))) |>
  mutate(proteins = pmap_int(list(parent_id, child_id, arm), \(parent_id, child_id, arm) {
    sum(by_arm$gene[by_arm$set_id == parent_id & by_arm$arm == arm] %in% sets[[child_id]])
  }))

page_caption <- \(text) str_wrap(text, 115)
page_theme <- theme(
  plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
  plot.caption.position = "plot"
)

overlap_page <- leading_edge |>
  mutate(label = factor(str_trunc(label, 60), unique(str_trunc(label, 60)))) |>
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
    caption = page_caption(paste(
      "S20 Figure. Leading-edge overlap between the legs. Every pathway that fgsea calls at FDR",
      sprintf("< 0.05 in both training contrasts, in the same direction (%d,", nrow(leading_edge)),
      "each a main set after collapse in at least one leg). Bars split each pathway's",
      "leading-edge proteins into those shared by both legs and those in one leg only; lowest",
      "overlap at the bottom of each collection. Leading edges move with small ranking changes,",
      "so this describes and does not test. Data: leading_edge sheet of 07_leading_edge.xlsx."
    ))
  ) +
  theme_minimal(base_size = 9) +
  theme(
    axis.text.y = element_text(size = 6), strip.text.y = element_text(angle = 0),
    legend.position = "top", panel.grid.major.y = element_blank()
  ) +
  page_theme

# Sam's Venns: one per pathway, with the shared proteins named beneath, six to a page.
venn <- function(id) {
  counts <- filter(leading_edge, set_id == id)
  fit <- eulerr::euler(c(BFR = counts$BFR, HLRT = counts$HLRT, "BFR&HLRT" = counts$both))
  shared <- sort(filter(proteins, set_id == id, in_arm == "both")$gene)
  wrap_elements(full = plot(
    fit,
    fills = list(fill = arm_colours, alpha = 0.6), quantities = list(fontsize = 8),
    labels = list(fontsize = 8)
  )) +
    labs(
      title = str_wrap(counts$label, 45), subtitle = str_wrap(paste(shared, collapse = " "), 70)
    ) +
    theme(
      plot.title = element_text(size = 8, face = "bold"),
      plot.subtitle = element_text(size = 5, lineheight = 1.1)
    )
}
venn_pages <- split(leading_edge$set_id, ceiling(seq_along(leading_edge$set_id) / 6)) |>
  imap(\(ids, page) {
    wrap_plots(map(ids, venn), ncol = 2) +
      plot_annotation(caption = page_caption(sprintf(paste(
        "S21 Figure (page %s of %d). Leading edges by leg, one Venn per pathway in the order of",
        "S20. Areas follow the counts of proteins in each leg's leading edge only and in both;",
        "the proteins shared by both legs are listed under each title. Data: proteins sheet of",
        "07_leading_edge.xlsx."
      ), page, ceiling(nrow(leading_edge) / 6))), theme = page_theme)
  })

child_pages <- go_children |>
  mutate(
    proteins = if_else(arm == "BFR", -proteins, proteins),
    parent = str_wrap(parent, 45), child = str_trunc(child, 50)
  ) |>
  (\(rows) split(rows, ceiling(match(rows$parent, unique(rows$parent)) / 8)))() |>
  imap(\(rows, page) {
    ggplot(rows, aes(proteins, child, fill = arm)) +
      geom_col(width = 0.7) +
      geom_vline(xintercept = 0, linewidth = 0.3) +
      facet_wrap(~parent, ncol = 2, scales = "free_y") +
      scale_fill_manual(values = arm_colours, name = NULL) +
      scale_x_continuous(labels = abs) +
      labs(
        x = "parent's leading-edge proteins in the child term (BFR left, HLRT right)", y = NULL,
        caption = page_caption(sprintf(paste(
          "S22 Figure (page %s of %d). GO:BP pathways both legs move, broken into their tested",
          "child terms. Bars count the parent's leading-edge proteins in each child, BFR to the",
          "left and HLRT to the right. Child terms overlap, so the counts are not tested. Data:",
          "go_children sheet of 07_leading_edge.xlsx."
        ), page, ceiling(n_distinct(go_children$parent) / 8)))
      ) +
      theme_minimal(base_size = 8) +
      theme(legend.position = "top", strip.text = element_text(size = 6.5)) +
      page_theme
  })

pdf(file.path(stage, "b_reports", "07_leading_edge.pdf"), width = 8.27, height = 11.69)
walk(c(list(overlap_page), venn_pages, child_pages), print)
invisible(dev.off())

sheets <- list(leading_edge = leading_edge, proteins = proteins, go_children = go_children)
overview <- data.frame(
  sheet = names(sheets), rows = sapply(sheets, nrow),
  description = c(
    "One row per pathway: leading-edge proteins in both legs or one, Jaccard index, NES per leg.",
    "One row per pathway and leading-edge protein, with the legs it appears in.",
    "GO:BP pathways both legs move, their tested child terms, and each leg's count there."
  )
)
write_xlsx(c(list(overview = overview), sheets), file.path(stage, "c_data", "07_leading_edge.xlsx"))
sessionInfo()
