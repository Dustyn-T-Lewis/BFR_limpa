# BFR against HLRT NES per set: do both modalities move the same biology? No new test.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, patchwork, readxl, writexl)

stage <- here("03_Pathway_Enrichment", "03_enrich_scatter_fgsea")
set_tests <- read_excel(
  here("03_Pathway_Enrichment", "01_run_fgsea_and_fry", "c_data", "01_run_fgsea_and_fry.xlsx"),
  "set_tests"
)

x_contrast <- "BFR_Post-Pre"
y_contrast <- "HLRT_Post-Pre"
paired <- set_tests |>
  filter(method == "fgsea", contrast %in% c(x_contrast, y_contrast)) |>
  mutate(contrast = if_else(contrast == x_contrast, "x", "y")) |>
  pivot_wider(
    id_cols = c(set_id, database, pathway, n),
    names_from = contrast, values_from = c(NES, padj, main)
  ) |>
  mutate(
    significance = case_when(
      padj_x < 0.05 & padj_y < 0.05 ~ "Both",
      padj_x < 0.05 ~ x_contrast,
      padj_y < 0.05 ~ y_contrast,
      .default = "NS"
    ) |> factor(c("Both", x_contrast, y_contrast, "NS")),
    # Discordant: opposite sides of zero. Each discordant set here is significant in one arm
    # only, so the opposite sign rests on the other arm's null.
    discordant = sign(NES_x) != sign(NES_y),
    survivor = main_x | main_y,
    label = enrichVolcano::clean_label(pathway)
  )

concordance <- function(data) {
  test <- suppressWarnings(cor.test(data$NES_x, data$NES_y, method = "spearman"))
  hits <- filter(data, significance != "NS")
  tibble(
    sets = nrow(data), significant = nrow(hits), discordant = sum(hits$discordant),
    rho = round(unname(test$estimate), 3), p = test$p.value,
    same_sign = round(mean(!hits$discordant), 3)
  )
}

# plot_scatter() hides a set under collapse = TRUE unless collapsePathways kept it in a contrast.
enrichment <- set_tests |>
  filter(method == "fgsea", contrast %in% c(x_contrast, y_contrast)) |>
  transmute(
    contrast, database, pathway, NES, pval = p, padj, size = n, leadingEdge,
    dedup_status = if_else(main, "kept", "redundant")
  ) |>
  enrichVolcano::as_enrichment(enrichment_test = "fgsea")
scatter <- function(...) {
  enrichVolcano::plot_scatter(
    enrichment, x_contrast, y_contrast, ..., label_n = 10,
    theme = enrichVolcano::plot_theme(base_size = 9, base_family = "sans")
  )
}

supplement <- function(number, title, text, tags = "A") {
  plot_annotation(
    caption = str_wrap(sprintf("S%d Figure. %s. %s", number, title, text), 100),
    tag_levels = tags,
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
      plot.caption.position = "plot"
    )
  )
}
scatter_text <- paste(
  "Points are fgsea sets; grey reach neither threshold, coloured are significant (BH < 0.05",
  "within contrast) in one or both. Shaded quadrants move the same way in both contrasts, with",
  "their counts in the corners; dashed line is identity. The subtitle gives Spearman's rho",
  "across the plotted sets. Data: nes_scatter sheet of 03_enrich_scatter_fgsea.xlsx."
)

survivors <- filter(paired, survivor)
curated <- filter(paired, database %in% c("Hallmark", "GO_Slim"))
figures <- list(
  wrap_plots(
    scatter(collapse = FALSE) + labs(title = "All collections"),
    scatter(collapse = TRUE) + labs(title = "Collapse survivors"),
    ncol = 1
  ) +
    supplement(10, "NES concordance, all collections", paste(
      sprintf(
        "fgsea NES for each of %d sets in %s (x) against %s (y). (A) Every set. (B) Sets that",
        nrow(paired), x_contrast, y_contrast
      ),
      "survived collapsePathways in either contrast.", scatter_text
    )),
  scatter(databases = c("Hallmark", "GO_Slim"), collapse = FALSE) +
    labs(title = "Hallmark and GO Slim") +
    supplement(11, "NES concordance, Hallmark and GO Slim", paste(
      sprintf("All %d Hallmark and GO Slim sets, whose members do not nest.", nrow(curated)),
      scatter_text
    ), tags = NULL)
)

pdf(
  file.path(stage, "b_reports", "03_enrich_scatter_fgsea_figures.pdf"),
  width = 8.27, height = 11.69
)
walk(figures, print)
invisible(dev.off())

summary_table <- list(
  "all collections" = paired, "collapse survivors" = survivors, "Hallmark and GO Slim" = curated
) |>
  map(concordance) |>
  list_rbind(names_to = "population")

export <- paired |>
  transmute(
    set_id, database, pathway, label,
    genes = n,
    nes_x = round(NES_x, 3), nes_y = round(NES_y, 3),
    padj_x = signif(padj_x, 4), padj_y = signif(padj_y, 4),
    significance = as.character(significance), discordant, survivor
  ) |>
  arrange(significance, desc(abs(nes_x) + abs(nes_y)))
sheets <- list(
  concordance = summary_table,
  discordant = filter(export, discordant, significance != "NS"),
  nes_scatter = export
)
overview <- data.frame(
  sheet = names(sheets), rows = map_int(sheets, nrow), columns = map_int(sheets, ncol),
  description = c(
    "NES agreement between training contrasts",
    "Sets with opposite NES signs",
    "NES and adjusted p, both contrasts"
  )
)
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "03_enrich_scatter_fgsea.xlsx")
)
sessionInfo()
