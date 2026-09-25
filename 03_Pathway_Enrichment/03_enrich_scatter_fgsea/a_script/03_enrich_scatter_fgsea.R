# BFR against HLRT NES per set: do both modalities move the same biology? No new test.

pacman::p_load(
  here, dplyr, tidyr, tibble, purrr, stringr, ggplot2, ggrepel, patchwork, readxl, writexl
)

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

set_colours <- set_names(
  c("#6A3D9A", "#D7301F", "#2B6CB0"), c("Both", x_contrast, y_contrast)
)

# Builds all six panels. `labelled` is the subset that gets names, so a dense cloud and a
# zoomed handful differ only in what is passed in.
nes_panel <- function(data, title, labelled = data[0, ], pad = 0.12) {
  span <- range(c(data$NES_x, data$NES_y))
  limit <- span + c(-1, 1) * diff(span) * pad
  shown <- filter(data, significance != "NS")
  stats <- concordance(data)
  # A rank correlation over a handful of points is noise, so rho appears only from 30 sets up.
  rho <- if (stats$sets >= 30) sprintf(" | rho %.2f", stats$rho) else ""
  caption <- sprintf(
    "%d sets | %d significant%s | %d discordant",
    stats$sets, stats$significant, rho, stats$discordant
  )
  ggplot(data, aes(NES_x, NES_y)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_abline(slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_point(
      data = filter(data, significance == "NS"),
      colour = "grey82", size = 0.45, alpha = 0.3
    ) +
    geom_point(aes(colour = significance, size = n), data = shown, alpha = 0.85) +
    geom_text_repel(
      data = labelled, aes(label = label), size = 2.3, colour = "grey15",
      segment.colour = "grey60", segment.size = 0.25, min.segment.length = 0,
      max.overlaps = Inf, seed = 1, force = 6, lineheight = 0.85
    ) +
    scale_colour_manual(values = set_colours, name = "significant in", drop = FALSE) +
    scale_size_continuous(
      range = c(1, 4.5), name = "genes",
      limits = range(paired$n), breaks = c(50, 150, 300)
    ) +
    coord_fixed(xlim = limit, ylim = limit) +
    labs(
      x = paste("NES,", x_contrast), y = paste("NES,", y_contrast),
      title = title, subtitle = caption
    ) +
    theme_minimal(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      plot.subtitle = element_text(size = 7.5, colour = "grey30")
    )
}

supplement <- function(number, title, text) {
  plot_annotation(
    caption = str_wrap(sprintf("S%d Figure. %s. %s", number, title, text), 115),
    tag_levels = "A",
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
      plot.caption.position = "plot"
    )
  )
}

# Figure one: all collections, collapse survivors, and the discordant sets on their own axes.
discordant_sets <- filter(paired, discordant, significance != "NS")
survivors <- filter(paired, survivor)
figure_all <- wrap_plots(
  nes_panel(paired, "All collections"),
  nes_panel(survivors, "Collapse survivors"),
  # Too few points for a size key, and patchwork will not merge guide sets that differ.
  nes_panel(discordant_sets, "Discordant", labelled = discordant_sets, pad = 0.35) +
    guides(colour = "none", size = "none"),
  ncol = 2, guides = "collect"
) +
  supplement(10, "NES agreement, BFR against HLRT", paste(
    sprintf(
      "fgsea NES for each of %d sets in %s (x) against %s (y), BH within contrast.",
      nrow(paired), x_contrast, y_contrast
    ),
    "(A) Every set. (B) Sets that survived collapsePathways in either contrast.",
    sprintf(
      "(C) The %d discordant sets, on opposite sides of zero and significant in one arm only,",
      nrow(discordant_sets)
    ),
    "rescaled so each can be named. Dashed line is identity; grey points reach neither",
    "threshold; point size is gene count. Data: nes_scatter sheet of 03_enrich_scatter_fgsea.xlsx."
  )) &
  theme(legend.position = "bottom")

# Figure two: the two collections whose members do not nest, then each concordant quadrant
# scaled to its own points so every set can be named.
curated <- filter(paired, database %in% c("Hallmark", "GO_Slim"))
quadrant <- function(direction) {
  rows <- filter(curated, significance != "NS", (NES_x > 0) == direction)
  nes_panel(rows, if (direction) "Up in both" else "Down in both", labelled = rows, pad = 0.22)
}
figure_curated <- (
  nes_panel(
    curated, "Hallmark and GO Slim",
    labelled = slice_min(filter(curated, significance != "NS"), padj_x + padj_y, n = 8)
  ) / (quadrant(TRUE) | quadrant(FALSE))
) +
  plot_layout(guides = "collect", heights = c(1, 1)) +
  supplement(11, "NES agreement, Hallmark and GO Slim", paste(
    sprintf("(A) All %d Hallmark and GO Slim sets, whose members do not nest.", nrow(curated)),
    "(B, C) The significant ones rescaled by direction so each can be named. Axes, colours",
    "and point size as in S10 Figure. Data: nes_scatter sheet of 03_enrich_scatter_fgsea.xlsx."
  )) &
  theme(legend.position = "bottom")

pdf(
  file.path(stage, "b_reports", "03_enrich_scatter_fgsea_figures.pdf"),
  width = 8.27, height = 11.69
)
print(figure_all)
print(figure_curated)
invisible(dev.off())

summary_table <- bind_rows(
  mutate(concordance(paired), population = "all collections"),
  mutate(concordance(survivors), population = "collapse survivors"),
  mutate(concordance(curated), population = "Hallmark and GO Slim")
) |>
  relocate(population)
print(as.data.frame(summary_table))

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
overview <- tibble(
  sheet = names(sheets),
  rows = map_int(sheets, nrow),
  columns = map_int(sheets, ncol),
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
message("wrote 03_enrich_scatter_fgsea.xlsx and a 2-page figure PDF")
sessionInfo()
