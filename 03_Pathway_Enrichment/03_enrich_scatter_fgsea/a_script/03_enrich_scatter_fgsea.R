# NES in the two training contrasts on one pair of axes: fgsea scores each set once per contrast,
# so the scatter asks whether both modalities move the same biology. No new test; reshapes the
# 01_run_fgsea_and_fry output.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
  library(patchwork)
})

stage <- here("03_Pathway_Enrichment", "03_enrich_scatter_fgsea")
figure_dir <- file.path(stage, "b_reports")
out <- file.path(stage, "c_data")
walk(c(figure_dir, out), dir.create, recursive = TRUE, showWarnings = FALSE)

inputs <- c(set_tests = "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) stop("Run 01_run_fgsea_and_fry first.")
fg <- readRDS(paths[["set_tests"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

x_contrast <- "BFR_Post-Pre"
y_contrast <- "HLRT_Post-Pre"
paired <- fg$set_tests |>
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
    label = enrichVolcano::ev_clean_label(pathway)
  )
stopifnot(nrow(paired) > 0, !anyNA(paired$NES_x), !anyNA(paired$NES_y))

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
  caption <- sprintf(
    "%d sets | %d significant | %d discordant",
    stats$sets, stats$significant, stats$discordant
  )
  if (stats$sets >= 30) {
    caption <- sprintf(
      "%d sets | %d significant | rho %.2f | %d discordant",
      stats$sets, stats$significant, stats$rho, stats$discordant
    )
  }
  ggplot(data, aes(NES_x, NES_y)) +
    geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
    geom_abline(slope = 1, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
    geom_point(
      data = filter(data, significance == "NS"),
      colour = "grey82", size = 0.45, alpha = 0.3
    ) +
    geom_point(aes(colour = significance, size = n), data = shown, alpha = 0.85) +
    ggrepel::geom_text_repel(
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

save_composite <- function(figure, name, width, height) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)), figure,
      width = width, height = height, dpi = 300, bg = "white"
    )
  })
  message("wrote ", name)
}

page_labels <- function(title, subtitle, caption) {
  plot_annotation(
    title = title, subtitle = subtitle, caption = caption, tag_levels = "A",
    theme = theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 8.5, colour = "grey30"),
      plot.caption = element_text(size = 7.5, colour = "grey45", hjust = 0)
    )
  )
}

# Composite one: all collections, collapse survivors, and the discordant sets on their own axes.
discordant_sets <- filter(paired, discordant, significance != "NS")
survivors <- filter(paired, survivor)
composite_all <- wrap_plots(
  nes_panel(paired, "All collections"),
  nes_panel(survivors, "Collapse survivors"),
  # Too few points for a size key, and patchwork will not merge guide sets that differ.
  nes_panel(discordant_sets, "Discordant", labelled = discordant_sets, pad = 0.35) +
    guides(colour = "none", size = "none"),
  nrow = 1
) +
  page_labels(
    "NES concordance, all collections",
    sprintf("fgsea NES per contrast, %d sets, BH within contrast", nrow(paired)),
    paste(
      "One point per gene set, scored in both training contrasts. Dashed line is identity;",
      "grey points reach neither threshold. Discordant sets sit on opposite sides of zero and",
      "are significant in one arm only. Panel C rescales those eight.",
      "Table: c_data/nes_scatter.csv."
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_composite(composite_all, "01_nes_concordance_all", 13, 6)

# Composite two: the two collections whose members do not nest, then each concordant quadrant
# scaled to its own points so every set can be named.
curated <- filter(paired, database %in% c("Hallmark", "GO_Slim"))
quadrant <- function(direction) {
  rows <- filter(curated, significance != "NS", (NES_x > 0) == direction)
  nes_panel(rows, if (direction) "Up in both" else "Down in both", labelled = rows, pad = 0.22)
}
composite_curated <- (
  nes_panel(
    curated, "Hallmark and GO Slim",
    labelled = slice_min(filter(curated, significance != "NS"), padj_x + padj_y, n = 8)
  ) | (quadrant(TRUE) / quadrant(FALSE))
) +
  page_labels(
    "NES concordance, Hallmark and GO Slim",
    sprintf("fgsea NES per contrast, %d non-nesting sets, BH within contrast", nrow(curated)),
    paste(
      "Panel A is every Hallmark and GO Slim set; B and C rescale the significant ones by",
      "direction so each can be named. Point size is gene count, colour the contrast a set",
      "reached FDR 0.05 in. Table: c_data/nes_scatter.csv."
    )
  ) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_composite(composite_curated, "02_nes_concordance_curated", 12, 7)

summary_table <- bind_rows(
  mutate(concordance(paired), population = "all collections"),
  mutate(concordance(survivors), population = "collapse survivors"),
  mutate(concordance(curated), population = "Hallmark and GO Slim")
) |>
  relocate(population)
print(as.data.frame(summary_table))
stopifnot(
  nrow(discordant_sets) == sum(paired$discordant & paired$significance != "NS"),
  sum(curated$significance != "NS") == nrow(filter(curated, significance != "NS"))
)

packages <- c("here", "fgsea", "dplyr", "ggplot2", "ggrepel", "patchwork", "enrichVolcano")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
export <- paired |>
  transmute(
    set_id, database, pathway, label,
    genes = n,
    nes_x = round(NES_x, 3), nes_y = round(NES_y, 3),
    padj_x = signif(padj_x, 4), padj_y = signif(padj_y, 4),
    significance = as.character(significance), discordant, survivor
  ) |>
  arrange(significance, desc(abs(nes_x) + abs(nes_y)))
readr::write_csv(export, file.path(out, "nes_scatter.csv"))
writexl::write_xlsx(
  list(
    concordance = summary_table,
    discordant = filter(export, discordant, significance != "NS"),
    nes_scatter = export, input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "03_enrich_scatter_fgsea.xlsx")
)
combined <- file.path(figure_dir, "03_enrich_scatter_fgsea_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message("wrote 03_enrich_scatter_fgsea.xlsx and a ", length(pages), "-page figure PDF")
