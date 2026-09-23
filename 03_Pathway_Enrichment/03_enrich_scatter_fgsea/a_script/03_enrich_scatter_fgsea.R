# Put the two training contrasts on one pair of axes. fgsea scores each set once per contrast,
# so plotting one NES against the other asks directly whether the two modalities move the same
# biology. Computes no new test; it reads 01_run_fgsea_and_fry and reshapes.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
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

# Hallmark and GO Slim only, the two collections whose members do not nest inside one another.
# GO:BP's 1,366 sets would bury the plot and, being largely restatements of each other, would
# pull the rank correlation toward whichever branch happens to be densest.
x_contrast <- "BFR_Post-Pre"
y_contrast <- "HLRT_Post-Pre"
paired <- fg$set_tests |>
  filter(
    method == "fgsea", database %in% c("Hallmark", "GO_Slim"),
    contrast %in% c(x_contrast, y_contrast)
  ) |>
  mutate(contrast = if_else(contrast == x_contrast, "x", "y")) |>
  pivot_wider(
    id_cols = c(set_id, database, pathway, n),
    names_from = contrast, values_from = c(NES, padj)
  ) |>
  mutate(
    significance = case_when(
      padj_x < 0.05 & padj_y < 0.05 ~ "Both",
      padj_x < 0.05 ~ x_contrast,
      padj_y < 0.05 ~ y_contrast,
      .default = "NS"
    ) |> factor(c("Both", x_contrast, y_contrast, "NS")),
    label = enrichVolcano::ev_clean_label(pathway)
  )
stopifnot(nrow(paired) > 0, !anyNA(paired$NES_x), !anyNA(paired$NES_y))

# Spearman over every set, so the number describes the experiment rather than a selection.
concordance <- cor.test(paired$NES_x, paired$NES_y, method = "spearman", exact = FALSE)
signif_sets <- filter(paired, significance != "NS")
agreement <- mean(sign(signif_sets$NES_x) == sign(signif_sets$NES_y))
quadrants <- count(signif_sets, up_x = NES_x > 0, up_y = NES_y > 0)
summary_row <- tibble(
  sets = nrow(paired), significant = nrow(signif_sets),
  both = sum(paired$significance == "Both"),
  x_only = sum(paired$significance == x_contrast),
  y_only = sum(paired$significance == y_contrast),
  rho = round(concordance$estimate, 3), p = concordance$p.value,
  same_sign = round(agreement, 3)
)
print(as.data.frame(summary_row))

limit <- max(abs(c(paired$NES_x, paired$NES_y))) * 1.08
figure <- ggplot(paired, aes(NES_x, NES_y)) +
  geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_point(
    data = filter(paired, significance == "NS"),
    colour = "grey78", size = 0.7, alpha = 0.5
  ) +
  geom_point(data = signif_sets, aes(colour = significance, size = n), alpha = 0.8) +
  ggrepel::geom_text_repel(
    data = slice_min(mutate(signif_sets, lobe = NES_x > 0), padj_x + padj_y,
      n = 6, by = lobe, with_ties = FALSE
    ),
    aes(label = label), size = 2.4, colour = "grey15", segment.colour = "grey60",
    segment.size = 0.25, min.segment.length = 0, max.overlaps = Inf, seed = 1, force = 6
  ) +
  scale_colour_manual(
    values = set_names(c("#6A3D9A", "#D7301F", "#2B6CB0"), c("Both", x_contrast, y_contrast)),
    name = "significant in", drop = FALSE
  ) +
  scale_size_continuous(range = c(1.2, 5), name = "genes") +
  coord_fixed(xlim = c(-limit, limit), ylim = c(-limit, limit)) +
  labs(
    x = paste("NES,", x_contrast), y = paste("NES,", y_contrast),
    title = "The two training modalities move the same pathways",
    subtitle = sprintf(
      paste(
        "Hallmark and GO Slim, %d sets (%d significant) | Spearman rho = %.2f",
        "| %.0f%% of significant sets agree in sign"
      ),
      summary_row$sets, summary_row$significant, summary_row$rho, 100 * summary_row$same_sign
    )
  ) +
  theme_minimal(base_size = 9) +
  theme(legend.position = "bottom", plot.subtitle = element_text(size = 7.5))

walk(c("png", "pdf"), \(extension) {
  ggsave(
    file.path(figure_dir, paste0("01_nes_scatter.", extension)), figure,
    width = 7, height = 7, dpi = 300, bg = "white"
  )
})

packages <- c("here", "fgsea", "dplyr", "ggplot2", "ggrepel", "enrichVolcano")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
export <- paired |>
  transmute(set_id, database, pathway,
    genes = n,
    nes_x = round(NES_x, 3), nes_y = round(NES_y, 3),
    padj_x = signif(padj_x, 4), padj_y = signif(padj_y, 4),
    significance = as.character(significance)
  ) |>
  arrange(significance, desc(abs(nes_x) + abs(nes_y)))
readr::write_csv(export, file.path(out, "nes_scatter.csv"))
writexl::write_xlsx(
  list(
    concordance = summary_row, quadrants = quadrants, nes_scatter = export,
    input_manifest = manifest, package_versions = versions
  ),
  file.path(out, "03_enrich_scatter_fgsea.xlsx")
)
message("wrote 01_nes_scatter and 03_enrich_scatter_fgsea.xlsx")
