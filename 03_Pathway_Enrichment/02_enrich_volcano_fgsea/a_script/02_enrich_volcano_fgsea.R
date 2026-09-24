# Protein volcanoes with collapse-surviving fgsea sets ringed. Computes nothing. Two claims share
# a panel: colour and count badges read protein-level BH FDR, rings read set-level fgsea FDR.
# The pre-training control is not drawn; 01_run_fgsea_and_fry reports it.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})

stage <- here("03_Pathway_Enrichment", "02_enrich_volcano_fgsea")
figure_dir <- file.path(stage, "b_reports")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(set_tests = "03_Pathway_Enrichment/01_run_fgsea_and_fry/c_data/set_tests.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run 01_run_fgsea_and_fry first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
fg <- readRDS(paths[["set_tests"]])

# Keep exact p-values in the export; bound only the logarithm for plotting.
protein_results <- mutate(fg$protein_results, plot_p = pmax(P.Value, .Machine$double.xmin))
fgsea_results <- filter(fg$set_tests, method == "fgsea", main)
stopifnot(nrow(fgsea_results) > 0, !is.null(protein_results$label))

# volcano_ring() matches leading-edge genes against point labels, which carry an accession when
# a symbol sits on more than one protein. Untranslated edges make the tick lines draw nothing.
gene_to_label <- with(
  filter(protein_results, contrast == contrast[1], !is.na(gene)),
  set_names(label, gene)
)
fgsea_results$leadingEdge <- map(fgsea_results$leadingEdge, \(genes) {
  unname(gene_to_label[genes[genes %in% names(gene_to_label)]])
})

# Default palette: red up, blue down, dark blue to dark red NES ramp.
plot_theme <- enrichVolcano::volcano_ring_theme(
  base_size = 12, base_family = "sans", palette = "default", ns = "#9a9a9a"
)
# Contrast names carry their own algebra and are the panel titles. Only the interaction needs
# expanding: a difference of differences has no one-line name.
contrast_subtitle <- c(
  Modality_x_Time_Interaction = "(BFR_Post - BFR_Pre) - (HLRT_Post - HLRT_Pre)"
)

make_volcano <- function(contrast, rank_by = "fdr") {
  points <- filter(protein_results, .data$contrast == .env$contrast)
  ring <- fgsea_results |>
    filter(.data$contrast == .env$contrast, padj < 0.05) |>
    slice_min(padj, n = 8) |>
    # With five overlapping collections two ringed sets can share a name once volcano_ring()
    # strips the prefix: GOBP_MUSCLE_CONTRACTION and REACTOME_MUSCLE_CONTRACTION would label
    # two arcs identically.
    mutate(
      stem = sub("^[A-Z0-9]+_", "", pathway),
      pathway = if_else(
        duplicated(stem) | duplicated(stem, fromLast = TRUE),
        paste0(pathway, " ", database), pathway
      )
    )
  if (rank_by == "pi") {
    labels <- points |>
      arrange(pi_score, protein) |>
      slice_head(n = 5) |>
      pull(label)
    subtitle <- paste0(
      contrast_subtitle[contrast] %||% "",
      "\nLabels ranked by pi-score; colours and counts use BH FDR"
    )
  } else {
    labels <- points |>
      filter(adj.P.Val < 0.05) |>
      arrange(adj.P.Val, P.Value, protein) |>
      slice_head(n = 5) |>
      pull(label)
    subtitle <- paste0(
      contrast_subtitle[contrast] %||% "",
      "\nProtein significance: BH FDR < 0.05"
    )
  }
  enrichVolcano::volcano_ring(
    volc_df = select(points, label, logFC, plot_p, adj.P.Val),
    enrich_df = ring,
    gene_col = "label", pval_col = "plot_p", padj_col = "padj",
    volc_sig_col = "adj.P.Val", genes_col = "leadingEdge",
    p_threshold = 0.05, logfc_threshold = 0,
    title = contrast, subtitle = subtitle,
    label_mode = if (length(labels)) "by_genes" else "none",
    label_genes = labels, label_n = 5,
    label_size = 3.1, axis_size = 3.1, count_size = 3.3,
    point_size = 1.2, theme = plot_theme
  ) +
    # volcano_ring() draws with clip = "off" on a square panel and puts the NES key on the right,
    # so a label on that edge lands on the colourbar. Below the plot nothing collides.
    guides(fill = guide_colorbar(direction = "horizontal", title.position = "top")) +
    theme(
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10, face = "plain"),
      plot.margin = margin(8, 8, 8, 8),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 9, hjust = 0.5),
      legend.key.height = unit(2.5, "mm"),
      legend.key.width = unit(22, "mm")
    )
}
save_volcano <- function(volcano, name) {
  walk(c("png", "pdf"), \(extension) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)),
      plot = volcano, width = 7, height = 6.5, units = "in", dpi = 300, bg = "white"
    )
  })
}

# Primary question first, then the two training responses and the between-treatment comparison.
plot_order <- c(
  "Modality_x_Time_Interaction", "BFR_Post-Pre", "HLRT_Post-Pre", "BFR_Post-HLRT_Post"
)
for (contrast in plot_order) {
  save_volcano(make_volcano(contrast), paste0("protein_volcano_fdr_", contrast))
  message("drew ", contrast)
}
# Same points, colours and rings; only the labels move. Pi ranks and selects nothing, so no
# protein named on these two panels is a discovery.
for (contrast in c("BFR_Post-Pre", "HLRT_Post-Pre")) {
  volcano <- make_volcano(contrast, rank_by = "pi")
  save_volcano(volcano, paste0("protein_volcano_pi_rank_", contrast))
  message("drew ", contrast, " ranked by pi")
}

drawn <- sort(basename(list.files(figure_dir, pattern = "[.]png$")))
print(tibble(file = drawn))
stopifnot(length(drawn) == length(plot_order) + 2)

combined <- file.path(figure_dir, "02_enrich_volcano_fgsea_figures.pdf")
pages <- setdiff(list.files(figure_dir, "[.]pdf$", full.names = TRUE), combined)
invisible(qpdf::pdf_combine(sort(pages), combined))
message("wrote a ", length(pages), "-page figure PDF")
