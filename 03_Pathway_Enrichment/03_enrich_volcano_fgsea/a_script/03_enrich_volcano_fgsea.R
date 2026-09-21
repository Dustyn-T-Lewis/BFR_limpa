# Draw the protein volcanoes with the surviving fgsea pathways ringed. Computes nothing.
# Two separate claims share a panel: point colour and the count badges read protein-level BH
# FDR, the ring reads set-level fgsea FDR. Only sets that survived collapsePathways are ringed.
# The pre-training control is not drawn; 01_run_fgsea reports what it returns.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(ggplot2)
})
here::i_am("config.yml")

cfg <- yaml::read_yaml(here("config.yml"))
stage <- here("03_Pathway_Enrichment", "03_enrich_volcano_fgsea")
figure_dir <- file.path(stage, "b_reports", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  fgsea = "03_Pathway_Enrichment/01_run_fgsea/c_data/fgsea.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds"
)
paths <- map_chr(inputs, here::here)
if (!all(file.exists(paths))) {
  stop(
    "Run 01_run_fgsea first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
fg <- readRDS(paths[["fgsea"]])
d <- readRDS(paths[["design"]])

# Keep exact p-values in the export; bound only the logarithm for plotting.
protein_results <- mutate(fg$protein_results, plot_p = pmax(P.Value, .Machine$double.xmin))
fgsea_results <- filter(fg$fgsea_results, main)
stopifnot(nrow(fgsea_results) > 0, !is.null(protein_results$label))

# volcano_ring() matches leading-edge genes against the point labels, and a label carries its
# accession when a symbol sits on more than one protein. Translate the edges into label space
# or the tick lines silently draw nothing.
gene_to_label <- with(
  filter(protein_results, contrast == contrast[1], !is.na(gene)),
  set_names(label, gene)
)
fgsea_results$leadingEdge <- map(fgsea_results$leadingEdge, \(genes) {
  unname(gene_to_label[genes[genes %in% names(gene_to_label)]])
})

# The default palette is red for up, blue for down, with a dark blue to dark red NES ramp.
plot_theme <- enrichVolcano::volcano_ring_theme(
  base_size = 12, base_family = "sans", palette = "default", ns = "#9a9a9a"
)
contrast_titles <- c(
  BFR_post_vs_pre = "BFR: post vs pre",
  HLRT_post_vs_pre = "High load: post vs pre",
  BFR_vs_HLRT_at_T1 = "BFR vs high load: pre-training control",
  BFR_vs_HLRT_at_T2 = "BFR vs high load: post-training",
  interaction = "Interaction: difference in training responses"
)
# Written out rather than generated, because the nested form says what an interaction is and a
# flat signed sum does not. Each string is evaluated against the fitted contrast matrix below,
# so a notation that drifts from 01_Design stops the run instead of mislabelling a panel.
contrast_notation <- c(
  BFR_post_vs_pre = "BFR_T2 - BFR_T1",
  HLRT_post_vs_pre = "HLRT_T2 - HLRT_T1",
  BFR_vs_HLRT_at_T1 = "BFR_T1 - HLRT_T1",
  BFR_vs_HLRT_at_T2 = "BFR_T2 - HLRT_T2",
  interaction = "(BFR_T2 - BFR_T1) - (HLRT_T2 - HLRT_T1)"
)
levels_used <- rownames(d$contrasts)[rowSums(abs(d$contrasts)) > 0]
indicators <- set_names(asplit(diag(length(levels_used)), 1), levels_used)
stopifnot(map_lgl(set_names(colnames(d$contrasts)), \(name) {
  fitted <- d$contrasts[levels_used, name]
  all(eval(parse(text = contrast_notation[[name]]), indicators) == fitted)
}))

make_volcano <- function(contrast, rank_by = "fdr") {
  points <- filter(protein_results, .data$contrast == .env$contrast)
  ring <- fgsea_results |>
    filter(.data$contrast == .env$contrast, padj < cfg$set_fdr) |>
    slice_min(padj, n = cfg$volcano_ring_n) |>
    # Ten overlapping databases mean two ringed sets can carry the same name, and
    # volcano_ring() strips the database prefix before drawing. GOBP_MUSCLE_CONTRACTION and
    # REACTOME_MUSCLE_CONTRACTION then label two arcs identically.
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
      slice_head(n = cfg$volcano_label_n) |>
      pull(label)
    subtitle <- paste0(
      contrast_notation[contrast],
      "\nLabels ranked by pi-score; colours and counts use BH FDR"
    )
  } else {
    labels <- points |>
      filter(adj.P.Val < cfg$protein_fdr) |>
      arrange(adj.P.Val, P.Value, protein) |>
      slice_head(n = cfg$volcano_label_n) |>
      pull(label)
    subtitle <- paste0(
      contrast_notation[contrast],
      "\nProtein significance: BH FDR < ", cfg$protein_fdr
    )
  }
  enrichVolcano::volcano_ring(
    volc_df = select(points, label, logFC, plot_p, adj.P.Val),
    enrich_df = ring,
    gene_col = "label", pval_col = "plot_p", padj_col = "padj",
    volc_sig_col = "adj.P.Val", genes_col = "leadingEdge",
    p_threshold = cfg$protein_fdr, logfc_threshold = 0,
    title = unname(contrast_titles[contrast]), subtitle = subtitle,
    label_mode = if (length(labels)) "by_genes" else "none",
    label_genes = labels, label_n = cfg$volcano_label_n,
    label_size = 3.1, axis_size = 3.1, count_size = 3.3,
    point_size = 1.2, theme = plot_theme
  ) +
    # volcano_ring() draws with clip = "off" against a square panel and puts the NES key on
    # the right, so a label anchored on that edge lands on top of the colourbar. Below the
    # plot there is nothing to collide with.
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
  for (extension in c("png", "pdf")) {
    ggsave(file.path(figure_dir, paste0(name, ".", extension)),
      plot = volcano, width = 7, height = 6.5, units = "in", dpi = 300, bg = "white"
    )
  }
}

# The primary question first, then the two training responses and the between-treatment
# comparison. The pre-training control is fitted and reported in 01_run_fgsea, not drawn here.
plot_order <- c("interaction", "BFR_post_vs_pre", "HLRT_post_vs_pre", "BFR_vs_HLRT_at_T2")
for (contrast in plot_order) {
  save_volcano(make_volcano(contrast), paste0("protein_volcano_fdr_", contrast))
  message("drew ", contrast)
}
# Same points, colours, rings and notation; only the labels move. A pi label ranks and selects
# nothing, so nothing named on these two panels has been discovered.
for (contrast in c("BFR_post_vs_pre", "HLRT_post_vs_pre")) {
  volcano <- make_volcano(contrast, rank_by = "pi")
  save_volcano(volcano, paste0("protein_volcano_pi_rank_", contrast))
  message("drew ", contrast, " ranked by pi")
}

drawn <- sort(basename(list.files(figure_dir, pattern = "[.]png$")))
print(tibble(file = drawn))
stopifnot(length(drawn) == length(plot_order) + 2)
