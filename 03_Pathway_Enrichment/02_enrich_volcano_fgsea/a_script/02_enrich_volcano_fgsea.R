# Protein volcanoes with fgsea rings. Colour reads protein BH FDR, rings read set FDR. The
# control is not drawn.

pacman::p_load(here, dplyr, purrr, stringr, ggplot2, patchwork, readxl)

tests <- here(
  "03_Pathway_Enrichment", "01_run_fgsea_and_fry", "c_data", "01_run_fgsea_and_fry.xlsx"
)
protein_results <- read_excel(tests, "protein_results")
fgsea_results <- read_excel(tests, "set_tests") |>
  filter(method == "fgsea", main) |>
  mutate(leadingEdge = str_split(leadingEdge, ";"))

# plot_volcano_ring() matches leading-edge genes against point labels, which carry an accession when
# a symbol sits on more than one protein. The representative protein is the one fgsea ranked.
gene_to_label <- with(
  filter(protein_results, contrast == contrast[1], selected),
  set_names(label, gene)
)
fgsea_results$leadingEdge <- map(fgsea_results$leadingEdge, \(genes) unname(gene_to_label[genes]))

da <- enrichVolcano::as_da(
  select(protein_results, protein, gene = label, contrast, logFC, t, P.Value, adj.P.Val),
  species = NULL
)
# Default palette: red up, blue down, dark blue to dark red NES ramp.
ring_theme <- enrichVolcano::plot_theme(
  base_size = 10, base_family = "sans", ns = "#9a9a9a"
)
# Contrast names carry their own algebra and are the panel titles. Only the interaction needs
# expanding: a difference of differences has no one-line name.
contrast_subtitle <- c(
  Modality_x_Time_Interaction = "(BFR_Post - BFR_Pre) - (HLRT_Post - HLRT_Pre)"
)

make_volcano <- function(cn, rank_by = "fdr") {
  points <- filter(protein_results, contrast == cn)
  ring <- fgsea_results |>
    filter(contrast == cn) |>
    # Once the prefix goes, GOBP_MUSCLE_CONTRACTION and REACTOME_MUSCLE_CONTRACTION share a label.
    mutate(
      stem = sub("^[A-Z0-9]+_", "", pathway),
      pathway = if_else(
        duplicated(stem) | duplicated(stem, fromLast = TRUE),
        paste0(pathway, " ", database), pathway
      )
    )
  if (rank_by == "pi") {
    labels <- arrange(points, pi_score, protein)
    note <- "Labels ranked by pi-score; colours and counts use BH FDR"
  } else {
    labels <- arrange(filter(points, adj.P.Val < 0.05), adj.P.Val, P.Value, protein)
    note <- "Protein significance: BH FDR < 0.05"
  }
  labels <- head(labels$label, 5)
  enrichment <- ring |>
    select(contrast, database, pathway, NES, pval = p, padj, size = n, leadingEdge) |>
    enrichVolcano::as_enrichment(enrichment_test = "fgsea")
  # The ring draws the twelve survivors with the lowest adjusted p, in either direction.
  enrichVolcano::plot_volcano_ring(
    da, enrichment,
    contrast = cn, collapse = FALSE, title = cn,
    subtitle = paste(na.omit(c(contrast_subtitle[cn], note)), collapse = "\n"),
    label_mode = if (length(labels)) "by_genes" else "none",
    label_genes = labels,
    label_size = 2.6, axis_size = 2.6, count_size = 2.8,
    point_size = 0.9, theme = ring_theme
  ) +
    # The ring draws with clip = "off" on a square panel and puts the NES key on the right, so a
    # label on that edge lands on the colourbar. Below the plot nothing collides.
    guides(fill = guide_colorbar(direction = "horizontal", title.position = "top")) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      plot.subtitle = element_text(size = 8, face = "plain"),
      plot.margin = margin(6, 28, 6, 28),
      legend.position = "bottom",
      legend.justification = "center",
      legend.title = element_text(size = 8, hjust = 0.5),
      legend.key.height = unit(2, "mm"),
      legend.key.width = unit(16, "mm")
    )
}

supplement <- function(number, title, text) {
  plot_annotation(
    caption = str_wrap(sprintf("S%d Figure. %s. %s", number, title, text), 115),
    tag_levels = "A",
    theme = theme(plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2))
  )
}
panel_list <- function(contrasts) {
  paste(sprintf("(%s) %s.", LETTERS[seq_along(contrasts)], contrasts), collapse = " ")
}
ring_text <- paste(
  "The ring shows up to twelve fgsea sets that survived collapsePathways, lowest adjusted p",
  "first, in either direction, coloured by NES; ticks mark their leading-edge proteins. Protein",
  "colour and counts read protein-level BH FDR, the ring set-level fgsea FDR: two separate",
  "claims. Data: protein_results and set_tests sheets of 01_run_fgsea_and_fry.xlsx."
)

# Primary question first, then the two training responses and the between-treatment comparison.
fdr_order <- c(
  "Modality_x_Time_Interaction", "BFR_Post-Pre", "HLRT_Post-Pre", "BFR_Post-HLRT_Post"
)
# Same points, colours and rings on the second page; only the labels move. Pi ranks and selects
# nothing, so no protein named there is a discovery.
pi_order <- c("BFR_Post-Pre", "HLRT_Post-Pre")
figures <- list(
  wrap_plots(map(fdr_order, make_volcano), ncol = 2) +
    supplement(8, "Protein volcanoes with collapse-surviving fgsea sets ringed", paste(
      panel_list(fdr_order), "Points are proteins, coloured when BH FDR < 0.05 within the",
      "contrast, with the five lowest-FDR proteins named.", ring_text
    )),
  wrap_plots(map(pi_order, make_volcano, rank_by = "pi"), ncol = 1) +
    supplement(9, "Protein volcanoes labelled by pi-score", paste(
      panel_list(pi_order), "As S8 Figure, but the five proteins named are those with the",
      "smallest pi-score, P.Value^|logFC| (Xiao et al. 2014). Pi controls no error rate, so no",
      "protein named here is a discovery.", ring_text
    ))
)

figure_dir <- here("03_Pathway_Enrichment", "02_enrich_volcano_fgsea", "b_reports")
pdf(file.path(figure_dir, "02_enrich_volcano_fgsea_figures.pdf"), width = 8.27, height = 11.69)
walk(figures, print)
invisible(dev.off())
message("wrote a ", length(figures), "-page figure PDF")
sessionInfo()
