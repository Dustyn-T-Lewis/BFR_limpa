# Pathway maps of the two training responses: which pathways each leg moves, how they cluster by
# shared proteins, and which proteins link them. compareCluster runs fgsea on the same sets and
# rankings as 01_run_fgsea_and_fry and returns the object enrichplot draws. Figures only: the
# tests and tables stay in 01. Adapted from Sam's PR #26 (sc-norton).

pacman::p_load(here, dplyr, purrr, stringr, ggplot2, readxl, clusterProfiler, enrichplot)

stage <- here("03_Pathway_Enrichment", "08_pathway_maps")
book <- here(
  "03_Pathway_Enrichment", "01_run_fgsea_and_fry", "c_data", "01_run_fgsea_and_fry.xlsx"
)
protein_results <- read_excel(book, "protein_results")
set_tests <- read_excel(book, "set_tests")
sets <- readRDS(here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data", "gene_sets.rds"))
arms <- c(BFR = "BFR_Post-Pre", HLRT = "HLRT_Post-Pre")

ranked <- map(arms, \(cn) {
  with(filter(protein_results, contrast == cn, selected), sort(set_names(t, gene), TRUE))
})
labels <- set_tests |>
  distinct(set_id, database, pathway) |>
  mutate(label = enrichVolcano::clean_label(pathway, width = 1000))

# The main sets 01 kept, one per cluster of overlapping sets, significant in either leg. Only
# GO:BP has more than 30, so the 30 with the smallest padj.
drawn <- set_tests |>
  filter(method == "fgsea", contrast %in% arms, padj < 0.05, main) |>
  summarise(padj = min(padj), .by = c(set_id, database)) |>
  slice_min(padj, n = 30, by = database, with_ties = FALSE) |>
  left_join(select(labels, set_id, label), by = "set_id")

caption <- function(number, title, db, text) {
  labs(caption = str_wrap(sprintf("S%d Figure. %s, %s. %s", number, title, db, text), 115))
}
arm_fill <- scale_fill_manual(values = c(BFR = "#B2182B", HLRT = "#2166AC"), name = NULL)
caption_theme <- theme(
  plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
  plot.caption.position = "plot"
)
source_note <- paste(
  "Pathways are the collapse-surviving fgsea sets that 01_run_fgsea_and_fry calls at FDR < 0.05",
  "in either training contrast. Data: set_tests sheet of 01_run_fgsea_and_fry.xlsx."
)

maps <- function(db) {
  members <- sets[startsWith(names(sets), paste0(db, "|"))]
  term_gene <- stack(members)[, c("ind", "values")]
  term_name <- select(filter(labels, database == db), set_id, label)
  show <- filter(drawn, database == db)$label
  set.seed(1)
  # compareCluster does not pass eps on to fgsea, so p below 1e-10 is floored with a warning.
  # The floor only shades the dots; significance comes from 01.
  compared <- suppressWarnings(compareCluster(
    ranked,
    fun = "GSEA", TERM2GENE = term_gene, TERM2NAME = term_name,
    minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1, verbose = FALSE
  ))
  # Similarity only among the drawn pathways: GO:BP holds 1,366 sets.
  shown <- filter(compared, Description %in% show)
  similar <- pairwise_termsim(shown, showCategory = length(show))
  list(
    # dotplot logs its colour scale, built for p-values; NES is signed, so a diverging scale.
    dotplot(shown, showCategory = length(show), color = "NES", label_format = 45) +
      scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", name = "NES") +
      caption(23, "Pathways by leg", db, paste(
        "Colour is NES, dot size the share of the pathway in its leading edge.", source_note
      )),
    # On the two small collections (Hallmark, GO Slim) few pathways share proteins, and emapplot
    # warns from an internal max() over no edges; the map still draws.
    suppressWarnings(
      emapplot(similar, showCategory = length(show), label_format = 30, size_category = 1.6)
    ) + arm_fill +
      caption(24, "Enrichment map", db, paste(
        "Nodes are pathways, pies split by leg, joined when their proteins overlap.", source_note
      )),
    treeplot(
      similar,
      showCategory = length(show), nCluster = min(5, length(show) - 1), label_format = 32,
      fontsize_tiplab = 2.6, fontsize_cladelab = 2.6, hexpand = 0.5, cladelab_offset = 9
    ) + arm_fill +
      caption(25, "Pathway tree", db, paste(
        "Pathways clustered by shared proteins, each cluster named by its most frequent words.",
        source_note
      )),
    cnetplot(shown, showCategory = head(show, 5), node_label = "category") + arm_fill +
      caption(26, "Pathways and their proteins", db, paste(
        "The five strongest pathways and their leading-edge proteins; pies show which leg each",
        "protein leads in. Shared proteins are named in S21 of 07_leading_edge.", source_note
      ))
  ) |>
    map(\(p) p + caption_theme)
}

walk(unique(drawn$database), \(db) {
  pdf(file.path(stage, "b_reports", paste0("08_maps_", db, ".pdf")), width = 8.27, height = 11.69)
  walk(maps(db), print)
  invisible(dev.off())
})
sessionInfo()
