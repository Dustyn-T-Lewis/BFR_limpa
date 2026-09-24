# fgsea and fry, side by side; neither is simply better. fgsea is competitive on moderated t and
# assumes exchangeable proteins, which they are not. fry is self-contained and rotates residuals
# of the participant design, so correlation cannot inflate its null; under the global training
# effect it flags about a third of sets.

pacman::p_load(
  here, dplyr, tibble, tidyr, purrr, stringr, limma, ggplot2, patchwork, readxl, writexl
)

stage <- here("03_Pathway_Enrichment", "01_run_fgsea_and_fry")
gene_sets <- here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data")
sets <- readRDS(file.path(gene_sets, "gene_sets.rds"))
gene_set_book <- file.path(gene_sets, "00_build_gene_sets.xlsx")
set_catalog <- read_excel(gene_set_book, "set_catalog")
protein_map <- read_excel(gene_set_book, "protein_gene_map")
fit <- readRDS(here("02_Differential_Expression", "02_Differential", "c_data", "fit.rds"))
d <- readRDS(here("02_Differential_Expression", "01_Design", "c_data", "design.rds"))
proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))

contrast_names <- colnames(d$contrasts)
# A stale fit still lines up by count, so these compare names rather than lengths.
stopifnot(
  identical(rownames(fit$coefficients), protein_map$protein),
  identical(colnames(fit$coefficients), contrast_names),
  identical(rownames(proteins$E), protein_map$protein)
)

# All five contrasts, BH once per contrast as upstream applied it; re-adjusting on the mapped
# subset would change what every FDR here means. pi_score is the Xiao et al. (2014, PMID
# 22321699) score: it controls no error rate, so it orders a contrast and selects nothing.
protein_results <- map(set_names(contrast_names), function(contrast) {
  topTable(fit, coef = contrast, number = Inf, adjust.method = "BH", sort.by = "none") |>
    rownames_to_column("protein") |>
    select(protein, logFC, AveExpr, t, P.Value, adj.P.Val)
}) |>
  list_rbind(names_to = "contrast") |>
  left_join(protein_map, by = "protein", relationship = "many-to-one") |>
  mutate(pi_score = P.Value^abs(logFC))

protein_summary <- protein_results |>
  summarise(
    tested = n(),
    up = sum(adj.P.Val < 0.05 & logFC > 0),
    down = sum(adj.P.Val < 0.05 & logFC < 0),
    min_fdr = signif(min(adj.P.Val), 3), .by = contrast
  )
print(protein_summary)

# One vector per contrast, representative proteins only, keyed by gene symbol so fgsea's
# leadingEdge comes back in the namespace the volcanoes label with.
ranked <- protein_results |>
  filter(selected) |>
  split(~contrast) |>
  map(\(result) set_names(result$t, result$gene))
ranked <- ranked[contrast_names]

set.seed(1)
fgsea_raw <- map(ranked, \(stats) fgsea::fgsea(sets, stats, minSize = 15, maxSize = 500))

# fry indexes rows of proteins$E, so sets map back through the representative protein of each
# symbol. limma reads y$weights, so limpa's per-observation precision reaches fry: a protein
# rebuilt largely from missing precursors counts for less. fry takes no block argument, because
# 01_Design fixes participant in the design and the residual within-leg correlation is 0.027.
gene_map <- filter(protein_map, selected)
set_rows <- map(sets, \(genes) {
  match(gene_map$protein[match(genes, gene_map$gene)], protein_map$protein)
})
proteins$weights <- fit$EList$weights

# Both packages return their own BH column over the list they were given: fgsea's padj and
# fry's FDR. Neither is recomputed here.
set_tests <- bind_rows(
  map(contrast_names, \(contrast) {
    as_tibble(fgsea_raw[[contrast]]) |>
      transmute(
        set_id = pathway, contrast, method = "fgsea", n = size,
        direction = if_else(NES > 0, "Up", "Down"), NES, p = pval, padj,
        leadingEdge
      )
  }),
  map(contrast_names, \(contrast) {
    fry(proteins, set_rows, d$design, d$contrasts[, contrast], sort = "none") |>
      rownames_to_column("set_id") |>
      transmute(
        set_id, contrast,
        method = "fry", n = NGenes,
        direction = Direction, NES = NA_real_, p = PValue, padj = FDR,
        leadingEdge = list(NULL)
      )
  })
)

# Four collections overlap: glycolysis is tested in Hallmark, KEGG, Reactome and GO.
# collapsePathways re-runs each significant set conditioned on a stronger set's leading edge and
# keeps it only if it stands alone. It needs the results, so it runs after testing, and it
# prunes without re-adjusting, so the surviving FDR is conservative, not inflated.
main_sets <- map(set_names(contrast_names), function(contrast) {
  significant <- fgsea_raw[[contrast]][padj < 0.05][order(pval)]
  if (nrow(significant) < 2) {
    return(significant$pathway)
  }
  fgsea::collapsePathways(significant, sets, ranked[[contrast]], pval.threshold = 0.05)$mainPathways
})

main_lookup <- imap(main_sets, \(ids, cn) tibble(contrast = cn, set_id = ids, kept = TRUE)) |>
  list_rbind()
set_tests <- set_tests |>
  left_join(main_lookup, by = c("contrast", "set_id")) |>
  mutate(main = if_else(method == "fgsea", coalesce(kept, FALSE), NA), kept = NULL) |>
  left_join(
    select(set_catalog, set_id, database, pathway, source_size, description),
    by = "set_id"
  ) |>
  relocate(contrast, method, set_id, database, pathway)

# Sets called per test and contrast, and how many survived collapse. The negative control sits
# in the table unflagged, for comparison.
set_summary <- set_tests |>
  summarise(
    sets = n_distinct(set_id),
    fgsea = sum(method == "fgsea" & padj < 0.05),
    collapsed = sum(method == "fgsea" & padj < 0.05 & main),
    fry = sum(method == "fry" & padj < 0.05),
    .by = contrast
  )
print(as.data.frame(set_summary))


# ---- figures: S1 to S7, one A4 page each ---------------------------------------------------

supplement <- function(number, title, text, tags = "A") {
  plot_annotation(
    caption = str_wrap(sprintf("S%d Figure. %s. %s", number, title, text), 115),
    tag_levels = tags,
    theme = theme(plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2))
  )
}

shown <- c("BFR_Post-Pre", "HLRT_Post-Pre", "Modality_x_Time_Interaction", "BFR_Post-HLRT_Post")
collections <- unique(set_catalog$database[set_catalog$qualifies])
survivors <- set_tests |>
  filter(method == "fgsea", padj < 0.05, main, contrast %in% shown) |>
  mutate(contrast = factor(contrast, shown))

dotplot <- function(rows, by_database, sizes, fdr) {
  top <- rows |>
    slice_min(padj, n = 10, with_ties = FALSE) |>
    mutate(label = enrichVolcano::clean_label(pathway, width = 45))
  colour <- if (by_database) {
    list(
      scale_colour_brewer(palette = "Dark2", limits = collections, name = NULL),
      guides(colour = guide_legend(override.aes = list(size = 3)))
    )
  } else {
    scale_colour_viridis_c(
      option = "rocket", direction = -1, end = 0.9, limits = fdr,
      name = expression(-log[10] ~ FDR)
    )
  }
  ggplot(top, aes(NES, reorder(label, NES), size = n)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
    geom_point(
      aes(colour = if (by_database) database else -log10(padj)),
      alpha = 0.9, show.legend = TRUE
    ) +
    colour +
    scale_size_continuous(range = c(1.5, 5), limits = sizes, name = "genes") +
    labs(
      x = "normalised enrichment score", y = NULL,
      title = sprintf("%s, %d of %d survivors", rows$contrast[1], nrow(top), nrow(rows))
    ) +
    theme_minimal(base_size = 9) +
    theme(panel.grid.major.y = element_blank(), plot.title = element_text(size = 9))
}

dotplot_figure <- function(rows, number, name, by_database = FALSE) {
  by_contrast <- split(rows, rows$contrast, drop = TRUE)
  wrap_plots(
    map(by_contrast, dotplot,
      by_database = by_database, sizes = range(rows$n), fdr = range(-log10(rows$padj))
    ),
    ncol = 1, guides = "collect"
  ) +
    supplement(
      number, paste("Strongest fgsea sets after collapse,", name),
      paste0(
        paste(sprintf("(%s) %s.", LETTERS[seq_along(by_contrast)], names(by_contrast)),
          collapse = " "
        ),
        " Up to ten collapsePathways survivors per contrast, lowest adjusted p first; FDR is",
        " BH within contrast. Point size is the number of measured genes in the set; colour is ",
        if (by_database) "the collection" else "-log10 FDR",
        ". Contrasts with no survivor are omitted. Data: set_tests sheet of",
        " 01_run_fgsea_and_fry.xlsx."
      )
    )
}

figures <- c(
  list(dotplot_figure(survivors, 1, "all collections", by_database = TRUE)),
  imap(unname(collections), \(db, i) {
    dotplot_figure(filter(survivors, database == db), i + 1, db)
  })
)

# Shows how much redundancy collapsePathways removed, rather than asserting it.
collapse_effect <- set_tests |>
  filter(method == "fgsea", padj < 0.05, contrast %in% shown) |>
  summarise(before = n(), after = sum(main), .by = c(contrast, database)) |>
  pivot_longer(c(before, after), names_to = "stage", values_to = "sets") |>
  mutate(stage = factor(stage, c("before", "after")), contrast = factor(contrast, shown))
figures <- c(figures, list(
  ggplot(collapse_effect, aes(stage, sets, fill = database)) +
    geom_col(position = "dodge") +
    facet_wrap(~contrast, scales = "free_y", ncol = 2) +
    scale_fill_brewer(palette = "Dark2", name = NULL) +
    labs(x = NULL, y = "significant sets") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom") +
    supplement(
      length(figures) + 1, "Significant fgsea sets before and after collapse",
      paste(
        "Sets at FDR 0.05 per collection and contrast, before and after collapsePathways, which",
        "keeps a set only if it stands alone given a stronger set's leading edge. Data: set_tests",
        "sheet of 01_run_fgsea_and_fry.xlsx."
      ),
      tags = NULL
    )
))

pdf(file.path(stage, "b_reports", "01_run_fgsea_and_fry_figures.pdf"), width = 8.27, height = 11.69)
walk(figures, print)
invisible(dev.off())

sheets <- list(
  set_summary = set_summary,
  protein_summary = protein_summary,
  set_tests = mutate(set_tests, leadingEdge = map_chr(leadingEdge, paste, collapse = ";")),
  protein_results = protein_results
)
overview <- tibble(
  sheet = names(sheets),
  rows = map_int(sheets, nrow),
  columns = map_int(sheets, ncol),
  description = c(
    "Sets per contrast significant by fgsea, by fgsea after collapse, and by fry, at FDR 0.05",
    "Proteins per contrast up and down at FDR 0.05",
    "Every set, contrast and method; main marks fgsea collapse survivors; leadingEdge is ;-joined",
    "Every protein and contrast from the saved fit, with its gene mapping and pi score"
  )
)
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "01_run_fgsea_and_fry.xlsx")
)
message("wrote 01_run_fgsea_and_fry.xlsx and a ", length(figures), "-page figure PDF")
sessionInfo()
