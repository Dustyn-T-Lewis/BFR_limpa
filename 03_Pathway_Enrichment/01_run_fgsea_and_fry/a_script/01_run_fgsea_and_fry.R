# Two set tests, reported side by side.
#
# fgsea is competitive: it ranks proteins by moderated t and asks whether a set piles up at one
# end, relative to every other protein. That assumes the ranked proteins are exchangeable, and
# they are not.
#
# fry is self-contained: it asks whether a set moved at all, rotating the residuals of the fitted
# model rather than shuffling gene labels. Correlation rotates with the data, so it cannot inflate
# the null, and the design's participant structure is built in.
#
# Neither is simply better. fry flags roughly a third of all sets on the training contrasts,
# because with a global effect most sets did move a little. Both ship as rows of one table, and
# set_summary puts their counts side by side, the negative control among them.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(limma)
  library(ggplot2)
  library(stringr)
})

out <- here("03_Pathway_Enrichment", "01_run_fgsea_and_fry", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  fit = "02_Differential_Expression/02_Differential/c_data/fit.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds",
  proteins = "01_Preprocess/02_Quantification/c_data/proteins.rds"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop(
    "Run the upstream stages first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
gs <- readRDS(paths[["gene_sets"]])
fit <- readRDS(paths[["fit"]])
d <- readRDS(paths[["design"]])
proteins <- readRDS(paths[["proteins"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

sets <- gs$sets
protein_map <- gs$protein_map
protein_ids <- protein_map$protein
contrast_names <- colnames(d$contrasts)
# A stale fit still lines up by count, so these compare names rather than lengths.
stopifnot(
  identical(rownames(fit$coefficients), protein_ids),
  identical(colnames(fit$coefficients), contrast_names),
  identical(rownames(proteins$E), protein_ids),
  d$negative_control %in% contrast_names
)

# topTable rebuilds all five contrasts from the saved fit. BH stays as upstream applied it,
# one adjustment per contrast; re-adjusting on the mapped subset would change what every FDR
# here means. pi_score is the Xiao et al. (2014, PMID 22321699) score: it controls no error
# rate, so it orders a contrast and selects nothing.
protein_results <- map(set_names(contrast_names), function(contrast) {
  topTable(fit,
    coef = contrast, number = Inf, adjust.method = "BH", sort.by = "none"
  ) |>
    rownames_to_column("protein") |>
    select(protein, logFC, AveExpr, t, P.Value, adj.P.Val)
}) |>
  list_rbind(names_to = "contrast") |>
  left_join(protein_map, by = "protein", relationship = "many-to-one") |>
  mutate(pi_score = P.Value^abs(logFC))
stopifnot(nrow(protein_results) == length(protein_ids) * length(contrast_names))

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
  filter(selected, !is.na(gene)) |>
  split(~contrast) |>
  map(\(result) set_names(result$t, result$gene))
ranked <- ranked[contrast_names]

# fry indexes rows of proteins$E, not genes, so the sets are mapped back through the
# representative protein chosen in 00_build_gene_sets. A missing index would silently test the
# wrong rows.
gene_map <- filter(protein_map, selected)
set_rows <- map(sets, \(genes) {
  match(gene_map$protein[match(genes, gene_map$gene)], protein_ids)
})
stopifnot(!any(map_lgl(set_rows, anyNA)))

# limma reads y$weights, so attaching limpa's per-observation precision to the EList carries it
# into fry without depending on an argument name. A protein rebuilt largely from missing
# precursors counts for less.
weights <- fit$EList$weights
stopifnot(identical(dim(weights), dim(proteins$E)), all(is.finite(weights)), all(weights > 0))
proteins$weights <- weights

# fry takes no block argument: 01_Design fixes participant in the design and the residual
# within-leg correlation is 0.027.
set.seed(1)
fgsea_raw <- map(ranked, \(stats) fgsea::fgsea(sets, stats, minSize = 15, maxSize = 500))

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

# Four collections overlap, so glycolysis is tested in Hallmark, KEGG, Reactome and again in
# GO. collapsePathways re-runs each significant set conditioned on a more significant one's
# leading edge and keeps it only if it still stands alone. It needs the results, so it can only
# run after testing, and it prunes rather than re-adjusting: the surviving list's FDR is
# conservative, not inflated.
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
    select(gs$set_catalog, set_id, database, pathway, theme, source_size, description),
    by = "set_id"
  ) |>
  relocate(contrast, method, set_id, database, pathway, theme)

# One row per contrast: how many sets each test called, and how many survived collapse. The
# control sits in the table unflagged, for the reader to compare against.
set_summary <- set_tests |>
  summarise(
    sets = n_distinct(set_id),
    fgsea = sum(method == "fgsea" & padj < 0.05),
    collapsed = sum(method == "fgsea" & padj < 0.05 & main),
    fry = sum(method == "fry" & padj < 0.05),
    .by = contrast
  )
print(as.data.frame(set_summary))

themed_hits <- set_tests |>
  filter(method == "fgsea", padj < 0.05, main) |>
  mutate(theme = coalesce(theme, paste0("(", database, ", no hierarchy)"))) |>
  count(contrast, theme, sort = TRUE, name = "sets")

# ---- figures -------------------------------------------------------------------------------

# One directory per collection plus all_db, one file per contrast, so a panel can be read at full
# size. Each shows the ten strongest collapse survivors by adjusted p, the rule the volcano rings
# use. Labels come from enrichVolcano::ev_clean_label, the same function the volcanoes use, with
# its line breaks flattened because these sit on an axis.
figure_root <- here("03_Pathway_Enrichment", "01_run_fgsea_and_fry", "b_reports")
flat_label <- function(x) gsub("\n", " ", enrichVolcano::ev_clean_label(x))
shown <- c("BFR_Post-Pre", "HLRT_Post-Pre", "Modality_x_Time_Interaction", "BFR_Post-HLRT_Post")
survivors <- set_tests |>
  filter(method == "fgsea", padj < 0.05, main, contrast %in% shown)

draw_dotplot <- function(rows, colour_by, file) {
  top <- rows |>
    slice_min(padj, n = 10, with_ties = FALSE) |>
    mutate(label = str_trunc(flat_label(pathway), 46))
  figure <- ggplot(top, aes(NES, reorder(label, NES), size = n, colour = .data[[colour_by]])) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey75") +
    geom_point(alpha = 0.9) +
    scale_size_continuous(range = c(2, 6), name = "genes") +
    labs(x = "normalised enrichment score", y = NULL, title = unique(top$contrast)) +
    theme_minimal(base_size = 10) +
    theme(panel.grid.major.y = element_blank())
  figure <- if (colour_by == "database") {
    figure + scale_colour_brewer(palette = "Dark2", name = NULL)
  } else {
    figure + scale_colour_viridis_c(
      option = "rocket", direction = -1, end = 0.9,
      name = expression(-log[10] ~ FDR)
    )
  }
  ggsave(file, figure, width = 7.5, height = 1.6 + 0.22 * nrow(top), dpi = 200, bg = "white")
}

collections <- unique(gs$set_catalog$database[gs$set_catalog$qualifies])
for (db in c(collections, "all_db")) {
  dir <- file.path(figure_root, db)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  rows <- if (db == "all_db") survivors else filter(survivors, database == db)
  drawn <- intersect(shown, unique(rows$contrast))
  for (cn in drawn) {
    draw_dotplot(
      mutate(filter(rows, contrast == cn), `-log10 FDR` = -log10(padj)),
      if (db == "all_db") "database" else "-log10 FDR",
      file.path(dir, paste0("01_dotplot_", cn, ".png"))
    )
  }
  message("drew ", db, ": ", length(drawn), " contrasts")
}

# How much redundancy collapsePathways removed, so the pruning is visible rather than asserted.
collapse_effect <- set_tests |>
  filter(method == "fgsea", padj < 0.05, contrast %in% shown) |>
  summarise(before = n(), after = sum(main), .by = c(contrast, database)) |>
  tidyr::pivot_longer(c(before, after), names_to = "stage", values_to = "sets") |>
  mutate(
    stage = factor(stage, c("before", "after")),
    contrast = factor(contrast, levels = shown)
  )
ggsave(
  file.path(figure_root, "all_db", "02_collapse_before_after.png"),
  ggplot(collapse_effect, aes(stage, sets, fill = database)) +
    geom_col(position = "dodge") +
    facet_wrap(~contrast, scales = "free_y", nrow = 1) +
    scale_fill_brewer(palette = "Dark2", name = NULL) +
    labs(x = NULL, y = "significant sets", title = "collapsePathways removes redundancy only") +
    theme_minimal(base_size = 9),
  width = 10, height = 3, dpi = 200, bg = "white"
)

packages <- c("here", "limma", "fgsea", "dplyr", "purrr", "enrichVolcano")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
flat <- mutate(set_tests, leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

saveRDS(
  list(
    protein_results = protein_results, set_tests = set_tests, set_summary = set_summary,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "set_tests.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    set_summary = set_summary,
    significant = filter(flat, padj < 0.05),
    themed_hits = themed_hits,
    protein_summary = protein_summary,
    protein_results = protein_results,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "01_run_fgsea_and_fry.xlsx")
)
readr::write_csv(flat, file.path(out, "set_tests.csv"))
message("wrote set_tests.rds, 01_run_fgsea_and_fry.xlsx and set_tests.csv")
