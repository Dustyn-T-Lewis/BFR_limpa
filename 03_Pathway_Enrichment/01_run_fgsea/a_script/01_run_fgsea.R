# Rank proteins by moderated t and ask fgsea whether a set piles up at one end, then let
# collapsePathways prune the significant sets whose signal a more significant one already
# carries. The negative-control count prints before any result, because fgsea assumes the
# ranked proteins are exchangeable and they are not.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
  library(data.table)
})
here::i_am("config.yml")

cfg <- yaml::read_yaml(here("config.yml"))
out <- here("03_Pathway_Enrichment", "01_run_fgsea", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  fit = "02_Differential_Expression/02_Differential/c_data/fit.rds",
  design = "02_Differential_Expression/01_Design/c_data/design.rds"
)
paths <- map_chr(inputs, here::here)
if (!all(file.exists(paths))) {
  stop(
    "Run the upstream stages first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
gs <- readRDS(paths[["gene_sets"]])
fit <- readRDS(paths[["fit"]])
d <- readRDS(paths[["design"]])
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
  d$negative_control %in% contrast_names
)

# topTable rebuilds all five contrasts from the saved fit. BH stays as upstream applied it,
# one adjustment per contrast; re-adjusting on the mapped subset would change what every FDR
# here means. pi_score is the Xiao et al. (2014, PMID 22321699) score: it controls no error
# rate, so it orders a contrast and selects nothing.
protein_results <- map(set_names(contrast_names), function(contrast) {
  limma::topTable(fit,
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
    up = sum(adj.P.Val < cfg$protein_fdr & logFC > 0),
    down = sum(adj.P.Val < cfg$protein_fdr & logFC < 0),
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

set.seed(cfg$seed)
fgsea_raw <- map(ranked, \(stats) {
  fgsea::fgsea(sets, stats, minSize = cfg$min_measured, maxSize = cfg$max_set_size)
})

# The control compares two legs of one person before either was trained, so it cannot contain
# signal. Whatever it calls significant is the cost of fgsea's independence assumption. Read
# every other count against this number, not against zero.
control_first <- map(fgsea_raw, \(res) {
  tibble(
    tested = nrow(res), significant = sum(res$padj < cfg$set_fdr),
    min_padj = signif(min(res$padj), 3)
  )
}) |>
  list_rbind(names_to = "contrast")
message("negative control is ", d$negative_control, ":")
print(control_first)

# Ten collections overlap, so glycolysis is tested in Hallmark, KEGG, Reactome, WikiPathways
# and several times over in GO. collapsePathways re-runs each significant set conditioned on a
# more significant one's leading edge and keeps it only if it still stands alone. It reads the
# data, not an overlap cutoff or a name, and touches only significant sets.
main_sets <- map(set_names(contrast_names), function(contrast) {
  significant <- fgsea_raw[[contrast]][padj < cfg$set_fdr][order(pval)]
  if (nrow(significant) < 2) {
    return(significant$pathway)
  }
  fgsea::collapsePathways(
    significant, sets, ranked[[contrast]],
    pval.threshold = cfg$collapse_pval
  )$mainPathways
})

fgsea_results <- imap(fgsea_raw, function(res, this_contrast) {
  as.data.frame(res) |>
    mutate(
      contrast = this_contrast, main = pathway %in% main_sets[[this_contrast]], .before = 1
    )
}) |>
  list_rbind() |>
  rename(set_id = pathway) |>
  left_join(
    select(gs$set_catalog, set_id, database, pathway, theme, source_size, description),
    by = "set_id"
  ) |>
  relocate(contrast, set_id, database, pathway, theme)

set_summary <- fgsea_results |>
  summarise(
    significant = sum(padj < cfg$set_fdr),
    after_collapse = sum(padj < cfg$set_fdr & main), .by = contrast
  )
print(set_summary)

themed_hits <- fgsea_results |>
  filter(padj < cfg$set_fdr, main) |>
  mutate(theme = coalesce(theme, paste0("(", database, ", no hierarchy)"))) |>
  count(contrast, theme, sort = TRUE, name = "sets")

packages <- c("here", "limma", "fgsea", "dplyr", "purrr", "data.table")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
flat <- mutate(fgsea_results, leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

saveRDS(
  list(
    protein_results = protein_results, fgsea_results = fgsea_results,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config = cfg,
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "fgsea.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    set_summary = set_summary,
    negative_control = control_first,
    fgsea_significant = filter(flat, padj < cfg$set_fdr),
    themed_hits = themed_hits,
    protein_summary = protein_summary,
    protein_results = protein_results,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "01_run_fgsea.xlsx")
)
readr::write_csv(flat, file.path(out, "fgsea_results.csv"))
message("wrote fgsea.rds, 01_run_fgsea.xlsx and fgsea_results.csv")
