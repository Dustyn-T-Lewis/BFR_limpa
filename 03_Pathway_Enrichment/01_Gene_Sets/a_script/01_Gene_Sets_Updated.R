# Parameters declared once at the top of the file
# Purpose of this script: Prevent reproducibility failure when running the analysis
# in the future. Databases that store gene sets tend to change frequently. These gene sets
# must be frozen by you so that when you re run the same gene sets and do statistics on it
# in the future, you are pulling out the same data from the same gene set and NOT the updated one.
params <- list(
  msigdb_version = "2026.1.Hs",
  collections    = c("Hallmark", "Reactome", "GOBP"), # The three collections
  min_set_size   = 15, # A pathway has to be a sensible size on paper — 15 to 500 genes
  max_set_size   = 500, # A pathway has to be a sensible size on paper — 15 to 500 genes
  min_measured   = 15, # we have to have measured at least 15 of its members
  jaccard_cutoff = 0.5,# Because gene ontologies are hierarchical and lots of pathways are near-copies of each other, we drop any pathway that overlaps more than half its genes with one we've already kept.
  protein_fdr    = 0.05,
  label_n        = 5
)

required <- c(
  "here", "limma", "limpa", "dplyr", "tibble", "tidyr", "purrr", "readr",
  "writexl", "digest", "Matrix", "ggplot2", "knitr"
)

library(here)
library(limpa)
library(limma)
library(dplyr)
library(tibble)
library(readr)
library(msigdbr)
library(Matrix)
library(purrr)
library(ggplot2)

# Pulling out data to use for downstream analysis 
here::i_am("03_Pathway_Enrichment/01_Gene_Sets/a_script/01_Gene_Sets_Updated.R")
stage <- here("03_Pathway_Enrichment", "01_Gene_Sets")
out <- file.path(stage, "c_data")
cache_dir <- file.path(out, "cache")
figure_dir <- file.path(stage, "b_reports", "figures")
for (path in c(out, cache_dir, figure_dir)) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
}

# Loading the three objects written by three earlier stages
input_paths <- c(
  proteins = here("01_Preprocess/02_Quantification/c_data/proteins.rds"), #proteins is an EList: $E is a 3,042 × 131 matrix of log₂ protein abundances
  fit = here("02_Differential_Expression/02_Differential/c_data/fit.rds"), # fit is the fitted limma model, its coefficients are estimated contrast effects, one row per protein, one column per contrast.
  design = here("02_Differential_Expression/01_Design/c_data/design.rds") # d holds the design matrix and the contrast matrix.
)
if (!all(file.exists(input_paths))) {
  stop(
    "Run the upstream stages first. Missing: ",
    paste(input_paths[!file.exists(input_paths)], collapse = ", ")
  )
}
# Recall built-in R function used to restore a single R object saved in RDS
proteins <- readRDS(input_paths[["proteins"]])
fit <- readRDS(input_paths[["fit"]])
d <- readRDS(input_paths[["design"]])
protein_ids <- rownames(proteins$E) # Pull out the protein_ids from the rows proteins$E and name it protein_ids tp store it as a vector
contrast_names <- colnames(d$contrasts) # Pull out the contrast names from columns of d$contrasts and name it protein_ids tp store it as a vector
training_contrasts <- c("BFR_post_vs_pre", "HLRT_post_vs_pre")
# limma compares only counts
# If design has 131 rows and E has 131 columns, limma proceeds — even if row 7 of the design is participant 12's BFR pre-sample while column 7 of the matrix
# is participant 4's HLRT post-sample. Every downstream number would be computed without error, and every one would be wrong.
# These code lines safeguard the downstream analyses
stopifnot(
  !anyDuplicated(protein_ids), !anyDuplicated(colnames(proteins$E)),
  identical(rownames(proteins$genes), protein_ids),
  identical(rownames(fit$coefficients), protein_ids),
  identical(colnames(proteins$E), rownames(d$design)),
  identical(colnames(proteins$E), proteins$targets$sample_id),
  identical(colnames(fit$coefficients), contrast_names),
  identical(fit$design, d$design), identical(fit$contrasts, d$contrasts),
  identical(fit$EList$E, proteins$E),
  identical(fit$EList$other, proteins$other),
  identical(fit$EList$genes, proteins$genes),
  all(is.finite(proteins$E)),
  identical(dimnames(proteins$other$standard.error), dimnames(proteins$E)),
  identical(dimnames(proteins$other$n.observations), dimnames(proteins$E)),
  all(is.finite(proteins$other$standard.error)),
  all(proteins$other$standard.error > 0),
  d$negative_control %in% contrast_names,
  all(training_contrasts %in% contrast_names)
)

stopifnot(
  is.numeric(params$min_set_size), params$min_set_size >= 1,
  params$min_set_size == as.integer(params$min_set_size),
  params$max_set_size >= params$min_set_size,
  params$min_measured >= 1,
  params$min_measured == as.integer(params$min_measured),
  params$min_measured <= params$max_set_size,
  params$jaccard_cutoff > 0, params$jaccard_cutoff <= 1,
  params$protein_fdr > 0, params$protein_fdr < 1,
  params$label_n >= 0, params$label_n == as.integer(params$label_n))

tibble(
  n_proteins = nrow(proteins$E), n_samples = ncol(proteins$E),
  n_participants = n_distinct(proteins$targets$participant),
  n_contrasts = length(contrast_names)
)

# Freezing Gene Sets
# Defines a named list of configuration parameters specifying five distinct gene set collections for enrichment analysis.
# Basically ordering what pathways will be looked at initially.
collection_specs <- list(
  Hallmark = list(collection = "H", subcollection = NULL), # Hallmaek: 50 curated gene sets representing well-defined biological states or processes.
  Reactome = list(collection = "C2", subcollection = "CP:REACTOME"), # CP:Reactome = Curated canonical pathways sourced specifically from the Reactome database.
  GOBP = list(collection = "C5", subcollection = "GO:BP"), # BO:BP = Gene Ontology Biological Process ontology gene sets.
  GOCC = list(collection = "C5", subcollection = "GO:CC"), # GO:CC = Gene Ontology Cellular Component ontology gene sets.
  GOMF = list(collection = "C5", subcollection = "GO:MF") # GO:MF = Gene Ontology Molecular Function ontology gene sets.
)
# Place the order
collections <- unlist(params$collections, use.names = FALSE)
stopifnot(
  length(collections) > 0, !anyDuplicated(collections),
  all(collections %in% names(collection_specs)),
  grepl("^[0-9]{4}\\.[0-9]+\\.Hs$", params$msigdb_version)
)
# Write the order down
cache_request <- list(
  schema = 1L, species = "Homo sapiens", db_species = "HS",
  db_version = params$msigdb_version,
  collections = collection_specs[collections]
)
# Fingerprint the order and rename the file after it
cache_key <- substr(digest::digest(cache_request, algo = "sha256"), 1, 12)
cache_file <- file.path(cache_dir, paste0(
  "msigdb_", params$msigdb_version, "_", cache_key, ".rds"
))
checksum_file <- paste0(cache_file, ".md5")
# Download
if (!file.exists(cache_file)) {
  if (!requireNamespace("msigdbr", quietly = TRUE)) {
    stop(
      "The frozen cache is absent. Install msigdbr or restore ",
      "the matching RDS and .md5 files under c_data/cache/."
    )
  }
  fetched <- imap(collection_specs[collections], function(spec, database) {
    df <- msigdbr::msigdbr(
      db_species = "HS", species = "Homo sapiens",
      collection = spec$collection, subcollection = spec$subcollection
    )
    if (!identical(unique(df$db_version), params$msigdb_version)) {
      stop(
        "Requested MSigDB ", params$msigdb_version,
        " but msigdbr returned ", paste(unique(df$db_version), collapse = ", "),
        ". Restore the requested snapshot or explicitly change the version."
      )
    }
    df |>
      transmute(
        database,
        pathway = gs_name, set_id = paste(database, gs_name, sep = "|"),
        gene = gene_symbol, msigdb_id = gs_id, source_id = gs_exact_source,
        description = gs_description, url = gs_url, db_version
      ) |>
      filter(!is.na(gene), nzchar(gene)) |>
      distinct()
  })
  # Sealing the list
  frozen <- list(
    request = cache_request,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    msigdbr_version = as.character(packageVersion("msigdbr")),
    source = "https://www.gsea-msigdb.org/gsea/msigdb/",
    membership = list_rbind(fetched) # list_rbind(fetched) stacks the three per-collection tables into one.
  )
  # Write to a scratch name first, then rename.
  temporary <- tempfile(tmpdir = cache_dir, fileext = ".rds")
  saveRDS(frozen, temporary, compress = "xz")
  checksum <- unname(tools::md5sum(temporary))
  if (!file.rename(temporary, cache_file)) stop("Could not save MSigDB cache.")
  writeLines(checksum, checksum_file)
}
# Check the seal every single run
if (!file.exists(checksum_file) || !identical(
  unname(tools::md5sum(cache_file)), readLines(checksum_file, warn = FALSE)
)) {
  stop(
    "MSigDB cache checksum is missing or mismatched. Restore the RDS ",
    "and its matching .md5; do not reuse an unverified snapshot."
  )
}
# Load and confirm you got what you ordered
frozen <- readRDS(cache_file)
stopifnot(identical(frozen$request, cache_request))
membership <- frozen$membership
stopifnot(
  identical(unique(membership$db_version), params$msigdb_version),
  setequal(unique(membership$database), collections),
  !anyNA(membership$gene), all(nzchar(membership$gene))
)
tibble(
  database_release = frozen$request$db_version,
  msigdbr_at_freeze = frozen$msigdbr_version,
  frozen_utc = frozen$created_utc,
  md5 = unname(tools::md5sum(cache_file)))

# Section 3: One protein One Gene
# Mass spec produces accessions (a protein identifier) whereas the pathway datasets provide gene symbols
# This section works to connect the two data sets on a single base.
# Gene mapping, one protein per gene
# Set databases are keyed by gene symbol. Your matrix is keyed by protein-group accession. Something must bridge them, and the bridge is not one-to-one in either direction.
protein_map <- proteins$genes |> # a bookkeeping table with one row per protein, carrying a mapping_status label saying which of the four situations that protein was in
  rownames_to_column("protein") |>
  mutate(
    mean_obs = rowMeans(proteins$other$n.observations),
    gene = trimws(Genes),
    mapping_status = case_when(
      is.na(gene) | !nzchar(gene) ~ "missing_symbol",
      grepl("[;,|]", gene) ~ "multiple_symbols",
      TRUE ~ "candidate"
    ),
    gene = if_else(mapping_status == "candidate", gene, NA_character_)
  )
representatives <- protein_map |> # Where several proteins belong and claim one gene, pick one: based on mean_obs. If two tie on mean_obs, the one with more precursors wins. If they still tie, choose based off alphabetical accession.
  filter(mapping_status == "candidate") |>
  arrange(gene, desc(mean_obs), desc(NPrec), protein) |>
  distinct(gene, .keep_all = TRUE) |> # keeps the first row per gene and drops the rest
  pull(protein)
protein_map <- protein_map |> # write down what happened to every protein
  mutate(
    selected = protein %in% representatives,
    mapping_status = case_when(
      selected ~ "representative", # representative	= clean symbol, best evidence for that gene
      mapping_status == "candidate" ~ "duplicate_gene", # duplicate_gene = clean symbol, but another row represents this gene
      TRUE ~ mapping_status
    ),
    annotated = !is.na(gene) & gene %in% membership$gene,
    # Accessions distinguish duplicate symbols in protein-level plot labels.
    label = if_else(
      is.na(gene), protein,
      if_else(duplicated(gene) | duplicated(gene, fromLast = TRUE),
        paste0(gene, " (", protein, ")"), gene
      )
    )
  )
# Build the gene-level matrix
# 3,042 proteins collapse to about 3,040 symbols
gene_map <- filter(protein_map, selected)
gene_universe <- gene_map$gene # gene_universe is the headline output your list of measured gene symbols, each appearing exactly once.
stopifnot(!anyDuplicated(gene_universe), !anyDuplicated(protein_map$label))
gene_matrix <- proteins$E[gene_map$protein, , drop = FALSE] # gene_matrix takes your original abundance matrix, keeps only the winning rows, and relabels them with symbols instead of accessions
rownames(gene_matrix) <- gene_universe
count(protein_map, mapping_status, annotated, name = "proteins")

# ----- Section 4: Filtering ------
set_members <- membership |> distinct(set_id, gene) # One row per pathway, with three sizes and a verdict.
sets_full <- split(set_members$gene, set_members$set_id) |> map(sort) #split() reshapes the data from a long list of 'membership' into to a named list, which is the natural shape for "a pathway is a bag of genes"
sets_measured <- map(sets_full, intersect, y = gene_universe) # intersect: to a named list, which is the natural shape for "a pathway is a bag of genes"

set_catalog <- membership |>
  distinct(set_id, database, pathway, msigdb_id, source_id, description, url) |> # lengths() (note the "s") counts how many items are in each bag of the list
  mutate(
    source_size = lengths(sets_full)[set_id],
    measured_size = lengths(sets_measured)[set_id],
    measured_fraction = measured_size / source_size, # 74/200 = 0.87 means you saw most of the pathway; 11/300 = 0.037
    status = case_when(
      source_size < params$min_set_size ~ "source_too_small",
      source_size > params$max_set_size ~ "source_too_large",
      measured_size < params$min_measured ~ "too_few_measured",
      TRUE ~ "eligible"
    ),
    representative = NA_character_, jaccard = NA_real_
  )

stopifnot(!anyDuplicated(set_catalog$set_id))
eligible_ids <- set_catalog$set_id[set_catalog$status == "eligible"]
if (!length(eligible_ids)) stop("No gene sets passed the size/coverage filters.")

reduce_sets <- function(catalog) {
  catalog <- arrange(catalog, desc(measured_size), desc(source_size), set_id)
  sets <- sets_measured[catalog$set_id]
  incidence <- Matrix::sparseMatrix(
    i = match(unlist(sets, use.names = FALSE), gene_universe),
    j = rep(seq_along(sets), lengths(sets)), x = 1,
    dims = c(length(gene_universe), length(sets))
  )
  intersections <- Matrix::crossprod(incidence)
  kept <- integer()
  for (i in seq_along(sets)) {
    similarity <- numeric()
    if (length(kept)) {
      shared <- as.numeric(intersections[i, kept])
      similarity <- shared / (lengths(sets)[i] + lengths(sets)[kept] - shared)
    }
    if (length(similarity) && max(similarity) >= params$jaccard_cutoff) {
      j <- kept[which.max(similarity)]
      catalog$status[i] <- "redundant"
      catalog$representative[i] <- catalog$set_id[j]
      catalog$jaccard[i] <- max(similarity)
    } else {
      kept <- c(kept, i)
      catalog$status[i] <- "retained"
    }
  }
  catalog
}
reduced_catalog <- set_catalog |>
  filter(status == "eligible") |>
  split(~database) |>
  map(reduce_sets) |>
  list_rbind()
set_catalog <- bind_rows(
  filter(set_catalog, status != "eligible"), reduced_catalog
) |>
  arrange(match(database, collections), set_id)
sets_eligible <- sets_measured[eligible_ids]
sets_retained <- sets_measured[set_catalog$set_id[set_catalog$status == "retained"]]
set_indices <- map(sets_retained, function(genes) {
  match(gene_map$protein[match(genes, gene_map$gene)], protein_ids)
})
annotated_universes <- setNames(map(collections, function(database) {
  intersect(gene_universe, membership$gene[membership$database == database])
}), collections)
set_summary <- set_catalog |>
  summarise(
    downloaded = n(),
    source_size_pass = sum(
      source_size >= params$min_set_size & source_size <= params$max_set_size
    ),
    measured_size_pass = sum(status %in% c("retained", "redundant")),
    retained = sum(status == "retained"), .by = database
  )
set_summary

# --- Protein Statistics ---------
protein_results <- map(contrast_names, function(contrast) {
  limma::topTable(fit,
    coef = contrast, number = Inf,
    adjust.method = "BH", sort.by = "none"
  ) |>
    rownames_to_column("protein") |>
    select(protein, logFC, AveExpr, t, P.Value, adj.P.Val, B) |>
    mutate(contrast = contrast, .before = 1)
}) |>
  list_rbind() |>
  left_join(protein_map, by = "protein", relationship = "many-to-one") |>
  mutate(
    padj = adj.P.Val,
    # Keep exact p-values in the export; bound only the logarithm for plotting.
    plot_p = pmax(P.Value, .Machine$double.xmin),
    pi_score = P.Value^abs(logFC),
    pi_value = if_else(logFC == 0, 0, abs(logFC) * -log10(plot_p)),
    signed_pi = sign(logFC) * pi_value,
    fdr_significant = padj < params$protein_fdr,
    direction = case_when(
      fdr_significant & logFC > 0 ~ "up",
      fdr_significant & logFC < 0 ~ "down",
      TRUE ~ "not_significant"
    )
  )
stopifnot(
  nrow(protein_results) == length(protein_ids) * length(contrast_names),
  !anyDuplicated(protein_results[c("contrast", "protein")]),
  all(is.finite(protein_results$logFC)),
  all(is.finite(protein_results$t)),
  all(is.finite(protein_results$P.Value)),
  all(protein_results$P.Value >= 0 & protein_results$P.Value <= 1),
  all(is.finite(protein_results$padj)),
  all(protein_results$padj >= 0 & protein_results$padj <= 1)
)
for (contrast in contrast_names) {
  result <- filter(protein_results, .data$contrast == .env$contrast)
  stopifnot(isTRUE(all.equal(result$padj, p.adjust(result$P.Value, "BH"))))
}
gene_results <- filter(protein_results, selected)
ranked_t <- setNames(map(contrast_names, function(contrast) {
  result <- filter(gene_results, .data$contrast == .env$contrast)
  sort(setNames(result$t, result$gene), decreasing = TRUE)
}), contrast_names)
ranked_pi <- setNames(map(training_contrasts, function(contrast) {
  result <- filter(gene_results, .data$contrast == .env$contrast)
  sort(setNames(result$signed_pi, result$gene), decreasing = TRUE)
}), training_contrasts)
contrast_summary <- protein_results |>
  summarise(
    tested = n(), fdr_up = sum(direction == "up"),
    fdr_down = sum(direction == "down"), min_fdr = min(padj),
    zero_p_values = sum(P.Value == 0), .by = contrast
  )
contrast_summary |> mutate(min_fdr = signif(min_fdr, 3))

input_manifest <- tibble(
  input = names(input_paths),
  path = substring(unname(input_paths), nchar(here()) + 2),
  md5 = unname(tools::md5sum(input_paths))
)
package_versions <- tibble(
  package = required,
  version = vapply(required, function(package) {
    as.character(packageVersion(package))
  }, character(1))
)
provenance <- list(
  schema = 1L, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  inputs = input_manifest, parameters = params,
  script_md5 = unname(tools::md5sum(file.path(
    stage, "a_script", "01_gene_sets.qmd"
  ))),
  cache_file = file.path("cache", basename(cache_file)),
  cache_md5 = unname(tools::md5sum(cache_file)),
  cache_request = frozen$request, packages = package_versions,
  mapping = paste(
    "Exact single symbols; highest mean observed precursors,",
    "then NPrec, then accession; independent of differential statistics."
  ),
  fdr = "BH across all tested proteins separately within each contrast",
  pi_score = "P.Value^abs(logFC); smaller ranks higher; no error control",
  pi_value = paste(
    "abs(logFC) * -log10(max(P.Value, double.xmin));",
    "zero logFC gives zero; larger ranks higher; plotting bound only"
  ),
  signed_pi = "sign(logFC) * pi_value; exploratory training ranks only",
  reduction = "Within database, measured-member Jaccard; deterministic greedy",
  pathway_tests_run = FALSE,
  plotting_adapter = "Private empty-enrichment validation adapter; namespace unchanged"
)
gene_sets <- list(
  full = sets_full, eligible = sets_eligible, retained = sets_retained,
  catalog = set_catalog, measured_universe = gene_universe,
  annotated_universes = annotated_universes, provenance = provenance
)
pathway_inputs <- list(
  proteins = proteins, fit = fit, design = d,
  protein_results = protein_results, gene_results = gene_results,
  protein_map = protein_map, gene_matrix = gene_matrix,
  gene_sets = sets_retained, set_indices = set_indices,
  ranked_t = ranked_t, ranked_pi = ranked_pi,
  annotated_universes = annotated_universes, provenance = provenance
)
stopifnot(
  all(map_lgl(sets_retained, ~ all(.x %in% gene_universe))),
  all(lengths(sets_retained) >= params$min_measured),
  all(map_lgl(set_indices, ~ !anyNA(.x) && !anyDuplicated(.x))),
  all(map2_lgl(sets_retained, set_indices, function(genes, indices) {
    identical(genes, protein_map$gene[indices])
  })),
  identical(colnames(gene_matrix), rownames(d$design))
)
saveRDS(gene_sets, file.path(out, "gene_sets.rds"), compress = "xz")
saveRDS(pathway_inputs, file.path(out, "pathway_inputs.rds"), compress = "xz")
saveRDS(proteins$E, file.path(out, "protein_matrix.rds"), compress = "xz")
readr::write_csv(protein_results, file.path(out, "protein_statistics.csv"))
readr::write_csv(set_catalog, file.path(out, "gene_set_catalog.csv"))
readr::write_csv(input_manifest, file.path(out, "input_manifest.csv"))
matrix_table <- as.data.frame(proteins$E) |> rownames_to_column("protein")
writexl::write_xlsx(list(
  protein_matrix = matrix_table,
  samples = proteins$targets,
  protein_statistics = protein_results,
  gene_statistics = gene_results,
  protein_gene_map = protein_map,
  gene_set_catalog = set_catalog,
  set_summary = set_summary,
  contrast_summary = contrast_summary,
  input_manifest = input_manifest,
  package_versions = package_versions
), file.path(out, "01_gene_sets.xlsx"))
tibble(
  output = c(
    "01_gene_sets.xlsx", "protein_matrix.rds", "gene_sets.rds",
    "pathway_inputs.rds", "protein_statistics.csv", "gene_set_catalog.csv"
  ),
  size_mb = round(file.info(file.path(out, output))$size / 1024^2, 2))
      
