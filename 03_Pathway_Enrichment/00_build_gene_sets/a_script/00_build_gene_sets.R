# Freeze one MSigDB release, map the protein matrix onto gene symbols, and keep the sets
# large enough to test. Writes the set list that both 01_run_fgsea_and_fry and 04_run_singscore
# read, so membership is decided once.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
})


out <- here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data")
cache_dir <- file.path(out, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(proteins = "01_Preprocess/02_Quantification/c_data/proteins.rds")
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_Preprocess first. Missing: ", paste(inputs[!file.exists(paths)], collapse = ", "))
}
proteins <- readRDS(paths[["proteins"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

# Membership shifts between MSigDB releases, so one release is pinned. The first run fetches
# it and writes a snapshot with an md5; later runs read and verify that snapshot, needing
# neither the network nor msigdbr. The filename carries the release and the collections.
msigdb_release <- "2026.1.Hs"
collection_specs <- list(
  Hallmark = list("H", NULL),
  Reactome = list("C2", "CP:REACTOME"),
  KEGG_Legacy = list("C2", "CP:KEGG_LEGACY"),
  GOBP = list("C5", "GO:BP")
)
collections <- names(collection_specs)
cache_file <- file.path(cache_dir, paste0(
  "msigdb_", msigdb_release, "_", paste(collections, collapse = "-"), ".rds"
))
checksum_file <- paste0(cache_file, ".md5")

if (!file.exists(cache_file)) {
  message("fetching ", length(collections), " collections from msigdbr")
  fetched <- imap(collection_specs[collections], function(spec, database) {
    df <- msigdbr::msigdbr(
      db_species = "HS", species = "Homo sapiens",
      collection = spec[[1]], subcollection = spec[[2]]
    )
    if (!identical(unique(df$db_version), msigdb_release)) {
      stop(
        "Requested MSigDB ", msigdb_release, " but msigdbr returned ",
        paste(unique(df$db_version), collapse = ", "), "."
      )
    }
    df |>
      transmute(
        database,
        pathway = gs_name, set_id = paste(database, gs_name, sep = "|"),
        gene = gene_symbol, source_id = gs_exact_source, description = gs_description
      ) |>
      filter(!is.na(gene), nzchar(gene)) |>
      distinct()
  })
  frozen <- list(
    db_version = msigdb_release, collections = collections,
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    msigdbr_version = as.character(packageVersion("msigdbr")),
    membership = list_rbind(fetched)
  )
  # Write a complete snapshot before exposing the final cache filename.
  temporary <- tempfile(tmpdir = cache_dir, fileext = ".rds")
  saveRDS(frozen, temporary, compress = "xz")
  writeLines(unname(tools::md5sum(temporary)), checksum_file)
  if (!file.rename(temporary, cache_file)) stop("Could not save MSigDB cache.")
}
if (!file.exists(checksum_file) || !identical(
  unname(tools::md5sum(cache_file)), readLines(checksum_file, warn = FALSE)
)) {
  stop("MSigDB cache checksum missing or mismatched. Restore the RDS and its .md5 together.")
}
frozen <- readRDS(cache_file)
membership <- frozen$membership
stopifnot(identical(frozen$collections, collections))
message(
  "frozen: ", frozen$db_version, ", ", n_distinct(membership$set_id), " sets, ",
  "cached ", frozen$created_utc
)

# Set tests need one row per gene. A row with no symbol, or several, is dropped rather than
# split, because splitting invents measurements nobody made. Where rows share a symbol, the
# one with the most observed precursors represents it. Decided once, reading no fold change.
protein_map <- proteins$genes |>
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
representatives <- protein_map |>
  filter(mapping_status == "candidate") |>
  arrange(gene, desc(mean_obs), desc(NPrec), protein) |>
  distinct(gene, .keep_all = TRUE) |>
  pull(protein)
protein_map <- protein_map |>
  mutate(
    selected = protein %in% representatives,
    mapping_status = case_when(
      selected ~ "representative",
      mapping_status == "candidate" ~ "duplicate_gene",
      TRUE ~ mapping_status
    ),
    # Accessions distinguish duplicate symbols in protein-level plot labels.
    label = if_else(
      is.na(gene), protein,
      if_else(duplicated(gene) | duplicated(gene, fromLast = TRUE),
        paste0(gene, " (", protein, ")"), gene
      )
    )
  )
gene_map <- filter(protein_map, selected)
gene_universe <- gene_map$gene
# Labels are what the volcanoes print, and a collision would put two proteins on one point.
stopifnot(!anyDuplicated(protein_map$label))
mapping_summary <- count(protein_map, mapping_status, name = "proteins")
print(mapping_summary)
message("measured gene universe: ", length(gene_universe))

# Size is the only pre-test filter. The measured bar counts proteins this experiment detected,
# which is what governs power. Nothing is dropped for overlapping another set or for its name.
set_members <- distinct(membership, set_id, gene)
sets_full <- split(set_members$gene, set_members$set_id)
sets_measured <- map(sets_full, intersect, y = gene_universe)

set_catalog <- membership |>
  distinct(set_id, database, pathway, source_id, description) |>
  mutate(
    source_size = lengths(sets_full)[set_id],
    measured_size = lengths(sets_measured)[set_id],
    qualifies = source_size >= 15 & source_size <= 500 &
      measured_size >= 15
  )
stopifnot(!anyDuplicated(set_catalog$set_id))
sets <- sets_measured[set_catalog$set_id[set_catalog$qualifies]]
if (!length(sets)) stop("No gene sets passed the size filters.")

# GO Slim sets come from the GO Consortium's generic slim, a published 140-term list, frozen with
# an md5 like the MSigDB snapshot.
slim_file <- file.path(cache_dir, "goslim_generic.obo")
stopifnot(
  file.exists(slim_file),
  identical(unname(tools::md5sum(slim_file)), readLines(paste0(slim_file, ".md5"), warn = FALSE))
)
slim_offspring <- AnnotationDbi::mget(
  GSEABase::ids(GSEABase::getOBOCollection(slim_file)), GO.db::GOBPOFFSPRING,
  ifnotfound = NA
)
slim_offspring <- slim_offspring[!is.na(slim_offspring)]

# Each slim term is also a set: every measured gene annotated to it or to any term beneath it.
# Built from the frozen membership rather than from the qualifying subset, because a gene must
# not be lost when the GO:BP set carrying it falls outside the size filter. A slim term is broad
# by design, so the 15-to-500 rule is read on measured size, which is the size that governs
# whether the set is testable here.
go_genes <- membership |>
  filter(database == "GOBP") |>
  with(split(gene, source_id))
slim_sets <- imap(slim_offspring, function(descendants, slim_id) {
  covered <- intersect(c(slim_id, descendants), names(go_genes))
  sort(intersect(unique(unlist(go_genes[covered], use.names = FALSE)), gene_universe))
})
slim_catalog <- tibble(
  theme_id = names(slim_sets),
  pathway = unname(AnnotationDbi::Term(GO.db::GOTERM[theme_id])),
  measured_size = lengths(slim_sets)
) |>
  transmute(
    set_id = paste("GO_Slim", toupper(gsub("[^A-Za-z0-9]+", "_", pathway)), sep = "|"),
    database = "GO_Slim", pathway, source_id = theme_id,
    description = "GO Slim term: every measured gene under it in the GO:BP hierarchy",
    source_size = measured_size, measured_size,
    qualifies = measured_size >= 15 & measured_size <= 500
  )
names(slim_sets) <- slim_catalog$set_id
set_catalog <- bind_rows(set_catalog, slim_catalog)
sets <- c(sets, slim_sets[slim_catalog$set_id[slim_catalog$qualifies]])
stopifnot(!anyDuplicated(set_catalog$set_id), !anyDuplicated(names(sets)))
message("GO Slim sets: ", sum(slim_catalog$qualifies), " of ", nrow(slim_catalog), " testable")

collection_summary <- set_catalog |>
  summarise(
    in_msigdb = n(), qualifying = sum(qualifies),
    median_measured = median(measured_size[qualifies]), .by = database
  ) |>
  arrange(match(database, c(collections, "GO_Slim")))
print(collection_summary)
message("qualifying sets: ", length(sets))

packages <- c("here", "limpa", "msigdbr", "GO.db", "GSEABase", "dplyr", "purrr")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    sets = sets, set_catalog = set_catalog, protein_map = protein_map,
    gene_universe = gene_universe,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, packages = versions,
      msigdb_cache_md5 = unname(tools::md5sum(cache_file)),
      goslim_md5 = unname(tools::md5sum(slim_file))
    )
  ),
  file.path(out, "gene_sets.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    collection_summary = collection_summary,
    set_catalog = set_catalog,
    protein_gene_map = protein_map,
    mapping_summary = mapping_summary,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "00_build_gene_sets.xlsx")
)
message("wrote gene_sets.rds and 00_build_gene_sets.xlsx")
