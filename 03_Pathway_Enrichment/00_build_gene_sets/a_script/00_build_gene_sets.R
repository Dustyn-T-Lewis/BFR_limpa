# Freeze one MSigDB release, map the protein matrix onto gene symbols, and keep the sets
# large enough to test. Writes the set list that 01_run_fgsea and 02_run_singscore both read,
# so membership is decided once.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
})
here::i_am("config.yml")

cfg <- yaml::read_yaml(here("config.yml"))
cfg$collections <- unlist(cfg$collections)

out <- here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data")
cache_dir <- file.path(out, "cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(proteins = "01_Preprocess/02_Quantification/c_data/proteins.rds")
paths <- map_chr(inputs, here::here)
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
collection_specs <- list(
  Hallmark = list("H", NULL),
  KEGG_Legacy = list("C2", "CP:KEGG_LEGACY"),
  KEGG_Medicus = list("C2", "CP:KEGG_MEDICUS"),
  Reactome = list("C2", "CP:REACTOME"),
  WikiPathways = list("C2", "CP:WIKIPATHWAYS"),
  PID = list("C2", "CP:PID"),
  BioCarta = list("C2", "CP:BIOCARTA"),
  GOBP = list("C5", "GO:BP"),
  GOCC = list("C5", "GO:CC"),
  GOMF = list("C5", "GO:MF")
)
stopifnot(all(cfg$collections %in% names(collection_specs)))
cache_file <- file.path(cache_dir, paste0(
  "msigdb_", cfg$msigdb_version, "_", paste(cfg$collections, collapse = "-"), ".rds"
))
checksum_file <- paste0(cache_file, ".md5")

if (!file.exists(cache_file)) {
  message("fetching ", length(cfg$collections), " collections from msigdbr")
  fetched <- imap(collection_specs[cfg$collections], function(spec, database) {
    df <- msigdbr::msigdbr(
      db_species = "HS", species = "Homo sapiens",
      collection = spec[[1]], subcollection = spec[[2]]
    )
    if (!identical(unique(df$db_version), cfg$msigdb_version)) {
      stop(
        "Requested MSigDB ", cfg$msigdb_version, " but msigdbr returned ",
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
    db_version = cfg$msigdb_version, collections = cfg$collections,
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
stopifnot(identical(frozen$collections, cfg$collections))
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
    qualifies = source_size >= cfg$min_set_size & source_size <= cfg$max_set_size &
      measured_size >= cfg$min_measured
  )
stopifnot(!anyDuplicated(set_catalog$set_id))
sets <- sets_measured[set_catalog$set_id[set_catalog$qualifies]]
if (!length(sets)) stop("No gene sets passed the size filters.")

# A significant GO result is usually a dozen nested terms saying one thing. Each GO set takes
# its nearest GO.db ancestor holding at least theme_min_size genes as a theme. Grouping only:
# every set is still tested alone, and non-GO collections get no theme.
go_ancestors <- list(
  GOBP = GO.db::GOBPANCESTOR, GOCC = GO.db::GOCCANCESTOR, GOMF = GO.db::GOMFANCESTOR
)
go_size <- membership |>
  filter(database %in% names(go_ancestors)) |>
  distinct(source_id, gene) |>
  count(source_id) |>
  with(set_names(n, source_id))

themes <- imap(
  go_ancestors[intersect(names(go_ancestors), cfg$collections)],
  function(ontology, database) {
    ids <- unique(set_catalog$source_id[set_catalog$database == database])
    ids <- ids[!is.na(ids) & ids %in% AnnotationDbi::keys(GO.db::GO.db)]
    tibble(
      database = database, source_id = ids,
      theme_id = map_chr(AnnotationDbi::mget(ids, ontology, ifnotfound = NA), function(a) {
        above <- intersect(a, names(go_size))
        above <- above[go_size[above] >= cfg$theme_min_size]
        if (length(above)) names(which.min(go_size[above])) else NA_character_
      })
    )
  }
) |>
  list_rbind() |>
  mutate(theme = map_chr(theme_id, \(id) {
    if (is.na(id)) NA_character_ else AnnotationDbi::Term(GO.db::GOTERM[[id]])
  }))
set_catalog <- left_join(set_catalog, themes, by = c("database", "source_id"))

collection_summary <- set_catalog |>
  summarise(
    in_msigdb = n(), qualifying = sum(qualifies),
    median_measured = median(measured_size[qualifies]), .by = database
  ) |>
  arrange(match(database, cfg$collections))
theme_summary <- set_catalog |>
  filter(qualifies, !is.na(theme)) |>
  count(theme, sort = TRUE, name = "qualifying_sets")
print(collection_summary)
message("qualifying sets: ", length(sets))

packages <- c("here", "limpa", "msigdbr", "GO.db", "dplyr", "purrr")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
saveRDS(
  list(
    sets = sets, set_catalog = set_catalog, protein_map = protein_map,
    gene_universe = gene_universe,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config = cfg,
      inputs = manifest, packages = versions,
      msigdb_cache_md5 = unname(tools::md5sum(cache_file))
    )
  ),
  file.path(out, "gene_sets.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    collection_summary = collection_summary,
    set_catalog = set_catalog,
    theme_summary = theme_summary,
    protein_gene_map = protein_map,
    mapping_summary = mapping_summary,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "00_build_gene_sets.xlsx")
)
message("wrote gene_sets.rds and 00_build_gene_sets.xlsx")
