# Sourced by the three stages of 03_Pathway_Enrichment. Every other stage in this repo is
# self-contained, and this file is the one exception: the gene-set loader and the overlap
# reduction are needed by all three, and three copies would drift. Nothing here reads a file
# or writes one, so sourcing it has no side effect beyond defining these four functions.
# The stages source it by path; see this directory's README for the call.

library(dplyr)
library(purrr)

# Gene-set names that describe a disease, an infection or a tumour. Kept here as a default
# rather than hidden in a notebook, because excluding a set by name is a judgement call: a
# set called INFLUENZA_INFECTION is largely ribosome and translation machinery relabelled,
# so dropping it removes real biology under a misleading name. drop_disease_terms() reports
# everything it removes, and passing NULL turns the whole thing off.
DISEASE_PATTERNS <- c(
  "INFECTION", "INFECTIOUS", "DISEASE", "VIRUS", "VIRAL", "INFLUENZA", "HIV",
  "HEPATITIS", "SARS", "COVID", "CANCER", "CARCINOMA", "TUMOR", "TUMOUR",
  "LEUKEMIA", "LYMPHOMA", "MELANOMA", "DIABETES", "ALZHEIMER", "PARKINSON",
  "HUNTINGTON", "AMYOTROPHIC", "PRION"
)

#' Drop sets whose Jaccard overlap with an already-kept set reaches the cutoff.
#'
#' Walks the catalogue from largest measured membership down, keeping a set unless it
#' overlaps something already kept. `catalog` needs `set_id`, `measured_size`, `source_size`,
#' `status`, `representative` and `jaccard`; `sets` is a named list of measured members keyed
#' by `set_id`. Returns the catalogue with those last three columns filled in.
reduce_sets <- function(catalog, sets, universe, cutoff) {
  catalog <- arrange(catalog, desc(measured_size), desc(source_size), set_id)
  sets <- sets[catalog$set_id]
  sizes <- lengths(sets)
  # Column j marks the measured genes in set j, so crossprod gives every pairwise
  # intersection count in one pass. Indexing an empty `kept` yields an empty Jaccard
  # vector, which is how the first set takes the else branch.
  incidence <- vapply(sets, \(genes) universe %in% genes, logical(length(universe)))
  shared <- crossprod(incidence)
  kept <- integer()
  for (i in seq_along(sets)) {
    jaccard <- shared[i, kept] / (sizes[i] + sizes[kept] - shared[i, kept])
    if (length(jaccard) && max(jaccard) >= cutoff) {
      catalog$status[i] <- "redundant"
      catalog$representative[i] <- catalog$set_id[kept[which.max(jaccard)]]
      catalog$jaccard[i] <- max(jaccard)
    } else {
      kept <- c(kept, i)
      catalog$status[i] <- "retained"
    }
  }
  catalog
}

#' Load one or more databases from gene_sets.rds as a single named list.
#'
#' `scope` decides what the overlap reduction compares. "within" reduces each database
#' separately, so a GO term cannot displace a Hallmark set covering the same biology; that is
#' what you want when reporting per database. "across" pools the requested databases first, so
#' a Hallmark set and its GO twin compete and the duplicate story survives once; that is what
#' you want when several databases feed one figure. `dedup = FALSE` returns the raw union.
#'
#' Returns a list of `sets` (named list of character vectors) and `catalog` (one row per set,
#' recording which set displaced each casualty and by how much).
load_sets <- function(gene_sets, dbs, dedup = TRUE, scope = c("across", "within"),
                      cutoff = 0.5) {
  scope <- match.arg(scope)
  stopifnot(all(dbs %in% unique(gene_sets$catalog$database)))

  catalog <- gene_sets$catalog |>
    filter(database %in% dbs, status %in% c("retained", "redundant")) |>
    mutate(status = "eligible", representative = NA_character_, jaccard = NA_real_)
  sets <- gene_sets$eligible[catalog$set_id]

  if (!dedup) {
    return(list(sets = sets, catalog = catalog))
  }
  reduced <- if (scope == "within") {
    catalog |>
      split(~database) |>
      map(\(part) reduce_sets(part, sets, gene_sets$measured_universe, cutoff)) |>
      list_rbind()
  } else {
    reduce_sets(catalog, sets, gene_sets$measured_universe, cutoff)
  }
  # reduce_sets() returns the walk order, which depends on set size. Sort back to database
  # then set_id so the same request always yields the same list in the same order.
  reduced <- arrange(reduced, match(database, dbs), set_id)
  list(
    sets = gene_sets$eligible[reduced$set_id[reduced$status == "retained"]],
    catalog = reduced
  )
}

#' Split a results table into the rows kept and the rows a name pattern caught.
#'
#' `patterns` is matched case-insensitively against `column`. Returns both halves, because the
#' excluded rows are evidence: a reader has to be able to see what a name filter removed and
#' judge whether the removal was fair. `patterns = NULL` keeps everything.
drop_disease_terms <- function(results, patterns = DISEASE_PATTERNS, column = "pathway") {
  if (is.null(patterns) || !length(patterns)) {
    return(list(kept = results, excluded = results[0, , drop = FALSE]))
  }
  matched <- map_chr(results[[column]], \(name) {
    hit <- patterns[vapply(patterns, grepl, logical(1), x = name, ignore.case = TRUE)]
    if (length(hit)) hit[1] else NA_character_
  })
  list(
    kept = results[is.na(matched), , drop = FALSE],
    excluded = mutate(results[!is.na(matched), , drop = FALSE],
      matched_pattern = matched[!is.na(matched)]
    )
  )
}

#' Map gene-symbol sets onto row numbers of the protein matrix, which is what fry() indexes.
#'
#' The route is symbol -> representative accession -> row number, so the result indexes
#' `proteins$E` as 01_Preprocess wrote it, never anything built downstream.
set_row_indices <- function(sets, protein_map, protein_ids) {
  gene_map <- filter(protein_map, selected)
  map(sets, \(genes) {
    match(gene_map$protein[match(genes, gene_map$gene)], protein_ids)
  })
}

#' Resolve a stage's inputs, stop if any is absent, and read them all.
#'
#' Takes repo-relative paths, returns the loaded objects by name plus a `manifest` tibble of
#' path and md5 for provenance. Keeping the relative strings means the manifest needs no path
#' arithmetic later.
stage_inputs <- function(...) {
  files <- c(...)
  paths <- map_chr(files, here::here)
  if (!all(file.exists(paths))) {
    stop(
      "Run the upstream stages first. Missing: ",
      paste(files[!file.exists(paths)], collapse = ", ")
    )
  }
  objects <- map(set_names(paths, names(files)), \(path) {
    if (grepl("\\.csv$", path)) readr::read_csv(path, show_col_types = FALSE) else readRDS(path)
  })
  c(objects, list(
    manifest = tibble::tibble(
      input = names(files), path = unname(files),
      md5 = unname(tools::md5sum(paths))
    )
  ))
}

#' Write a stage's workbook, RDS files and CSVs, and stamp them with one provenance record.
#'
#' `name` is the workbook stem, matching the stage's script. `notes` carries whatever the
#' stage needs a later reader not to misunderstand — what a score means, which test ran, what
#' was excluded. Returns a tibble of what was written.
write_stage_outputs <- function(out, name, sheets, rds = list(), csv = list(),
                                manifest, parameters, packages, notes = list()) {
  dir.create(out, recursive = TRUE, showWarnings = FALSE)
  versions <- tibble::tibble(
    package = packages,
    version = map_chr(packages, \(p) as.character(utils::packageVersion(p)))
  )
  provenance <- c(
    list(
      schema = 2L, created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
      inputs = manifest, parameters = parameters, packages = versions
    ),
    notes
  )
  # Every workbook ends with the same two sheets, so a reader never has to open the RDS to
  # find out which inputs and which package versions produced the numbers above them.
  sheets <- c(sheets, list(input_manifest = manifest, package_versions = versions))
  iwalk(rds, \(object, name) {
    saveRDS(c(object, list(provenance = provenance)), file.path(out, paste0(name, ".rds")),
      compress = "xz"
    )
  })
  iwalk(csv, \(table, name) readr::write_csv(table, file.path(out, paste0(name, ".csv"))))
  writexl::write_xlsx(sheets, file.path(out, paste0(name, ".xlsx")))
  written <- c(paste0(name, ".xlsx"), paste0(names(rds), ".rds"), paste0(names(csv), ".csv"))
  tibble::tibble(
    output = written,
    size_mb = round(file.info(file.path(out, written))$size / 1024^2, 2)
  )
}

#' Run fgsea for every set collection against every contrast.
#'
#' `ranked` is a named list of named numeric vectors, one per contrast; `collections` a named
#' list of set lists. fgseaMultilevel samples internally, so `seed` is what makes a rerun
#' reproduce. Returns one long data.frame with `leadingEdge` still a list column, and
#' `as.data.frame()` applied per run because volcano_ring() cannot take a data.table.
run_fgsea_all <- function(ranked, collections, min_size, max_size, seed = 1) {
  set.seed(seed)
  map(collections, \(sets) {
    map(ranked, \(stats) {
      as.data.frame(fgsea::fgsea(
        pathways = sets, stats = stats, minSize = min_size, maxSize = max_size
      ))
    }) |>
      list_rbind(names_to = "contrast")
  }) |>
    list_rbind(names_to = "collection")
}

#' Run fry for every contrast.
#'
#' fry takes one contrast vector at a time, never a matrix. `weights` carries limpa's
#' per-observation precision into the set test. Pass `block` and `correlation` through `...`
#' only when the residual within-block correlation is far enough from zero to matter.
run_fry_all <- function(expression, index, design, contrasts, weights = NULL, ...) {
  map(set_names(colnames(contrasts)), \(contrast) {
    limma::fry(
      y = expression, index = index, design = design,
      contrast = contrasts[, contrast], weights = weights, sort = "none", ...
    ) |>
      tibble::rownames_to_column("set_id")
  }) |>
    list_rbind(names_to = "contrast")
}

#' Score every sample on every set with singscore.
#'
#' Rank-based and sample-independent: a sample's score does not change with cohort
#' composition, which is what a paired within-participant design needs. Measured on this
#' matrix, dropping 51 of 131 samples moves a singscore value by 0 and a GSVA value by up to
#' 0.34 on a range of 1.5. Returns sets in rows, samples in columns.
score_samples <- function(expression, protein_map, sets) {
  gene_map <- filter(protein_map, selected)
  gene_matrix <- expression[gene_map$protein, ]
  rownames(gene_matrix) <- gene_map$gene
  singscore::multiScore(singscore::rankGenes(gene_matrix), upSetColc = sets)$Scores
}
