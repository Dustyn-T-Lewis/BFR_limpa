# Score every sample on every set. singscore ranks each sample's proteins and scores sets
# against those ranks, so a score is rank-based and sample-independent: it does not move when
# the cohort changes, which is what a paired within-participant design needs. No p-value, and
# it never sees the contrast. The matrix is what 04_singscore_pheno_associations consumes.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(tibble)
  library(purrr)
})
here::i_am("config.yml")

cfg <- yaml::read_yaml(here("config.yml"))
out <- here("03_Pathway_Enrichment", "02_run_singscore", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  gene_sets = "03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds",
  proteins = "01_Preprocess/02_Quantification/c_data/proteins.rds"
)
paths <- map_chr(inputs, here::here)
if (!all(file.exists(paths))) {
  stop(
    "Run 00_build_gene_sets first. Missing: ",
    paste(inputs[!file.exists(paths)], collapse = ", ")
  )
}
gs <- readRDS(paths[["gene_sets"]])
proteins <- readRDS(paths[["proteins"]])
manifest <- tibble(
  input = names(inputs), path = unname(inputs), md5 = unname(tools::md5sum(paths))
)

gene_map <- filter(gs$protein_map, selected)
stopifnot(
  identical(gene_map$protein, intersect(gene_map$protein, rownames(proteins$E))),
  !anyDuplicated(gene_map$gene)
)

# Rows become gene symbols because that is the namespace the sets are keyed on.
gene_matrix <- proteins$E[gene_map$protein, ]
rownames(gene_matrix) <- gene_map$gene
scores <- singscore::multiScore(
  singscore::rankGenes(gene_matrix),
  upSetColc = gs$sets
)$Scores
stopifnot(
  identical(rownames(scores), names(gs$sets)),
  identical(colnames(scores), colnames(proteins$E))
)
message("scores: ", nrow(scores), " sets x ", ncol(scores), " samples")

# Spread across sets within a sample is set composition; spread across samples within a set is
# what the phenotype analysis has to work with, so both are recorded.
score_summary <- tibble(
  sets = nrow(scores), samples = ncol(scores),
  min = round(min(scores), 4), max = round(max(scores), 4),
  mean = round(mean(scores), 4),
  mean_sd_across_sets = round(mean(apply(scores, 2, sd)), 4),
  mean_sd_across_samples = round(mean(apply(scores, 1, sd)), 4)
)
print(score_summary)

# Participant identity dominates the raw scores, which is why the phenotype analysis works on
# the within-leg change and never on the raw value.
components <- prcomp(t(scores), scale. = FALSE)
variance <- summary(components)$importance[2, 1:4]
targets <- proteins$targets[match(rownames(components$x), proteins$targets$sample_id), ]
participant_share <- map_dbl(1:2, function(i) {
  terms <- summary(aov(components$x[, i] ~ targets$participant))[[1]]
  terms[1, 2] / sum(terms[, 2])
})
structure_check <- tibble(
  component = paste0("PC", 1:4),
  variance_explained = round(as.numeric(variance), 3),
  participant_share = c(round(participant_share, 3), NA, NA)
)
print(structure_check)

packages <- c("here", "limpa", "singscore", "dplyr", "purrr")
versions <- tibble(
  package = packages, version = map_chr(packages, \(p) as.character(packageVersion(p)))
)
score_table <- rownames_to_column(as.data.frame(scores), "set_id")

saveRDS(
  list(
    scores = scores, score_summary = score_summary, structure_check = structure_check,
    provenance = list(
      created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), config = cfg,
      inputs = manifest, packages = versions
    )
  ),
  file.path(out, "singscore.rds"),
  compress = "xz"
)
writexl::write_xlsx(
  list(
    score_summary = score_summary,
    structure_check = structure_check,
    set_scores = score_table,
    input_manifest = manifest,
    package_versions = versions
  ),
  file.path(out, "02_run_singscore.xlsx")
)
readr::write_csv(score_table, file.path(out, "set_scores.csv"))
message("wrote singscore.rds, 02_run_singscore.xlsx and set_scores.csv")
