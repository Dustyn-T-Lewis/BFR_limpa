# Score every sample on every set with singscore. Scores rest on within-sample ranks, so they do
# not move when the cohort changes, as a paired within-participant design needs. No p-value, no
# contrast. 05_classify_and_associate_sets reads set_scores.rds.

pacman::p_load(here, dplyr, tibble, purrr, stringr, ggplot2, patchwork, readxl, writexl)

stage <- here("03_Pathway_Enrichment", "04_run_singscore")
gene_sets <- here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data")
sets <- readRDS(file.path(gene_sets, "gene_sets.rds"))
set_catalog <- read_excel(file.path(gene_sets, "00_build_gene_sets.xlsx"), "set_catalog")
protein_map <- read_excel(file.path(gene_sets, "00_build_gene_sets.xlsx"), "protein_gene_map")
proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins.rds"))

# Sets are keyed on gene symbols.
gene_map <- filter(protein_map, selected)
gene_matrix <- proteins$E[gene_map$protein, ]
rownames(gene_matrix) <- gene_map$gene
scored <- singscore::multiScore(singscore::rankGenes(gene_matrix), upSetColc = sets)
scores <- scored$Scores
message("scores: ", nrow(scores), " sets x ", ncol(scores), " samples")

# Dispersion is the spread of a set's member ranks within one sample: low means the members sit
# together, high that they scatter. multiScore returns it beside the scores at no extra cost.
set_spread <- set_catalog |>
  filter(qualifies) |>
  transmute(set_id, database) |>
  mutate(
    score = rowMeans(scores[set_id, ]),
    dispersion = rowMeans(scored$Dispersions[set_id, ])
  )
collection_spread <- set_spread |>
  summarise(
    sets = n(), median_score = round(median(score), 4),
    median_dispersion = round(median(dispersion), 1), .by = database
  )
print(as.data.frame(collection_spread))

# Participant identity dominates the raw scores, so the phenotype analysis uses the within-leg
# change, never the raw value. targets rows follow the matrix columns.
components <- prcomp(t(scores), scale. = FALSE)
structure_check <- tibble(
  component = paste0("PC", 1:4),
  variance_explained = round(summary(components)$importance[2, 1:4], 3),
  participant_share = c(round(map_dbl(1:2, \(i) {
    summary(lm(components$x[, i] ~ proteins$targets$participant))$r.squared
  }), 3), NA, NA)
)
print(structure_check)

# singscore's plotDispersion and plotRankDensity draw one signature at a time. These cohort-scale
# views are plain ggplots over the score and dispersion matrices multiScore already returned.
group_scores <- tibble(
  group = rep(proteins$targets$group, each = nrow(scores)),
  score = as.vector(scores)
)
figure <- wrap_plots(
  ggplot(set_spread, aes(score, dispersion, colour = database)) +
    geom_vline(xintercept = 0, linewidth = 0.3, colour = "grey80") +
    geom_point(alpha = 0.4, size = 0.8) +
    scale_colour_brewer(palette = "Dark2", name = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "mean score across samples", y = "mean dispersion across samples") +
    theme_minimal(base_size = 10),
  ggplot(group_scores, aes(score, group)) +
    geom_violin(fill = "grey85", colour = NA) +
    geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.3) +
    labs(x = "singscore", y = NULL) +
    theme_minimal(base_size = 10),
  ncol = 1, heights = c(1.4, 1)
) +
  plot_annotation(
    caption = str_wrap(paste(
      "S12 Figure. singscore across the cohort.",
      sprintf("(A) One point per set, %d sets, score and dispersion averaged", nrow(scores)),
      sprintf("over %d samples.", ncol(scores)),
      "Dispersion is the spread of a set's member ranks within a sample: low means members sit",
      "together. (B) Every set in every sample, pooled within study group; the box is the",
      "interquartile range. Data: set_scores sheet of 04_run_singscore.xlsx."
    ), 115),
    tag_levels = "A",
    theme = theme(
      plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
      plot.caption.position = "plot"
    )
  )
pdf(file.path(stage, "b_reports", "04_run_singscore_figures.pdf"), width = 8.27, height = 11.69)
print(figure)
invisible(dev.off())

sheets <- list(
  collection_spread = collection_spread,
  structure_check = structure_check,
  set_scores = rownames_to_column(as.data.frame(scores), "set_id")
)
overview <- tibble(
  sheet = names(sheets),
  rows = map_int(sheets, nrow),
  columns = map_int(sheets, ncol),
  description = c(
    "Per collection: sets, median score and median dispersion",
    "Variance each leading component explains, and the share of it participant explains",
    "The set by sample singscore matrix; 05_classify_and_associate_sets reads set_scores.rds"
  )
)
# 05 reads the matrix from rds: singscore values carry exact rank ties, and a round trip through
# the workbook changes the last bit of some. A broken tie switches the Wilcoxon and Spearman tests
# from the approximate p to the exact one.
saveRDS(scores, file.path(stage, "c_data", "set_scores.rds"))
write_xlsx(
  c(list(overview = overview), sheets),
  file.path(stage, "c_data", "04_run_singscore.xlsx")
)
message("wrote set_scores.rds, 04_run_singscore.xlsx and a 1-page figure PDF")
sessionInfo()
