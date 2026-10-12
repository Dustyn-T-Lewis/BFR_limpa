# How much mitochondria each leg gained (MPE), and whether any part of the mitochondrion changed
# beyond that (mitochondrial normalisation). Both follow Granata et al. 2021 (Nat Commun 12:7056)
# and Reisman et al. 2024 (Nat Commun 15:7677), applied here to whole muscle, not isolates.

pacman::p_load(here, dplyr, tibble, tidyr, purrr, stringr, limma, limpa, ggplot2, patchwork,
               readxl, writexl)

stage <- here("05_Mitochondria", "02_test_mito")
quant <- here("01_Preprocess", "02_Quantification", "c_data")
proteins <- readRDS(file.path(quant, "proteins.rds"))
raw <- readRDS(file.path(quant, "quant_runs", "proteins_slope_0.7.rds"))[rownames(proteins), ]
d <- readRDS(here("02_Differential_Expression", "01_Design", "c_data", "design.rds"))
fit <- readRDS(here("02_Differential_Expression", "02_Differential", "c_data", "fit.rds"))
protein_map <- read_excel(
  here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data", "00_build_gene_sets.xlsx"),
  "protein_gene_map"
)
membership <- read_excel(
  here("05_Mitochondria", "01_build_mito_sets", "c_data", "01_build_mito_sets.xlsx"),
  "membership"
)

# proteins.rds is this raw run after cyclic loess and nothing else.
stopifnot(
  identical(protein_map$protein, rownames(proteins)),
  identical(colnames(raw), colnames(proteins)),
  identical(normalizeBetweenArrays(raw$E, method = "cyclicloess"), proteins$E)
)
contrast_names <- colnames(d$contrasts)
targets <- proteins$targets
# MitoCarta keeps small pathways because some functions are one or two proteins, so there is no
# size floor. The 20% coverage floor is ours, not a published rule: below it the measured part
# says little about the pathway.
coverage <- membership |>
  summarise(listed = n(), measured = sum(measured), .by = c(collection, set)) |>
  mutate(coverage = measured / listed)
sets <- with(filter(membership, measured), split(gene, set))
complexes <- paste(c("CI", "CII", "CIII", "CIV", "CV"), "subunits")
pathways <- sets[filter(coverage, collection == "MitoPathways", coverage >= 0.2)$set]

# MPE is a share of the raw signal: loess assumes the typical protein holds still, which is what a
# mitochondrial expansion breaks. Every protein group counts toward the total.
intensity <- 2^raw$E
is_mito <- protein_map$gene %in% sets$MitoCarta
mpe <- tibble(
  sample_id = colnames(raw), participant = targets$participant, leg = targets$leg,
  treatment = targets$treatment, timepoint = targets$timepoint,
  mpe = colSums(intensity[is_mito, ]) / colSums(intensity)
)

# One outcome row, so there is nothing for eBayes to borrow from: the plain least-squares t on
# the same design and contrasts as the protein fit.
mpe_fit <- contrasts.fit(lmFit(rbind(log2_mpe = log2(mpe$mpe)), d$design), d$contrasts)
mpe_test <- tibble(
  contrast = contrast_names, log2FC = mpe_fit$coefficients[1, ],
  se = mpe_fit$stdev.unscaled[1, ] * mpe_fit$sigma, df = mpe_fit$df.residual
) |>
  mutate(
    t = log2FC / se, p = 2 * pt(-abs(t), df),
    conf.low = log2FC - qt(0.975, df) * se, conf.high = log2FC + qt(0.975, df) * se
  )

# Citrate synthase activity is the literature's content marker; this is its protein, from the
# stage 02 fit, for comparison only.
cs <- protein_map$protein[protein_map$gene %in% "CS" & protein_map$selected]
citrate_synthase <- map(set_names(contrast_names), \(cn) {
  topTable(fit, coef = cn, number = Inf, sort.by = "none")[cs, c("logFC", "P.Value", "adj.P.Val")]
}) |>
  list_rbind(names_to = "contrast")

# Mitochondrial normalisation as both papers did it: subset the mitochondrial proteins from the
# unnormalised data, normalise the subset alone, refit.
mito <- raw[is_mito, ]
mito$E <- normalizeBetweenArrays(mito$E, method = "cyclicloess")
mito_fit <- dpcDE(mito, d$design, plot = FALSE, sample.weights = TRUE) |>
  contrasts.fit(d$contrasts) |>
  eBayes(robust = TRUE)
mito_proteins <- map(set_names(contrast_names), \(cn) {
  topTable(mito_fit, coef = cn, number = Inf, sort.by = "none") |>
    rownames_to_column("protein") |>
    select(protein, logFC, AveExpr, t, P.Value, adj.P.Val)
}) |>
  list_rbind(names_to = "contrast") |>
  left_join(select(protein_map, protein, gene, label, selected), by = "protein")

# fry reads limpa's precision weights and indexes each symbol's representative protein, as in
# 01_run_fgsea_and_fry. It returns its own BH FDR per contrast.
representative <- with(protein_map, set_names(if_else(selected, gene, NA), protein))
fry_sets <- function(y, weights, sets) {
  y$weights <- weights
  index <- ids2indices(sets, representative[rownames(y)])
  map(set_names(contrast_names), \(cn) {
    fry(y, index, d$design, d$contrasts[, cn], sort = "none") |> rownames_to_column("set")
  }) |>
    list_rbind(names_to = "contrast")
}
complex_tests <- bind_rows(
  standard = fry_sets(proteins, fit$EList$weights, sets[complexes]),
  mitochondria_normalised = fry_sets(mito, mito_fit$EList$weights, sets[complexes]),
  .id = "fit"
) |>
  select(fit, contrast, set, NGenes, Direction, PValue, FDR)

# fry is the confirmatory test, as in 01_run_fgsea_and_fry. After normalisation the shared
# expansion is gone, so it already reads change relative to the rest of the mitochondrion and a
# competitive test (camera) would ask the same question again. fgsea draws the rings. It shuffles
# proteins as if independent, so on correlated sets such as complex I its p runs small: between
# legs it calls pathways fry does not, and those are leads.
ranked <- map(set_names(contrast_names), \(cn) {
  with(filter(mito_proteins, contrast == cn, selected), set_names(t, gene))
})
set.seed(1)
fgsea_raw <- map(ranked, \(stats) fgsea::fgsea(pathways, stats)) |>
  list_rbind(names_to = "contrast")
pathway_tests <- bind_rows(
  fry = fry_sets(mito, mito_fit$EList$weights, pathways),
  fgsea = transmute(
    fgsea_raw, contrast, set = pathway, NGenes = size,
    Direction = if_else(NES > 0, "Up", "Down"), NES, PValue = pval, FDR = padj,
    leading_edge = map_chr(leadingEdge, paste, collapse = ";")
  ),
  .id = "method"
) |>
  select(method, contrast, set, NGenes, Direction, NES, PValue, FDR, leading_edge) |>
  left_join(select(coverage, set, listed, coverage), by = "set")

theme_set(theme_minimal(base_size = 9))
theme_update(
  plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
  plot.caption.position = "plot"
)
caption <- \(title, text) str_wrap(paste0(title, ". ", text), 115)
arm_colours <- c(BFR = "#B2182B", HLRT = "#2166AC")

mpe_page <- ggplot(mpe, aes(timepoint, 100 * mpe, colour = treatment)) +
  geom_line(aes(group = paste(participant, leg)), alpha = 0.35) +
  geom_boxplot(aes(group = timepoint), width = 0.3, outlier.shape = NA, fill = NA) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~treatment) +
  scale_colour_manual(values = arm_colours, guide = "none") +
  labs(
    x = NULL, y = "MPE, % of raw protein intensity",
    caption = caption("Mitochondrial protein enrichment by arm and time", paste(
      "MPE is the summed raw intensity of the", sum(is_mito), "measured MitoCarta 3.0 proteins",
      "over all", format(nrow(raw), big.mark = ","), "proteins, per biopsy. Lines join each leg",
      "before and after training. Tests: mpe_test sheet of 02_test_mito.xlsx."
    ))
  )

complex_page <- complex_tests |>
  mutate(
    score = if_else(Direction == "Up", 1, -1) * -log10(FDR),
    contrast = factor(contrast, contrast_names), set = factor(set, rev(complexes)),
    fit = factor(
      fit, c("standard", "mitochondria_normalised"),
      c("standard fit", "mitochondria-normalised fit")
    )
  ) |>
  ggplot(aes(contrast, set, fill = score)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = signif(FDR, 2)), size = 2.6) +
  facet_wrap(~fit, ncol = 1) +
  scale_fill_gradient2(low = "#2166AC", high = "#B2182B", name = "signed -log10 FDR") +
  labs(
    x = NULL, y = NULL,
    caption = caption("Complexes I to V, standard and mitochondria-normalised fits", paste(
      "fry on the structural subunits of each complex; labels give fry FDR within the contrast,",
      "colour its direction. The standard fit asks whether a complex changed; the normalised fit",
      "whether it changed beyond overall mitochondrial content. Data: complex_tests sheet."
    ))
  ) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "top")

da <- enrichVolcano::as_da(
  select(mito_proteins, protein, gene = label, contrast, logFC, t, P.Value, adj.P.Val),
  species = NULL
)
# The ring draws the twelve lowest-padj fgsea pathways, one per overlapping cluster.
enrichment <- fgsea_raw |>
  mutate(leadingEdge = map(leadingEdge, \(genes) {
    unname(set_names(protein_map$label, protein_map$gene)[genes])
  })) |>
  enrichVolcano::as_enrichment(enrichment_test = "fgsea", database = "MitoPathways") |>
  enrichVolcano::dedup_terms(pathways)
volcano <- function(cn) {
  enrichVolcano::plot_volcano_ring(
    da, enrichment,
    contrast = cn, title = cn,
    label_mode = "top_per_direction", label_size = 2.6, axis_size = 2.6, point_size = 0.9
  ) +
    # As in 02_enrich_volcano_fgsea: a horizontal NES key below the ring collides with nothing.
    guides(fill = guide_colorbar(direction = "horizontal", title.position = "top")) +
    theme(
      legend.position = "bottom", legend.key.height = unit(2, "mm"),
      legend.key.width = unit(16, "mm")
    )
}
drawn <- setdiff(contrast_names, d$negative_control)
volcano_page <- wrap_plots(map(drawn, volcano), ncol = 2) +
  plot_annotation(caption = caption("Mitochondria-normalised proteins", paste(
    "Points are the", nrow(mito), "MitoCarta proteins after normalising them alone, coloured at",
    "BH FDR < 0.05, five named per direction. The ring shows up to twelve MitoPathways at fgsea",
    "FDR < 0.05, one per overlapping cluster, coloured by NES; ticks mark leading-edge proteins.",
    "fgsea treats proteins as independent, so ringed pathways that fry does not call are leads.",
    "Data: mito_proteins and pathway_tests sheets."
  )))

pdf(file.path(stage, "b_reports", "02_test_mito.pdf"), width = 8.27, height = 11.69)
walk(list(mpe_page, complex_page, volcano_page), print)
invisible(dev.off())

sheets <- list(
  mpe = mpe, mpe_test = mpe_test, citrate_synthase = citrate_synthase,
  complex_tests = complex_tests,
  pathway_tests = pathway_tests, mito_proteins = mito_proteins
)
overview <- data.frame(
  sheet = names(sheets), rows = sapply(sheets, nrow),
  description = c(
    "MPE per biopsy, from raw intensities.",
    "log2 MPE through the saved design and contrasts; plain least-squares t, 95% CI.",
    "Citrate synthase protein from the stage 02 fit, for comparison with CS activity studies.",
    "fry on complexes I to V, standard and mitochondria-normalised fits; FDR within contrast.",
    "fry and fgsea on MitoPathways at 20% coverage, mitochondria-normalised; FDR within contrast.",
    "Mitochondria-normalised dpcDE, every contrast; adj.P.Val is BH within contrast."
  )
)
write_xlsx(c(list(overview = overview), sheets), file.path(stage, "c_data", "02_test_mito.xlsx"))
sessionInfo()
