# 03_Pathway_Enrichment/03_Enrichment/a_script/04_fgsea.R

if (!requireNamespace("fgsea", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("fgsea")
}

# Bioconductor first, tidyverse second  (your Part C rule)
library(fgsea)
library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(readr)

for (fn in c("select", "filter", "rename")) {
  if (!identical(environment(get(fn)), asNamespace("dplyr")))
    stop("dplyr::", fn, "() is masked by ", environmentName(environment(get(fn))))
}

set.seed(1)   # fgseaMultilevel is stochastic

stage <- here("03_Pathway_Enrichment", "01_Gene_Sets", "c_data")
out   <- here("03_Pathway_Enrichment", "03_Enrichment", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

pi_in <- readRDS(file.path(stage, "pathway_inputs.rds"))

ranked        <- pi_in$ranked_t          # named list: one ranked vector per contrast
sets          <- pi_in$gene_sets         # named list of gene symbols (retained sets only)
gene_universe <- pi_in$protein_map$gene[pi_in$protein_map$selected]
CONTROL       <- pi_in$design$negative_control
contrasts     <- names(ranked)

stopifnot(
  length(sets) > 0,
  CONTROL %in% contrasts,
  all(map_lgl(sets, ~ all(.x %in% gene_universe))),
  all(map_lgl(ranked, ~ !is.null(names(.x)) && !anyDuplicated(names(.x))))
)

# The handoff
fit_md5 <- unname(tools::md5sum(here("02_Differential_Expression", "02_Differential",
                                     "c_data", "fit.rds")))
stopifnot(identical(
  fit_md5,
  pi_in$provenance$inputs$md5[pi_in$provenance$inputs$input == "fit"]
))

cat(sprintf("%d sets | %d contrasts | %d genes | control = %s\n",
            length(sets), length(contrasts), length(gene_universe), CONTROL))

tie_report <- map(contrasts, \(cn) {
  s <- ranked[[cn]]
  tibble(contrast = cn, n = length(s),
         n_unique = length(unique(s)),
         pct_tied = round(100 * (1 - length(unique(s)) / length(s)), 2))
}) |> list_rbind()
tie_report

ctl <- fgsea(pathways = sets,
             stats    = ranked[[CONTROL]],
             minSize  = 15,
             maxSize  = 500,
             eps      = 0,
             scoreType = "std",
             nproc    = 1) |>
  as_tibble() |>
  arrange(pval)

ctl_summary <- tibble(
  n_sets       = nrow(ctl),
  sig_p01      = sum(ctl$pval < 0.01),
  sig_p05      = sum(ctl$pval < 0.05),
  expected_p05 = round(0.05 * nrow(ctl)),
  sig_padj05   = sum(ctl$padj < 0.05),
  min_p        = min(ctl$pval)
)
ctl_summary

ctl |> slice_head(n = 10) |>
  select(pathway, size, NES, pval, padj)

gsea_raw <- map(contrasts, \(cn) {
  set.seed(1)
  fgsea(pathways = sets, stats = ranked[[cn]],
        minSize = 15, maxSize = 500, eps = 0,
        scoreType = "std", nproc = 1) |>
    as_tibble() |>
    mutate(contrast = cn, .before = 1)
}) |> list_rbind()

gsea <- gsea_raw |>
  separate_wider_delim(pathway, delim = "|",
                       names = c("database", "set_name"), cols_remove = FALSE) |>
  mutate(padj_within = p.adjust(pval, method = "BH"),
         .by = c(contrast, database)) |>
  mutate(leading_edge_n = lengths(leadingEdge),
         leading_edge_frac = round(leading_edge_n / size, 3))

gsea |>
  summarise(n_sets = n(),
            sig = sum(padj_within < 0.05),
            up = sum(padj_within < 0.05 & NES > 0),
            down = sum(padj_within < 0.05 & NES < 0),
            min_p = min(pval),
            .by = c(contrast, database)) |>
  arrange(match(contrast, c(CONTROL, "interaction",
                            "BFR_post_vs_pre", "HLRT_post_vs_pre", "BFR_vs_HLRT_at_T2")))

concentration <- gsea |>
  filter(padj_within < 0.05) |>
  reframe(gene = unlist(leadingEdge), .by = c(contrast, database, pathway)) |>
  summarise(n_sets = n_distinct(pathway),
            n_genes = n_distinct(gene),
            genes_per_set = round(n_distinct(gene) / n_distinct(pathway), 1),
            .by = c(contrast, database))
concentration

top_drivers <- gsea |>
  filter(padj_within < 0.05) |>
  reframe(gene = unlist(leadingEdge), .by = c(contrast, pathway)) |>
  count(contrast, gene, sort = TRUE, name = "n_sets_driven") |>
  slice_head(n = 15, by = contrast)
top_drivers

sens <- map(c("BFR_post_vs_pre", "HLRT_post_vs_pre"), \(cn) {
  set.seed(1)
  pi_rank <- fgsea(pathways = sets, stats = pi_in$ranked_pi[[cn]],
                   minSize = 15, maxSize = 500, eps = 0, nproc = 1) |>
    as_tibble() |> select(pathway, NES_pi = NES, pval_pi = pval)
  gsea |>
    filter(contrast == cn) |>
    select(pathway, database, NES_t = NES, pval_t = pval, padj_within) |>
    inner_join(pi_rank, by = "pathway") |>
    mutate(contrast = cn, .before = 1)
}) |> list_rbind()

sens |>
  summarise(spearman_NES = round(cor(NES_t, NES_pi, method = "spearman"), 3),
            sig_t   = sum(padj_within < 0.05),
            sig_both = sum(padj_within < 0.05 & pval_pi < 0.05),
            .by = contrast)

gsea_out <- gsea |>
  mutate(leadingEdge = map_chr(leadingEdge, paste, collapse = ";"))

write_csv(gsea_out,      file.path(out, "fgsea_all.csv"))
write_csv(filter(gsea_out, padj_within < 0.05), file.path(out, "fgsea_significant.csv"))
write_csv(mutate(ctl, leadingEdge = map_chr(leadingEdge, paste, collapse = ";")),
          file.path(out, "fgsea_negative_control.csv"))
write_csv(concentration, file.path(out, "fgsea_leading_edge_concentration.csv"))
write_csv(sens,          file.path(out, "fgsea_metric_sensitivity.csv"))

saveRDS(list(
  results     = gsea,
  control     = ctl,
  control_summary = ctl_summary,
  sensitivity = sens,
  provenance  = list(
    created_utc  = format(Sys.time(), tz = "UTC", usetz = TRUE),
    fit_md5      = fit_md5,
    handoff      = pi_in$provenance,
    metric       = "limma moderated t; signed pi as sensitivity",
    background   = "the measured gene universe from 01_Gene_Sets; not the genome",
    correction   = "BH within contrast x database; fgsea's own padj not used",
    eps          = 0,
    seed         = 1,
    fgsea_version = as.character(packageVersion("fgsea")),
    caveat = paste(
      "fgsea permutes gene labels and therefore assumes independent genes.",
      "Muscle proteins are correlated. Read every contrast against the",
      "negative control, never against 0.05.")
  )
), file.path(out, "fgsea_objects.rds"))

# 05_Figures/a_script/06_fgsea_figures.R
# Dot plots + volcano ring plots after Lewis et al. (bioRxiv 2026.04.27.721183) Fig. 3

library(enrichVolcano)
library(ggplot2)
library(here)
library(dplyr)
library(tibble)
library(readr)
library(purrr)
library(stringr)
library(patchwork)

fig_dir <- here("05_Figures", "c_data", "Plots")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Decisions, fixed here -------------------------------------------------
SIG_MODE  <- "fdr"   # "fdr" = BH adj.P.Val (your pipeline). "pi" = Pi-score, as the paper.
SIG_CUT   <- 0.05
RING_N    <- 7       # pathways per direction in the ring
DOT_N     <- 12      # pathways per direction in the dot plot
PATHWAY_Q <- 0.05    # ring/dot plots show sets below this padj_within

CONTRASTS <- c(
  BFR_post_vs_pre   = "BFR: post vs pre",
  HLRT_post_vs_pre  = "High load: post vs pre",
  interaction       = "Interaction: BFR vs high load",
  BFR_vs_HLRT_at_T1 = "Negative control: untrained legs"
)

volc_all <- read_csv(
  here("03_Pathway_Enrichment", "01_Gene_Sets", "c_data", "protein_statistics.csv"),
  show_col_types = FALSE
)

enr_all <- read_csv(
  here("03_Pathway_Enrichment", "03_Enrichment", "c_data", "fgsea_all.csv"),
  show_col_types = FALSE
)

stopifnot(
  all(names(CONTRASTS) %in% volc_all$contrast),
  all(names(CONTRASTS) %in% enr_all$contrast),
  all(c("label", "gene", "logFC", "plot_p", "padj", "pi_score") %in% names(volc_all)),
  all(c("set_name", "database", "NES", "size", "padj_within", "leadingEdge") %in% names(enr_all))
)

volc_all <- volc_all |>
  mutate(
    # every point needs a non-NA id; representatives keep their true symbol so
    # the leading-edge ticks can match them
    gene_plot = if_else(is.na(gene), label, gene),
    sig_value = if (SIG_MODE == "pi") pi_score else padj
  )

# *** THE GOTCHA ***
# volcano_ring() uses ONE padj_col argument for BOTH frames. fgsea's own `padj`
# column is pooled across all three databases, which is not the correction you
# chose. Overwrite it with padj_within so the ring reads the right one.
enr_all <- enr_all |> mutate(padj = padj_within)

ABBREV <- c(
  "Oxidative Phosphorylation"              = "OXPHOS",
  "Aerobic Respiration And Respiratory Electron Transport" = "Aerobic Respiration",
  "Generation Of Precursor Metabolites And Energy" = "Gen. of Precursor\nMetabolites",
  "Mitochondrion Organization"             = "Mito. Org.",
  "Mitochondrial Gene Expression"          = "Mito. Gene Expr.",
  "Mitochondrial Translation"              = "Mito. Translation",
  "Mitochondrial Protein Degradation"      = "Mito. Protein Degr.",
  "Striated Muscle Contraction"            = "Striated Muscle\nContraction",
  "Muscle Cell Differentiation"            = "Muscle Cell Diff.",
  "Myotube Differentiation"                = "Myotube Diff.",
  "Regulation Of Muscle System Process"    = "Reg. Muscle\nSystem Process",
  "Tnfa Signaling Via Nfkb"                = "TNFA/NF-kB",
  "Epithelial Mesenchymal Transition"      = "EMT",
  "Small Molecule Catabolic Process"       = "Small Molecule\nCatabolism",
  "Lipid Catabolic Process"                = "Lipid Catabolism",
  "Heme Metabolism"                        = "Heme Metabolism"
)

clean_label <- function(x) {
  out <- x |>
    str_remove("^(HALLMARK|REACTOME|GOBP|GOCC|GOMF)_") |>
    str_replace_all("_", " ") |>
    str_to_title()
  ifelse(out %in% names(ABBREV), ABBREV[out], out)
}

prep_enrich <- function(cn, n_per_dir) {
  e <- enr_all |> filter(contrast == cn, padj_within < PATHWAY_Q)
  if (nrow(e) == 0) return(e[0, ])
  bind_rows(
    e |> filter(NES > 0) |> slice_min(padj_within, n = n_per_dir, with_ties = FALSE),
    e |> filter(NES < 0) |> slice_min(padj_within, n = n_per_dir, with_ties = FALSE)
  ) |>
    mutate(term_label = make.unique(clean_label(set_name))) |>
    arrange(desc(NES))
}

ring_data <- map(names(CONTRASTS), \(cn) prep_enrich(cn, RING_N)) |> set_names(CONTRASTS)
dot_data  <- map(names(CONTRASTS), \(cn) prep_enrich(cn, DOT_N))  |> set_names(names(CONTRASTS))
volc_data <- map(names(CONTRASTS), \(cn) filter(volc_all, contrast == cn)) |> set_names(CONTRASTS)

map_int(ring_data, nrow)   # check: control should read 3, interaction 10

NES_COLS <- c("#2166AC", "#67A9CF", "#F7F7F7", "#EF8A62", "#B2182B")

dotplot_one <- function(cn) {
  d <- dot_data[[cn]]
  ttl <- unname(CONTRASTS[cn])
  if (nrow(d) == 0) {
    return(ggplot() + theme_void() +
      labs(title = ttl, subtitle = sprintf("no gene set reaches BH q < %.2f", PATHWAY_Q)))
  }
  d <- d |> mutate(term_label = factor(term_label, levels = rev(term_label)))

  ggplot(d, aes(x = NES, y = term_label)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey55") +
    geom_segment(aes(x = 0, xend = NES, yend = term_label),
                 linewidth = 0.35, colour = "grey75") +
    geom_point(aes(size = size, fill = NES), shape = 21, colour = "grey25", stroke = 0.3) +
    scale_fill_gradientn(colours = NES_COLS, limits = c(-3, 3),
                         oob = scales::squish, name = "NES") +
    scale_size_continuous(range = c(2, 8), name = "set size") +
    facet_grid(database ~ ., scales = "free_y", space = "free_y") +
    labs(
      title = ttl,
      subtitle = sprintf(
        "top %d per direction at BH q < %.2f within contrast x database  |  %d sets significant",
        DOT_N, PATHWAY_Q, sum(enr_all$contrast == cn & enr_all$padj_within < PATHWAY_Q)),
      x = "Normalised enrichment score (NES)", y = NULL) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey92"),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey95", colour = NA),
      strip.text.y = element_text(face = "bold", size = 8),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = 8.5, colour = "grey35"),
      axis.text.y = element_text(size = 8))
}

for (cn in names(CONTRASTS)) {
  p <- dotplot_one(cn)
  h <- max(4, 1.6 + 0.22 * max(1, nrow(dot_data[[cn]])))
  ggsave(file.path(fig_dir, paste0("fgsea_dotplot_", cn, ".png")),
         p, width = 8.5, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0("fgsea_dotplot_", cn, ".pdf")),
         p, width = 8.5, height = h, bg = "white")
}

draw_volcano <- local({
  draw <- enrichVolcano::volcano_ring
  penv <- environment(draw)
  validate <- get("validate_enrich_df", envir = penv)
  aenv <- new.env(parent = penv)
  aenv$validate_enrich_df <- function(enrich_df, cols, ...) {
    if (is.data.frame(enrich_df) && nrow(enrich_df) == 0L) return(invisible(enrich_df))
    validate(enrich_df, cols, ...)
  }
  environment(draw) <- aenv
  draw
})

ring_theme <- volcano_ring_theme(base_size = 12, base_family = "sans",
                                 palette = "okabe", ns = "#9a9a95")

ring_subtitle <- function(cn) {
  v <- volc_data[[unname(CONTRASTS[cn])]]
  n_sig  <- sum(v$sig_value < SIG_CUT, na.rm = TRUE)
  n_path <- sum(enr_all$contrast == cn & enr_all$padj_within < PATHWAY_Q)
  n_tot  <- sum(enr_all$contrast == cn)
  sprintf("%d proteins (%s < %.2f)  |  %d / %d pathways (BH q < %.2f)",
          n_sig, if (SIG_MODE == "pi") "Pi" else "FDR", SIG_CUT,
          n_path, n_tot, PATHWAY_Q)
}

ring_one <- function(cn) {
  ttl <- unname(CONTRASTS[cn])
  draw_volcano(
    volc_df   = volc_data[[ttl]],
    enrich_df = ring_data[[ttl]],
    gene_col  = "gene_plot",
    logfc_col = "logFC",
    pval_col  = "plot_p",
    padj_col  = "padj",
    volc_sig_col = "sig_value",
    term_col  = "term_label",
    nes_col   = "NES",
    size_col  = "size",
    genes_col = "leadingEdge",
    genes_sep = ";",
    p_threshold     = SIG_CUT,
    logfc_threshold = 0,
    magnitude  = "neg_log_padj",
    arc_order  = "padj",
    nes_limits = c(-3, 3),
    title = ttl, subtitle = ring_subtitle(cn),
    label_mode = "none",
    show_counts = TRUE,
    point_size = 1.1, label_size = 2.7, axis_size = 2.3, count_size = 3.0,
    theme = ring_theme
  )
}

for (cn in names(CONTRASTS)) {
  p <- ring_one(cn)
  ggsave(file.path(fig_dir, paste0("volcano_ring_", cn, ".png")),
         p, width = 7.2, height = 6.8, dpi = 300, bg = "white")
  ggsave(file.path(fig_dir, paste0("volcano_ring_", cn, ".pdf")),
         p, width = 7.2, height = 6.8, bg = "white")
}

grid <- volcano_ring_grid(
  volc_dfs   = volc_data,
  enrich_dfs = ring_data,
  contrasts  = unname(CONTRASTS),
  subtitles  = map_chr(names(CONTRASTS), ring_subtitle),
  nrow = 2, ncol = 2,
  tag_levels = "A",
  guides = "collect",
  legend_position = "bottom", legend_width = 34,
  label_headroom = 1.2,
  gene_col = "gene_plot", logfc_col = "logFC", pval_col = "plot_p",
  padj_col = "padj", volc_sig_col = "sig_value",
  term_col = "term_label", nes_col = "NES", size_col = "size",
  genes_col = "leadingEdge", genes_sep = ";",
  p_threshold = SIG_CUT, logfc_threshold = 0,
  magnitude = "neg_log_padj", arc_order = "padj", nes_limits = c(-3, 3),
  label_mode = "none", show_counts = TRUE,
  point_size = 1.0, label_size = 2.5, axis_size = 2.1, count_size = 2.8,
  theme = ring_theme
)

ggsave(file.path(fig_dir, "volcano_ring_grid_4contrasts.pdf"),
       grid$plot, width = 13.5, height = 13.0, bg = "white")
ggsave(file.path(fig_dir, "volcano_ring_grid_4contrasts.png"),
       grid$plot, width = 13.5, height = 13.0, dpi = 300, bg = "white")

cat("wrote:\n"); cat(paste0("  ", list.files(fig_dir, pattern = "^(fgsea_dotplot|volcano_ring)"), "\n"))
