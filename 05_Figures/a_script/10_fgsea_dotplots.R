# Dot plots from fgsea results, D's catalog (1,040 sets).
# same abbreviations. Only the x-axis quantity differs (NES vs mean member t).

library(ggplot2); library(here); library(dplyr); library(tibble)
library(tidyr); library(readr); library(purrr); library(stringr); library(patchwork)

fig <- here("05_Figures", "c_data", "Plots")
tab <- here("05_Figures", "c_data")
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

FDR     <- 0.05
TOP_N   <- 12     # sets per direction
NULL_N  <- 10     # shown when nothing is significant
RANK_BY <- "q"    # "q" (evidence) or "NES" (magnitude)

CONTRASTS <- c(
  BFR_post_vs_pre   = "BFR: post vs pre",
  HLRT_post_vs_pre  = "High load: post vs pre",
  interaction       = "Interaction: BFR vs high load",
  BFR_vs_HLRT_at_T2 = "BFR vs high load at T2",
  BFR_vs_HLRT_at_T1 = "Negative control: untrained legs"
)

fg <- read_csv(here("03_Pathway_Enrichment","02_Set_Tests","c_data",
                    "three_way_all.csv"), show_col_types = FALSE) |>
  filter(method == "fgsea") |>
  select(contrast, set, database, n, direction, p, NES, q)

stopifnot(all(names(CONTRASTS) %in% fg$contrast),
          !anyNA(fg$NES), nrow(fg) == 1040 * n_distinct(fg$contrast))

cat(sprintf("fgsea: %d sets x %d contrasts | NES range %.2f to %.2f\n",
            1040, n_distinct(fg$contrast), min(fg$NES), max(fg$NES)))

ABBREV <- c(
  "Oxidative Phosphorylation" = "OXPHOS",
  "Aerobic Respiration And Respiratory Electron Transport" = "Aerobic respiration",
  "Generation Of Precursor Metabolites And Energy" = "Gen. precursor metabolites",
  "Mitochondrial Protein Degradation" = "Mito. protein degradation",
  "Citric Acid Cycle Tca Cycle" = "TCA cycle (Reactome)",
  "Tricarboxylic Acid Cycle" = "TCA cycle (GO)",
  "Mitochondrial Rna Metabolic Process" = "Mito. RNA metabolism",
  "Mitochondrial Gene Expression" = "Mito. gene expression",
  "Mitochondrial Translation" = "Mito. translation",
  "Mitochondrial Respiratory Chain Complex Assembly" = "Resp. chain assembly",
  "Amide Biosynthetic Process" = "Peptide/amide biosynthesis (translation)",
  "Striated Muscle Contraction" = "Striated muscle contraction",
  "Muscle Cell Differentiation" = "Muscle cell differentiation",
  "Myotube Differentiation" = "Myotube differentiation",
  "Small Molecule Catabolic Process" = "Small molecule catabolism",
  "Regulation Of Muscle System Process" = "Reg. muscle system process",
  "Positive Regulation Of Autophagy" = "Pos. reg. autophagy",
  "Calcium Ion Transport" = "Ca2+ transport"
)

clean_label <- function(x) {
  out <- x |> str_remove("^(HALLMARK|REACTOME|GOBP)_") |>
    str_replace_all("_", " ") |> str_to_sentence()
  out <- str_replace(out, "^(.)", toupper)
  ifelse(out %in% names(ABBREV), ABBREV[out], out)
}

MITO <- regex(paste0("MITOCHOND|OXIDATIVE_PHOS|TRICARBOX|CITRIC_ACID|RESPIRAT|",
                     "AEROBIC|CYTOCHROME|NADH|ATP_SYNTH|ELECTRON_TRANSPORT|",
                     "TCA_CYCLE|PEROXISOM|FATTY_ACID|BETA_OXIDATION"))

fg <- fg |>
  mutate(set_name = sub("^.*\\|", "", set),
         theme = if_else(str_detect(set_name, MITO),
                         "Mitochondrial / oxidative", "Other"))

pick <- function(cn) {
  d   <- filter(fg, contrast == cn)
  sig <- filter(d, q < FDR)

  sel <- if (nrow(sig) == 0) {
    slice_min(d, p, n = NULL_N, with_ties = FALSE) |> mutate(is_sig = FALSE)
  } else {
    top <- function(x, up) {
      if (!nrow(x)) return(x)
      if (RANK_BY == "q") slice_min(x, q, n = TOP_N, with_ties = FALSE)
      else if (up)        slice_max(x, NES, n = TOP_N, with_ties = FALSE)
      else                slice_min(x, NES, n = TOP_N, with_ties = FALSE)
    }
    bind_rows(top(filter(sig, NES > 0), TRUE),
              top(filter(sig, NES < 0), FALSE)) |> mutate(is_sig = TRUE)
  }

  sel |>
    mutate(label = make.unique(clean_label(set_name)),
           nlq   = -log10(pmax(q, .Machine$double.xmin))) |>
    arrange(database, NES)
}

plot_data <- map(names(CONTRASTS), pick) |> set_names(names(CONTRASTS))
map_int(plot_data, nrow)      # expect 24, 24, 10, 24, 3


pick <- function(cn) {
  d   <- filter(fg, contrast == cn)
  sig <- filter(d, q < FDR)

  sel <- if (nrow(sig) == 0) {
    slice_min(d, p, n = NULL_N, with_ties = FALSE) |> mutate(is_sig = FALSE)
  } else {
    top <- function(x, up) {
      if (!nrow(x)) return(x)
      if (RANK_BY == "q") slice_min(x, q, n = TOP_N, with_ties = FALSE)
      else if (up)        slice_max(x, NES, n = TOP_N, with_ties = FALSE)
      else                slice_min(x, NES, n = TOP_N, with_ties = FALSE)
    }
    bind_rows(top(filter(sig, NES > 0), TRUE),
              top(filter(sig, NES < 0), FALSE)) |> mutate(is_sig = TRUE)
  }

  sel |>
    mutate(label = make.unique(clean_label(set_name)),
           nlq   = -log10(pmax(q, .Machine$double.xmin))) |>
    arrange(database, NES)
}

plot_data <- map(names(CONTRASTS), pick) |> set_names(names(CONTRASTS))
map_int(plot_data, nrow)      # expect 24, 24, 10, 24, 3

XLIM <- max(abs(unlist(map(plot_data, "NES")))) * 1.15

dot_one <- function(cn) {
  d     <- plot_data[[cn]]
  ttl   <- unname(CONTRASTS[cn])
  n_sig <- sum(fg$contrast == cn & fg$q < FDR)

  sub <- if (n_sig == 0) {
    sprintf("0 / 1,040 sets at BH q < %.2f \u2014 %d strongest by raw p shown", FDR, nrow(d))
  } else {
    sprintf("%d / 1,040 sets at BH q < %.2f \u2014 top %d per direction by %s",
            n_sig, FDR, TOP_N, if (RANK_BY == "q") "evidence" else "NES")
  }

  d <- mutate(d, label = factor(label, levels = unique(d$label)))

  ggplot(d, aes(NES, label)) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey45") +
    geom_segment(aes(x = 0, xend = NES, yend = label),
                 linewidth = 0.3, colour = "grey78") +
    geom_point(aes(size = n, fill = nlq, shape = theme),
               colour = "grey25", stroke = 0.3) +
    scale_shape_manual(values = c("Mitochondrial / oxidative" = 21, "Other" = 24),
                       name = NULL, drop = FALSE) +
    scale_fill_viridis_c(option = "mako", direction = -1,
                         name = expression(-log[10]~(q))) +
    scale_size_continuous(range = c(2, 8), name = "measured\nmembers") +
    scale_x_continuous(limits = c(-XLIM, XLIM)) +
    facet_grid(database ~ ., scales = "free_y", space = "free_y") +
    labs(title = ttl, subtitle = sub,
         x = "Normalised enrichment score (NES)", y = NULL) +
    guides(shape = guide_legend(override.aes = list(size = 3.5, fill = "grey70"))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.major.y = element_line(linewidth = 0.25, colour = "grey92"),
          panel.grid.minor    = element_blank(),
          strip.background    = element_rect(fill = "grey95", colour = NA),
          strip.text.y        = element_text(face = "bold", size = 8),
          plot.title          = element_text(face = "bold"),
          plot.subtitle       = element_text(size = 8.5, colour = "grey35"),
          axis.text.y         = element_text(size = 8),
          legend.position     = "right", legend.box = "vertical")
}

for (cn in names(CONTRASTS)) {
  p <- dot_one(cn)
  h <- max(4, 1.9 + 0.23 * nrow(plot_data[[cn]]))
  ggsave(file.path(fig, paste0("fgsea_dotplot_", cn, ".png")), p,
         width = 9.5, height = h, dpi = 300, bg = "white")
  ggsave(file.path(fig, paste0("fgsea_dotplot_", cn, ".pdf")), p,
         width = 9.5, height = h, bg = "white")
}

combined <- wrap_plots(map(names(CONTRASTS), dot_one), ncol = 2, guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom", legend.box = "horizontal")

ggsave(file.path(fig, "fgsea_dotplot_grid_all_contrasts.pdf"), combined,
       width = 18, height = 18, bg = "white")
ggsave(file.path(fig, "fgsea_dotplot_grid_all_contrasts.png"), combined,
       width = 18, height = 18, dpi = 300, bg = "white")

list_rbind(map(names(CONTRASTS), \(cn) mutate(plot_data[[cn]], contrast = cn))) |>
  select(contrast, database, set_name, theme, n, direction, NES, p, q, is_sig, label) |>
  arrange(contrast, database, desc(NES)) |>
  write_csv(file.path(tab, "fgsea_dotplot_values.csv"))

cat("wrote:\n"); cat(paste0("  ", list.files(fig, "^fgsea_dotplot"), "\n"))



