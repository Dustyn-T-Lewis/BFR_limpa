# ETC Complex (I-V) logFC Comparison: BFR_post_vs_pre vs HLRT_post_vs_pre
#
# PURPOSE:
#   Plot 1: Classify proteins into mitochondrial ETC complexes I-V
#   (structural subunits) using MitoCarta3.0, then compare per-complex
#   logFC between the BFR_post_vs_pre and HLRT_post_vs_pre contrasts.
#   Plot 2: Same comparison for the components NOT included in Plot 1 --
#   per-complex assembly factors (CI-CV assembly factors), plus an
#   "Other" category for OXPHOS-related proteins that aren't tied to a
#   single complex (Cytochrome C carriers, general/respirasome assembly
#   factors).
#
#   Both plots: bar = mean logFC (dodged by group), SE error bars,
#   jittered points for every classified protein, and a paired
#   (protein-matched) t-test bracket + p-value label above each
#   complex/category's bars (bold when p < 0.05).
#
# DATA LAYOUT (02_differential.xlsx, sheet 1):
#   One row per protein, wide format:
#     protein | Genes | NPrec | <contrast>_logFC | <contrast>_P.Value | <contrast>_adj.P.Val ...
#   for each of: BFR_post_vs_pre, HLRT_post_vs_pre, BFR_vs_HLRT_at_T1,
#   BFR_vs_HLRT_at_T2, interaction. This script only pulls the two
#   post-vs-pre contrasts.
#
# ASSUMPTION (please confirm):
#   `protein` is treated as a single UniProt accession that matches
#   MitoCarta's `UniProt` column directly. If `protein` is actually a
#   DIA-NN-style Protein.Group with multiple semicolon-separated IDs
#   (e.g. "P12345;Q67890"), uncomment the `separate_rows()` line marked
#   below to explode it into one row per ID before joining -- otherwise
#   those rows will silently fail to match MitoCarta and get dropped.

# --- 0: Setup ------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, dplyr, tidyr, stringr, ggplot2)

# --- Configuration -----------------------------------------------------------
config <- list(
  diff_data = "02_Differential_Expression/02_Differential/c_data/02_differential.xlsx",
  mitocarta  = "00_input/Human.MitoCarta3.0.xls",
  out_dir    = "02_Differential_Expression/05_MitoCarta/b_reports"
)
dir.create(config$out_dir, showWarnings = FALSE, recursive = TRUE)

contrasts <- c("BFR_post_vs_pre", "HLRT_post_vs_pre")

# Bar/point colors (ColorBrewer Dark2-style, consistent with existing
# figure palette; distinct from the GO:BP/CC/MF hues already in use)
group_colors <- c(BFR_post_vs_pre = "darkred", HLRT_post_vs_pre = "darkblue")
group_labels <- c(BFR_post_vs_pre = "BFR", HLRT_post_vs_pre = "HLRT")

complex_order <- c("I", "II", "III", "IV", "V")
dodge_width   <- 0.7

cat("=== ETC Complex logFC Comparison: BFR vs HLRT ===\n\n")

mito <- read_excel(config$mitocarta, sheet = "A Human MitoCarta3.0")

# --- Read + reshape the single-sheet differential results ------------------
# Wide -> long: one row per protein x contrast, with logFC / P.Value /
# adj.P.Val as columns. This replaces the old per-contrast-sheet read loop.

diff_wide <- read_excel(config$diff_data, sheet = 1)

diff_long <- diff_wide %>%
  select(protein, gene = Genes, matches(paste0("^(", paste(contrasts, collapse = "|"), ")_"))) %>%
  pivot_longer(
    cols = -c(protein, gene),
    names_to = c("contrast", ".value"),
    names_pattern = paste0("^(", paste(contrasts, collapse = "|"), ")_(logFC|P\\.Value|adj\\.P\\.Val)$")
  ) %>%
  rename(P.Value = `P.Value`, adj.P.Val = `adj.P.Val`) %>%
  mutate(contrast = factor(contrast, levels = contrasts))

# If `protein` is a semicolon-delimited Protein.Group rather than a single
# UniProt accession, explode it before joining to MitoCarta:
# diff_long <- diff_long %>% separate_rows(protein, sep = ";")

cat(sprintf("Reshaped to %d protein x contrast rows (%d unique proteins)\n\n",
            nrow(diff_long), n_distinct(diff_long$protein)))

# --- Helper functions ------------------------------------------------------

# Classify proteins by matching MitoCarta's MitoPathways string against a
# capturing regex (capture group 1 = category label). Drops any protein
# matched to more than one category (ambiguous assignment), reporting
# what was dropped.
classify_by_regex <- function(mito, regex, valid_categories, label) {
  out <- mito %>%
    select(UniProt, Symbol, MitoCarta3.0_MitoPathways) %>%
    filter(!is.na(MitoCarta3.0_MitoPathways)) %>%
    mutate(category = str_match(MitoCarta3.0_MitoPathways, regex)[, 2]) %>%
    filter(!is.na(category), category %in% valid_categories) %>%
    distinct(UniProt, Symbol, category)
  
  dual <- out %>% count(UniProt) %>% filter(n > 1) %>% pull(UniProt)
  if (length(dual) > 0) {
    cat(sprintf(
      "  [%s] Dropping %d protein(s) with ambiguous multi-category assignment: %s\n",
      label, length(dual),
      paste(out$Symbol[out$UniProt %in% dual], collapse = ", ")
    ))
    out <- out %>% filter(!UniProt %in% dual)
  }
  out
}

# Merge a category-classification table with the reshaped differential
# results for both contrasts, then print detection counts.
build_category_data <- function(classification, category_order, label) {
  data <- diff_long %>%
    inner_join(classification, by = c("protein" = "UniProt")) %>%
    mutate(
      category = factor(category, levels = category_order),
      group    = factor(contrast, levels = contrasts)
    )
  
  cat(sprintf("  [%s] %d classified protein x contrast rows matched\n", label, nrow(data)))
  print(table(data$category, data$group))
  cat("\n")
  data
}

# Detected-in-dataset vs total-in-MitoCarta counts per category, for
# "(n = detected/total)" x-axis labels.
build_n_labels <- function(data, classification, category_order, prefix) {
  n_detected <- data %>% distinct(category, protein) %>% count(category, name = "n_detected")
  n_total    <- classification %>% count(category, name = "n_total")
  
  n_total %>%
    left_join(n_detected, by = "category") %>%
    mutate(
      n_detected = tidyr::replace_na(n_detected, 0),
      category   = factor(category, levels = category_order)
    ) %>%
    arrange(category) %>%
    mutate(label = sprintf("%s %s\n(n = %d/%d)", prefix, category, n_detected, n_total))
}

# Mean logFC + SE per category x group.
summarize_data <- function(data) {
  data %>%
    group_by(category, group) %>%
    summarise(mean_logFC = mean(logFC), se_logFC = sd(logFC) / sqrt(n()), n = n(), .groups = "drop")
}

# Paired (protein-matched) t-test per category comparing BFR_post_vs_pre
# vs HLRT_post_vs_pre logFC, plus bracket geometry for plotting.
test_categories <- function(data, summary, category_order) {
  stats <- lapply(category_order, function(cat_lvl) {
    subset <- data %>% filter(category == cat_lvl)
    
    if (nrow(subset) == 0) {
      return(data.frame(category = cat_lvl, p_value = NA_real_, n_paired = 0))
    }
    
    wide <- subset %>%
      select(protein, group, logFC) %>%
      pivot_wider(names_from = group, values_from = logFC)
    
    if (!all(contrasts %in% colnames(wide))) {
      return(data.frame(category = cat_lvl, p_value = NA_real_, n_paired = 0))
    }
    
    wide <- wide %>% filter(!is.na(.data[[contrasts[1]]]), !is.na(.data[[contrasts[2]]]))
    
    if (nrow(wide) < 2) {
      return(data.frame(category = cat_lvl, p_value = NA_real_, n_paired = nrow(wide)))
    }
    tt <- t.test(wide[[contrasts[1]]], wide[[contrasts[2]]], paired = TRUE)
    data.frame(category = cat_lvl, p_value = tt$p.value, n_paired = nrow(wide))
  }) %>%
    bind_rows() %>%
    mutate(
      category = factor(category, levels = category_order),
      p_label = ifelse(
        is.na(p_value), "NA",
        ifelse(p_value < 0.001, "p < 0.001", sprintf("p = %.3f", p_value))
      ),
      p_fontface = ifelse(!is.na(p_value) & p_value < 0.05, "bold", "plain")
    )
  
  y_range <- range(c(data$logFC, summary$mean_logFC + summary$se_logFC))
  y_pad   <- diff(y_range) * 0.08
  
  # Categories with fewer than 2 paired proteins (e.g. none detected in
  # this dataset) get no bracket -- nothing to draw a comparison over.
  bracket_df <- data %>%
    group_by(category) %>%
    summarise(y_top = max(logFC), .groups = "drop") %>%
    right_join(stats, by = "category") %>%
    filter(!is.na(y_top), !is.na(p_value)) %>%
    mutate(
      x_num     = as.numeric(category),
      x_left    = x_num - dodge_width / 4,
      x_right   = x_num + dodge_width / 4,
      y_bracket = y_top + y_pad,
      y_text    = y_bracket + y_pad * 0.5
    )
  
  list(stats = stats, bracket_df = bracket_df, y_range = y_range, y_pad = y_pad)
}

# Bar + SE + jitter + significance-bracket plot, shared by both figures.
build_bar_plot <- function(summary, data, bracket_df, y_range, y_pad, n_labels, title, xlab) {
  ggplot() +
    geom_col(
      data = summary, aes(x = category, y = mean_logFC, fill = group),
      position = position_dodge(width = dodge_width), width = 0.6,
      color = "black", linewidth = 0.3
    ) +
    geom_errorbar(
      data = summary,
      aes(x = category, ymin = mean_logFC - se_logFC, ymax = mean_logFC + se_logFC, group = group),
      position = position_dodge(width = dodge_width), width = 0.2, linewidth = 0.4
    ) +
    geom_point(
      data = data, aes(x = category, y = logFC, group = group),
      position = position_jitterdodge(jitter.width = 0.15, dodge.width = dodge_width),
      size = 1.4, alpha = 0.5, color = "black"
    ) +
    scale_fill_manual(
      values = setNames(unname(group_colors[contrasts]), contrasts),
      labels = group_labels[contrasts], name = "Group"
    ) +
    scale_x_discrete(labels = setNames(n_labels$label, n_labels$category), drop = FALSE) +
    geom_text(
      data = n_labels %>% filter(n_detected == 0),
      aes(x = category, y = 0, label = "n.d."),
      size = 3.2, fontface = "italic", color = "grey40", vjust = -1, inherit.aes = FALSE
    ) +
    geom_hline(yintercept = 0, linewidth = 0.3, color = "grey40") +
    geom_segment(
      data = bracket_df, aes(x = x_left, xend = x_right, y = y_bracket, yend = y_bracket),
      linewidth = 0.4, inherit.aes = FALSE
    ) +
    geom_text(
      data = bracket_df, aes(x = x_num, y = y_text, label = p_label, fontface = p_fontface),
      size = 3.2, inherit.aes = FALSE
    ) +
    coord_cartesian(ylim = c(y_range[1] - y_pad * 0.5, y_range[2] + y_pad * 3)) +
    labs(x = xlab, y = expression(log[2] ~ "Fold Change (Post vs Pre)"), title = title) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "top"
    )
}

# ============================================================================
# PLOT 1: ETC Complex I-V structural subunits
# ============================================================================
cat("--- Plot 1: Classifying Complex I-V structural subunits ---\n")

subunit_regex <- "OXPHOS > Complex (I{1,3}|IV|V) > C[IVX]+ subunits"
mito_complex  <- classify_by_regex(mito, subunit_regex, complex_order, "Complex subunits")

cat(sprintf("  Classified %d proteins across Complexes I-V\n", nrow(mito_complex)))
print(table(factor(mito_complex$category, levels = complex_order)))
cat("\n")

cat("--- Plot 1: Merging with differential results ---\n")
etc_data    <- build_category_data(mito_complex, complex_order, "Complex subunits")
etc_summary <- summarize_data(etc_data)
etc_nlabels <- build_n_labels(etc_data, mito_complex, complex_order, "Complex")

cat("--- Plot 1: Testing BFR vs HLRT logFC per complex ---\n")
etc_test <- test_categories(etc_data, etc_summary, complex_order)
print(etc_test$stats)
cat("\n")

p_subunits <- build_bar_plot(
  summary = etc_summary, data = etc_data, bracket_df = etc_test$bracket_df,
  y_range = etc_test$y_range, y_pad = etc_test$y_pad, n_labels = etc_nlabels,
  title = "ETC Complex Subunit logFC: BFR vs HLRT", xlab = "ETC Complex"
)
print(p_subunits)

# ============================================================================
# PLOT 2: Assembly factors + other OXPHOS components (not in Plot 1)
# ============================================================================
cat("--- Plot 2: Classifying assembly factors + other OXPHOS components ---\n")

# Per-complex assembly factors (CI-CV assembly factors)
af_regex   <- "OXPHOS > Complex (I{1,3}|IV|V) > C[IVX]+ assembly factors"
mito_af    <- classify_by_regex(mito, af_regex, complex_order, "Assembly factors")

# "Other": OXPHOS-related proteins not tied to a single complex --
# Cytochrome C carriers (CYCS, HCCS) and general/respirasome assembly
# factors (HIGD2A, RAB5IF).
mito_other_raw <- mito %>%
  select(UniProt, Symbol, MitoCarta3.0_MitoPathways) %>%
  filter(!is.na(MitoCarta3.0_MitoPathways)) %>%
  filter(str_detect(MitoCarta3.0_MitoPathways, "Cytochrome C|Respirasome assembly")) %>%
  filter(!str_detect(MitoCarta3.0_MitoPathways, "Complex ")) %>%
  mutate(category = "Other") %>%
  distinct(UniProt, Symbol, category)

other_order  <- c(complex_order, "Other")
mito_other   <- bind_rows(mito_af, mito_other_raw)

cat(sprintf("  Classified %d proteins (assembly factors + other)\n", nrow(mito_other)))
print(table(factor(mito_other$category, levels = other_order)))
cat("\n")

cat("--- Plot 2: Merging with differential results ---\n")
other_data    <- build_category_data(mito_other, other_order, "Assembly factors + other")
other_summary <- summarize_data(other_data)
other_nlabels <- build_n_labels(other_data, mito_other, other_order, "")
other_nlabels$label <- sprintf(
  "%s\n(n = %d/%d)",
  ifelse(other_nlabels$category == "Other", "Other", paste0("Complex ", other_nlabels$category, " AF")),
  other_nlabels$n_detected, other_nlabels$n_total
)

cat("--- Plot 2: Testing BFR vs HLRT logFC per category ---\n")
other_test <- test_categories(other_data, other_summary, other_order)
print(other_test$stats)
cat("\n")

p_other <- build_bar_plot(
  summary = other_summary, data = other_data, bracket_df = other_test$bracket_df,
  y_range = other_test$y_range, y_pad = other_test$y_pad, n_labels = other_nlabels,
  title = "OXPHOS Assembly Factors & Other Components logFC: BFR vs HLRT",
  xlab = "Category (per-complex assembly factors, plus non-complex-specific OXPHOS components)"
)
print(p_other)

# --- Save outputs ------------------------
cat("--- Saving outputs ---\n")

ggsave(file.path(config$out_dir, "ETC_complex_subunit_logFC_BFR_vs_HLRT.pdf"), p_subunits, width = 7, height = 5.5)
ggsave(file.path(config$out_dir, "ETC_complex_subunit_logFC_BFR_vs_HLRT.png"), p_subunits, width = 7, height = 5.5, dpi = 300)
ggsave(file.path(config$out_dir, "ETC_assembly_other_logFC_BFR_vs_HLRT.pdf"), p_other, width = 8, height = 5.5)
ggsave(file.path(config$out_dir, "ETC_assembly_other_logFC_BFR_vs_HLRT.png"), p_other, width = 8, height = 5.5, dpi = 300)

write.csv(etc_data, file.path(config$out_dir, "ETC_complex_subunit_logFC_protein_level.csv"), row.names = FALSE)
write.csv(etc_summary, file.path(config$out_dir, "ETC_complex_subunit_logFC_summary.csv"), row.names = FALSE)
write.csv(other_data, file.path(config$out_dir, "ETC_assembly_other_logFC_protein_level.csv"), row.names = FALSE)
write.csv(other_summary, file.path(config$out_dir, "ETC_assembly_other_logFC_summary.csv"), row.names = FALSE)

cat(sprintf("  Saved plots + data to %s\n", config$out_dir))
cat("\n=== Done ===\n")