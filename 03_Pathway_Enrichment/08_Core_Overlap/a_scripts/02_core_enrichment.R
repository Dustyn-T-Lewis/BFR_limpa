# Core-enrichment (leading-edge) overlap between contrasts (BFR)
#
# For pathways significant in the specified contrasts, compares which
# core-enrichment genes drive the enrichment in each contrast.
# Sections 11-13 use the GO hierarchy to ask which child processes of each
# significant GO term carry its leading edge, and whether that differs
# between contrasts.
# Input: gsea_results_<label>.rds from 03_gsea_BFR.R

# --- 0: Setup --------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(dplyr, tidyr, purrr, stringr, readxl, writexl, ggplot2,
               ggvenn, patchwork, clusterProfiler, GO.db, AnnotationDbi,
               org.Hs.eg.db)
# Load (not attach) the package that defines the gseaResult class, so
# Bioconductor generics (S4Vectors::first, etc.) don't mask dplyr.
loadNamespace("DOSE")

config <- list(
  gsea_dir   = "03_Pathway_Enrichment/07_GSEA/c_data",
  id_map     = "03_Pathway_Enrichment/07_GSEA/b_reports/protein_id_mapping.xlsx",
  report_dir = "03_Pathway_Enrichment/08_Core_Overlap/b_reports",
  data_dir   = "03_Pathway_Enrichment/08_Core_Overlap/c_data",
  
  # Contrasts to compare (labels used in gsea_results_<label>.rds); 2 or more
  contrasts  = c("Training_BFR", "Training_HLRT"),
  
  # Which pathways to keep:
  #   "all" = significant in every listed contrast
  #   "any" = significant in at least one (core genes still compared in all,
  #           since GSEA was run with pvalueCutoff = 1)
  require    = "all",
  padj_cut   = 0.05,
  
  # Database label = slot name in the GSEA .rds
  databases  = c(Hallmark = "gsea_hall",
                 Reactome = "gsea_reactome",
                 KEGG     = "gsea_kegg",
                 GO_BP    = "gsea_bp_clean",
                 GO_CC    = "gsea_cc_clean",
                 GO_MF    = "gsea_mf_clean"),
  
  plot_top_n = 30,  # max pathways per database in the heatmaps / Venns
  
  # Venn diagrams (2-4 contrasts)
  venn_direction = "same",   # "same" = only pathways with the same NES sign in
  # every contrast; "all" = include opposite-direction
  venn_per_page  = 6,
  venn_ncol      = 3,
  venn_cols      = c("#E69F00", "#56B4E9", "#009E73", "#CC79A7"),
  
  # NES-vs-NES scatter (every tested pathway, one plot per contrast pair)
  scatter_databases       = c("Hallmark", "Reactome", "KEGG", "GO_BP", "GO_CC", "GO_MF"),
  scatter_include_neither = TRUE,  # FALSE = drop pathways not significant in either
  
  # Dotplot: top N significant pathways per contrast (per database)
  dot_top_n = 20,
  
  # GO hierarchy (sections 11-13). Names = ontology (BP/CC/MF), values = GSEA
  # slot to read. If the *_clean slots were made with simplify(), point this at
  # the unsimplified slot if you saved one: simplify() drops redundant child
  # terms, which makes parents look less "explained" by their children.
  hier_slots     = c(BP = "gsea_bp_clean"),
  hier_min_child = 3,    # min detected proteins for a child term to be tested
  hier_top_n     = 25,   # parents plotted in the decomposition PDF
  hier_max_kids  = 15,   # max child terms shown per parent in the plot
  
  # Venn-unique core proteins -> child terms (section 14)
  uniq_min      = 3,     # min unique core proteins required in EACH contrast
  uniq_perm     = 2000,  # permutations for the per-parent divergence test
  uniq_top_n    = 25,    # parents plotted
  uniq_max_kids = 15     # child terms shown per parent
)

stopifnot(length(config$contrasts) >= 2, config$require %in% c("all", "any"))

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

# Fails with a clear message if the config list in memory is older than the
# section being run (e.g. a section was added without re-running section 0).
check_config <- function(keys) {
  missing <- setdiff(keys, names(config))
  if (length(missing)) {
    stop("config is missing: ", paste(missing, collapse = ", "),
         "\nRe-run the `config <- list(...)` block in section 0.", call. = FALSE)
  }
}

tag <- paste(config$contrasts, collapse = "_vs_")
k   <- length(config$contrasts)


# --- 1: Load GSEA results -----------------------------------------------------
gsea <- set_names(config$contrasts) %>%
  map(~ readRDS(file.path(config$gsea_dir, sprintf("gsea_results_%s.rds", .x))))

tidy_gsea <- function(res, contrast) {
  imap_dfr(config$databases, function(slot, db) {
    obj <- res[[slot]]
    if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
    as_tibble(obj@result) %>%
      transmute(database = db, contrast = contrast, ID, Description,
                setSize, NES, p.adjust, core_enrichment)
  })
}

all_res <- imap_dfr(gsea, tidy_gsea)

# --- 2: Core-enrichment IDs -> gene symbols ------------------------------------
# GO / Hallmark / Reactome core genes are Ensembl; KEGG core genes are Entrez.
id_map <- read_excel(config$id_map)
id2sym <- bind_rows(
  id_map %>% transmute(gene_id = Ensembl,         symbol = gene),
  id_map %>% transmute(gene_id = as.character(Entrez), symbol = gene)
) %>%
  filter(!is.na(gene_id)) %>%
  distinct(gene_id, .keep_all = TRUE)


# --- 3: Select pathways ---------------------------------------------------------
sel <- all_res %>%
  group_by(database, ID) %>%
  summarise(n_sig = sum(p.adjust < config$padj_cut), .groups = "drop") %>%
  filter(if (config$require == "all") n_sig == k else n_sig >= 1)

cat(sprintf("\nPathways passing '%s' filter (padj < %.2f):\n",
            config$require, config$padj_cut))
print(count(sel, database))
if (nrow(sel) == 0) stop("No pathways meet the significance filter for these contrasts.")

core_long <- all_res %>%
  semi_join(sel, by = c("database", "ID")) %>%
  mutate(gene_id = str_split(core_enrichment, "/")) %>%
  dplyr::select(-core_enrichment) %>%
  unnest(gene_id) %>%
  left_join(id2sym, by = "gene_id") %>%
  mutate(symbol = coalesce(symbol, gene_id))

core_sets <- core_long %>%
  group_by(database, ID, Description, contrast) %>%
  summarise(setSize  = dplyr::first(setSize),
            NES      = dplyr::first(NES),
            p.adjust = dplyr::first(p.adjust),
            genes    = list(unique(symbol)),
            .groups  = "drop")


# --- 4: Pairwise overlap ---------------------------------------------------------
collapse_genes <- function(x) paste(sort(x), collapse = "/")

pair_overlap <- function(cA, cB) {
  a <- core_sets %>% filter(contrast == cA) %>% dplyr::select(-contrast)
  b <- core_sets %>% filter(contrast == cB) %>% dplyr::select(-contrast, -setSize)
  
  inner_join(a, b, by = c("database", "ID", "Description"),
             suffix = c("_A", "_B")) %>%
    mutate(
      contrast_A      = cA,
      contrast_B      = cB,
      direction       = if_else(sign(NES_A) == sign(NES_B), "same", "opposite"),
      n_A             = lengths(genes_A),
      n_B             = lengths(genes_B),
      n_shared        = map2_int(genes_A, genes_B, ~ length(intersect(.x, .y))),
      jaccard         = n_shared / (n_A + n_B - n_shared),
      overlap_coef    = n_shared / pmin(n_A, n_B),
      # expected shared genes if the two leading edges were random draws
      # from the pathway's detected members (descriptive, not a test)
      expected_shared = n_A * n_B / setSize,
      shared_genes    = map2_chr(genes_A, genes_B, ~ collapse_genes(intersect(.x, .y))),
      only_A          = map2_chr(genes_A, genes_B, ~ collapse_genes(setdiff(.x, .y))),
      only_B          = map2_chr(genes_A, genes_B, ~ collapse_genes(setdiff(.y, .x)))
    ) %>%
    dplyr::select(database, ID, Description, contrast_A, contrast_B, setSize,
                  NES_A, NES_B, p.adjust_A, p.adjust_B, direction,
                  n_A, n_B, n_shared, expected_shared, jaccard, overlap_coef,
                  shared_genes, only_A, only_B)
}

pairwise <- combn(config$contrasts, 2, simplify = FALSE) %>%
  map_dfr(~ pair_overlap(.x[1], .x[2])) %>%
  arrange(database, dplyr::desc(jaccard))


# --- 5: Gene membership across all contrasts ---------------------------------------
membership <- core_long %>%
  distinct(database, ID, Description, contrast, symbol) %>%
  mutate(in_core = TRUE) %>%
  pivot_wider(names_from = contrast, values_from = in_core, values_fill = FALSE)

for (cc in setdiff(config$contrasts, names(membership))) membership[[cc]] <- FALSE

membership <- membership %>%
  mutate(n_contrasts = rowSums(across(all_of(config$contrasts)))) %>%
  arrange(database, ID, dplyr::desc(n_contrasts), symbol)

shared_all <- membership %>%
  group_by(database, ID, Description) %>%
  summarise(n_core_union  = n(),
            n_core_in_all = sum(n_contrasts == k),
            prop_in_all   = n_core_in_all / n_core_union,
            core_in_all   = collapse_genes(symbol[n_contrasts == k]),
            .groups = "drop") %>%
  arrange(database, dplyr::desc(prop_in_all))

# Wide NES / padj per pathway for reference
pathway_summary <- core_sets %>%
  dplyr::select(database, ID, Description, setSize, contrast, NES, p.adjust) %>%
  pivot_wider(names_from = contrast, values_from = c(NES, p.adjust),
              names_glue = "{contrast}_{.value}")


# --- 6: Save tables --------------------------------------------------------------
write_xlsx(
  list(pathway_summary = pathway_summary,
       pairwise_overlap = pairwise,
       shared_all       = shared_all,
       gene_membership  = membership),
  file.path(config$data_dir, sprintf("core_overlap_%s.xlsx", tag))
)
saveRDS(list(config = config, pairwise = pairwise, shared_all = shared_all,
             membership = membership, core_sets = core_sets),
        file.path(config$data_dir, sprintf("core_overlap_%s.rds", tag)))


# --- 7: Heatmap per database (Jaccard; label = n shared, * = opposite NES) -------
plot_db <- function(db) {
  d <- pairwise %>% filter(database == db)
  if (nrow(d) == 0) return(invisible(NULL))
  
  keep <- d %>%
    group_by(ID) %>%
    summarise(best = min(pmax(p.adjust_A, p.adjust_B)), .groups = "drop") %>%
    slice_min(best, n = config$plot_top_n, with_ties = FALSE) %>%
    pull(ID)
  
  d <- d %>%
    filter(ID %in% keep) %>%
    mutate(pair  = paste(contrast_A, "vs", contrast_B),
           label = paste0(n_shared, if_else(direction == "opposite", "*", "")),
           Description = str_trunc(Description, 60))
  
  p <- ggplot(d, aes(pair, reorder(Description, jaccard), fill = jaccard)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = label), size = 3) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1)) +
    labs(title = sprintf("%s: core-enrichment overlap", db),
         x = NULL, y = NULL, fill = "Jaccard",
         caption = "Label = shared core genes; * = opposite NES direction") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid = element_blank())
  
  ggsave(file.path(config$report_dir, sprintf("core_overlap_%s_%s.pdf", tag, db)),
         p, width = 4 + 1.2 * n_distinct(d$pair),
         height = 2 + 0.22 * n_distinct(d$ID), limitsize = FALSE)
}

walk(unique(pairwise$database), plot_db)


# --- 8: Venn diagrams of core-enrichment genes per pathway ------------------------
# One PDF per database, several pathways per page. Counts are core genes;
# the caption lists the genes shared by every contrast.
check_config(c("venn_direction", "venn_per_page", "venn_ncol", "venn_cols"))

venn_sets <- core_sets %>%
  group_by(database, ID, Description) %>%
  filter(n() == k) %>%                                  # present in every contrast
  summarise(
    worst_padj = max(p.adjust),
    same_dir   = n_distinct(sign(NES)) == 1,
    sets       = list(set_names(genes, contrast)[config$contrasts]),
    stats      = paste(sprintf("%s: NES %.2f, padj %.2g", contrast, NES, p.adjust),
                       collapse = "\n"),
    .groups    = "drop"
  )

make_venn <- function(sets, title, subtitle) {
  shared     <- Reduce(intersect, sets)
  shared_txt <- if (length(shared)) paste(sort(shared), collapse = ", ") else "none"
  
  ggvenn(sets,
         fill_color      = config$venn_cols[seq_along(sets)],
         fill_alpha      = 0.4,
         stroke_size     = 0.4,
         set_name_size   = 3,
         text_size       = 3.2,
         show_percentage = FALSE) +
    labs(title    = str_wrap(title, 45),
         subtitle = subtitle,
         caption  = str_wrap(paste("Shared:", str_trunc(shared_txt, 300)), 70)) +
    theme(plot.title    = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 7.5),
          plot.caption  = element_text(size = 7, hjust = 0))
}

plot_venns_db <- function(db) {
  d <- venn_sets %>% filter(database == db)
  if (config$venn_direction == "same") d <- d %>% filter(same_dir)
  if (nrow(d) == 0) return(invisible(NULL))
  
  d <- d %>% slice_min(worst_padj, n = config$plot_top_n, with_ties = FALSE)
  
  plots <- pmap(list(d$sets, d$Description, d$stats), make_venn)
  pages <- split(plots, ceiling(seq_along(plots) / config$venn_per_page))
  nrow  <- ceiling(config$venn_per_page / config$venn_ncol)
  
  pdf(file.path(config$report_dir, sprintf("core_venn_%s_%s.pdf", tag, db)),
      width = 4 * config$venn_ncol, height = 4.5 * nrow)
  walk(pages, ~ print(wrap_plots(.x, ncol = config$venn_ncol, nrow = nrow)))
  dev.off()
}

if (k <= 4) {
  walk(unique(venn_sets$database), plot_venns_db)
} else {
  message("Venn diagrams skipped: ggvenn supports at most 4 contrasts ",
          "(use the gene_membership sheet for an UpSet plot instead).")
}

# --- 9: NES vs NES scatter of all tested pathways (per contrast pair) --------------
# x = NES in contrast A, y = NES in contrast B; colour = database;
# shape = significant (padj < padj_cut) in both / one / neither contrast.
# Quadrant labels = number and % of all plotted pathways in that quadrant.
check_config(c("scatter_databases", "scatter_include_neither"))

nes_pair_df <- function(cA, cB) {
  a <- all_res %>%
    filter(contrast == cA, database %in% config$scatter_databases) %>%
    dplyr::select(database, ID, Description, NES_A = NES, padj_A = p.adjust)
  b <- all_res %>%
    filter(contrast == cB, database %in% config$scatter_databases) %>%
    dplyr::select(database, ID, NES_B = NES, padj_B = p.adjust)
  
  sig_levels <- c("Both", paste(cA, "only"), paste(cB, "only"), "Neither")
  
  inner_join(a, b, by = c("database", "ID")) %>%
    filter(!is.na(NES_A), !is.na(NES_B)) %>%
    mutate(
      contrast_A = cA,
      contrast_B = cB,
      sig = case_when(
        padj_A < config$padj_cut & padj_B < config$padj_cut ~ "Both",
        padj_A < config$padj_cut                            ~ paste(cA, "only"),
        padj_B < config$padj_cut                            ~ paste(cB, "only"),
        TRUE                                                ~ "Neither"),
      sig = factor(sig, levels = sig_levels),
      quadrant = factor(case_when(
        NES_A > 0 & NES_B > 0 ~ "Up / Up",
        NES_A < 0 & NES_B > 0 ~ "Down / Up",
        NES_A < 0 & NES_B < 0 ~ "Down / Down",
        TRUE                  ~ "Up / Down"),
        levels = c("Up / Up", "Down / Up", "Down / Down", "Up / Down")),
      database = factor(database, levels = config$scatter_databases)
    ) %>%
    { if (config$scatter_include_neither) . else filter(., sig != "Neither") }
}

quadrant_counts <- function(d) {
  d %>%
    dplyr::count(quadrant, .drop = FALSE) %>%
    mutate(pct   = 100 * n / sum(n),
           label = sprintf("%d (%.1f%%)", n, pct),
           x     = if_else(quadrant %in% c("Up / Up", "Up / Down"),   Inf, -Inf),
           y     = if_else(quadrant %in% c("Up / Up", "Down / Up"),   Inf, -Inf),
           hjust = if_else(x > 0, 1.1, -0.1),
           vjust = if_else(y > 0, 1.5, -0.5))
}

plot_nes_scatter <- function(cA, cB) {
  d <- nes_pair_df(cA, cB)
  if (nrow(d) == 0) return(NULL)
  
  q    <- quadrant_counts(d)
  lim  <- max(abs(c(d$NES_A, d$NES_B)), na.rm = TRUE) * 1.05
  rho  <- cor(d$NES_A, d$NES_B, method = "spearman")
  lvls <- levels(d$sig)
  
  sig_title <- sprintf("Significant\n(padj < %.2f)", config$padj_cut)
  
  p <- ggplot(d, aes(NES_A, NES_B)) +
    geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.4) +
    geom_abline(slope = 1, intercept = 0, colour = "grey80", linetype = "dashed") +
    geom_point(aes(colour = database, shape = sig, size = sig, alpha = sig)) +
    geom_text(data = q, aes(x = x, y = y, label = label, hjust = hjust, vjust = vjust),
              inherit.aes = FALSE, size = 3.6, fontface = "bold") +
    scale_colour_brewer(palette = "Dark2", drop = FALSE) +
    scale_shape_manual(values = setNames(c(17, 16, 15, 1),       lvls), drop = FALSE, name = sig_title) +
    scale_size_manual( values = setNames(c(2.6, 2, 2, 1.2),      lvls), drop = FALSE, name = sig_title) +
    scale_alpha_manual(values = setNames(c(0.95, 0.8, 0.8, 0.3), lvls), drop = FALSE, name = sig_title) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
    labs(x = paste("NES:", cA), y = paste("NES:", cB), colour = "Database",
         title = sprintf("%s vs %s: pathway NES", cA, cB),
         subtitle = sprintf("%d pathways | Spearman rho = %.2f", nrow(d), rho)) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())
  
  ggsave(file.path(config$report_dir, sprintf("nes_scatter_%s_vs_%s.pdf", cA, cB)),
         p, width = 9, height = 7.5)
  
  list(data = d, quadrants = q %>% dplyr::select(quadrant, n, pct) %>%
         mutate(contrast_A = cA, contrast_B = cB, .before = 1))
}

scatter_out <- combn(config$contrasts, 2, simplify = FALSE) %>%
  map(~ plot_nes_scatter(.x[1], .x[2])) %>%
  compact()

write_xlsx(
  list(quadrant_counts = map_dfr(scatter_out, "quadrants"),
       pathway_nes     = map_dfr(scatter_out, "data") %>%
         mutate(across(c(sig, quadrant, database), as.character))),
  file.path(config$data_dir, sprintf("nes_scatter_%s.xlsx", tag))
)

# --- 10: Dotplot of top pathways per contrast --------------------------------------
# One PDF per database. Takes the top `dot_top_n` significant pathways (lowest
# padj) from EACH contrast, pools them, and shows every contrast's result for
# the pooled set. Fill = NES, size = -log10(padj), faded = not significant.
check_config(c("dot_top_n"))

dot_ids <- all_res %>%
  filter(p.adjust < config$padj_cut) %>%
  group_by(database, contrast) %>%
  slice_min(p.adjust, n = config$dot_top_n, with_ties = FALSE) %>%
  ungroup() %>%
  distinct(database, ID)

dot_df <- all_res %>%
  semi_join(dot_ids, by = c("database", "ID")) %>%
  mutate(contrast      = factor(contrast, levels = config$contrasts),
         sig           = factor(p.adjust < config$padj_cut, levels = c(TRUE, FALSE)),
         neglog10_padj = -log10(p.adjust))

plot_dot_db <- function(db) {
  d <- dot_df %>% filter(database == db)
  if (nrow(d) == 0) return(invisible(NULL))
  
  # Order pathways by mean NES across contrasts
  ord <- d %>%
    group_by(ID) %>%
    summarise(m = mean(NES), .groups = "drop") %>%
    arrange(m) %>%
    pull(ID)
  labs_y <- d %>%
    distinct(ID, Description) %>%
    { setNames(str_trunc(.$Description, 60), .$ID) }
  lim <- max(abs(d$NES))
  
  p <- ggplot(d, aes(contrast, factor(ID, levels = ord))) +
    geom_point(aes(size = neglog10_padj, fill = NES, alpha = sig),
               shape = 21, colour = "grey20", stroke = 0.3) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-lim, lim)) +
    scale_size_continuous(range = c(1.5, 6), name = "-log10(padj)") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.3), drop = FALSE,
                       labels = c(sprintf("padj < %.2f", config$padj_cut), "n.s."),
                       name = NULL) +
    scale_y_discrete(labels = labs_y) +
    labs(title = sprintf("%s: top %d pathways per contrast", db, config$dot_top_n),
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.minor = element_blank())
  
  ggsave(file.path(config$report_dir, sprintf("dotplot_%s_%s.pdf", tag, db)),
         p, width = 6 + 0.8 * k, height = 2 + 0.2 * length(ord), limitsize = FALSE)
}

walk(unique(dot_df$database), plot_dot_db)


# --- 11: GO hierarchy: are parents explained by their significant descendants? ----
# For each selected GO term (from section 3) and each contrast, finds the
# significant descendant terms (via GO*OFFSPRING) and asks what fraction of the
# parent's core genes they cover. Low coverage = the parent carries signal that
# no more specific significant term captures ("unexplained" core genes).
check_config(c("hier_slots", "hier_min_child", "hier_top_n", "hier_max_kids"))

go_env <- list(
  BP = list(offspring = GOBPOFFSPRING, children = GOBPCHILDREN),
  CC = list(offspring = GOCCOFFSPRING, children = GOCCCHILDREN),
  MF = list(offspring = GOMFOFFSPRING, children = GOMFCHILDREN)
)
stopifnot(all(names(config$hier_slots) %in% names(go_env)))

to_sym   <- function(ids) coalesce(id2sym$symbol[match(ids, id2sym$gene_id)], ids)
go_terms <- function(ids) if (length(ids)) unname(AnnotationDbi::Term(ids)) else character(0)

# Core IDs kept as Ensembl here so they match GO annotation; symbols for output
hier_res <- imap_dfr(config$hier_slots, function(slot, ont) {
  imap_dfr(gsea, function(res, cc) {
    obj <- res[[slot]]
    if (is.null(obj) || nrow(obj@result) == 0) return(NULL)
    as_tibble(obj@result) %>%
      transmute(ont = ont, contrast = cc, ID, Description, setSize, NES, p.adjust,
                core = str_split(core_enrichment, "/"))
  })
})

bad_ids <- unique(hier_res$ID[!str_detect(hier_res$ID, "^GO:\\d{7}$")])
if (length(bad_ids)) {
  stop("GO slots must use GO IDs (e.g. GO:0006936); found: ",
       paste(head(bad_ids, 3), collapse = ", "), call. = FALSE)
}

hier_sig <- hier_res %>% filter(p.adjust < config$padj_cut)

# Parents = GO pathways selected in section 3
sel_go <- sel %>%
  mutate(ont = str_remove(database, "^GO_")) %>%
  filter(ont %in% names(config$hier_slots))

# Descendant lookup, cached once per ontology (takes a few seconds)
off_all <- map(set_names(names(config$hier_slots)), ~ as.list(go_env[[.x]]$offspring))

coverage <- hier_sig %>%
  semi_join(sel_go, by = c("ont", "ID")) %>%
  group_by(ont, contrast) %>%
  group_modify(function(par, key) {
    sig_here  <- hier_sig %>% filter(ont == key$ont, contrast == key$contrast)
    core_here <- set_names(sig_here$core, sig_here$ID)
    pmap_dfr(par, function(ID, Description, core, ...) {
      kids    <- intersect(off_all[[key$ont]][[ID]], sig_here$ID)
      covered <- unique(unlist(core_here[kids]))
      tibble(ID, Description,
             n_core          = length(core),
             n_sig_desc      = length(kids),
             frac_explained  = mean(core %in% covered),
             sig_descendants = paste(go_terms(kids), collapse = " | "),
             unexplained     = collapse_genes(to_sym(setdiff(core, covered))))
    })
  }) %>%
  ungroup() %>%
  arrange(ont, contrast, frac_explained)

coverage_wide <- coverage %>%
  dplyr::select(ont, ID, Description, contrast, n_core, n_sig_desc, frac_explained) %>%
  pivot_wider(names_from = contrast,
              values_from = c(n_core, n_sig_desc, frac_explained),
              names_glue = "{contrast}_{.value}")


# --- 12: GO hierarchy: decompose each parent's core genes into its children ----
# For each selected parent and contrast: universe = the parent's proteins
# detected in that contrast's ranked list; for each direct child (is_a,
# part_of, regulates ...), a hypergeometric test asks whether the parent's
# core genes are concentrated in that child. Child gene sets include
# annotations propagated from the child's own descendants.
# Descriptive: nested, non-independent tests on small sets.

goall_keys <- AnnotationDbi::keys(org.Hs.eg.db, keytype = "GOALL")  # cache once

annot_sets <- function(ids, ont) {
  # Children with no human annotations (or newer than the installed
  # org.Hs.eg.db) aren't valid keys, and select() errors if none are.
  ids <- intersect(ids, goall_keys)
  if (!length(ids)) return(list())
  suppressMessages(
    AnnotationDbi::select(org.Hs.eg.db, keys = ids, keytype = "GOALL",
                          columns = "ENSEMBL")
  ) %>%
    filter(ONTOLOGYALL == ont, !is.na(ENSEMBL)) %>%
    { split(.$ENSEMBL, .$GOALL) }
}

# Prefer the gene sets stored in the gseaResult (same annotation the GSEA used);
# fall back to org.Hs.eg.db for terms not stored there.
term_sets <- function(ids, obj, ont) {
  gs   <- obj@geneSets[intersect(ids, names(obj@geneSets))]
  miss <- setdiff(ids, names(gs))
  c(gs, annot_sets(miss, ont))
}

decompose_term <- function(parent_id, contrast, ont,
                           slot      = config$hier_slots[[ont]],
                           min_child = config$hier_min_child) {
  obj <- gsea[[contrast]][[slot]]
  row <- obj@result[obj@result$ID == parent_id, ]
  if (nrow(row) == 0) return(NULL)
  
  core     <- str_split(row$core_enrichment, "/")[[1]]
  universe <- intersect(term_sets(parent_id, obj, ont)[[parent_id]], names(obj@geneList))
  core     <- intersect(core, universe)
  
  kids <- AnnotationDbi::mget(parent_id, go_env[[ont]]$children, ifnotfound = NA)[[1]]
  kids <- kids[!is.na(kids)]
  if (!length(kids) || !length(core)) return(NULL)
  # Named character vector (child ID -> relation type[s]); a child can be linked
  # by more than one relation, e.g. is_a and part_of
  relation <- vapply(split(names(kids), unname(kids)),
                     function(r) paste(unique(r), collapse = ","), character(1))
  
  sets <- term_sets(unique(unname(kids)), obj, ont) %>%
    map(~ intersect(.x, universe)) %>%
    keep(~ length(.x) >= min_child)
  if (!length(sets)) return(NULL)
  
  U <- length(universe); n <- length(core)
  imap_dfr(sets, function(g, child) {
    hits <- intersect(core, g)
    m <- length(g); x <- length(hits)
    tibble(parent_ID        = parent_id,
           parent           = go_terms(parent_id),
           contrast         = contrast,
           child_ID         = child,
           child            = go_terms(child),
           relation         = unname(relation[child]),
           parent_detected  = U,
           parent_core      = n,
           child_detected   = m,
           core_in_child    = x,
           frac_parent_core = x / n,              # share of parent's leading edge
           fold_enrich      = (x / n) / (m / U),
           p                = phyper(x - 1, m, U - m, n, lower.tail = FALSE),
           child_core_genes = collapse_genes(to_sym(hits)))
  }) %>%
    mutate(p.adjust = stats::p.adjust(p, method = "BH")) %>%
    arrange(p)
}

# Usable interactively too, e.g. to drill down a level:
#   decompose_term("GO:0006936", "Training_BFR", "BP")
decomp <- sel_go %>%
  dplyr::select(ont, ID) %>%
  tidyr::crossing(contrast = config$contrasts) %>%
  pmap_dfr(function(ont, ID, contrast) decompose_term(ID, contrast, ont))

# Side-by-side: does each contrast's leading edge sit in the same children?
decomp_wide <- decomp %>%
  dplyr::select(parent_ID, parent, child_ID, child, relation, contrast,
                child_detected, core_in_child, frac_parent_core, fold_enrich,
                p.adjust, child_core_genes) %>%
  pivot_wider(names_from = contrast,
              values_from = c(child_detected, core_in_child, frac_parent_core,
                              fold_enrich, p.adjust, child_core_genes),
              names_glue = "{contrast}_{.value}")

write_xlsx(
  list(coverage_by_contrast = coverage,
       coverage_wide        = coverage_wide,
       decomposition        = decomp,
       decomposition_wide   = decomp_wide),
  file.path(config$data_dir, sprintf("go_hierarchy_%s.xlsx", tag))
)
saveRDS(list(coverage = coverage, decomp = decomp),
        file.path(config$data_dir, sprintf("go_hierarchy_%s.rds", tag)))


# --- 13: Plot parent -> child decomposition ------------------------------------------
# One page per parent (top `hier_top_n` by worst padj across contrasts).
# Rows = child terms [relation]; fill = share of the parent's core genes in that
# child; size = number of core genes; black outline = child padj < padj_cut.

plot_decomp <- function(pid) {
  d <- decomp %>%
    filter(parent_ID == pid) %>%
    group_by(child_ID) %>%
    filter(max(core_in_child) > 0) %>%
    ungroup()
  if (nrow(d) == 0) return(NULL)
  
  keep <- d %>%
    group_by(child_ID) %>%
    summarise(m = max(frac_parent_core), .groups = "drop") %>%
    slice_max(m, n = config$hier_max_kids, with_ties = FALSE) %>%
    pull(child_ID)
  
  d <- d %>%
    filter(child_ID %in% keep) %>%
    mutate(contrast = factor(contrast, levels = config$contrasts),
           sig      = factor(p.adjust < config$padj_cut, levels = c(TRUE, FALSE)),
           label    = str_trunc(sprintf("%s [%s]", child, relation), 70))
  
  ggplot(d, aes(contrast, reorder(label, frac_parent_core))) +
    geom_point(aes(size = core_in_child, fill = frac_parent_core, colour = sig),
               shape = 21, stroke = 0.8) +
    scale_fill_gradient(low = "#f7fbff", high = "#08519c", limits = c(0, 1),
                        name = "Share of parent\ncore genes") +
    scale_colour_manual(values = c(`TRUE` = "black", `FALSE` = "grey75"),
                        labels = c(sprintf("padj < %.2f", config$padj_cut), "n.s."),
                        drop = FALSE, name = "Child test") +
    scale_size_area(max_size = 8, name = "Core genes\nin child") +
    labs(title    = str_wrap(dplyr::first(d$parent), 60),
         subtitle = sprintf("%s | core genes: %s", pid,
                            paste(sprintf("%s %d", config$contrasts,
                                          tapply(d$parent_core, d$contrast, dplyr::first)),
                                  collapse = ", ")),
         x = NULL, y = NULL) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          panel.grid.minor = element_blank())
}

if (nrow(decomp) > 0) {
  top_parents <- core_sets %>%
    filter(ID %in% unique(decomp$parent_ID)) %>%
    group_by(ID) %>%
    summarise(worst = max(p.adjust), .groups = "drop") %>%
    slice_min(worst, n = config$hier_top_n, with_ties = FALSE) %>%
    pull(ID)
  
  pdf(file.path(config$report_dir, sprintf("go_decomposition_%s.pdf", tag)),
      width = 8.5, height = 6)
  walk(top_parents, function(pid) {
    p <- plot_decomp(pid)
    if (!is.null(p)) print(p)
  })
  dev.off()
} else {
  message("No GO parents had testable child terms; decomposition plot skipped.")
}

# --- 14: Venn-unique core proteins -> child terms ---------------------------------
# For each GO parent in the Venn diagrams (same direction filter as section 8)
# and each contrast pair, splits the core genes into the Venn regions
# (only A / shared / only B) and asks whether the proteins unique to each
# contrast fall in different child terms.
#   Per child:  Fisher's exact test, only-A vs only-B, in child vs not.
#               OR > 1 = child over-represented among A-unique proteins.
#   Per parent: divergence = sum over children of |frac(only A in child) -
#               frac(only B in child)|, incl. a "(no child term)" bin; p-value
#               from permuting the A/B labels among the unique proteins
#               (valid even though children share proteins).
#   Threshold check: median position of each contrast's unique proteins in the
#               OTHER contrast's ranked list (0 = at its enriched end). Values
#               near 0 mean they just missed the other leading edge.
check_config(c("uniq_min", "uniq_perm", "uniq_top_n", "uniq_max_kids"))

child_sets_for <- function(parent_id, obj, ont, universe, min_child) {
  kids <- AnnotationDbi::mget(parent_id, go_env[[ont]]$children, ifnotfound = NA)[[1]]
  kids <- kids[!is.na(kids)]
  if (!length(kids)) return(list(sets = list(), relation = character(0)))
  relation <- vapply(split(names(kids), unname(kids)),
                     function(r) paste(unique(r), collapse = ","), character(1))
  sets <- term_sets(unique(unname(kids)), obj, ont) %>%
    map(~ intersect(.x, universe)) %>%
    keep(~ length(.x) >= min_child)
  list(sets = sets, relation = relation)
}

# Median percentile of `ids` in obj's ranked list, measured from the end the
# pathway is enriched at (top if NES > 0, bottom if NES < 0)
edge_pct <- function(ids, obj, nes) {
  gl  <- names(sort(obj@geneList, decreasing = TRUE))
  pct <- match(ids, gl) / length(gl)
  if (nes < 0) pct <- 1 - pct
  median(pct, na.rm = TRUE)
}

uniq_decompose <- function(parent_id, ont, cA, cB) {
  slot <- config$hier_slots[[ont]]
  oA <- gsea[[cA]][[slot]]; oB <- gsea[[cB]][[slot]]
  rA <- oA@result[oA@result$ID == parent_id, ]
  rB <- oB@result[oB@result$ID == parent_id, ]
  if (!nrow(rA) || !nrow(rB)) return(NULL)
  
  universe <- intersect(term_sets(parent_id, oA, ont)[[parent_id]],
                        union(names(oA@geneList), names(oB@geneList)))
  coreA  <- intersect(str_split(rA$core_enrichment, "/")[[1]], universe)
  coreB  <- intersect(str_split(rB$core_enrichment, "/")[[1]], universe)
  onlyA  <- setdiff(coreA, coreB)
  onlyB  <- setdiff(coreB, coreA)
  shared <- intersect(coreA, coreB)
  if (length(onlyA) < config$uniq_min || length(onlyB) < config$uniq_min) return(NULL)
  
  cs <- child_sets_for(parent_id, oA, ont, universe, config$hier_min_child)
  if (!length(cs$sets)) return(NULL)
  
  genes <- c(onlyA, onlyB)
  isA   <- seq_along(genes) <= length(onlyA)
  M <- vapply(cs$sets, function(g) genes %in% g, logical(length(genes)))
  M <- matrix(M, nrow = length(genes), dimnames = list(genes, names(cs$sets)))
  M <- cbind(M, `(no child term)` = rowSums(M) == 0)
  
  div_stat <- function(a) sum(abs(colMeans(M[a, , drop = FALSE]) -
                                    colMeans(M[!a, , drop = FALSE])))
  obs      <- div_stat(isA)
  perm     <- replicate(config$uniq_perm, div_stat(sample(isA)))
  p_global <- (1 + sum(perm >= obs - 1e-12)) / (1 + config$uniq_perm)
  
  in_any <- unique(unlist(cs$sets))
  map_dfr(colnames(M), function(ch) {
    is_none  <- ch == "(no child term)"
    inA      <- sum(M[isA, ch]); inB <- sum(M[!isA, ch])
    sh_genes <- if (is_none) setdiff(shared, in_any) else intersect(shared, cs$sets[[ch]])
    ft <- fisher.test(matrix(c(inA, length(onlyA) - inA,
                               inB, length(onlyB) - inB), nrow = 2))
    tibble(contrast_A      = cA, contrast_B = cB,
           parent_ID       = parent_id, parent = go_terms(parent_id),
           child_ID        = ch,
           child           = if (is_none) ch else go_terms(ch),
           relation        = if (is_none) NA_character_ else unname(cs$relation[ch]),
           child_detected  = if (is_none) NA_integer_ else length(cs$sets[[ch]]),
           n_only_A        = length(onlyA), n_only_B = length(onlyB),
           n_shared        = length(shared),
           only_A_in_child = inA, only_B_in_child = inB,
           shared_in_child = length(sh_genes),
           frac_only_A     = inA / length(onlyA),
           frac_only_B     = inB / length(onlyB),
           odds_ratio      = unname(ft$estimate),
           p               = ft$p.value,
           genes_only_A    = collapse_genes(to_sym(genes[isA & M[, ch]])),
           genes_only_B    = collapse_genes(to_sym(genes[!isA & M[, ch]])),
           genes_shared    = collapse_genes(to_sym(sh_genes)),
           divergence      = obs,
           p_global        = p_global,
           only_A_pct_in_B = edge_pct(onlyA, oB, rB$NES),
           only_B_pct_in_A = edge_pct(onlyB, oA, rA$NES))
  }) %>%
    mutate(p.adjust = stats::p.adjust(p, method = "BH"),
           favours  = case_when(frac_only_A > frac_only_B ~ cA,
                                frac_only_B > frac_only_A ~ cB,
                                TRUE ~ "neither"))
}

uniq_parents <- venn_sets %>%
  filter(database %in% paste0("GO_", names(config$hier_slots)))
if (config$venn_direction == "same") uniq_parents <- uniq_parents %>% filter(same_dir)

contrast_pairs <- combn(config$contrasts, 2, simplify = FALSE) %>%
  map_dfr(~ tibble(cA = .x[1], cB = .x[2]))

set.seed(1)
uniq_child <- uniq_parents %>%
  transmute(ont = str_remove(database, "^GO_"), ID) %>%
  tidyr::crossing(contrast_pairs) %>%
  pmap_dfr(function(ont, ID, cA, cB) uniq_decompose(ID, ont, cA, cB))

if (nrow(uniq_child) > 0) {
  # One row per parent x pair; BH across parents for the divergence test
  uniq_summary <- uniq_child %>%
    group_by(contrast_A, contrast_B, parent_ID, parent) %>%
    summarise(
      n_only_A        = dplyr::first(n_only_A),
      n_only_B        = dplyr::first(n_only_B),
      n_shared        = dplyr::first(n_shared),
      n_children      = sum(child_ID != "(no child term)"),
      divergence      = dplyr::first(divergence),
      p_global        = dplyr::first(p_global),
      only_A_pct_in_B = dplyr::first(only_A_pct_in_B),
      only_B_pct_in_A = dplyr::first(only_B_pct_in_A),
      top_child_A     = child[which.max(frac_only_A - frac_only_B)],
      top_child_B     = child[which.max(frac_only_B - frac_only_A)],
      n_child_padj    = sum(p.adjust < config$padj_cut),
      .groups = "drop") %>%
    group_by(contrast_A, contrast_B) %>%
    mutate(p_global_adj = stats::p.adjust(p_global, method = "BH")) %>%
    ungroup() %>%
    arrange(p_global)
  
  write_xlsx(
    list(summary_by_parent = uniq_summary,
         child_detail      = uniq_child %>% arrange(parent_ID, p)),
    file.path(config$data_dir, sprintf("go_unique_children_%s.xlsx", tag))
  )
  
  # Diverging bars: left = proteins unique to A, right = unique to B,
  # coloured to match the Venns; label gives shared count; * = child padj < cut
  plot_uniq <- function(pid, cA, cB) {
    d <- uniq_child %>%
      filter(parent_ID == pid, contrast_A == cA, contrast_B == cB,
             only_A_in_child + only_B_in_child > 0) %>%
      slice_max(only_A_in_child + only_B_in_child, n = config$uniq_max_kids,
                with_ties = FALSE) %>%
      mutate(label = str_trunc(sprintf("%s%s (shared %d)",
                                       if_else(p.adjust < config$padj_cut, "* ", ""),
                                       child, shared_in_child), 70),
             label = reorder(label, frac_only_B - frac_only_A))
    if (nrow(d) == 0) return(NULL)
    s <- uniq_summary %>% filter(parent_ID == pid, contrast_A == cA, contrast_B == cB)
    
    bars <- bind_rows(
      d %>% transmute(label, contrast = cA, n = -only_A_in_child),
      d %>% transmute(label, contrast = cB, n =  only_B_in_child)) %>%
      mutate(contrast = factor(contrast, levels = c(cA, cB)))
    cols <- setNames(config$venn_cols[match(c(cA, cB), config$contrasts)], c(cA, cB))
    
    ggplot(bars, aes(n, label, fill = contrast)) +
      geom_col(width = 0.7, alpha = 0.8) +
      geom_vline(xintercept = 0, colour = "grey40") +
      scale_x_continuous(labels = abs) +
      scale_fill_manual(values = cols, name = "Unique core\nproteins in") +
      labs(title = str_wrap(dplyr::first(d$parent), 60),
           subtitle = sprintf(
             "%s | only %s: %d, only %s: %d, shared: %d\ndivergence %.2f, permutation p = %.3g (BH %.3g)",
             pid, cA, s$n_only_A, cB, s$n_only_B, s$n_shared,
             s$divergence, s$p_global, s$p_global_adj),
           x = "Unique core proteins in child term", y = NULL,
           caption = "* = child Fisher padj < cutoff; bins overlap (a protein can sit in several children)") +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank(),
            plot.subtitle = element_text(size = 8))
  }
  
  top_uniq <- uniq_summary %>% slice_min(p_global, n = config$uniq_top_n, with_ties = FALSE)
  pdf(file.path(config$report_dir, sprintf("go_unique_children_%s.pdf", tag)),
      width = 9, height = 6)
  pwalk(top_uniq %>% dplyr::select(parent_ID, contrast_A, contrast_B),
        function(parent_ID, contrast_A, contrast_B) {
          p <- plot_uniq(parent_ID, contrast_A, contrast_B)
          if (!is.null(p)) print(p)
        })
  dev.off()
} else {
  message("No GO parents had >= ", config$uniq_min,
          " unique core proteins in both contrasts; section 14 skipped.")
}

cat(sprintf("\nDone. Tables in %s, plots in %s\n",
            config$data_dir, config$report_dir))