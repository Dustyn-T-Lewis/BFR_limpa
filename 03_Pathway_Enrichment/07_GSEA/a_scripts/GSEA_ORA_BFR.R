# GSEA / ORA Script (BFR) -- limpa input

# --- 0: Setup --------------------------------------------------------------
setwd(rprojroot::find_rstudio_root_file())
getwd()

if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(readxl, readr, dplyr, tidyr, stringr, purrr, limma,
               writexl, clusterProfiler, org.Hs.eg.db, enrichplot,
               msigdbr, ReactomePA, GO.db, AnnotationDbi)

config <- list(
  # Directories
  report_dir = "03_Pathway_Enrichment/02_GSEA/b_reports",
  data_dir   = "03_Pathway_Enrichment/02_GSEA/c_data",
  # Data files
  diff_data  = "02_Differential_Expression/02_Differential/c_data/02_differential.xlsx",
  diff_sheet = "DEP_matrix",
  diff_rds   = "02_Differential_Expression/02_Differential/c_data/fit.rds",
  # ORA hit definition (the old pipeline used Pi.Val, which limpa doesn't output)
  ora_p_col  = "adj.P.Val",   # "P.Value" or "adj.P.Val"
  ora_p_cut  = 0.05
)

dir.create(config$report_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(config$data_dir,   showWarnings = FALSE, recursive = TRUE)

organism <- org.Hs.eg.db

# Analysis label -> limpa contrast name (column prefix in the xlsx / colname in fit$t)
contrast_map <- c(
  Baseline      = "BFR_vs_HLRT_at_T1",
  Training_BFR  = "BFR_post_vs_pre",
  Training_HLRT = "HLRT_post_vs_pre",
  Post_training = "BFR_vs_HLRT_at_T2",
  Interaction   = "interaction"
)


# --- 1: Load limpa results + fit -------------------------------------------
limpa <- read_excel(config$diff_data, sheet = config$diff_sheet)
fit   <- readRDS(config$diff_rds)

stopifnot(
  all(contrast_map %in% colnames(fit$t)),
  all(limpa$protein %in% rownames(fit$t)),
  !anyDuplicated(limpa$protein)
)


# --- 2: Gene symbol -> Entrez / Ensembl mapping ----------------------------
# Tiered lookup to Entrez (org.Hs.eg.db's central key), then Entrez -> Ensembl,
# so the Ensembl (GSEA/GO) and Entrez (KEGG/ORA) IDs always refer to the same gene.
#   1. DIA-NN gene symbol (first symbol if ";"-separated, e.g. "HARS1;HARS")
#   2. UniProt accession (rescues outdated/renamed symbols)
#   3. Symbol as alias (last resort; aliases can be ambiguous)
safe_map <- function(ids, keytype, column) {
  u     <- unique(na.omit(ids))
  valid <- u[u %in% AnnotationDbi::keys(organism, keytype = keytype)]
  out   <- setNames(rep(NA_character_, length(u)), u)
  if (length(valid)) {
    out[valid] <- suppressMessages(
      AnnotationDbi::mapIds(organism, keys = valid, keytype = keytype,
                            column = column, multiVals = "first")
    )
  }
  unname(out[ids])
}

annot <- limpa %>%
  dplyr::select(protein, Genes) %>%
  mutate(
    gene         = sub(";.*$", "", Genes),
    entrez_sym   = safe_map(gene,    "SYMBOL",  "ENTREZID"),
    entrez_uni   = safe_map(protein, "UNIPROT", "ENTREZID"),
    entrez_alias = safe_map(gene,    "ALIAS",   "ENTREZID"),
    Entrez       = coalesce(entrez_sym, entrez_uni, entrez_alias),
    map_source   = case_when(!is.na(entrez_sym)   ~ "SYMBOL",
                             !is.na(entrez_uni)   ~ "UNIPROT",
                             !is.na(entrez_alias) ~ "ALIAS",
                             TRUE                 ~ NA_character_),
    Ensembl      = safe_map(Entrez, "ENTREZID", "ENSEMBL")
  ) %>%
  dplyr::select(protein, Genes, gene, Entrez, Ensembl, map_source)

cat("\nID mapping summary:\n")
print(table(annot$map_source, useNA = "ifany"))
cat(sprintf("Mapped to Ensembl: %d / %d proteins\n",
            sum(!is.na(annot$Ensembl)), nrow(annot)))

write_xlsx(annot, file.path(config$report_dir, "protein_id_mapping.xlsx"))


# --- 3: Per-contrast tables (logFC/P from xlsx, t from fit) -----------------
build_contrast_df <- function(contrast) {
  df <- tibble(
    protein   = limpa$protein,
    logFC     = limpa[[paste0(contrast, "_logFC")]],
    P.Value   = limpa[[paste0(contrast, "_P.Value")]],
    adj.P.Val = limpa[[paste0(contrast, "_adj.P.Val")]],
    t         = unname(fit$t[limpa$protein, contrast])
  ) %>%
    left_join(annot %>% dplyr::select(protein, gene, Entrez, Ensembl),
              by = "protein")
  
  # Guard against row/contrast misalignment between the xlsx and fit
  # (xlsx values are rounded, hence the tolerance)
  stopifnot(isTRUE(all.equal(unname(fit$coefficients[limpa$protein, contrast]),
                             df$logFC, tolerance = 1e-4)))
  df
}

limpa_results <- purrr::map(contrast_map, build_contrast_df)

# Collapse proteins sharing a gene ID (e.g. MOCS2, DUSP13) -> keep max |t|.
# GSEA needs unique names in the ranked list.
dedup_by <- function(df, key) {
  df %>%
    filter(!is.na(.data[[key]]), !is.na(t)) %>%
    group_by(.data[[key]]) %>%
    slice_max(abs(t), n = 1, with_ties = FALSE) %>%
    ungroup()
}


# --- 4: Gene sets (loaded once, reused for every contrast) ------------------
hallmark_t2g <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::distinct(gs_name, ensembl_gene)
hallmark_t2n <- msigdbr(species = "Homo sapiens", collection = "H") %>%
  dplyr::distinct(gs_name, gs_description)

react_t2g <- msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:REACTOME") %>%
  dplyr::distinct(gs_name, ensembl_gene)
react_t2n <- msigdbr(species = "Homo sapiens", collection = "C2",
                     subcollection = "CP:REACTOME") %>%
  dplyr::distinct(gs_name, gs_description)


# --- 5: Pathway / Enrichment functions --------------------------------------
run_gsea_pipeline <- function(df, label) {
  
  cat(sprintf("\n>> Running GSEA pipeline for: %s\n", label))
  set.seed(26)
  
  # --- Ranked lists (limma moderated t) ---
  df_ens    <- dedup_by(df, "Ensembl")
  gene_list <- sort(setNames(df_ens$t, df_ens$Ensembl), decreasing = TRUE)
  
  df_ent           <- dedup_by(df, "Entrez")
  entrez_gene_list <- sort(setNames(df_ent$t, df_ent$Entrez), decreasing = TRUE)
  
  cat(sprintf("   Ranked list: %d Ensembl | %d Entrez\n",
              length(gene_list), length(entrez_gene_list)))
  
  # --- Hallmark GSEA ---
  gsea_hall <- GSEA(
    geneList = gene_list, TERM2GENE = hallmark_t2g, TERM2NAME = hallmark_t2n,
    eps = 0, minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
    pAdjustMethod = "BH", seed = 26, verbose = FALSE
  )
  
  # --- Reactome GSEA ---
  gsea_reactome <- GSEA(
    geneList = gene_list, TERM2GENE = react_t2g, TERM2NAME = react_t2n,
    eps = 0, minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
    pAdjustMethod = "BH", seed = 26, verbose = FALSE
  )
  
  # --- GO GSEA (BP / CC / MF) ---
  run_gse_go <- function(ont) {
    res <- gseGO(
      geneList = gene_list, ont = ont, keyType = "ENSEMBL", OrgDb = organism,
      seed = 26, eps = 0, minGSSize = 15, maxGSSize = 500,
      pvalueCutoff = 1, pAdjustMethod = "BH", verbose = FALSE
    )
    clean        <- res
    clean@result <- res@result %>% filter(!is.na(NES), !is.na(pvalue))
    simple       <- simplify(clean, cutoff = 0.7, by = "p.adjust", select_fun = min)
    list(clean = clean, simple = simple)
  }
  go_bp <- run_gse_go("BP")
  go_cc <- run_gse_go("CC")
  go_mf <- run_gse_go("MF")
  
  # --- KEGG GSEA ---
  kegg <- gseKEGG(
    geneList = entrez_gene_list, organism = "hsa", keyType = "ncbi-geneid",
    eps = 0, minGSSize = 15, maxGSSize = 500, pvalueCutoff = 1,
    seed = 26, pAdjustMethod = "BH"
  )
  
  # --- GO ORA (background = all detected proteins with an Entrez ID) ---
  sig      <- df_ent %>% filter(.data[[config$ora_p_col]] < config$ora_p_cut)
  sig_up   <- sig$Entrez[sig$logFC > 0]
  sig_down <- sig$Entrez[sig$logFC < 0]
  universe <- df_ent$Entrez
  
  run_go_ora <- function(ids, ont_type, direction) {
    if (length(ids) == 0) return(NULL)
    tryCatch(
      enrichGO(gene = ids, universe = universe, OrgDb = organism,
               ont = ont_type, pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n",
                    ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  results <- list(
    label          = label,
    gsea_bp_clean  = go_bp$clean, gsea_bp_simple = go_bp$simple,
    gsea_cc_clean  = go_cc$clean, gsea_cc_simple = go_cc$simple,
    gsea_mf_clean  = go_mf$clean, gsea_mf_simple = go_mf$simple,
    gsea_kegg      = kegg,
    gsea_hall      = gsea_hall,
    gsea_reactome  = gsea_reactome,
    ora_bp_up      = run_go_ora(sig_up,   "BP", "UP"),
    ora_cc_up      = run_go_ora(sig_up,   "CC", "UP"),
    ora_mf_up      = run_go_ora(sig_up,   "MF", "UP"),
    ora_bp_down    = run_go_ora(sig_down, "BP", "DOWN"),
    ora_cc_down    = run_go_ora(sig_down, "CC", "DOWN"),
    ora_mf_down    = run_go_ora(sig_down, "MF", "DOWN")
  )
  
  saveRDS(results, file.path(config$data_dir, sprintf("gsea_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  invisible(results)
}


# ORA only (for predefined protein sets). df needs Entrez + logFC;
# bg_df needs Entrez.
run_ora_pipeline <- function(df, label, bg_df) {
  
  cat(sprintf("\n>> Running ORA pipeline for: %s\n", label))
  
  universe <- unique(na.omit(bg_df$Entrez))
  df       <- df %>% filter(!is.na(Entrez)) %>% distinct(Entrez, .keep_all = TRUE)
  sig_up   <- df$Entrez[df$logFC > 0]
  sig_down <- df$Entrez[df$logFC < 0]
  
  run_go_ora <- function(ids, ont_type, direction) {
    if (length(ids) == 0) return(NULL)
    tryCatch(
      enrichGO(gene = ids, universe = universe, OrgDb = organism,
               ont = ont_type, pAdjustMethod = "BH",
               pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
      error = function(e) {
        cat(sprintf("   GO ORA (%s, %s) failed for %s: %s\n",
                    ont_type, direction, label, e$message))
        NULL
      }
    )
  }
  
  results <- list(
    label       = label,
    ora_bp_up   = run_go_ora(sig_up,   "BP", "UP"),
    ora_cc_up   = run_go_ora(sig_up,   "CC", "UP"),
    ora_mf_up   = run_go_ora(sig_up,   "MF", "UP"),
    ora_bp_down = run_go_ora(sig_down, "BP", "DOWN"),
    ora_cc_down = run_go_ora(sig_down, "CC", "DOWN"),
    ora_mf_down = run_go_ora(sig_down, "MF", "DOWN")
  )
  
  saveRDS(results, file.path(config$data_dir, sprintf("gsea_results_%s.rds", label)))
  cat(sprintf("   Saved .rds for: %s\n", label))
  invisible(results)
}


# --- 6a: Run GSEA for every contrast ----------------------------------------
purrr::iwalk(limpa_results, run_gsea_pipeline)

# --- 6b: CONCORDANT / DIVERGENT ORA (Training_BFR vs Training_HLRT) --------
# CONCORDANT: significant in both contrasts, same direction
# DIVERGENT:  significant in both contrasts, opposite direction
# Background: full detected proteome (same tested universe for all contrasts).
bg_full <- limpa_results$Training_BFR

df_bfr_ora <- limpa_results$Training_BFR %>%
  dplyr::select(protein, gene, Entrez, logFC_BFR = logFC,
                p_BFR = all_of(config$ora_p_col))
df_hlrt_ora <- limpa_results$Training_HLRT %>%
  dplyr::select(protein, logFC_HLRT = logFC,
                p_HLRT = all_of(config$ora_p_col))

df_train_ora <- inner_join(df_bfr_ora, df_hlrt_ora, by = "protein") %>%
  mutate(sig = case_when(
    p_BFR < config$ora_p_cut & p_HLRT < config$ora_p_cut ~ "Both",
    p_BFR < config$ora_p_cut                             ~ "BFR only",
    p_HLRT < config$ora_p_cut                            ~ "HLRT only",
    TRUE                                                 ~ "Neither"
  ))

# logFC = logFC_HLRT is the up/down split key (HLRT is the y-axis in
# fig_05_panel_A_protein_BFR.R's scatter)
concordant_df <- df_train_ora %>%
  filter(sig == "Both", (logFC_BFR > 0) == (logFC_HLRT > 0)) %>%
  mutate(logFC = logFC_HLRT)

divergent_df <- df_train_ora %>%
  filter(sig == "Both", (logFC_BFR > 0) != (logFC_HLRT > 0)) %>%
  mutate(logFC = logFC_HLRT)

cat(sprintf("\nCONCORDANT proteins: %d | DIVERGENT proteins: %d\n",
            nrow(concordant_df), nrow(divergent_df)))

run_ora_pipeline(concordant_df, "CONCORDANT", bg_df = bg_full)
run_ora_pipeline(divergent_df,  "DIVERGENT",  bg_df = bg_full)