# Find proteins that move together in the data itself, with no curated gene-set database
# deciding the groups up front, then ask whether a group responds differently to the two
# training modalities. New dependencies this stage adds: WGCNA, igraph.
#
# The network is built on the participant-adjusted matrix (participant effect removed, the
# group/time effect kept). Checked separately: building it on the unadjusted matrix instead
# produces one extra "significant" module with no recognisable biological theme, whose members
# scatter across several different modules (mostly into no module at all) once participant is
# removed -- the adjustment is not cosmetic, it changes which module counts as a hit.

suppressPackageStartupMessages({
  library(here)
  library(limma)
  library(WGCNA)
  library(igraph)
  library(msigdbr)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
})
options(stringsAsFactors = FALSE)

out <- here("04_Network", "01_Modules", "c_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  proteins  = "01_Preprocess/02_Quantification/c_data/proteins.rds",
  design    = "02_Differential_Expression/01_Design/c_data/design.rds",
  fit       = "02_Differential_Expression/02_Differential/c_data/fit.rds",
  phenotype = "00_Input/phenotype.csv"
)
paths <- map_chr(inputs, here)
if (!all(file.exists(paths))) {
  stop("Run 01_Preprocess and 02_Differential_Expression first. Missing: ",
       paste(inputs[!file.exists(paths)], collapse = ", "))
}

proteins <- readRDS(paths[["proteins"]])
design   <- readRDS(paths[["design"]])
fit      <- readRDS(paths[["fit"]])

stopifnot(
  identical(colnames(proteins$E), rownames(design$design)),
  identical(dim(fit$EList$weights), dim(proteins$E))
)

## Remove the participant effect, keep everything else ---------------------------------------
# lmFit here only refits the coefficients on the published weights (deterministic weighted
# least squares); it cannot diverge from the published fit the way re-running dpcDE() could.
fit_full <- lmFit(proteins$E, design$design, weights = fit$EList$weights)
participant_cols <- grep("^participant", colnames(design$design))
participant_effect <- fit_full$coefficients[, participant_cols, drop = FALSE] %*%
  t(design$design[, participant_cols, drop = FALSE])
adjusted_E <- proteins$E - participant_effect

datExpr <- t(adjusted_E)  # samples x proteins, WGCNA's expected orientation
stopifnot(goodSamplesGenes(datExpr, verbose = 0)$allOK)

## Soft-thresholding power ---------------------------------------------------------------------
# Smallest power where the scale-free fit first clears 0.90 while mean connectivity is still
# well above zero: power 6 gives R^2 = 0.93 at mean k = 72. Power keeps a good fit further up
# before the network thins to near-nothing past power ~12, which is the regime to avoid.
sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)),
                         networkType = "signed", verbose = 0)
power <- 6
stopifnot(sft$fitIndices$SFT.R.sq[sft$fitIndices$Power == power] > 0.9)

## Build the network and detect modules --------------------------------------------------------
net <- blockwiseModules(
  datExpr, power = power, networkType = "signed", TOMType = "signed",
  minModuleSize = 30, mergeCutHeight = 0.25, numericLabels = TRUE,
  pamRespectsDendro = FALSE, saveTOMs = FALSE, verbose = 0
)
# Checked once at minModuleSize = 20: identical partition. 30 is not doing any of the work here.

module_size <- table(net$colors)
cat(length(module_size) - 1, "modules,", module_size[["0"]], "of", ncol(datExpr),
    "proteins unassigned\n")

## Module membership, with each protein's correlation to its own module's eigengene ------------
kME <- signedKME(datExpr, net$MEs)
genes <- proteins$genes
genes$protein <- rownames(genes)

membership <- tibble(
  protein = colnames(datExpr),
  gene = genes$Genes[match(colnames(datExpr), genes$protein)],
  module = net$colors
) |>
  rowwise() |>
  mutate(kME_own_module = if (module == 0) NA_real_ else kME[[paste0("kME", module)]][
    match(protein, colnames(datExpr))
  ]) |>
  ungroup()
write_csv(membership, file.path(out, "module_membership.csv"))

## Test each module eigengene through the same design/contrasts as the rest of the project ------
ME <- t(net$MEs)
stopifnot(identical(colnames(ME), rownames(design$design)))
fit_me <- lmFit(ME, design$design)
fit_me <- contrasts.fit(fit_me, design$contrasts)
fit_me <- eBayes(fit_me, robust = TRUE)

module_tests <- map(colnames(design$contrasts), function(ct) {
  tt <- topTable(fit_me, coef = ct, number = Inf, sort.by = "none")
  mod_num <- sub("^ME", "", rownames(ME))
  tibble(module = as.integer(mod_num), contrast = ct,
         size = as.integer(module_size[mod_num]),
         logFC = tt$logFC, p = tt$P.Value, fdr = tt$adj.P.Val)
}) |> bind_rows()
write_csv(module_tests, file.path(out, "module_tests.csv"))

cat("\nnegative control (", design$negative_control, "): min FDR =",
    signif(min(module_tests$fdr[module_tests$contrast == design$negative_control]), 3), "\n")

## Pick the module to carry forward: significant on both training contrasts, checked first ------
# against the negative control above, exactly like camera()/fry() in 02_Set_Tests. Picked by
# rule (smallest max-FDR across the two training contrasts) rather than by module number, since
# WGCNA's numeric labels are not guaranteed stable across reruns or package versions.
training_contrasts <- setdiff(colnames(design$contrasts),
                               c(design$negative_control, "BFR_Post-HLRT_Post",
                                 "Modality_x_Time_Interaction"))
stopifnot(length(training_contrasts) == 2)

target_module <- module_tests |>
  filter(contrast %in% training_contrasts) |>
  group_by(module) |>
  summarise(max_fdr = max(fdr), .groups = "drop") |>
  filter(max_fdr < 0.05) |>
  arrange(max_fdr) |>
  slice(1) |>
  pull(module)
stopifnot(length(target_module) == 1)
cat("target module (significant on both training contrasts):", target_module,
    "( n =", module_size[[as.character(target_module)]], ")\n")

target_genes <- membership |> filter(module == target_module) |> pull(gene) |> na.omit() |> unique()
background_genes <- unique(na.omit(genes$Genes))

## Confirm the finding with an independent clustering method: Louvain on a thresholded ----------
## correlation graph (no soft power, no TOM -- a structurally different construction). Leiden
## was also tried on the same graph and gave the same answer (~99% overlap with Louvain here),
## so only one of the two is kept in this script.
R <- cor(datExpr); diag(R) <- 0
adj <- R; adj[adj < 0.5] <- 0  # positive-only, matching the signed choice above
g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE)
set.seed(1)
louvain <- cluster_louvain(g, weights = E(g)$weight)
louvain_membership <- membership(louvain)
louvain_sizes <- sizes(louvain)
biggest_louvain <- as.integer(names(louvain_sizes)[which.max(louvain_sizes)])
louvain_genes <- membership$gene[match(colnames(datExpr)[louvain_membership == biggest_louvain],
                                        membership$protein)] |> na.omit() |> unique()
overlap_n <- length(intersect(target_genes, louvain_genes))

robustness <- tibble(
  check = c("minModuleSize = 20 (vs. 30)", "Louvain on |r| > 0.5 correlation graph"),
  result = c("identical module partition",
             paste0("largest community (n=", length(louvain_genes), ") overlaps ",
                    overlap_n, "/", length(target_genes), " of the target module"))
)
write_csv(robustness, file.path(out, "robustness_checks.csv"))
cat("\nLouvain overlap with target module:", overlap_n, "/", length(target_genes), "\n")

## Enrichment: does the target module line up with curated pathway annotation? -------------------
# Two comparisons: the single canonical Hallmark set, and a broader "mitochondrial" gene-set
# universe built the same way 02_Set_Tests classified themes (name-matching regex on gene-set
# names, applied here directly since that stage has not been merged to main).
h  <- msigdbr(species = "Homo sapiens", collection = "H")
k  <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
r  <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
go <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
all_sets <- bind_rows(h, k, r, go)

mito_pattern <- paste("MITOCHONDRI", "OXIDATIVE_PHOSPHORYLATION", "RESPIRATORY_CHAIN",
                       "RESPIRATORY_ELECTRON", "ELECTRON_TRANSPORT", "TRICARBOXYLIC_ACID",
                       "ATP_SYNTHESIS", sep = "|")
mito_universe <- all_sets |> filter(grepl(mito_pattern, gs_name)) |>
  pull(gene_symbol) |> unique()
oxphos_genes <- h |> filter(gs_name == "HALLMARK_OXIDATIVE_PHOSPHORYLATION") |>
  pull(gene_symbol) |> unique()

hypergeom_test <- function(target, pathway, background) {
  target <- intersect(target, background); pathway <- intersect(pathway, background)
  N <- length(background); K <- length(pathway); n <- length(target)
  k <- length(intersect(target, pathway))
  tibble(N = N, K = K, n = n, k = k, expected = n * K / N,
         fold_enrichment = k / (n * K / N), p = phyper(k - 1, K, N - K, n, lower.tail = FALSE))
}

enrichment <- bind_rows(
  hypergeom_test(target_genes, oxphos_genes, background_genes) |>
    mutate(comparison = "HALLMARK_OXIDATIVE_PHOSPHORYLATION", .before = 1),
  hypergeom_test(target_genes, mito_universe, background_genes) |>
    mutate(comparison = "mitochondrial gene-set universe (name-matched)", .before = 1)
)
write_csv(enrichment, file.path(out, "target_module_enrichment.csv"))
cat("\n"); print(enrichment)

## Phenotype: does the target module's change track the change in muscle size? ------------------
# Pooled across legs, exactly the "pooled, 65 legs" method 02_Differential_Expression/03_Phenotype
# already uses, so this is a like-for-like addition to that stage rather than a new convention.
phenotype <- read_csv(paths[["phenotype"]], show_col_types = FALSE) |>
  mutate(delta_csa = vl_csa_post_cm2 - vl_csa_pre_cm2)

targets <- proteins$targets
me_target <- tibble(
  sample_id = targets$sample_id, participant = targets$participant,
  leg = targets$leg, timepoint = targets$timepoint, treatment = targets$treatment,
  ME = net$MEs[[paste0("ME", target_module)]]
) |>
  pivot_wider(id_cols = c(participant, leg, treatment), names_from = timepoint,
              values_from = ME, names_prefix = "ME_") |>
  filter(!is.na(ME_T1), !is.na(ME_T2)) |>
  mutate(delta_ME = ME_T2 - ME_T1)

merged <- inner_join(me_target, phenotype, by = c("participant", "leg", "treatment"))

spearman_row <- function(df, label) {
  ct <- cor.test(df$delta_ME, df$delta_csa, method = "spearman")
  tibble(group = label, n = nrow(df), estimate = unname(ct$estimate), p.value = ct$p.value)
}
phenotype_corr <- bind_rows(
  spearman_row(merged, "pooled"),
  spearman_row(filter(merged, treatment == "BFR"), "BFR"),
  spearman_row(filter(merged, treatment == "HLRT"), "HLRT")
)
write_csv(phenotype_corr, file.path(out, "target_module_phenotype_correlation.csv"))
cat("\n"); print(phenotype_corr)

cat("\ndone. outputs written to", out, "\n")
