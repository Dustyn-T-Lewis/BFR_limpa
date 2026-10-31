# Three-way set-test calibration: fry, camera, fgsea on IDENTICAL inputs.
# Purpose: separate METHOD from DATA and GENE SETS. Everything below runs on
# one fit, one catalog, one correction. The negative control is the headline.

library(limma)
library(fgsea)
library(here)
library(dplyr)
library(tibble)
library(tidyr)
library(purrr)
library(readr)
library(ggplot2)

for (fn in c("select", "filter", "rename")) {
  if (!identical(environment(get(fn)), asNamespace("dplyr")))
    stop("dplyr::", fn, "() is masked by ", environmentName(environment(get(fn))))
}

set.seed(1)

out <- here("03_Pathway_Enrichment", "02_Set_Tests", "c_data")
fig <- here("05_Figures", "c_data", "Plots")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

pi_in <- readRDS(here("03_Pathway_Enrichment", "01_Gene_Sets", "c_data",
                      "pathway_inputs.rds"))

E       <- pi_in$fit$EList$E          # 3,042 x 131 log2 abundances
W       <- pi_in$fit$EList$weights    # limpa precision weights
design  <- pi_in$design$design
contr   <- pi_in$design$contrasts
CONTROL <- pi_in$design$negative_control
idx     <- pi_in$set_indices          # row numbers into E, per set
sets    <- pi_in$gene_sets            # gene symbols, same sets, same order
ranked  <- pi_in$ranked_t             # moderated t, named by gene symbol

CONTRASTS <- colnames(contr)

stopifnot(
  identical(dim(E), dim(W)),
  identical(rownames(design), colnames(E)),
  CONTROL %in% CONTRASTS,
  # fry/camera sets and fgsea sets are the same sets, in the same order
  identical(names(idx), names(sets)),
  # every index is valid, unique, and points at the gene the set names
  all(map_lgl(idx, ~ !anyNA(.x) && !anyDuplicated(.x) && max(.x) <= nrow(E))),
  all(map2_lgl(sets, idx, \(g, i) identical(g, pi_in$protein_map$gene[i]))),
  # fgsea's ranked vector covers the same genes the sets are drawn from
  all(map_lgl(sets, ~ all(.x %in% names(ranked[[1]]))))
)

# The handoff must descend from the fit currently on disk
fit_md5 <- unname(tools::md5sum(here("02_Differential_Expression", "02_Differential",
                                     "c_data", "fit.rds")))
stopifnot(identical(
  fit_md5, pi_in$provenance$inputs$md5[pi_in$provenance$inputs$input == "fit"]))

db_of <- function(x) sub("\\|.*$", "", x)          # "Hallmark|HALLMARK_X" -> "Hallmark"

cat(sprintf("%d sets | %d contrasts | %d proteins | control = %s\n",
            length(idx), length(CONTRASTS), nrow(E), CONTROL))
print(table(db_of(names(idx))))

tidy_limma <- function(x, method, cn) {
  as.data.frame(x) |>
    rownames_to_column("set") |>
    transmute(contrast = cn, method = method, set,
              database = db_of(set),
              n = NGenes,
              direction = Direction,
              p = PValue)
}

# BH within contrast x database x method — identical family for every method
add_bh <- function(df) {
  df |> mutate(q = p.adjust(p, "BH"), .by = c(contrast, method, database))
}

RHOS <- c(0.00, 0.01, 0.02, 0.05, 0.10, 0.15)

cam_at <- function(cn, rho) {
  camera(y = E, index = idx, design = design, contrast = contr[, cn],
         weights = W, inter.gene.cor = rho, use.ranks = FALSE, sort = FALSE) |>
    tidy_limma(paste0("camera_rho", rho), cn)
}

sweep <- map(RHOS, \(r) cam_at(CONTROL, r)) |> list_rbind() |> add_bh() |>
  summarise(n_sets = n(),
            raw_p05 = sum(p < 0.05),
            expected = round(0.05 * n()),
            bh_05 = sum(q < 0.05),
            min_p = min(p),
            .by = c(method, database)) |>
  arrange(database, method)
sweep

# smallest rho at which the control rejects nothing in ANY database
clean <- sweep |> summarise(bad = sum(bh_05), .by = method) |> filter(bad == 0)
RHO <- if (nrow(clean)) {
  min(as.numeric(sub("camera_rho", "", clean$method)))
} else {
  warning("No rho gives a clean control. Reporting at 0.02; camera is uncalibrable here.")
  0.02
}
cat("camera inter.gene.cor =", RHO, "\n")

run_fry <- function(cn) {
  set.seed(1)
  fry(y = E, index = idx, design = design, contrast = contr[, cn],
      weights = W, sort = "none") |>
    tidy_limma("fry", cn)
}

run_camera <- function(cn) {
  camera(y = E, index = idx, design = design, contrast = contr[, cn],
         weights = W, inter.gene.cor = RHO, use.ranks = FALSE, sort = FALSE) |>
    tidy_limma("camera", cn)
}

run_fgsea <- function(cn) {
  set.seed(1)
  fgsea(pathways = sets, stats = ranked[[cn]],
        minSize = 1, maxSize = Inf,       # the catalog already bounded size
        eps = 0, scoreType = "std", nproc = 1) |>
    as_tibble() |>
    transmute(contrast = cn, method = "fgsea", set = pathway,
              database = db_of(pathway), n = size,
              direction = if_else(NES > 0, "Up", "Down"),
              p = pval, NES = NES)
}

three_way <- map(CONTRASTS, \(cn) {
  bind_rows(run_fry(cn), run_camera(cn), run_fgsea(cn)) |> mutate(contrast = cn)
}) |> list_rbind() |> add_bh()

stopifnot(
  nrow(three_way) == length(idx) * length(CONTRASTS) * 3,
  !anyNA(three_way$p)
)

calibration <- three_way |>
  filter(contrast == CONTROL) |>
  summarise(n_sets      = n(),
            raw_p_05    = sum(p < 0.05),
            expected_05 = round(0.05 * n()),
            ratio       = round(sum(p < 0.05) / (0.05 * n()), 2),
            bh_05       = sum(q < 0.05),
            min_p       = min(p),
            .by = c(method, database)) |>
  arrange(database, match(method, c("fry", "camera", "fgsea")))
calibration

# and pooled across databases
calibration_overall <- three_way |>
  filter(contrast == CONTROL) |>
  summarise(n_sets = n(), raw_p_05 = sum(p < 0.05),
            expected_05 = round(0.05 * n()),
            ratio = round(sum(p < 0.05) / (0.05 * n()), 2),
            bh_05 = sum(q < 0.05), min_p = min(p),
            .by = method)
calibration_overall

interaction_wide <- three_way |>
  filter(contrast == "interaction") |>
  select(set, database, n, method, p, q) |>
  pivot_wider(names_from = method, values_from = c(p, q),
              names_glue = "{method}_{.value}") |>
  arrange(fgsea_p)

# every set fgsea calls significant, with what the other two say
interaction_wide |>
  filter(fgsea_q < 0.05) |>
  select(set, database, n, fgsea_p, fgsea_q, camera_p, camera_q, fry_p, fry_q) |>
  mutate(across(ends_with(c("_p", "_q")), \(x) signif(x, 3))) |>
  print(n = Inf)

# how many sets does each method call, in each contrast
headline <- three_way |>
  summarise(sig = sum(q < 0.05), min_p = signif(min(p), 3),
            .by = c(contrast, method)) |>
  pivot_wider(names_from = method, values_from = c(sig, min_p)) |>
  arrange(match(contrast, c(CONTROL, "interaction", "BFR_post_vs_pre",
                            "HLRT_post_vs_pre", "BFR_vs_HLRT_at_T2")))
headline

qq <- three_way |>
  filter(contrast == CONTROL) |>
  arrange(p, .by_group = TRUE) |>
  mutate(expected = -log10(ppoints(n())),
         observed = -log10(p),
         .by = method) |>
  mutate(method = factor(method, levels = c("fry", "camera", "fgsea")))

METH_COLS <- c(fry = "#2e8b6b", camera = "#c9a227", fgsea = "#b2182b")

p_qq <- ggplot(qq, aes(expected, observed, colour = method)) +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.5, colour = "grey55") +
  geom_point(size = 0.9, alpha = 0.75) +
  scale_colour_manual(values = METH_COLS, name = NULL) +
  labs(
    title = "Negative control: two untrained legs of the same participant",
    subtitle = sprintf(
      "%s | %d gene sets | BH within database | camera at inter.gene.cor = %.2f",
      CONTROL, length(idx), RHO),
    x = expression(Expected~-log[10](p)),
    y = expression(Observed~-log[10](p))) +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        legend.position = c(0.15, 0.85),
        legend.background = element_rect(fill = "white", colour = "grey80"),
        plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"))

ggsave(file.path(fig, "three_way_calibration_qq.png"), p_qq,
       width = 6.5, height = 5.5, dpi = 300, bg = "white")
ggsave(file.path(fig, "three_way_calibration_qq.pdf"), p_qq,
       width = 6.5, height = 5.5, bg = "white")

write_csv(three_way,            file.path(out, "three_way_all.csv"))
write_csv(calibration,          file.path(out, "three_way_calibration.csv"))
write_csv(calibration_overall,  file.path(out, "three_way_calibration_overall.csv"))
write_csv(interaction_wide,     file.path(out, "three_way_interaction.csv"))
write_csv(headline,             file.path(out, "three_way_headline.csv"))
write_csv(sweep,                file.path(out, "three_way_rho_sweep.csv"))

saveRDS(list(
  results = three_way, calibration = calibration, headline = headline,
  rho_sweep = sweep, rho_used = RHO,
  provenance = list(
    created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    fit_md5 = fit_md5,
    handoff = pi_in$provenance,
    n_sets = length(idx),
    correction = "BH within contrast x database, applied identically to all three methods",
    camera_rho = RHO,
    camera_rho_basis = "smallest value giving zero BH rejections on the negative control",
    fgsea_size_filter = "none applied in fgsea; catalog already bounded at 15-500 measured",
    seed = 1,
    versions = c(limma = as.character(packageVersion("limma")),
                 fgsea = as.character(packageVersion("fgsea"))),
    note = paste(
      "All three tests run on one fit, one catalog, one correction.",
      "Any difference between methods is attributable to the method.")
  )
), file.path(out, "three_way_objects.rds"))

cat("\nwrote:\n"); cat(paste0("  ", list.files(out, "^three_way"), "\n"))
