# Each subunit and assembly factor of complexes I to V, its training change in each leg, side by
# side. Adapted from Sam's PR #26 (sc-norton). Descriptive only: the inference is fry on the
# complexes in 02_test_mito. A paired t across subunits would count correlated proteins as
# replicates.

pacman::p_load(here, dplyr, tidyr, purrr, stringr, ggplot2, readxl, writexl)

stage <- here("05_Mitochondria", "03_complex_subunits")
de <- read_excel(
  here("02_Differential_Expression", "02_Differential", "c_data", "02_differential.xlsx"),
  "DEP_matrix"
)
membership <- read_excel(
  here("05_Mitochondria", "01_build_mito_sets", "c_data", "01_build_mito_sets.xlsx"),
  "membership"
)
complex_tests <- read_excel(
  here("05_Mitochondria", "02_test_mito", "c_data", "02_test_mito.xlsx"),
  "complex_tests"
)

complexes <- c("CI", "CII", "CIII", "CIV", "CV")
arms <- c(BFR = "BFR_Post-Pre", HLRT = "HLRT_Post-Pre")
subunits <- membership |>
  filter(measured, str_detect(set, "^C(I|II|III|IV|V) (subunits|assembly factors)$")) |>
  separate_wider_delim(set, " ", names = c("complex", "role"), too_many = "merge") |>
  inner_join(de, by = c(gene = "Genes")) |>
  transmute(
    complex = factor(complex, complexes), role, gene, protein,
    BFR_logFC = .data[[paste0(arms["BFR"], "_logFC")]],
    BFR_fdr = .data[[paste0(arms["BFR"], "_adj.P.Val")]],
    HLRT_logFC = .data[[paste0(arms["HLRT"], "_logFC")]],
    HLRT_fdr = .data[[paste0(arms["HLRT"], "_adj.P.Val")]],
    difference = `Modality_x_Time_Interaction_logFC`,
    difference_fdr = `Modality_x_Time_Interaction_adj.P.Val`
  ) |>
  arrange(complex, role, desc(difference))

# Panel titles carry the one valid set-level test: fry on the subunits, standard fit, interaction.
interaction_fry <- complex_tests |>
  filter(fit == "standard", contrast == "Modality_x_Time_Interaction") |>
  transmute(complex = str_remove(set, " subunits"), fry_fdr = FDR)

theme_set(theme_minimal(base_size = 9))
theme_update(
  plot.caption = element_text(hjust = 0, size = 9, lineheight = 1.2),
  plot.caption.position = "plot", legend.position = "top", panel.grid.major.x = element_blank()
)
role_page <- function(which_role, title_text) {
  rows <- filter(subunits, role == which_role)
  headers <- rows |>
    count(complex) |>
    left_join(interaction_fry, by = "complex") |>
    mutate(panel = if_else(
      is.na(fry_fdr), sprintf("%s (%d)", complex, n),
      sprintf("%s (%d), fry interaction FDR %s", complex, n, signif(fry_fdr, 2))
    ))
  rows |>
    pivot_longer(c(BFR_logFC, HLRT_logFC), names_to = "arm", values_to = "logFC") |>
    mutate(
      arm = str_remove(arm, "_logFC"), gene = factor(gene, unique(rows$gene)),
      panel = headers$panel[match(complex, headers$complex)]
    ) |>
    ggplot(aes(gene, logFC, fill = arm)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey50") +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    facet_wrap(~panel, ncol = 1, scales = "free_x") +
    scale_fill_manual(values = c(BFR = "#B2182B", HLRT = "#2166AC"), name = NULL) +
    labs(
      x = NULL, y = "log2 change, post - pre",
      caption = str_wrap(paste(
        title_text, "Bars are each measured protein's training change in the BFR and the HLRT",
        "leg from the stage 02 fit, ordered within complex by BFR minus HLRT. Bars are not",
        "tested against each other: proteins of a complex are correlated, so the set-level test",
        "is fry in 02_test_mito, given in the panel title where it exists. Data: subunits sheet",
        "of 03_complex_subunits.xlsx."
      ), 115)
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6))
}

pdf(file.path(stage, "b_reports", "03_complex_subunits.pdf"), width = 8.27, height = 11.69)
walk(
  list(
    role_page("subunits", "Complex subunits by leg."),
    role_page("assembly factors", "Complex assembly factors by leg.")
  ),
  print
)
invisible(dev.off())

overview <- data.frame(
  sheet = "subunits", rows = nrow(subunits),
  description = paste(
    "Measured subunits and assembly factors of complexes I to V: training logFC and BH FDR per",
    "leg, and the interaction (BFR minus HLRT) from the stage 02 fit."
  )
)
write_xlsx(
  list(overview = overview, subunits = subunits),
  file.path(stage, "c_data", "03_complex_subunits.xlsx")
)
sessionInfo()
