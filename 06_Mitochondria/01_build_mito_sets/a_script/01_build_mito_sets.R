# MitoCarta 3.0 lists for 02_test_mito: every MitoCarta protein and the 149 MitoPathways. The
# complex readouts are the pathways "CI subunits" to "CV subunits": structural subunits with the
# mtDNA-encoded ones, without the assembly factors that the parent "Complex I" etc. also hold.

pacman::p_load(here, dplyr, tibble, tidyr, readxl, writexl)

out <- here("06_Mitochondria", "01_build_mito_sets", "c_data")
# Downloaded 2026-09-28 from personal.broadinstitute.org/scalvo/MitoCarta3.0/.
mitocarta <- here("06_Mitochondria", "00_Input", "Human.MitoCarta3.0.xls")
protein_map <- read_excel(
  here("03_Pathway_Enrichment", "00_build_gene_sets", "c_data", "00_build_gene_sets.xlsx"),
  "protein_gene_map"
)

# col_types = "text" stops readxl guessing types column by column, which warns on this file.
# The 149 MitoPathways nest, so a parent pathway holds its children's genes.
pathways <- read_excel(mitocarta, "C MitoPathways", col_types = "text") |>
  filter(!is.na(MitoPathway)) |>
  transmute(collection = "MitoPathways", set = MitoPathway, gene = Genes) |>
  separate_longer_delim(gene, ", ")
all_mito <- tibble(
  collection = "MitoCarta", set = "MitoCarta",
  gene = read_excel(mitocarta, "A Human MitoCarta3.0", col_types = "text")$Symbol
)

# Measured means a representative protein carries the symbol, the rows 01_run_fgsea_and_fry tests.
membership <- bind_rows(all_mito, pathways) |>
  mutate(measured = gene %in% protein_map$gene[protein_map$selected])

overview <- data.frame(
  sheet = "membership", rows = nrow(membership),
  description = "One row per set and gene, flagged if measured. 02_test_mito reads this sheet."
)
write_xlsx(
  list(overview = overview, membership = membership),
  file.path(out, "01_build_mito_sets.xlsx")
)
sessionInfo()
