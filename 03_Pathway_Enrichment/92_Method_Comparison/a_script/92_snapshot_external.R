# Copies the set-test results of the other pathway branches into c_data/external/, pinned to a
# commit, so 92_method_comparison.qmd can compare against them without checking those branches
# out. Run from the repository root, once, whenever a source branch moves:
#
#   git fetch origin pull/11/head pull/13/head pull/15/head
#   Rscript 03_Pathway_Enrichment/92_Method_Comparison/a_script/92_snapshot_external.R
#
# Each snapshot is the file exactly as committed, minus the long text columns (leadingEdge,
# description) that nothing here reads. external_manifest.csv records the commit, the path and
# the md5 of the original file, so anyone can check a snapshot against the branch it came from
# with `git show <commit>:<path>`.

library(readr)
library(dplyr)

root <- here::here()
out <- file.path(root, "03_Pathway_Enrichment", "92_Method_Comparison", "c_data", "external")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

sources <- tribble(
  ~source,              ~owner,    ~pull_request, ~commit,    ~path,
  "dustyn_pr11_fgsea",  "Dustyn",  11L,           "b16fc4f",  "03_Pathway_Enrichment/01_Gene_Sets/c_data/fgsea_results.csv",
  "dustyn_pr11_fry",    "Dustyn",  11L,           "b16fc4f",  "03_Pathway_Enrichment/02_Set_Tests/c_data/fry_results.csv",
  "dustyn_pr15_fgsea",  "Dustyn",  15L,           "aa66aa0",  "03_Pathway_Enrichment/01_run_fgsea/c_data/fgsea_results.csv",
  "matheus_pr13",       "Matheus", 13L,           "62a2404",  "03_Pathway_Enrichment/02_Set_Tests/c_data/geneset_tests_long.csv"
)

manifest <- lapply(seq_len(nrow(sources)), function(i) {
  s <- sources[i, ]
  raw <- system2("git", c("-C", shQuote(root), "show", paste0(s$commit, ":", s$path)),
                 stdout = TRUE)
  if (!is.null(attr(raw, "status"))) {
    stop("Could not read ", s$commit, ":", s$path, ". Run the git fetch above first.")
  }
  tmp <- tempfile(fileext = ".csv")
  writeLines(raw, tmp, useBytes = TRUE)
  table <- read_csv(tmp, show_col_types = FALSE, guess_max = 1e5) |>
    select(-any_of(c("leadingEdge", "description")))
  write_csv(table, file.path(out, paste0(s$source, ".csv.gz")))
  mutate(s, full_commit = system2("git", c("-C", shQuote(root), "rev-parse", s$commit),
                                  stdout = TRUE),
         original_md5 = unname(tools::md5sum(tmp)), rows = nrow(table),
         snapshot_utc = format(Sys.time(), tz = "UTC", usetz = TRUE))
}) |>
  bind_rows()

write_csv(manifest, file.path(out, "external_manifest.csv"))
print(manifest)
