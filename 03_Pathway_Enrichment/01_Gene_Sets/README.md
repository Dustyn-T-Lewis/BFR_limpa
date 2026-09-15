# 03_Pathway_Enrichment / 01_Gene_Sets

Freezes human gene sets, prepares the protein matrix and fitted-model handoff,
and draws protein volcanoes with enrichVolcano. It runs no pathway tests.

| | |
|---|---|
| **Script** | `a_script/01_gene_sets.qmd` |
| **Reads** | Quantification `proteins.rds`, differential `fit.rds`, design `design.rds` |
| **Report** | `b_reports/01_gene_sets.html` |
| **Figures** | `b_reports/figures/`: five FDR volcanoes and two pi-ranked training views, PNG and PDF |
| **Data** | `c_data/01_gene_sets.xlsx`, `protein_matrix.rds`, `gene_sets.rds`, `pathway_inputs.rds`, CSV tables |
| **Frozen source** | `c_data/cache/msigdb_2026.1.Hs_<key>.rds` and matching `.md5` |

From the project root:

```sh
quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd --output-dir ../b_reports
```

## Cache and gene sets

The default is human MSigDB **2026.1.Hs**, Hallmark, Reactome and GO Biological
Process. `msigdbr` supplies the snapshot on the first run, using its existing
machine cache when present. The project snapshot is frozen before any
data-dependent filtering. Later renders verify its checksum and need no
network access or installed `msigdbr`.

The cache key includes species, database release, collection selection, and
schema. A requested release that differs from what `msigdbr` supplies stops
the render rather than silently substituting a newer release. Restore both
cache files together when a checksum fails. To intentionally rebuild a
snapshot, move its RDS and checksum aside first. To adopt a new release,
change `msigdb_version`; the previous snapshot remains on disk.

Both size bars and the overlap cutoff are Quarto parameters. By default a set
needs 15 to 500 source genes and at least 15 measured representatives, and
sets overlapping an already retained set at Jaccard 0.5 or more are dropped
within each database. Ties go to the larger measured set, then the larger
source set, then the set ID. The catalog records every size exclusion, every
redundancy exclusion, and which set displaced each casualty. `gene_sets.rds`
keeps all three collections: full, eligible and retained. This reduction is
structural and fixed in advance. It filters nothing on significance and makes
no claim that the survivors are independent.

For example, to change the measured-membership floor:

```sh
quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd --output-dir ../b_reports -P min_measured:20
```

`collections` also accepts `GOCC` and `GOMF`; specify an R/Quarto parameter
list when changing multiple collections. Derived outputs are rebuilt on
every render; changing the matrix or parameters cannot reuse stale overlap
filters or stale scores. Input checksums, settings, package versions and the
source-cache checksum are embedded in the RDS provenance.

## Matrix and identifiers

The protein matrix contains **every upstream protein**, unchanged and in the
same sample order. The workbook contains it as the `protein_matrix` sheet;
`protein_matrix.rds` is the numeric matrix. The complete EList, standard errors,
observation counts, sample metadata, weighted differential fit, and participant
design remain in `pathway_inputs.rds`. The stage checks that the saved fit
actually contains the same matrix, uncertainty, annotations, and design.

Gene-set inputs match on a single symbol, exactly. A protein annotated with no
symbol, or with several, is left out rather than split into separate genes, and
no alias is inferred. Where several proteins share a symbol, the one with the
highest mean observed precursor count represents it, with `NPrec` and then
accession breaking ties. That choice holds across all five contrasts and reads
no p-value and no fold change. The workbook's `protein_gene_map` sheet records
every decision. A representative no database annotates still counts as measured,
and one annotated background per collection is exported beside it.

## FDR and pi-score volcanoes

[enrichVolcano](https://github.com/Dustyn-T-Lewis/enrichVolcano) is a plotting
package. This stage uses `volcano_ring()` with an empty enrichment table:
the protein volcano is drawn, and no untested pathway is presented as enriched.

Version 0.3.0.9000 rejects empty enrichment tables in its validator even though
its renderer handles an empty ring. A small adapter in the notebook allows
that empty-table case in a private function environment. It uses the package's
actual drawing function, preserves validation for nonempty enrichment tables,
does not modify the installed namespace, and introduces no dummy pathways.
The package version and renderer-body checksum are recorded in provenance.

All five contrasts use raw-p height, log2-fold-change x, and **BH FDR < 0.05**
for colours and counts. No fold-change floor is added. The two training
contrasts also have a view whose labels are ranked by pi-value. Those views
keep the same FDR colours and counts, and a pi-ranked label can land on a
protein that never passed FDR. The package rescales coordinates for each panel,
so read exact values from the tables rather than comparing cloud size between
contrasts.

- `padj` / `adj.P.Val`: BH adjustment across all tested proteins within each
  contrast, preserved after gene mapping.
- `pi_score = P.Value^abs(logFC)`: the upstream inverted score; smaller ranks higher.
- `pi_value = abs(logFC) * -log10(P.Value)`: equivalent published scale;
  larger ranks higher. Zero p-values are bounded only for the logarithm.
- `signed_pi = sign(logFC) * pi_value`: a directional ranking, without error control.

Pi selects no hit list, for discovery or for ORA. Every protein score is
exported, while the pi-ranked plots, top tables and signed rank vectors cover
only the two training contrasts. Detection counts travel with the scores. The
negative control and the interaction stay visible in the FDR plots and in the
full exports.

## Next-stage contract

```r
x <- readRDS("03_Pathway_Enrichment/01_Gene_Sets/c_data/pathway_inputs.rds")
x$proteins        # Original EList: E, standard errors, observations, targets
x$fit             # Existing limpa/limma fit, including precision/sample weights
x$design          # Participant design, contrasts, negative-control name
x$protein_map     # Gene representative chosen per protein; what set_indices indexes into
x$gene_matrix     # One fixed representative per gene; normalised log2 values
x$protein_results # All tested proteins, all contrasts, FDR and pi scores
x$gene_results    # Representative rows, retaining original protein-level FDR
x$gene_sets       # Retained sets of gene symbols
x$set_indices     # Matching row indices into x$proteins$E (not x$gene_matrix)
x$ranked_t        # Named moderated-t vectors for every contrast
x$ranked_pi       # Named signed-pi vectors for the two training contrasts
x$annotated_universes # Measured background per collection
x$provenance      # Input checksums, settings, package versions, source-cache checksum
```

A whole-set test reading this handoff has to keep the participant design and
limpa's uncertainty. The rank vectors are here because some tools want them,
not as a licence to run a method that assumes independent proteins. Report the
baseline control before anyone reads the interaction. Per-sample pathway
scoring and ORA belong to the later stages. The approximation documented
upstream, where cyclic loess moves the abundances and leaves the standard
errors behind, carries forward into all of this.

Dependencies are checked, never auto-installed. Use an enrichVolcano version
with `volcano_ring()` and its `volc_sig_col` argument (verified here with
0.3.0.9000). Other packages: `here`, `limma`, `limpa`, `dplyr`, `tibble`,
`tidyr`, `purrr`, `readr`, `writexl`, `digest`, `Matrix`, `ggplot2`, `knitr`;
`msigdbr` only to create a missing source snapshot.
