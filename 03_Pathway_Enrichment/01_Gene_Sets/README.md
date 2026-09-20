# 03_Pathway_Enrichment / 01_Gene_Sets

Freezes human gene sets, maps proteins onto genes, and runs fgsea on each database.

| | |
|---|---|
| **Script** | `a_script/01_gene_sets.qmd` |
| **Reads** | `proteins.rds`, `fit.rds`, `design.rds` |
| **Writes** | `c_data/gene_sets.rds`, `pathway_inputs.rds`, `01_gene_sets.xlsx`, four CSVs |
| **Frozen source** | `c_data/cache/msigdb_2026.1.Hs_Hallmark-Reactome-GOBP.rds` and its `.md5` |

```sh
quarto render 03_Pathway_Enrichment/01_Gene_Sets/a_script/01_gene_sets.qmd --output-dir ../b_reports
```

## The freeze

Set membership changes between MSigDB releases, so one release is pinned and kept. The first
render fetches human MSigDB 2026.1.Hs through `msigdbr` and writes a snapshot with an md5 beside
it; every render after that reads the snapshot and checks the md5, so it needs neither the network
nor `msigdbr` installed. The filename carries the release and the collections, which are what make
one snapshot different from another.

Do not rebuild it casually. Restore the RDS and its checksum together when a check fails, and move
both aside first if a rebuild is genuinely wanted.

## Which sets are kept

A set needs 15 to 500 source genes and at least 15 of them measured here. The second bar counts
proteins this experiment actually detected, which is what governs power.

Many surviving sets then say almost the same thing. Within each database the stage walks the
qualifying sets from largest measured membership down and drops any whose Jaccard overlap with a
set already kept reaches 0.5. Databases reduce separately, so a GO term cannot erase a Hallmark set
covering the same biology. Every dropped set records which set displaced it and by how much.

1,837 sets qualify and 1,040 survive reduction. That cut is structural and fixed in advance: it
reads no p-value, and it does not make what survives independent.

## One protein per gene

Set testing needs one row per gene and the protein matrix does not supply that cleanly. A protein
with no gene symbol, or with several, is left out rather than split, since splitting would invent
measurements nobody made. Where several proteins share a symbol, the one with the most observed
precursors per sample represents it. That choice is made once, for all five contrasts, and reads no
fold change and no p-value.

Every protein stays in the matrix and in the volcano data regardless. Only the gene-level view
drops them, and the `protein_gene_map` sheet records every decision.

## fgsea, and why its control matters

fgsea walks the proteins ordered by moderated t and asks whether a set clusters at one end. It is
fast, needs no per-sample data, and is what `enrichVolcano` consumes, which is why it runs here.
It also assumes the ranked proteins are exchangeable, and they are not.

That costs something measurable. On `BFR_vs_HLRT_at_T1` — two legs of one person before either
was trained, so no signal is possible — fgsea calls one set significant on the combined
collection and two on Reactome. It calls 59 on `BFR_vs_HLRT_at_T2`, where `02_Set_Tests` finds
none. The negative-control count therefore prints before any other result in the report, and
every other count should be read against it rather than against zero.

`02_Set_Tests` runs `fry()`, which takes the design and the participant structure and returns
nothing on both null contrasts. That is the test to quote.

## Excluding sets by name

`exclude_disease: true` drops sets whose names mention a disease, infection or tumour. This is a
judgement, not a statistic: `REACTOME_INFLUENZA_INFECTION` is largely ribosome and translation
machinery under a misleading label, so the filter removes real biology along with the label.

It is therefore built to be reversible and visible. The unfiltered table is exported beside the
filtered one, `excluded_terms` lists every removal with the pattern that caught it, and setting
the parameter to `false` reproduces the unfiltered result exactly.

## Next-stage contract

`pathway_inputs.rds` holds what this stage derives, and nothing it merely read:

```r
x <- readRDS("03_Pathway_Enrichment/01_Gene_Sets/c_data/pathway_inputs.rds")
x$protein_results # All tested proteins, all contrasts, FDR and pi scores
x$protein_map     # Gene representative chosen per protein
x$gene_sets       # Retained sets of gene symbols
x$set_indices     # Row indices into proteins$E, not into anything built here
x$set_lists       # Per database and pooled, as load_sets() returned them
x$fgsea_results   # All collections, all contrasts, leadingEdge as a list column
x$provenance      # Input checksums, settings, package versions, cache checksum
```

The matrix and the fitted model stay where the upstream stages wrote them. A set test opens
`01_Preprocess/02_Quantification/c_data/proteins.rds` and
`02_Differential_Expression/02_Differential/c_data/fit.rds` itself, the way every stage in this
repo reads its upstream. Copying them here would only let the copies go stale.

`fgsea_results` is keyed on gene symbols, which is what lets `03_Volcanoes` match `leadingEdge`
against its point labels. Rank on anything else and the tick lines vanish without a warning.

Whatever runs next has to keep the participant design and limpa's uncertainty, and it should report
the baseline control before anyone reads the interaction. The upstream approximation carries
forward: cyclic loess moved the abundances and left the standard errors behind.

`msigdbr` is needed only to create a missing snapshot.
