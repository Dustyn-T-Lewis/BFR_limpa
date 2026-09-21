# 03 · Pathway Enrichment

Testing protein by protein asks three thousand separate questions. Pathway analysis asks whether
a group of proteins that work together moved together, which is the right question when the
protein-level result is thin.

Five steps, each named for what it runs, each one job. Data passes through disk, so any step can
re-run alone.

| Step | Runs | Writes |
|---|---|---|
| [`00_build_gene_sets`](00_build_gene_sets/README.md) | freeze MSigDB, map proteins to genes, filter on size | `gene_sets.rds` |
| [`01_run_fgsea`](01_run_fgsea/README.md) | topTable, fgsea, collapsePathways, GO themes | `fgsea.rds` |
| [`02_run_singscore`](02_run_singscore/README.md) | score every sample on every set | `singscore.rds` |
| [`03_enrich_volcano_fgsea`](03_enrich_volcano_fgsea/README.md) | volcano figures with fgsea rings | 6 PNG + 6 PDF |
| [`04_classify_and_associate_themes`](04_classify_and_associate_themes/README.md) | group sets into GO themes, classify them, associate them with the phenotype | one workbook + 2 figures |

```sh
Rscript 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.R
Rscript 03_Pathway_Enrichment/01_run_fgsea/a_script/01_run_fgsea.R
Rscript 03_Pathway_Enrichment/02_run_singscore/a_script/02_run_singscore.R
Rscript 03_Pathway_Enrichment/03_enrich_volcano_fgsea/a_script/03_enrich_volcano_fgsea.R
Rscript 03_Pathway_Enrichment/04_classify_and_associate_themes/a_script/04_classify_and_associate_themes.R
```

About a minute end to end. Settings live in `config.yml` at the repo root, so a threshold is
written once and cannot drift between steps.

## Two questions, two methods

**Did a pathway move?** `fgsea`, in `01_run_fgsea`. Ranks proteins by moderated t and asks whether
a set piles up at one end. Fast, needs no per-sample data, and it is what `enrichVolcano` consumes.

**Where does each sample sit?** `singscore`, in `02_run_singscore`. One number per set per sample,
rank-based and sample-independent. No p-value, never sees the contrast. That 2,732 x 131 matrix is
what the phenotype analysis consumes.

Neither replaces the other. A rotation test (`limma::fry`) that takes the participant design was
tried and removed by an explicit decision to run a single set test.

## What fgsea's assumption costs

`fgsea` assumes the ranked proteins are exchangeable. They are not — mitochondrial and ribosomal
proteins move together, and every participant contributed four samples.

**On `BFR_vs_HLRT_at_T1`, two legs of one person before either was trained, fgsea calls 14 sets
significant and 6 survive collapse**, at adjusted p down to 5e-6, and they are actin, thin-filament
and contractile sets. They look exactly like a training effect. Read every count against that six.

| Contrast | Significant | After collapse | Shared with the control |
|---|---:|---:|---:|
| BFR_post_vs_pre | 466 | 114 | 1 |
| HLRT_post_vs_pre | 375 | 88 | 1 |
| BFR_vs_HLRT_at_T1 *(control)* | 14 | 6 | 6 |
| BFR_vs_HLRT_at_T2 | 202 | 60 | 4 |
| interaction | 39 | 15 | 0 |

The artifact is specific to the between-leg contrasts. The two training contrasts share one set
each with the control out of a hundred-odd, so their results are mostly not this.

## Redundancy, and why no set is excluded by name

Ten collections overlap, so glycolysis is tested in Hallmark, KEGG, Reactome, WikiPathways and
several times over in GO. Nothing is removed before testing except on size: all 2,732 qualifying
sets are tested and BH corrects across all of them.

`collapsePathways` then prunes among the significant hits, re-running each conditioned on a more
significant set's leading edge. It reads the data, not an overlap cutoff and not a name. The `main`
column marks survivors; nothing is deleted, so both counts stay visible.

**No set is excluded by name.** `REACTOME_INFLUENZA_INFECTION` is largely ribosome and translation
machinery under a misleading label. An earlier disease-name regex removed real biology for a label
and has been retired.

## The phenotype result

Every participant contributed one BFR leg and one HLRT leg, so `04_` asks the differential question
within participant rather than between two groups of legs. That roughly doubles the power and
removes every between-participant confound.

**Training separates, treatment does not.** At theme level, pre-to-post classification clears
chance by 4.3x under BFR and 2.9x under high load. The pre-training control sits *below* chance,
as it must. BFR against high load reaches 2.1x post-training and 0.7x on the change score.

No phenotype association clears its chance expectation in either direction. With 32 pairs the
study resolves a paired correlation of about 0.5, so that is a null at this resolution, not
evidence that no difference exists.
