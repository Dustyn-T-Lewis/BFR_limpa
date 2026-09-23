# 03 · Pathway Enrichment

Testing protein by protein asks three thousand separate questions. Pathway analysis asks whether
a group of proteins that work together moved together, which is the right question when the
protein-level result is thin.

Six steps, each named for what it runs, each one job. Data passes through disk, so any step can
re-run alone.

| Step | Runs | Writes |
|---|---|---|
| [`00_build_gene_sets`](00_build_gene_sets/README.md) | freeze MSigDB, map proteins to genes, filter on size, build GO Slim sets | `gene_sets.rds` |
| [`01_run_fgsea_and_fry`](01_run_fgsea_and_fry/README.md) | topTable, fgsea and fry, collapsePathways | `set_tests.rds` + 24 figures |
| [`02_enrich_volcano_fgsea`](02_enrich_volcano_fgsea/README.md) | volcano figures with fgsea rings | six volcanoes |
| [`03_enrich_scatter_fgsea`](03_enrich_scatter_fgsea/README.md) | BFR against high load on one pair of axes | the concordance scatter |
| [`04_run_singscore`](04_run_singscore/README.md) | score every sample on every set | `singscore.rds` + two figures |
| [`05_classify_and_associate_sets`](05_classify_and_associate_sets/README.md) | classify every set, associate it with the phenotype | one workbook + seven figures |

```sh
Rscript 03_Pathway_Enrichment/00_build_gene_sets/a_script/00_build_gene_sets.R
Rscript 03_Pathway_Enrichment/01_run_fgsea_and_fry/a_script/01_run_fgsea_and_fry.R
Rscript 03_Pathway_Enrichment/02_enrich_volcano_fgsea/a_script/02_enrich_volcano_fgsea.R
Rscript 03_Pathway_Enrichment/03_enrich_scatter_fgsea/a_script/03_enrich_scatter_fgsea.R
Rscript 03_Pathway_Enrichment/04_run_singscore/a_script/04_run_singscore.R
Rscript 03_Pathway_Enrichment/05_classify_and_associate_sets/a_script/05_classify_and_associate_sets.R
```

About 80 seconds end to end. Thresholds are written where they are used, with a comment saying
what each one is for.

## Two questions, two methods

**Did a pathway move?** Two tests answer this, in `01_run_fgsea_and_fry`, and both columns ship.
`fgsea` is **competitive**: it ranks proteins by moderated t and asks whether a set piles up at one
end relative to every other protein. `limma::fry` is **self-contained**: it asks whether a set moved
at all given the fitted design, rotating the model residuals rather than shuffling gene labels, so
correlation cannot inflate the null and the participant structure is built in.

**Where does each sample sit?** `singscore`, in `04_run_singscore`. One number per set per sample,
rank-based and sample-independent. No p-value, never sees the contrast. That 1,990 x 131 matrix is
what the classification and phenotype analysis consumes.

Neither replaces the other.

## What the two tests say, and where they differ

`fgsea` assumes the ranked proteins are exchangeable. They are not — mitochondrial and ribosomal
proteins move together, and every participant contributed four samples. `fry` makes no such
assumption.

| Contrast | fgsea sig | after collapse | **fry sig** |
|---|---:|---:|---:|
| BFR_Post-Pre | 376 | 88 | 656 |
| HLRT_Post-Pre | 281 | 68 | 589 |
| BFR_Post-HLRT_Post | 159 | 52 | **0** |
| Modality_x_Time_Interaction | 26 | 10 | **0** |
| **BFR_Pre-HLRT_Pre** *(control)* | 1 | 1 | **0** |

**Neither is simply better.** fry flags 656 of 1,990 sets on BFR training, a third, because with a
strong global effect most sets genuinely did move a little. It is liberal where fgsea is
conservative. Where they disagree is what matters: fry returns nothing on all three between-leg
contrasts, including the interaction's 10 fgsea pathways, so that result is a direction to follow
up rather than a finding.

The calibration difference is not specific to this implementation. Compared on 996 shared sets
across three independent implementations of fry and four of fgsea, on the same negative control,
fry returned 0, 0, 0 significant at a p < 0.05 rate of 0.016 to 0.019, and fgsea returned 3, 3, 1
and 14 at a rate of 0.057 to 0.068. Under a true null that rate should be 0.05.

**The control is down from 14 significant to 1.** Six of the ten MSigDB collections were dropped
after measuring what each contributed, GO:CC among them. GO:CC is the ontology of physical
complexes, whose members are co-regulated by stoichiometry, which is the between-gene correlation
a competitive test assumes away; three of the control's original six collapse survivors were
GO:CC sets. `00_build_gene_sets/README.md` records the accounting.

One set remains on the control, `REACTOME_STRIATED_MUSCLE_CONTRACTION`. It is contractile and
looks like a training effect, and it is the residue of the same artifact. Read every count against
that one.

## Redundancy, and why no set is excluded by name

Five collections overlap, so glycolysis is tested in Hallmark, KEGG and Reactome and several times
over in GO:BP. Nothing is removed before testing except on size: all 1,990 qualifying sets are
tested and BH corrects across all of them, one family per contrast.

The order is deliberate, and it was measured rather than assumed. Deduplicating by Jaccard overlap
before testing shrank the correction by 1.75x, to 1,102 sets from the 1,931 tested at the time,
and lost two thirds of the discoveries: 363 down to 133 on BFR training. Bourgon, Gentleman and
Huber (PNAS 2010) is the reason: a pre-test filter only buys
power when it is independent of the test statistic under the null, and it works by removing
*hopeless* tests. A redundancy filter removes *redundant* ones, which are frequently the
significant ones. `collapsePathways` also needs the results to run at all, so it can only ever
come after.

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

**Training separates, treatment does not, and it replicates across five curations.** Every
collection is read against its own chance expectation rather than grouped, so a 41-set collection
is not judged by a 1,366-set one's denominator.

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| **pre vs post, BFR** | 1.4x | 3.8x | 3.4x | 4.0x | **7.3x** |
| **pre vs post, HLRT** | 4.8x | 4.0x | 4.2x | 3.7x | **4.7x** |
| post BFR vs HL | 0.95x | 0.64x | 0.37x | 1.16x | 0.67x |
| change, BFR vs HL | 0.95x | 0.00x | 0.14x | 0.53x | 0.67x |
| baseline *(control)* | 0.95x | 0.43x | 0.47x | 0.48x | **0.00x** |

Under BH within collection and task, 9 sets survive on BFR training and 49 on high load. Nothing
survives on any between-leg task, and no phenotype association survives at set level or at protein
level (see `02_Differential_Expression/03_Phenotype`). With 32 pairs the study resolves a paired
correlation of about 0.5, so that is a null at this resolution, not evidence that no difference
exists.

The scatter in `03_enrich_scatter_fgsea` says why the treatment comparison is empty. On 100
non-nesting sets the two arms' NES correlate at rho 0.82, and all 25 significant sets agree in
sign: both modalities moved the same biology by the same amount.
