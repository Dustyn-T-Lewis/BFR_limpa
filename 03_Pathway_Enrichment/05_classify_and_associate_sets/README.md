# 05_classify_and_associate_sets

Answers two questions about every tested set: how well it separates the study's groups, and
whether it tracks the phenotype. The unit is the set, with no collapse and no grouping.

| | |
|---|---|
| **Script** | `a_script/05_classify_and_associate_sets.R` |
| **Reads** | `gene_sets.rds`, `set_tests.rds`, `singscore.rds`, `proteins.rds`, `phenotype.csv` |
| **Writes** | one workbook, `set_results.rds`, seven figures as PNG and PDF |

## Each collection carries its own chance expectation

Reading 1,990 sets against one denominator is unreadable. Grouping them into themes forces a
grouping choice onto the primary result. Splitting by collection avoids both: a 41-set collection
is read against 2.1 expected hits, a 1,366-set one against 68.3.

| Task | Hallmark | KEGG | Reactome | GO:BP | GO Slim |
|---|---:|---:|---:|---:|---:|
| **pre vs post, BFR** | 1.4x | 3.8x | 3.4x | 4.0x | **7.3x** |
| **pre vs post, HLRT** | 4.8x | 4.0x | 4.2x | 3.7x | **4.7x** |
| post BFR vs HL | 0.95x | 0.64x | 0.37x | 1.16x | 0.67x |
| change, BFR vs HL | 0.95x | 0.00x | 0.14x | 0.53x | 0.67x |
| baseline *(control)* | 0.95x | 0.43x | 0.47x | 0.48x | **0.00x** |

Training clears chance in all five collections. Every between-leg task sits at or below it, and
the pre-training control sits below it in all five, returning nothing at all from GO Slim. GO Slim
gives the strongest separation because its sets are broad and do not nest, so each is an
independent look rather than another slice of the same branch.

BH within collection and task survives only on the training contrasts: 9 sets under BFR and 49
under high load. It returns nothing on the three between-leg tasks and nothing on any phenotype
association.

The 32,000 Spearman tests behind those tables read `estimate` and `p.value` straight off each
`htest`. Passing them through `broom::tidy` first, which builds one tibble per test, cost 18 of
the stage's 27 seconds and returned the same numbers to the bit.

## Every comparison here is paired

Pre versus post is the same leg twice; BFR versus high load is two legs of one person. A standard
Mann-Whitney AUC treats them as independent and returns a p-value for a design this study did not
run. So **AUC comes from `pROC` as the effect-size descriptor and p from the paired Wilcoxon
signed-rank test.** `pROC::roc()` auto-orients unless `direction` is fixed, which would invert
every below-chance set and flip the red and blue on the figures, so the call pins it.
Correlations use `cor.test`, which computes the exact Spearman p at these sample sizes.

## The figures

A task can reach nominal p in hundreds of sets, so each figure draws the strongest two from each
collection. Reading the collections side by side is the point.

- `roc_<task>.png`, one per task except the baseline control, which is a number to read rather
  than a panel to present.
- `chance_by_database.png`, the table above as bars.
- `association_pooled.png` and `association_differential.png`.

## Handoff

```r
sr <- readRDS("03_Pathway_Enrichment/05_classify_and_associate_sets/c_data/set_results.rds")
sr$chance_expectation # read this first: nominal against chance, per collection
sr$set_auc            # every set x every task, AUC, paired p, BH within collection
sr$set_association    # every set x every outcome, pooled and differential
sr$by_arm             # the same correlation inside BFR and inside high load, descriptive
```
