# 03_Pathway_Enrichment / 02_Set_Tests

Runs three families of set test on the negative control, keeps the ones that stay calibrated, and
tests the remaining contrasts with those. Per-sample set scoring is not here; it is a separate
sub-stage and depends on this one choosing a method.

| | |
|---|---|
| **Script** | `a_script/02_set_tests.qmd` |
| **Reads** | `03_Pathway_Enrichment/01_Gene_Sets/c_data/pathway_inputs.rds` |
| **Writes** | `c_data/set_test_results.csv`, `c_data/negative_control_summary.csv`, `c_data/set_tests.rds` |

```sh
quarto render 03_Pathway_Enrichment/02_Set_Tests/a_script/02_set_tests.qmd --output-dir ../b_reports
```

## The control decides the method

`BFR_vs_HLRT_at_T1` compares the two legs of one person before either was trained, with treatment
counterbalanced. The labels are exchangeable and the truth is that nothing separates them. A
method that returns pathways here has told us what its false-positive rate looks like on this
data, and its results on the interaction cannot be read.

So the control runs first, and a method above `control_tolerance` significant sets is dropped
before the real contrasts are computed. Its results on those contrasts are never printed rather
than printed with a caveat, because a caveat beside a table is read as a table.

The p-value histogram is drawn beside the counts. A method well below a flat null is conservative
rather than wrong, which is not disqualifying, but it bounds the claim: an interaction returning
nothing under an under-calling test says less than an interaction returning nothing under a
calibrated one. `02_Differential` already documented exactly that pattern at the protein level for
all three between-leg contrasts.

## Why three

| Family | Function | Assumes | Sees the design |
|---|---|---|---|
| Self-contained, rotation | `fry` | nothing about the other sets | yes, with weights |
| Competitive, correlation-aware | `camera` | set against the rest, inter-gene correlation estimated | yes, with weights |
| Preranked, competitive | `fgsea` | genes are exchangeable and independent | no, a vector only |

`fgsea` is a comparator, not a candidate. It receives a ranked vector and has no way to know that
participants contributed four samples each. The stage README names an independence-assuming test
that reported significant sets on the negative control; running it here is how that description
gets checked against this data instead of being taken on trust.

## Precision weights have to be recovered

`dpcDE()` calls `voomaLmFitWithImputation`, which turns the DPC-Quant standard errors into vooma
precision weights and fits with them. Those weights are the only channel by which quantification
uncertainty reaches any test downstream. `contrasts.fit()` drops them, so the fit saved by
`02_Differential` no longer carries them, and they cannot be rebuilt as `1 / standard.error^2`
because vooma combines the standard errors with a fitted mean-variance trend.

This stage therefore repeats the fit with the same call `02_Differential` used, including
`sample.weights = TRUE`, and asserts that the recomputed coefficients match the published ones
before any weight is used. Without that check the tests could silently run on a different model
than the protein-level results everyone has already read.

## Multiplicity

BH runs within one method, one contrast, and one database. The three databases are separate
families, and Hallmark is the primary one, declared before any test ran. A set significant only
in Reactome or GO:BP is a lead, not a result. Nothing is pooled across methods: the methods are
alternatives under comparison, not replicates.

## What it found

On the negative control, fry calls nothing (median p 0.583, 1.6% below 0.05), camera calls one set
and fgsea three. The sets are not random: `REACTOME_STRIATED_MUSCLE_CONTRACTION` (camera and
fgsea), `REACTOME_MUSCLE_CONTRACTION` and `HALLMARK_HEME_METABOLISM` (fgsea). Contractile proteins
are the most abundant and most tightly co-varying in the sample, and heme is blood.

The interaction then shows why that matters. fry and camera call nothing. fgsea calls eight sets,
and every one of them belongs to one of those two families: muscle cell differentiation, myotube
differentiation, muscle contraction, membrane repolarisation, and `HALLMARK_HEME_METABOLISM` again,
now with the opposite sign.

| Set | fgsea p, interaction | fry p, interaction |
|---|---|---|
| `GOBP_MUSCLE_CELL_DIFFERENTIATION` | 2.0e-07 | 0.26 |
| `GOBP_MYOTUBE_DIFFERENTIATION` | 1.7e-05 | 0.23 |
| `HALLMARK_HEME_METABOLISM` | 3.6e-04 | 0.24 |
| `GOBP_MUSCLE_CONTRACTION` | 2.9e-04 | 0.34 |

Same sets, same matrix. The only difference is whether the test knows that 131 samples come from
33 people. The heme set is significant under fgsea on both the control (FDR 0.029) and the
interaction (FDR 0.015); fry gives it p = 0.20 and 0.24.

On the two training contrasts all three methods agree on direction, and fry is the most powerful:
271 sets on `BFR_post_vs_pre`, 235 on `HLRT_post_vs_pre`. Oxidative phosphorylation is the
strongest Hallmark call in both legs.
