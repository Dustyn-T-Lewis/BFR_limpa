# 04_classify_and_associate_themes

Groups the tested sets into GO themes, then answers two questions about them: how well a theme
separates the study's groups, and whether it tracks the phenotype.

| | |
|---|---|
| **Script** | `a_script/04_classify_and_associate_themes.R` |
| **Reads** | `gene_sets.rds`, `fgsea.rds`, `singscore.rds`, `proteins.rds`, `phenotype.csv` |
| **Writes** | one workbook, `theme_results.rds`, two figures |

## A theme is this pipeline's module

2,732 sets is unreadable and most of them restate each other. A theme is a GO ancestor holding at
least `theme_min_sets` qualifying members, scored as the **mean singscore of those members**. At
the current threshold of 8 that is **27 themes covering 261 sets**.

## Every comparison here is paired

Pre versus post is the same leg twice. BFR versus high load is two legs of one person. A standard
Mann-Whitney AUC treats them as independent and returns a p-value for a design this study did not
run.

So **AUC comes from the rank sums as the effect-size descriptor, and p comes from the paired
Wilcoxon signed-rank test.** Both are reported. Doing it properly costs one second across all
13,660 set-level tests, so there is no reason to take the shortcut.

## Nominal p only

Nothing in this substage is corrected and nothing is gated on q. The honest substitute is
`chance_expectation`, the first data sheet in the workbook: observed nominal count against the
count a table that size returns under the null. At 27 themes and p < 0.05 that is 1.4 per task.

`baseline_BFR_vs_HLRT` is the empirical floor. It compares the two legs of one person before
either was trained, so whatever it returns is noise by construction.

## The five tasks

| Task | Compares | Pairs | Nominal of 27 | Ratio to chance |
|---|---|---:|---:|---:|
| `pre_vs_post_BFR` | T2 vs T1, BFR legs | 32 | 6 | **4.29** |
| `pre_vs_post_HLRT` | T2 vs T1, high-load legs | 33 | 4 | **2.86** |
| `baseline_BFR_vs_HLRT` | BFR vs HLRT leg at T1, control | 33 | 1 | 0.71 |
| `post_BFR_vs_HLRT` | BFR vs HLRT leg at T2 | 32 | 3 | 2.14 |
| `delta_BFR_vs_HLRT` | change, BFR vs HLRT leg | 32 | 1 | 0.71 |

Training separates pre from post well above chance in both arms. The baseline control sits below
chance, which is what it should do. BFR against high load barely clears it, and the change task
sits at the floor.

## The figures

**Only significant results are drawn.** The full map of every theme against every task and every
outcome lives in the workbook, not on a figure. One figure per comparison, so each answers one
question.

| Figure | Panels |
|---|---:|
| `roc_pre_vs_post_BFR` | 6 |
| `roc_pre_vs_post_HLRT` | 4 |
| `roc_post_BFR_vs_HLRT` | 3 |
| `roc_delta_BFR_vs_HLRT` | 1 |
| `roc_baseline_BFR_vs_HLRT` *(control)* | 1 |
| `association_pooled` | 2 |
| `association_differential` | 6 |

Each ROC figure carries one curve per theme reaching nominal p, with AUC and p in the strip.
**Red means the theme is higher in the favoured group, blue means lower** — an AUC below 0.5 is a
result about direction, not a failure, so it is coloured rather than hidden.

The association figures draw one point per leg, a fit within each arm in the same red and blue,
and each arm's own r and p inline. The header r is the correlation that selected the panel.

Nothing is capped. A `stopifnot` compares the panels drawn against the nominal counts, so "all
significant ones" stays literal.

## The workbook

One file, index first. `read_me` names every sheet and its row count. `theme_*` sheets are the
27-row analysis; `set_*` sheets carry the same computations across all 2,732 sets so a strong
single set inside a quiet theme is still findable; `theme_membership` is the drill-down path from
a theme to the sets behind it.
