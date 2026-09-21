# 05_classify_pathway_themes

Groups the tested sets into GO themes, scores each theme per sample, and lays out how themes move
with training and how they couple to the phenotype.

| | |
|---|---|
| **Script** | `a_script/05_classify_pathway_themes.R` |
| **Reads** | `gene_sets.rds`, `fgsea.rds`, `singscore.rds`, `proteins.rds`, `phenotype.csv` |
| **Writes** | `c_data/pathway_themes.rds`, `05_classify_pathway_themes.xlsx`, two figures |

## A theme is this pipeline's module

2,732 sets is too many to read and most of them restate each other. A theme is a GO ancestor
holding at least `theme_min_sets` qualifying members, and its score is the **mean singscore of
those members**, so one row replaces a dozen nested terms saying the same thing.

At the current threshold of 8: **27 themes covering 261 sets**. Raise `theme_min_sets` for fewer
and broader rows, lower it for more and narrower ones.

## Nominal p only

Nothing in this stage is corrected and nothing is gated on q. A **dashed cell border** marks
p < 0.05 and **bold text** marks |r| at or above 0.4. Both are reading aids, not claims.

`nominal_counts` is the first sheet in the workbook and gives the observed count against the
count expected by chance at that threshold. At 27 themes, one hit per outcome is what chance
produces.

## The figures

**`theme_classification.png`** — panel A is the association map, three blocks left to right in
increasing strength of claim: mean NES per theme on the training contrasts, then Spearman r
against four phenotype change scores pooled across all 65 legs, then the same r taken
**within participant** as BFR minus high load across 32 pairs. Rows are ordered by clustering the
NES profile, so themes that behave alike sit together. Two colour scales, because NES and r are
different quantities.

Panel B draws the six strongest differential couplings: one point per leg, a fit per arm, and an
inline box giving r and nominal p for each arm. The header carries the within-participant r that
selected the panel.

**`theme_members.png`** — the supplement. Every qualifying member set behind each drawn theme,
positioned by its own NES, sized by measured genes, coloured by whether it reaches nominal p. A
theme is only as coherent as its members, and this is where that gets checked.

## What it shows

The training block is the readable result: cellular respiration, mitochondrion organization and
carbohydrate metabolism move **up** under both BFR and high load; muscle contraction and muscle
system process move **down** under both. The interaction column sits near zero for most themes.

Both modalities move the same themes in the same direction, and no phenotype coupling clears its
chance expectation.
