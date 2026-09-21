# 03_Pathway_Enrichment / 02_Set_Tests

Tests whether whole gene sets moved, instead of asking three thousand separate
protein-by-protein questions.

| | |
|---|---|
| **Script** | `a_script/02_set_tests.qmd` |
| **Reads** | `01_Preprocess/02_Quantification/c_data/proteins.rds`, `02_Differential_Expression/01_Design/c_data/design.rds` |
| **Writes** | `c_data/geneset_tests_long.csv`, `theme_cut_significant.csv` |

## Why `camera()` and `fry()`

The parent README rules out two classes of test: one that treats proteins as independent (wrong
for co-regulated sets — mitochondrial and ribosomal proteins move together — and it reported hits
on the negative control when tried) and one that takes a repeated-measures argument and silently
ignores it. Neither `camera()` (Wu & Smyth 2012) nor `fry()` needs that argument here: participant
is a fixed column inside `design`, the same design `02_differential.qmd` fits on, so the
repeated-measures structure is inside the model itself, not passed in on the side where it could
be dropped.

`camera()` models inter-gene correlation directly and is competitive — it asks whether a set moved
more than the rest of the (also shifting) proteome. `fry()` is a fast rotation approximation to
`roast()`/`mroast()` and is self-contained — it asks whether a set's own genes moved, without
reference to the rest. They are not run as two parallel counts: see "Reading significance" below
for how the two combine into one number.

## Gene sets

Four full MSigDB collections (`msigdbr`, `db_version` printed in the render):

| Collection | Sets before matching | Sets tested (≥5 matched proteins) |
|---|---|---|
| Hallmark | 50 | 48 |
| KEGG (legacy) | 186 | 149 |
| Reactome | 1,839 | 991 |
| GO:BP | 7,538 | 3,439 |

A set below 5 matched proteins is dropped before testing, not after — both tests degrade to
unstable low-df estimates below that, and keeping them would only add noise to the correction.

A fifth, single targeted set (`WP_EXERCISE_AND_HYPERTROPHY_IN_SKELETAL_MUSCLE`) was pre-registered
and run in an earlier version of this script, before any result was seen. It was dropped on review,
before this version's numbers were produced: it was not one of the four collections the group had
actually asked for, and only 4 of its 20 genes were quantified here, so it could never have
returned an answer either direction. Dropping an under-powered, non-requested addition on review is
not the failure pre-registration guards against — that failure is dropping a test *because its
result was inconvenient*, and this set never reached the point of producing a result to be tempted
by. If a manuscript reviewer asks whether this set was tried: yes, and here is why it was removed.

## Contrasts

Three of `01_Design`'s five contrasts, plus the negative control, run here — not all five, for the
same reason `02_differential.qmd` treats them as one family rather than five independent ones:

| Contrast | Role |
|---|---|
| `interaction` | primary — the question the study was designed to answer |
| `BFR_post_vs_pre` | descriptive |
| `HLRT_post_vs_pre` | descriptive |
| `BFR_vs_HLRT_at_T1` | negative control — checked first, see below |

## Multiple testing

Benjamini-Hochberg within each contrast-by-collection cell, never pooled. Hallmark's 48 sets and
GO:BP's 3,439 are different families; correcting them together would let GO:BP's bulk set the
threshold for Hallmark.

## Reading significance: `camera()` decides, `fry()` confirms

`camera_fdr` and `fry_fdr` are not two parallel counts of "how many pathways matter" — reporting
them that way is not how gene-set tests are normally read in an applied paper, and it asks the
reader to referee a statistics disagreement instead of reading a result. The convention used here:
**`camera()` is the significance test; `fry()` is a confirmation filter.** A set counts as a `hit`
(a column in `geneset_tests_long.csv`) only when `camera_fdr < 0.05` **and** `fry_fdr < 0.05`. The
reverse — `fry()` alone flagging something `camera()` misses — is not a hit, and is not reported as
one: `fry()` is self-contained, so on its own it says nothing about whether a set is specifically
implicated or just swept up in a broad shift.

This rule is not arbitrary — the negative control below is the concrete case it exists to catch.

## Negative control

Checked before anything else, per the parent README's instruction. Two untrained legs of the same
person, across 4,627 tests: `camera()` alone returns 1 hit in GO:BP (of 3,439) and 1 in Reactome
(of 991) — at or below the 5% the correction targets on its own — but `fry()` does not confirm
either one. **`n_hit` is zero in every collection.** Both the false-positive rate `camera()` alone
would suggest, and the correction the `hit` rule applies to it, check out here.

## Results

| Contrast | Collection | Tested | `camera` alone | `fry` alone | **`hit` (both)** |
|---|---|---|---|---|---|
| `interaction` | all four | 4,627 | 0 | 0 | **0** |
| `BFR_post_vs_pre` | Hallmark / KEGG / Reactome / GO:BP | 48 / 149 / 991 / 3,439 | 4 / 18 / 151 / 126 | 10 / 36 / 270 / 777 | **4 / 17 / 144 / 117** |
| `HLRT_post_vs_pre` | Hallmark / KEGG / Reactome / GO:BP | 48 / 149 / 991 / 3,439 | 2 / 18 / 120 / 61 | 7 / 35 / 280 / 708 | **2 / 18 / 118 / 61** |

`hit` stays close to the `camera()`-alone count on the two training contrasts (most of what
`camera()` finds, `fry()` agrees with) and drops to zero on the negative control (where `camera()`'s
one finding in each of two collections is not confirmed) — the pattern the rule is meant to
produce, not a coincidence of this one run.

The interaction — the primary contrast — is null in every collection, on both tests. That matches
the protein-level result in `02_Differential` and the power calculation already on record: this
design cannot detect an interaction of the size actually observed, so a null here confirms rather
than adds information.

Both training arms shifted large parts of the measured proteome, distributed across hundreds of
sets in Reactome and GO:BP — consistent with the parent README's premise that pathway testing is
the right stage when the protein-level result is thin.

## For the manuscript

Report `hit` (the `camera()`+`fry()`-confirmed count) in the main text as the finding. Keep
`camera_p`/`camera_fdr`/`fry_p`/`fry_fdr` side by side in a supplementary table, including the
cases where they disagree (the negative control's two lines above) — that is what lets a reader
audit where the two tests agreed and where the confirmation rule did its job, without asking the
main text to carry two competing counts.

## Thematic classification

Every tested set is matched against keyword patterns chosen for the study's hypotheses — applied to
**all** sets before any significance filter, so the classification itself cannot be shaped by which
sets happened to come out significant. One theme per set; a set matching no pattern is left
unclassified (most of GO:BP's 3,439 sets, which cover far more than this study asks about).

| Theme | What it covers |
|---|---|
| `Myofibrillar` | sarcomere / contractile apparatus (skeletal, not cardiac/smooth) |
| `Sarcoplasmic_Glycolytic` | glycolysis, glycogen metabolism |
| `Sarcoplasmic_SR_Calcium` | sarcoplasmic reticulum, calcium transport |
| `Sarcoplasmic_Mitochondrial` | mitochondria, oxidative phosphorylation, TCA cycle |
| `Metabolic_Stress` | hypoxia, angiogenesis/vascular, AMPK, oxidative/ER stress, heat shock |
| `Ribosomal` | ribosome biogenesis and structure |
| `Protein_Synthesis` | translation (initiation/elongation/termination), mTORC1 |
| `Protein_Degradation` | ubiquitin-proteasome, autophagy-lysosome |

`Metabolic_Stress` folds hypoxia and angiogenesis/vascular terms together — both trace back to the
same BFR mechanism (local ischemia under occlusion), not two separate questions. `Ribosomal` and
`Protein_Synthesis` are kept apart deliberately: more ribosomes without more translation initiation
is a different finding from the reverse.

Three sets from the earlier pilot analysis (`REACTOME_STRIATED_MUSCLE_CONTRACTION`,
`KEGG_GLYCOLYSIS_GLUCONEOGENESIS`, `GOBP_SARCOPLASMIC_RETICULUM_CALCIUM_ION_TRANSPORT`) are flagged
with `is_pilot_set = TRUE` so they can be found in this larger run without a separate lookup file.
Only the last of the three is confirmed here (`hit = TRUE` in both training contrasts) — the pilot's
own Reactome myofibrillar set clears `camera()` alone but not `fry()`, which is exactly the failure
mode the `hit` rule exists to catch, not a contradiction of the pilot.

### On-theme, confirmed hits (`theme_cut_significant.csv`)

| Contrast | Theme | Confirmed hits |
|---|---|---|
| `BFR_post_vs_pre` | Myofibrillar | 12 |
| `BFR_post_vs_pre` | Sarcoplasmic_Mitochondrial | 33 |
| `BFR_post_vs_pre` | Sarcoplasmic_SR_Calcium | 5 |
| `BFR_post_vs_pre` | Protein_Degradation | 12 |
| `BFR_post_vs_pre` | Metabolic_Stress | 2 |
| `HLRT_post_vs_pre` | Myofibrillar | 5 |
| `HLRT_post_vs_pre` | Sarcoplasmic_Mitochondrial | 30 |
| `HLRT_post_vs_pre` | Sarcoplasmic_SR_Calcium | 2 |
| `HLRT_post_vs_pre` | Sarcoplasmic_Glycolytic | 1 |
| `HLRT_post_vs_pre` | Protein_Degradation | 12 |
| `HLRT_post_vs_pre` | Metabolic_Stress | 2 |
| any | Ribosomal, Protein_Synthesis | 0 |
| `interaction` | any theme | 0 |

Myofibrillar falls and mitochondrial rises in **both** arms, not one versus the other — see the
direction columns (`camera_dir`/`fry_dir`) in the file for the per-set detail; the interaction
being null everywhere is the same finding stated a different way. Protein degradation
(ubiquitin-proteasome + autophagy-lysosome) moves with training in both arms too, at an identical
count (12 and 12) — worth a look for direction and overlap with the mitochondrial/myofibrillar
sets before reading into it. Ribosomal and protein-synthesis themes have zero confirmed hits in
either arm: no evidence here of a translation-machinery-specific response, in either direction.

## Outputs

Two files. `geneset_tests_long.csv` is every set tested, with `theme`/`is_pilot_set` attached, so
any cut — significant or not, on-theme or not — is one filter away without re-running anything.
`theme_cut_significant.csv` is the narrow, curated shortlist: on-theme sets confirmed by both
tests. The negative-control summary and the contrast-by-collection counts in this README are not
written to disk separately — both are one line of code on `geneset_tests_long.csv` and already sit
in the render (`b_reports/02_set_tests.html`) and here, so a third small file would only duplicate
them. The protein-index object a per-sample scoring step will need is not written here either — the
parent README scopes per-sample scoring ("score every set in every sample") as this stage's own
second half, not `03_Enrichment` (a separate over-representation analysis on hit lists, described
below in that README's own words — see the note at the end of this section). It takes under a
second to rebuild from this script when that half is actually started.

**Correction (2026-09-18):** an earlier version of this README called the per-sample-scoring step
`03_Enrichment`. It isn't — the parent `03_Pathway_Enrichment/README.md` puts "score every set in
every sample" inside `02_Set_Tests` itself, and reserves `03_Enrichment` for over-representation
analysis on a hit list ("for comparisons with hits, ask which processes those hits belong to"), a
different method entirely. Per-sample scoring (ssGSEA) is unstarted work that belongs in *this*
folder, not a pointer to a stage that hasn't begun.
