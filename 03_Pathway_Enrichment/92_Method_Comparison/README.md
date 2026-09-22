# 03_Pathway_Enrichment / 92_Method_Comparison

Every pathway set test run on this dataset, from every branch, in one table with the same column
names, compared on the same sets, the same contrasts and the same negative control. This stage
runs no test of its own.

| | |
|---|---|
| **Scripts** | `a_script/92_method_comparison.qmd`, `a_script/92_snapshot_external.R` |
| **Reads** | `90_Set_Test_Calibration/c_data/set_test_results.csv`, `91_Sample_Scores_ssGSEA/c_data/score_contrast_results.csv`, `c_data/external/*.csv.gz` |
| **Writes** | `c_data/92_method_comparison.xlsx`, `comparison_long.csv.gz`, `comparison_shared_sets.csv.gz`, `control_calibration.csv`, `same_method_agreement.csv`, `call_counts_shared_sets.csv`, `hallmark_side_by_side.csv` |

```sh
quarto render 03_Pathway_Enrichment/92_Method_Comparison/a_script/92_method_comparison.qmd --output-dir ../b_reports
```

## What is compared

| Source | Branch / PR | Commit | Tests |
|---|---|---|---|
| `diego_90` | this branch, `90_Set_Test_Calibration` | this commit | fry, camera, fgsea |
| `diego_91` | this branch, `91_Sample_Scores_ssGSEA` | this commit | ssGSEA scores tested with limma |
| `dustyn_pr11_fgsea` | PR #11 | `b16fc4f` | fgsea, per database |
| `dustyn_pr11_fgsea_pooled` | PR #11 | `b16fc4f` | fgsea, databases pooled (long table only) |
| `dustyn_pr11_fry` | PR #11 | `b16fc4f` | fry |
| `dustyn_pr15_fgsea` | PR #15 | `aa66aa0` | fgsea, ten collections |
| `matheus_pr13` | PR #13 | `62a2404` | camera and fry |

The other branches' results are snapshots, copied by `92_snapshot_external.R` and pinned in
`c_data/external/external_manifest.csv` with the commit, the path and the md5 of the original file.
To check any snapshot against its branch:

```sh
git fetch origin pull/11/head pull/13/head pull/15/head
git show b16fc4f:03_Pathway_Enrichment/02_Set_Tests/c_data/fry_results.csv | md5sum
```

If a branch moves, rerun the snapshot script and then the notebook. Singscore is not in the table:
it gives per-sample scores without a test, so there is no p-value to line up.

## How to read it

- **One long table.** `comparison_long.csv.gz` has one row per test, contrast and set, with the
  same columns for every source: `test`, `source`, `method`, `contrast`, `set_id`, `database`,
  `pathway`, `n`, `direction`, `effect`, `p`, `padj_reported`, `fdr_family`. `set_id` is the
  MSigDB name with its database prefix, so rows from different branches join on it directly.
- **Same sets.** The set lists differ between branches (floors of 5, 10 or 15 measured members,
  different overlap filters, different collections). Counts are compared only on the 996 sets
  every test ran (41 Hallmark, 195 Reactome, 760 GO:BP).
- **Same FDR family.** Each branch corrected across a different number of sets, so the reported
  FDR is not comparable. `padj_common` is recomputed as BH within one test, one contrast and one
  database, over the shared sets. Raw p-values need no correction and are compared as they are.
- **Excel.** `92_method_comparison.xlsx` holds the small tables, with a `read_me` sheet.

## What it shows

**The same function gives the same answer across branches.** Our fry and PR #11's fry agree to
the last digit on all 5,140 shared set-contrast pairs (largest difference in log10 p: 7e-14):
same weights, same fit, same sets. Matheus's fry and camera agree in direction on every set and
rank sets identically (Spearman 1.0); their p-values differ slightly because PR #13 refits with
`dpcDE(sample.weights = TRUE)`. fgsea agrees across branches up to its sampling noise.

**The negative control separates the methods, in every branch.** On the shared sets:

| Test | False calls at FDR < 0.05 | Median p | p < 0.05 |
|---|---|---|---|
| fry (ours, PR #11, PR #13) | 0 | 0.58 | 1.6% |
| ssGSEA + limma (ours) | 0 | 0.51 | 3.0% |
| camera (ours, PR #13) | 1 | 0.54 | 4.2% |
| fgsea (ours, PR #11, PR #15) | 3 | 0.51 | 5.9% |

The false calls are the same sets wherever they appear: `REACTOME_STRIATED_MUSCLE_CONTRACTION`
(camera and fgsea), `REACTOME_MUSCLE_CONTRACTION` and `HALLMARK_HEME_METABOLISM` (fgsea). On its
full ten-collection list PR #15's fgsea calls 14, the number its own README reports.

**The interaction is null for every test that passes the control cleanly.** fry, camera and ssGSEA
call nothing. fgsea, in all three branches that run it, calls heme metabolism and the contractile
GO terms: the same families it calls on the control.

**The training effects converge.** Oxidative phosphorylation is up in both legs under all nine
tests. In BFR, fatty-acid metabolism up and myogenesis down are called by eight of nine. Calls on
`BFR_vs_HLRT_at_T2` come only from camera and fgsea, the two tests that also call sets on the
control; fry and ssGSEA call nothing there.

The numbers above come from the committed render; `control_calibration.csv`,
`same_method_agreement.csv` and `hallmark_side_by_side.csv` hold them in full.
