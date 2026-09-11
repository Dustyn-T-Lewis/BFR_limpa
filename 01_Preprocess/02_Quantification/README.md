# 01_Preprocess / 02_Quantification

Produces the analysis-ready protein matrix: fits the detection curve, rolls precursors up to
proteins, filters on detection, normalises. This is the slow stage.

| | |
|---|---|
| **Script** | `a_script/02_quantify.qmd` |
| **Reads** | `01_Filtering/c_data/precursors_filtered.rds` |
| **Writes** | `c_data/proteins.rds`, the curve parameters, a per-protein quality table |
| **Beside it** | `a_script/sensitivity_dpc_slope.qmd`, a sensitivity check |

## The detection curve

Faint precursors go missing more often than abundant ones. `dpc()` fits that, so `dpcQuant()` can
treat a missing value as evidence the protein was low rather than filling in a number.

Two estimators are reported, complete normal and observed normal. The complete-normal fit is the one
used; they disagree by about a factor of two on this data. Neither is replaced with a preset. limpa's
usable range is 0.1 to 1.0, and a shallow slope recovers less from missing values, so it errs toward
finding nothing. `sensitivity_dpc_slope.qmd` refits at a fixed 0.7 and reports how far the
estimates move.

The curve's plot is not the diagnostic; it compares a fitted curve against proportions computed
differently. Judge the slope.

## What comes out

`proteins.rds` carries three matrices: abundance, standard error, and how many precursors were
detected behind each value. The standard errors travel into `dpcDE()` as precision weights, so a
protein built mostly from missing precursors counts for less than one measured directly.

`NPrec` is what the installed version writes; the help page calls it `NPeptides`.

## Detection filter

Runs after quantification, which is the order limpa specifies. A protein must be detected in at
least as many samples as the smallest group holds.

## Normalisation

limpa places this after `dpcQuant()` and before `dpcDE()`, and performs none itself, because DIA-NN
has usually normalised the precursors already.

The diagnostic is the MA curve: each sample is compared against the average of all samples, and the
difference is fitted against abundance. A flat curve means a sample disagrees by a constant, which a
shift corrects. A curve that travels means the disagreement depends on abundance.

Here the curve travels 0.45 log2 in the median sample, against 0.13 for the spread in per-sample
medians, and the curvature is unrelated to treatment or timepoint. Cyclic loess is therefore applied:
it is the only available method that corrects an abundance-dependent difference. It also has the most
freedom of the available methods, so it can absorb real biology that tracks abundance.

The matrix is normalised before it is written, so `proteins.rds` is what the model is fitted to.

## Cost

About 100 minutes and 11 GB. The quantification chunk is cached and keyed to the md5 of its input
file. Delete `a_script/02_quantify_cache/` to force it.
