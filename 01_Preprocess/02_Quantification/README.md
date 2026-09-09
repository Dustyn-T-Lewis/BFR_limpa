# 01_Preprocess / 02_Quantification

Fits the detection curve, rolls peptides up to proteins, filters on detection. The slow step.

| | |
|---|---|
| **Script** | `a_script/02_quantify.qmd` |
| **Reads** | `01_Filtering/c_data/precursors_filtered.rds` |
| **Writes** | `c_data/proteins.rds`, the curve parameters, a per-protein quality table |

## The idea

Faint peptides go missing more often than abundant ones. `dpc()` fits that, so a missing value
becomes "this was below the detection limit" rather than a blank to fill in. `dpcQuant()` then
gives every protein a value in every sample.

This is not imputation. Each protein's abundance is estimated from the peptides that were seen
plus the probability that the missing ones were faint, and the answer comes with a standard
error.

## What comes out

`proteins.rds` carries three matrices: abundance, standard error, and how many peptides were
detected behind each value. The standard errors travel into the statistics as weights, so a
protein built mostly from missing peptides counts for less than one measured directly.

## The slope

Two estimators are reported. The complete-normal model is the one used, and the observed-normal
model is printed beside it as an independent check; they disagree by about a factor of two on this
data. Neither is replaced with a preset. limpa's usable range is 0.1 to 1.0, and a shallow slope
recovers less from missing values, so it errs toward finding nothing.

`a_script/02b_quantify_slope07.qmd` refits at a fixed slope of 0.7 as a sensitivity check. It
re-runs the slow step, so it keeps its own cache.

The curve's plot is not the diagnostic; it compares a fitted curve against proportions computed
differently. Judge the slope.

## Detection filter

Runs after quantification, which is the order limpa specifies. A protein must be detected in at
least as many samples as the smallest group holds.

The last line reports the range of per-sample medians. Any normalisation decision rests on that
number, and limpa places the correction in the next stage, between `dpcQuant()` and `dpcDE()`.

About 100 minutes and 11 GB. Cached and keyed to its input file. Delete
`a_script/02_quantify_cache/` to force it.
