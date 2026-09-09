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

Reported as fitted, never replaced with a preset. limpa's usable range is 0.1 to 1.0. A shallow
slope recovers less from missing values, so it errs toward finding nothing.

The curve's plot is not the diagnostic; it compares a fitted curve against proportions computed
differently. Judge the slope.

## Detection filter

Runs after quantification, which is the order limpa specifies. A protein must be detected in at
least as many samples as the smallest group holds.

About 100 minutes and 11 GB. Cached and keyed to its input file. Delete
`a_script/02_quantify_cache/` to force it.
