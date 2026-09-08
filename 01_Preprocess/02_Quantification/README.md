# 01_Preprocess / 02_Quantification

Learns how missingness relates to abundance, rolls peptides up to proteins, then drops proteins
nobody detected often enough. This is the slow step.

| | |
|---|---|
| **Script** | `a_script/02_quantify.qmd` |
| **Reads** | `01_Filtering/c_data/precursors_filtered.rds` |
| **Writes** | `c_data/proteins.rds`, plus the fitted curve parameters and a per-protein quality table |
| **Next** | `02_Differential_Expression` |

## The idea

Faint peptides go missing more often than abundant ones. limpa fits that as a curve, so a
missing value becomes the statement "this was below the detection limit" rather than a blank to
be filled in. Every protein then gets a value in every sample.

This is not imputation. limpa estimates each protein's abundance from the peptides it did see
plus the probability that the ones it did not see were faint, and reports how uncertain the
answer is.

## What comes out

`proteins.rds` carries three matrices:

- **abundance**, one number per protein per sample, no gaps
- **standard error**, how much to trust each number
- **detections**, how many peptides were actually seen behind each number

The standard error is why this project uses limpa. It travels into the statistics as a weight, so
a protein built mostly from missing peptides counts for less than one measured directly. Anything
that pulls the abundances out into a plain matrix throws that away.

## Reading the curve

The notebook reports the fitted slope and never substitutes a preset. A shallow slope means less
information is recovered from missing values, which errs toward finding nothing rather than
finding too much, so it is the safe direction to be wrong in. limpa's guidance puts the usable
range between 0.1 and 1.0.

The curve's plot is not the diagnostic. It compares a fitted curve against proportions computed
a different way, so some mismatch is expected. Judge the slope.

## The detection filter

Runs after quantification, which is the order limpa specifies. A protein has to be detected in at
least as many samples as the smallest group holds. Set it lower and proteins resting on almost
nothing reach the statistics; set it higher and a protein present in one condition and absent in
the other gets dropped, which in a training study may be the interesting case.

## Cost

About 100 minutes and 11 GB. The step is cached and keyed to its input file, so it re-runs when
filtering produces a new matrix and not when the surrounding text changes. Delete
`a_script/02_quantify_cache/` to force it.

One protein group holds over three thousand peptides, and cost grows steeply with that count, so
a handful of very large proteins account for most of the runtime.
