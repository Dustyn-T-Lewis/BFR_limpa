# 01_Preprocess / 02_Quantification

Fits the detection probability curve, rolls precursors up to protein groups, then filters on
detection. This is the slow stage.

| | |
|---|---|
| **Script** | `a_script/02_quantify.qmd` |
| **Reads** | `01_Filtering/c_data/precursors_filtered.rds` |
| **Writes** | `c_data/proteins.rds`, `c_data/dpc_parameters.csv`, `c_data/protein_quality.csv` |
| **Read by** | both sub-stages of `02_Differential_Expression` |

## Why the DPC is not its own stage

It is three lines and two plots, and limpa's case study runs it inline just before
`dpcQuant()`. A separate stage bought one checkpoint, and `#| cache: true` on the
`dpcQuant` chunk already gives you that: the curve and its plots render in seconds at the
top of the notebook, so you see a bad fit long before the slow step runs again.

## What a DPC is

Faint precursors go missing more often than bright ones. `dpc()` fits that relationship,
giving the probability a precursor is detected as a function of how much of it was there.
`dpcQuant()` then reads a missing value as "this was below the detection threshold" rather
than filling in a number and forgetting where it came from.

Two defaults are worth knowing. `subset = 2000` fits the curve on a systematic sample
stratified by missingness and mean intensity, not on every row. `robust` applies only when
`model = "on"`, and the default is `model = "cn"`, so outlier downweighting is off and the
subsample supplies the robustness. Under `model = "cn"` the help page says `dpc()` returns
`dpcCN()`, which the FAQ recommends for noisy data, so calling `dpcCN()` separately would
change nothing.

The slope is reported as fitted. The stage README says why we never substitute a preset.

## What dpcQuant() returns

An `EList` with one row per protein group, and two extra matrices in `$other`.
`standard.error` says how certain each value is. `n.observations` counts the precursors
actually detected behind it. It also adds `NPrec` and `PropObs` to `$genes`, though the help
page calls the first of those `NPeptides`.

Those standard errors are why this project uses limpa. A protein rebuilt mostly from missing
precursors still gets a value, but a wide one, and `dpcDE()` turns that width into a
precision weight. The FAQ says plainly that `lmFit`, `arrayWeights` and `vooma` must not run
on a limpa object, and that pulling the matrix out and handing it to ordinary limma costs
power and can cost error-rate control.

The `dpcQuant` chunk is cached, so re-rendering to fix prose does not run it again. Delete
`a_script/02_quantify_cache/` to force it.

## Detection filtering

limpa fixes this order explicitly. `filterByDetection()`'s help says it runs after
`dpcQuant` and before DE, and it returns a logical vector rather than an object, so you
write `y[keep, ]`.

We set `n.samples` to the smallest group size, which the help page suggests for a small
experiment. That group is `BFR_T2` with 32 samples, one short because nobody acquired S06's
left leg at T2.

## The MDS plot

`plotMDSUsingSEs()` rather than `plotMDS()`, because limpa's version accounts for the
standard errors. Do not expect samples to separate by group. Every participant contributes
all four cells, so differences between participants dominate the first dimensions, and the
participant term in the design is what removes them. Investigate a sample sitting far from
everything else. Do not read anything into the lack of group structure.
