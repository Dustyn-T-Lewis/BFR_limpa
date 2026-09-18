# 01_Preprocess / 02_Quantification

Produces the analysis-ready protein matrix: fits the detection curve, rolls precursors up to
proteins, filters on detection, normalises. This is the slow stage.

| | |
|---|---|
| **Scripts** | `a_script/02_quantify_run.R` computes, `a_script/02_quantify.qmd` reports |
| **Reads** | `01_Filtering/c_data/precursors_filtered.rds` |
| **Writes** | `c_data/proteins.rds`, `c_data/02_quantify.xlsx`, `c_data/quant_runs/` |

Two files, because `dpcQuant()` costs about 100 minutes per call and a render should not. The
notebook **only loads**: it reads `c_data/quant_runs/proteins_slope_0.7.rds`, stops with the command
below if that file is absent, and cannot start a long run whatever you pass it.

```sh
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
```

The long run is yours to start, by hand, when the precursor matrix changes:

```sh
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R
Rscript 01_Preprocess/02_Quantification/a_script/02_quantify_run.R --sensitivity
```

The first writes the primary run at the preset slope; the second adds the fitted-slope sensitivity
run, another 100 minutes. Each file records the md5 of the precursor matrix it came from and the
limpa version that produced it. Nothing checks either, so rerunning `01_Filtering` means rerunning
this too: the notebook will happily load a checkpoint built from older precursors.

`proteins_slope_0.7.rds` is committed, so a fresh clone renders without running anything. It
predates the version field, so it records only the input md5.

## Why a preset slope

`dpc()` fits the curve that lets `dpcQuant()` treat a missing value as evidence the protein was
low rather than filling in a number. On this data the fit lands below the 0.7-to-0.9 band limpa's
FAQ calls typical for DIA-NN searched with match-between-runs, though inside the 0.1-to-1.0 range
it calls usable.

Quantification therefore uses a preset **0.7**, the rule the FAQ gives for exactly this case, and
reports the fitted slope as a sensitivity run. The choice follows the documented rule and was made
before any downstream count was consulted. The notebook prints both slopes side by side.

One direction is worth recording because it is not the one the documentation predicts. The FAQ
says a slope set too low makes the analysis more conservative; here the lower fitted slope gave
the smaller standard errors and the larger hit counts.

## Numbers measured separately

Refitting the five contrasts under each slope is **not** part of this stage. Done separately, it
gave 60 and 107 proteins on the two training contrasts under the preset against 82 and 114 under
the fitted slope, with the negative control and the interaction empty under both.

The same applies to the detection filter. The proteins it drops carried a mean standard error of
1.56 against 0.66 for those kept, and 0.15 detected precursors per sample against 7.48. Relaxing
the threshold to limpa's default of 3 recovered 371 of them, of which exactly one reached
BH < 0.05 anywhere, while the extra tests cost 14 and 29 proteins on the two training contrasts.

## Order

`filterByDetection()`'s help page fixes it: after `dpcQuant()`, before normalisation and
`dpcDE()`. limpa performs no normalisation itself, because DIA-NN has usually normalised the
precursors already. The notebook measures whether a correction is needed before applying one, and
`proteins.rds` is written afterwards, so it is the matrix the model is fitted to.

## Caveat to carry forward

Cyclic loess remaps `$E` non-linearly while `$other$standard.error` is left untouched, so after
normalisation the stored errors describe the unnormalised scale and `dpcDE()` reads them as though
they described the new one. Smyth recommends quantile and cyclic loess for this data all the same,
and the mapping is close to linear, but it is an approximation. The notebook repeats this beside
the code.

`NPrec` is what the installed version writes; the help page calls it `NPeptides`.
