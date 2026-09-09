# 02 · Differential Expression

Takes `proteins.rds` and asks which proteins changed. Two sub-stages, so a design problem shows up
before an hour of model fitting does.

```
01_Preprocess/02_Quantification/c_data/proteins.rds
  01_Design       design matrix, five contrasts, one diagnostic  -> design.rds
  02_Differential normalise, fit, contrasts, results             -> protein_contrasts_long.csv
```

`02_Differential/a_script/02b_normalization_report.qmd` decides which normalisation the protein
matrix needs. Run it when that matrix changes, not on every render.

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd        --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

## The five comparisons

| Name | Compares | Reads as |
|---|---|---|
| `BFR_post_vs_pre` | one leg over time | what restricted training did |
| `HLRT_post_vs_pre` | the other leg over time | what conventional training did |
| `BFR_vs_HLRT_at_T1` | two legs before training | **control**, must find nothing |
| `BFR_vs_HLRT_at_T2` | two legs after training | leg difference at the end |
| `interaction` | difference of the two time effects | **the study question** |

## Three things that matter

**The control is the most important one.** Two untrained legs of the same person should differ in
nothing. If they do, the cause is upstream: a mislabelled sample, contamination, annotation. The
stage asserts it comes back empty and writes results to disk before that check runs.

**Participant is a fixed term in the design.** Every comparison happens inside one person, so this
makes pre-to-post paired. It also means sex cannot be tested, since it does not vary within a
participant.

**Adjusted within each comparison, never pooled.** The five share participants and the interaction
is built from two of the others.

Testing goes through `dpcDE()`, which reads the standard errors from `proteins.rds`. A plain
`lmFit()` would discard them.

## What comes out

`protein_contrasts_long.csv` has one row per protein per contrast, sorted by p-value within each
contrast. `contrast_summary.csv` counts hits per contrast. `fit.rds` is the fitted model.

Hit counts near the FDR boundary move with the normalisation method, so report the two time
contrasts as approximate. The interaction is empty under every method tried.
