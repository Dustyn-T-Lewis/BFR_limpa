# 02 · Differential Expression

Takes `proteins.rds` and asks which proteins changed. Two sub-stages, so a design problem surfaces
before anything is fitted.

```
01_Preprocess/02_Quantification/c_data/proteins.rds
  01_Design       design matrix, five contrasts, one diagnostic  -> design.rds
  02_Differential fit, contrasts, results tables                 -> fit.rds
```

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
```

## The five comparisons

| Name | Role | Compares | Reads as |
|---|---|---|---|
| `interaction` | **primary** | difference of the two time effects | the study question |
| `BFR_vs_HLRT_at_T1` | **control** | two legs before training | must find nothing |
| `BFR_post_vs_pre` | descriptive | one leg over time | what restricted training did |
| `HLRT_post_vs_pre` | descriptive | the other leg over time | what conventional training did |
| `BFR_vs_HLRT_at_T2` | descriptive | two legs after training | leg difference at the end |

`interaction` equals `BFR_vs_HLRT_at_T2` minus `BFR_vs_HLRT_at_T1` exactly, so these are not five
independent questions. Report the primary and the control as the result; the other three describe
what happened on the way.

## Three things that matter

**The control is the most important one.** Two untrained legs of the same person should differ in
nothing. If they do, the cause is upstream: a mislabelled sample, contamination, annotation.
`02_Differential` writes results to disk before it checks, and prints the count rather than
asserting on it, so a failure lands in the report instead of halting the render.

**Participant is a fixed term in the design.** Every comparison happens inside one person, so this
makes pre-to-post paired. It also means sex cannot be tested, since it does not vary within a
participant.

**Adjusted within each comparison, never pooled.** The five share participants and the interaction
is built from two of the others.

Testing goes through `dpcDE()`, which reads the standard errors from `proteins.rds`. A plain
`lmFit()` would discard them.

## What comes out

`02_differential.xlsx` holds two sheets: `DEP_matrix`, one row per protein with `logFC`,
`P.Value` and `adj.P.Val` for each contrast, and `contrast_summary`, the hit counts. `fit.rds`
is the complete fitted model. All carry full precision; rounding happens only in the report.

## Reading the counts

Hit counts near the FDR boundary move with the normalisation applied upstream, so treat the two
time contrasts as approximate. The interaction has returned no hits under every variant tried,
which is a statement about what this study could detect, not a demonstration that the two
training modes act alike. `02_Differential` prints the effect size it had 80% power to find, and
that number is what bounds the claim.
