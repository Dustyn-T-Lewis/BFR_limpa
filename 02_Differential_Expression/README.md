# 02 · Differential Expression

Takes `proteins.rds` and asks what happened to the proteins. Three sub-stages, so a design problem
surfaces before anything is fitted.

```
01_Preprocess/02_Quantification/c_data/proteins.rds
  01_Design       design matrix, five contrasts, one diagnostic  -> design.rds
  02_Differential fit, contrasts, results tables                 -> fit.rds
  03_Phenotype    every protein against the ultrasound outcomes  -> phenotype.rds
```

```sh
quarto render 02_Differential_Expression/01_Design/a_script/01_design.qmd             --output-dir ../b_reports
quarto render 02_Differential_Expression/02_Differential/a_script/02_differential.qmd --output-dir ../b_reports
quarto render 02_Differential_Expression/03_Phenotype/a_script/03_phenotype.qmd       --output-dir ../b_reports
```

`03_Phenotype` asks a question the contrasts cannot: not whether a protein changed, but whether
its change tracks how much the muscle changed. Nothing survives BH across 24,336 correlations and
every nominal count sits within what 3,042 tests return under the null, which is a clean negative
at this sample size.

## The five comparisons

| Name | Role | Compares | Reads as |
|---|---|---|---|
| `Modality_x_Time_Interaction` | **primary** | difference of the two time effects | the study question |
| `BFR_Pre-HLRT_Pre` | **control** | two legs before training | must find nothing |
| `BFR_Post-Pre` | descriptive | one leg over time | what restricted training did |
| `HLRT_Post-Pre` | descriptive | the other leg over time | what conventional training did |
| `BFR_Post-HLRT_Post` | descriptive | two legs after training | leg difference at the end |

`Modality_x_Time_Interaction` equals `BFR_Post-HLRT_Post` minus `BFR_Pre-HLRT_Pre` exactly, so these are not five
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

**The interaction has returned no hits under every variant tried.** That is a statement about what
this study could detect, not a demonstration that the two training modes act alike. `02_Differential`
prints the effect size it had 80% power to find, and the two numbers have to be quoted together.

`02_Differential/README.md` covers what comes out, how the null contrasts are read, and why
Benjamini-Hochberg runs within each contrast rather than across the five.
`03_Phenotype/README.md` covers the phenotype correlations and why classification is left to
`02_Differential`.
