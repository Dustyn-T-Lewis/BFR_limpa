# 02 · Differential Expression

**Planned.** This folder holds READMEs. Pilot code exists offline and is not ready to read.

Takes `proteins.rds` and asks which proteins changed. Two sub-stages, so a design problem shows
up before an hour of model fitting does.

| Sub-stage | Will do |
|---|---|
| `01_Design` | the design matrix, five contrasts, and the assumption underneath them |
| `02_Differential` | fit, apply contrasts, adjust, write the tables |

## The five comparisons

| Name | Compares | Reads as |
|---|---|---|
| BFR post vs pre | one leg over time | what restricted training did |
| HLRT post vs pre | the other leg over time | what conventional training did |
| BFR vs HLRT at T1 | two legs before training | **control**, must find nothing |
| BFR vs HLRT at T2 | two legs after training | leg difference at the end |
| Interaction | difference of the two time effects | **the study question** |

## Three things that matter

**The control is the most important one.** Two untrained legs of the same person should differ in
nothing. If they do, the cause is upstream: a mislabelled sample, contamination, annotation. The
stage asserts it comes back empty and writes results to disk before that check runs.

**Participant goes in the design as a fixed term.** Every comparison happens inside one person,
so this makes pre-to-post paired. It also means sex cannot be tested, since it does not vary
within a participant.

**Adjusted within each comparison, never pooled.** The five share participants and the
interaction is built from two of the others.

Testing goes through limpa's own function, which reads the standard errors. A plain linear model
would discard them.
