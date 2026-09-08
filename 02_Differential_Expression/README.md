# 02 · Differential Expression

**Planned.** This folder holds this README. Pilot code exists offline and is not ready to be
read or relied on.

Takes the protein table from `01_Preprocess` and asks which proteins changed. Two sub-stages, so
that a problem with the design shows up before an hour of model fitting does.

```
01_Preprocess/02_Quantification/c_data/proteins.rds
  01_Design       build the design and the contrasts   -> design.rds
  02_Differential fit the model, write the results     -> one table per protein per contrast
```

| Sub-stage | Will do |
|---|---|
| `01_Design` | write the design matrix, define the five comparisons, check the assumption underneath them |
| `02_Differential` | fit, apply the contrasts, adjust for multiple testing, write the tables |

## The five comparisons

| Name | Compares | Reads as |
|---|---|---|
| BFR post vs pre | one leg over time | what restricted training did |
| HLRT post vs pre | the other leg over time | what conventional training did |
| BFR vs HLRT at T1 | two legs before training | **control**, must find nothing |
| BFR vs HLRT at T2 | two legs after training | leg difference at the end |
| Interaction | the difference of the two time effects | **the study question** |

The interaction is what the trial was designed to answer. The rest are there to make it
interpretable.

## The control is the most important one

Comparing two untrained legs of the same person on the same day should find nothing. If it finds
something, the cause is upstream: a mislabelled sample, contamination, a problem in the
annotation. The stage will assert that it comes back empty and write the results to disk before
that check runs, so a failure leaves the evidence behind rather than losing it.

## Why participant goes in the design

Every comparison the study asks for happens inside one person. Putting participant in the design
as a fixed term means nobody's before is ever compared against somebody else's after, which is
what makes the pre-to-post tests paired. The alternative, treating participant as a random
effect, is for designs that also compare between people, and this one does not.

One consequence worth stating: sex does not vary within a participant, so this model absorbs it
and cannot test it. A different model would be needed for that question.

## Multiple testing

Adjusted within each comparison, never pooled across the five. They share participants and the
interaction is built from two of the others, so pooling would claim an independence the design
does not have.

## What has to be true

Differential testing goes through limpa's own function, which reads the standard errors from
`proteins.rds` and turns them into weights. Handing the abundances to a plain linear model would
discard the uncertainty and treat every protein as equally reliable, which is the behaviour limpa
exists to replace.

## Downstream

`03_Pathway_Enrichment` and `04_Network` both read the protein tables this stage writes.
`05_Figures` reads them too.

## Packages

`limpa` for the fit, `limma` for the contrasts, multiple-testing adjustment and result tables.
