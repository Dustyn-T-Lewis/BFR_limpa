# 02_Differential_Expression / 01_Design

**Planned.** This folder holds this README.

Builds the design matrix and the five comparisons, and checks the assumption they rest on. It is
separate from the fitting step so a design problem surfaces in seconds rather than after an hour
of modelling.

| | |
|---|---|
| **Will read** | the protein table from `01_Preprocess/02_Quantification` |
| **Will write** | the design and contrasts, plus a small diagnostic table |
| **Next** | `02_Differential` |

## The design

Four groups, one per treatment-and-timepoint combination, entered so each column is a group
average rather than a difference from a reference. Participant enters as a fixed term.

Fixed, not random, because every comparison the study asks for happens inside one person. With
participant in the design, treatments are only ever compared within the same person, which is what
makes pre-to-post paired. A random participant effect is for designs that also compare between
people, and this one does not.

## The assumption worth checking

Samples share a participant, and within that they share a leg. The participant term removes the
first exactly. The second stays behind in the noise. If it were large, the within-leg comparisons
would be tested too cautiously and the between-leg ones too generously.

This stage will measure that leftover correlation on the real data and report it, so the choice is
checked rather than assumed. If it turns out large, the fix is one argument, not a redesign.

## Two checks before anything is fitted

The design has to be full rank, or the comparisons are not estimable. And the group names have to
survive being read as R code, because the contrast function parses its arguments. A group named
like `2E-T1` would be read as subtraction.
