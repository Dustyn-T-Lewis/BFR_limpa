# 02_Differential_Expression / 01_Design

**Planned.** This folder holds this README.

Builds the design matrix and the five contrasts, and checks the assumption they rest on. Separate
from the fit so a design problem surfaces in seconds.

## The design

Four groups, one per treatment-and-timepoint combination, entered so each column is a group
average. Participant enters as a fixed term.

Fixed, not random, because every comparison happens inside one person. With participant in the
design, treatments are only compared within the same person, which is what makes pre-to-post
paired. A random participant effect is for designs that also compare between people.

## The assumption to check

Samples share a participant, and within that they share a leg. The participant term removes the
first exactly; the second stays in the noise. This stage measures that leftover correlation on
the real data rather than assuming it is small. If it turns out large, the fix is one argument.

## Two checks before anything is fitted

The design must be full rank, or the comparisons are not estimable. And the group names must
survive being read as R code, because the contrast function parses its arguments.
