# 05 · Figures

**Planned.** This folder holds this README and nothing else. The panels are still moving while
the scope of the first paper settles, so writing figure code now would mean rewriting it.

Every other stage keeps its own diagnostic plots in its own `b_reports/`. This stage holds only
panels meant for someone who will never run the pipeline.

## What we would build

**The design.** Participants, two legs each, two timepoints, treatment counterbalanced. A reader
needs this before any result means anything, and it is the one panel that cannot go out of date,
because the design is fixed.

**What filtering removed.** The steps from the raw peptide table to the final one, with the share
of signal each step took. This is where the search-database problem becomes visible, and it is
where a reviewer who suspects contaminant removal ate real muscle can check.

**The detection curve.** The fitted curve beside the raw relationship between abundance and
missingness. This panel justifies the method: it shows missing values carry information, which is
why the analysis models them instead of filling them in.

**Quantification uncertainty.** Standard error against abundance. It should slope down. This shows
the uncertainty is tracking how much evidence there is, rather than being decorative, and it tells
a reader which proteins the analysis trusts.

**The two pre-to-post comparisons.** One panel each, showing that training changed the muscle
proteome in each leg.

**The interaction.** Whatever it turns out to be. If it is null, the panel has to show effect
sizes and their intervals, so a reader can see how large a difference the study could have
detected. A blank plot reads as a failed experiment, which is a different claim.

**The control.** Two untrained legs of one person, finding nothing. Placing it next to the
interaction panel is what lets a reader tell "no effect" apart from "no sensitivity".

## What we would not draw

Anything from `03_` or `04_` until those stages exist and their results have survived their own
checks. A pathway figure built on an untested set list is the fastest route to a panel that has to
be withdrawn.

## Practical points

One script per panel, one image out. Panels are assembled into figures outside the pipeline, so
nothing here tries to lay out a multi-panel figure in code.

A shared theme and colour scheme for the four study groups will live with these scripts, so every
panel matches without each one redefining it.

## Cost

Cheap, once the upstream tables are final.
