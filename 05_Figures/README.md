# 05 · Figures

**Planned.** This folder holds this README. Panels are still moving while the first paper's scope
settles, so writing figure code now would mean rewriting it.

Every other stage keeps its diagnostics in its own `b_reports/`. This stage holds only panels
meant for someone who will never run the pipeline.

## What we would build

**The design.** Participants, two legs, two timepoints, treatment counterbalanced. A reader needs
this before any result means anything.

**What filtering removed.** The steps from the raw peptide table to the final one, with the share
of signal each took. This is where the search-database problem becomes visible.

**The detection curve.** The fitted curve beside the raw abundance-against-missingness cloud. It
shows missing values carry information, which justifies the method.

**Quantification uncertainty.** Standard error against abundance. It should slope down, showing
which proteins the analysis trusts.

**The two pre-to-post comparisons.** One panel each.

**The interaction.** If it is null, the panel needs effect sizes and intervals so a reader can see
how large a difference the study could have detected. A blank plot reads as a failed experiment,
which is a different claim.

**The control.** Two untrained legs finding nothing. Beside the interaction panel, it lets a
reader tell "no effect" from "no sensitivity".

Nothing from `03_` or `04_` until those stages exist and their results survive their own checks.

One script per panel, one image out. Figures are assembled outside the pipeline.
