# 04 · Network

**Planned.** This folder holds this README.

`03_` asks whether proteins somebody already grouped moved together. This stage asks which
proteins move together in our data, whether or not anyone has grouped them before.

| Sub-stage | Will do |
|---|---|
| `01_Modules` | cluster proteins by how they covary, summarise each cluster, test the summaries |
| `02_Mechanism` | check clusters against a protein-interaction database, then the literature |

## What each would tell us

**Modules** are groups this muscle regulates together, discovered rather than assumed. A module
that moves on the interaction would be the strongest result this dataset could produce.

One decision governs whether it works: modules must be built on abundance centred within each
participant. Built on raw abundance they encode who the person is, because differences between
people dominate the variation.

**The database check** asks whether a module's proteins interact more than chance allows, tested
against a null built by shuffling module membership rather than the database's own background.
These databases are denser around well-studied proteins, so testing against the raw background
rewards famous proteins regardless of biology. Pilot work made that mistake and had to withdraw
the claim.

A module that is not enriched is not refuted. It may be proteins nobody has studied together.

**The literature search** produces hypotheses, never evidence. Ask it for paths between random
protein pairs and it finds those too, so without that control it says nothing.
