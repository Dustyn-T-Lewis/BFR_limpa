# 04 · Network

**Planned.** This folder holds this README and nothing else.

`03_Pathway_Enrichment` asks whether proteins somebody already grouped moved together. This stage
asks the opposite: which proteins move together in our data, whether or not anyone has grouped
them before. The two can disagree, and the disagreement is the interesting part.

```
02_Differential_Expression results + the protein table
  01_Modules    find groups the data defines   -> module summaries and their contrasts
  02_Mechanism  ask whether anyone agrees      -> interaction test and literature paths
```

| Sub-stage | Will do |
|---|---|
| `01_Modules` | cluster proteins by how they covary across samples, summarise each cluster, test the summaries |
| `02_Mechanism` | check the clusters against a protein-interaction database, then search the literature for mechanisms |

## What each part would tell us

**Modules** are groups of proteins this muscle regulates together, discovered rather than assumed.
Each becomes a single summary value per sample, which can then be tested like a protein. A module
that moves on the interaction would be the strongest result this dataset could produce: a
coordinated group, defined by our own data, responding differently to restriction than to heavy
load.

One decision governs whether it works. Modules have to be built on abundance centred within each
participant. Built on raw abundance they encode who the person is rather than what training did,
because differences between people dominate the variation. Get that backwards and every module is
a participant fingerprint.

**The interaction-database check** asks whether the proteins in a module interact more often than
chance allows. It has to be tested against a null built by shuffling module membership, not
against the database's own background, because these databases are denser around well-studied
proteins. Test against the raw background and any module full of famous proteins looks enriched
regardless of biology. Pilot work made exactly that mistake and had to withdraw the claim.

A module that is enriched is one where covariation we measured lines up with coupling somebody
else established independently. A module that is not enriched is not refuted. It may be proteins
nobody has studied together, which in a discovery experiment is the point.

**The literature search** produces hypotheses, never evidence. It suggests why two proteins our
data associates might be connected. Its output belongs in a discussion as a lead for the next
experiment, and should never be counted or reported as an enrichment. Ask it for paths between
random protein pairs and it will find those too, so without that control it says nothing.

## What would make this convincing

A module that survives all three views: it moves on the interaction, its members interact more
than shuffling allows, and the literature offers a mechanism. That is enough to design a
follow-up experiment on.

Expect at least one of the three to come back null, and report it that way.

## Packages

WGCNA for the clustering, a protein-interaction source such as STRING, and a literature source
such as INDRA. Both external services, so both need their query date recorded.

## Cost

Around twenty minutes, mostly the shuffling test. Worth caching.
