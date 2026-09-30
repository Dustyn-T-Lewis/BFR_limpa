# 06 · Mitochondria

Both legs gain mitochondrial protein with training. This stage asks two narrower questions. Did
one leg gain more mitochondria than the other? And once that overall gain is set aside, did any
part of the mitochondrion change on its own? The BFR literature disagrees on the first: Trinks
et al. (2026) report a larger citrate synthase rise after BFR in 10 people per arm, Groennebaek et
al. (2018) similar gains in 34.

Both methods come from the Bishop group. Mitochondrial protein enrichment (MPE) is the share of
raw protein intensity that comes from mitochondrial proteins (Granata et al. 2021, Nat Commun
12:7056). Mitochondrial normalisation takes the mitochondrial proteins out of the unnormalised
data, normalises them alone and refits, so each protein is read against total mitochondrial
content (Granata et al. 2021; Reisman et al. 2024, Nat Commun 15:7677). Three things differ from
those papers: we use whole muscle rather than isolated mitochondria, MitoCarta 3.0 rather than
IMPI, and our own normaliser (cyclic loess) rather than VSN.

`00_Input/` holds `Human.MitoCarta3.0.xls`, downloaded from the Broad Institute on 2026-09-28.

| Step | Question |
|---|---|
| `01_build_mito_sets` | Which MitoCarta proteins and MitoPathways do we measure? |
| `02_test_mito` | Did MPE rise differently between legs, and did any complex or pathway change beyond it? |

MPE rose from 29.6% to 32.8% in the BFR leg and from 29.0% to 32.2% in the HLRT leg. The
difference between legs is −0.003 log2 (95% CI −0.123 to 0.117), so a gap larger than about 8%
either way is unlikely. After mitochondrial normalisation no protein, complex or pathway differs
between legs by fry or protein-level FDR.

MPE from whole muscle mixes mitochondrial content with fibre type, since type I fibres hold more
mitochondria. `07_Composition` reads MPE from `02_test_mito` for that reason.
