# 05 · Mitochondria

Both legs gain mitochondrial protein with training. This stage asks two narrower questions. Did
one leg gain more mitochondria than the other? And once that overall gain is set aside, did any
part of the mitochondrion change on its own? In earlier BFR training studies, Groennebaek et al.
(2018, Front Physiol 9:1796) found similar rises in mitochondrial respiration and protein synthesis
after low-load BFR and high-load training, while citrate synthase activity did not change with
training compared with a control group. Trinks et al. (2026, Cell Metab 38:812) report increased
muscle mitochondrial content after BFR resistance training in people with type 2 diabetes.

Both methods come from the Bishop group. Mitochondrial protein enrichment (MPE) is the share of a
sample's summed protein values that comes from mitochondrial proteins (Granata et al. 2021, Nat
Commun 12:7056). Here it is the share that MitoCarta proteins take of the summed dpcQuant protein
values before normalisation, each protein counted once at its typical precursor intensity.
Mitochondrial normalisation takes the mitochondrial proteins out of the unnormalised data,
normalises them alone and refits, so each protein is read against total mitochondrial content
(Granata et al. 2021; Reisman et al. 2024, Nat Commun 15:7677).

Granata et al. used isolated mitochondria, marked mitochondrial proteins with IMPI, summed MaxQuant
LFQ intensities for MPE and normalised with VSN. We use whole-muscle lysate, MitoCarta 3.0, summed
dpcQuant protein values and cyclic loess. Reisman et al. are the precedent for whole-cell lysate
with MitoCarta 3.0; they used fibre pools where we use whole muscle.

`00_Input/` holds `Human.MitoCarta3.0.xls`, downloaded from the Broad Institute on 2026-09-28.

The steps run in order:

| Step | Question |
|---|---|
| [`01_build_mito_sets`](01_build_mito_sets/README.md) | Which MitoCarta proteins and MitoPathways do we measure? |
| [`02_test_mito`](02_test_mito/README.md) | Did MPE rise differently between legs, and did any complex or pathway change beyond it? |
| [`03_complex_subunits`](03_complex_subunits/README.md) | Which subunits and assembly factors carry each complex's change in each leg? |

Mean MPE rose from 29.7% to 32.9% in the BFR leg and from 29.1% to 32.3% in the HLRT leg. The
interaction, BFR's rise minus HLRT's, is −0.004 log2, and its 95% CI (−0.124 to 0.117) excludes
a difference in the MPE rise beyond about ±8%. After mitochondrial normalisation no protein,
complex or pathway differs between legs by fry or protein-level FDR.

MPE from whole muscle mixes mitochondrial content with fibre type, since type I fibres hold more
mitochondria. `06_Composition` estimates each biopsy's fibre mix for that reason.
