# 05 · Mitochondria

Asks whether one leg gained more mitochondrial protein than the other, and whether any part of the
mitochondrion changed beyond that overall gain. Low-load BFR and high-load training raised
mitochondrial respiration and protein synthesis alike in Groennebaek et al. (2018, Front Physiol
9:1796), and BFR training raised mitochondrial content in type 2 diabetes in Trinks et al. (2026,
Cell Metab 38:812).

Both methods follow Granata et al. (2021, Nat Commun 12:7056) and Reisman et al. (2024, Nat Commun
15:7677): mitochondrial protein enrichment (MPE), the share of summed protein values from MitoCarta
proteins, and mitochondrial normalisation, which normalises the mitochondrial proteins alone and
refits. Granata et al. used isolated mitochondria, IMPI, MaxQuant LFQ and VSN; here it is whole
muscle, MitoCarta 3.0, dpcQuant protein values before normalisation and cyclic loess. Reisman et
al. are the precedent for lysate with MitoCarta 3.0, on fibre pools.

`00_Input/` holds `Human.MitoCarta3.0.xls`, downloaded from the Broad Institute on 2026-09-28.

- [`01_build_mito_sets`](01_build_mito_sets/README.md): which MitoCarta genes and MitoPathways we
  measure.
- [`02_test_mito`](02_test_mito/README.md): whether MPE rose differently between legs, and whether
  any complex or pathway changed beyond it.
- [`03_complex_subunits`](03_complex_subunits/README.md): which subunits and assembly factors
  carry each complex's change.

Mean MPE rose from 29.7% to 32.9% in the BFR leg and from 29.1% to 32.3% in the HLRT leg; the
interaction is 0.017 log2 (95% CI −0.099 to 0.133). Relative to that content, both legs shift
away from mitophagy, autophagy and redox proteins. After mitochondrial normalisation no protein,
complex or pathway differs between legs by fry or protein-level FDR; the one lead is complex I,
lower in BFR relative to content after training (fgsea; fry FDR 0.055). Whole-muscle MPE mixes
mitochondrial content with fibre type, which `06_Composition` estimates.
