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

Input: `00_Input/` holds `Human.MitoCarta3.0.xls`, downloaded from the Broad Institute on
2026-09-28.

- [`01_build_mito_sets`](01_build_mito_sets/README.md): which MitoCarta genes and MitoPathways we
  measure.
- [`02_test_mito`](02_test_mito/README.md): whether MPE rose differently between legs, and whether
  any complex or pathway changed beyond it.
- [`03_complex_subunits`](03_complex_subunits/README.md): which subunits and assembly factors
  carry each complex's change.

Each step README gives its inputs, outputs and results. Whole-muscle MPE mixes mitochondrial
content with fibre type, which `06_Composition` estimates.
