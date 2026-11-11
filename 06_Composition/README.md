# 06 · Composition

A biopsy mixes slow (type I) and fast (type II) fibres, and type I fibres hold more mitochondria, so
a shift in fibre mix between biopsies could pass for a training effect. This stage estimates each
biopsy's slow fibre fraction with the FibeRtypeR method (Van der Stede et al. 2026, J Physiol
604:2624) and asks whether training or arm moved it.

[`02_fibre_type`](02_fibre_type/README.md) rebuilds the method in R from the paper's markers and
sex-specific reference and tests the slow fraction with the stage 02 design. Fibre type did not
change with training in either arm, so it enters no model; participant already absorbs most of it.

MPE from `05_Mitochondria` is never a covariate: it is a training outcome.
