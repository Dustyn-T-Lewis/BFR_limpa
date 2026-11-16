# 06 · Composition

A shift in fibre mix between biopsies could pass for a training effect, because type I fibres hold
more mitochondria. [`02_fibre_type`](02_fibre_type/README.md) estimates each biopsy's slow fibre
fraction with the FibeRtypeR method (Van der Stede et al. 2026, J Physiol 604:2624) and tests it
with the stage 02 design.

Results: fibre type did not change with training in either arm, so it enters no model. MPE from
`05_Mitochondria` is never a covariate: it is a training outcome.
