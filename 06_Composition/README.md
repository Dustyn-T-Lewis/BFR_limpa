# 06 · Composition

A biopsy mixes slow (type I) and fast (type II) fibres, and type I fibres hold more mitochondria, so
a shift in fibre mix between biopsies could pass for a training effect. This stage estimates each
biopsy's slow and fast fibre fractions with FibeRtypeR (Van der Stede et al. 2026, J Physiol
604:2624), a web app that reads each sample as a mix of pure slow and pure fast reference profiles,
and checks whether the stage 02 results hold once composition is in the model.

Only the first step is built. `01_fibre_type_export` writes the upload file; its README covers the
upload and the 43 proteins the file lacks. Two steps wait on the app run: `02_fibre_type_import`
will read the app's fractions per biopsy, and a third will refit the stage 02 model with the slow
fraction as a covariate and compare the results.

The covariate is the FibeRtypeR slow fraction alone; the app returns slow and fast fractions only.
MPE from `05_Mitochondria` is never a covariate: it is a training outcome.
