# 01 · Modules

Groups proteins into co-expression modules with WGCNA and tests each module's eigengene for
training, the BFR against HLRT contrasts and association with muscle growth. Participant is removed
first with `removeBatchEffect()` (group protected). The network is signed, with bicor and the soft
power from `pickSoftThreshold()`, and modules come from `blockwiseModules()` with the tutorial
settings. Eigengenes go through the stage 02 design and contrasts and the stage 02 association
designs (plain least-squares t, BH within contrast or outcome across modules). Modules are
enriched with `fgsea::fora()` against the project set list, and their preservation from pre- to
post-training samples is checked with `modulePreservation()`.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`,
`02_Differential_Expression/03_Protein_Association/c_data/association_designs.rds`, and
`gene_sets.rds` and the `protein_gene_map` sheet from `03_Pathway_Enrichment/00_build_gene_sets/`.

```sh
quarto render 04_Network/01_modules/a_script/01_modules.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_modules.html`: the walkthrough.
- `b_reports/01_modules.pdf`: the scale-free fit by soft power; the protein dendrogram with module
  colours; a heatmap of module t statistics across the five contrasts and eight association
  slopes, labelled with the BH-adjusted p where it is below 0.05.
- `c_data/01_modules.xlsx`:
  - `soft_threshold`: one row per power: signed scale-free R² and mean connectivity;
  - `membership`: one row per protein: module (M0 is unassigned) and kME with its own module;
  - `eigengene_tests`: one row per module and contrast: estimate, t, p, 95% interval and BH q;
  - `eigengene_association`: one row per module, analysis and outcome: the growth slope, t, p,
    interval and BH q;
  - `enrichment`: one row per module and set: overlap and hypergeometric p and q;
  - `preservation`: one row per module: size, Zsummary (above 10 strong, 2 to 10 moderate) and
    median rank, pre-training modules tested in post-training samples.

The module labels come from all samples, so the preservation check is mildly optimistic.
