# 01 · Modules

Groups proteins into co-expression modules with WGCNA, four times over:

- main: all samples, participant removed with `removeBatchEffect()` (group protected);
- change: each leg's post- minus pre-training abundance;
- pre and post: one timepoint's samples, participant kept.

Each network is signed, with bicor and the soft power from `pickSoftThreshold()`. Modules come from
`blockwiseModules()` with at least 10 proteins and `deepSplit = 3`, finer than the tutorial's 30
and 2. Modules are enriched with `fgsea::fora()` against the project set list and named by their
top set when it passes BH, and by their hub protein. Main-network eigengenes go through the stage 02
design and contrasts. Change-network modules are tested on their mean member change: each arm
against zero, and BFR minus HLRT within participant. Least-squares t throughout, BH within each
contrast or test across modules. Preservation of the pre-network modules in the post-training
samples is checked with `modulePreservation()`.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`, and `gene_sets.rds` and the
`protein_gene_map` sheet from `03_Pathway_Enrichment/00_build_gene_sets/`.

```sh
quarto render 04_Network/01_modules/a_script/01_modules.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_modules.html`: the walkthrough.
- `b_reports/01_modules.pdf`: the scale-free fit by soft power for each network, with the chosen
  power circled, then one protein dendrogram with module colours per network.
- `c_data/01_modules.xlsx`, read by `02_module_traits`:
  - `networks`: one row per network: samples, soft power, modules, share of proteins unassigned
    and median kME of assigned proteins;
  - `soft_threshold`: one row per network and power: signed scale-free R² and mean connectivity;
  - `membership`: one row per network and protein: module (M0 is unassigned) and kME with its own
    module;
  - `modules`: one row per network and module: size, hub protein, and top set and its BH q where
    one passes;
  - `enrichment`: one row per network, module and set with nominal p < 0.05: overlap and
    hypergeometric p and q (BH over every set tested for the module);
  - `eigengene_tests`: one row per main-network module and contrast: estimate, t, p, 95% interval
    and BH q;
  - `change_tests`: one row per change-network module and test (BFR training, HLRT training,
    BFR - HLRT): mean log2 change, t, p, 95% interval and BH q;
  - `preservation`: one row per pre-network module: size, Zsummary (above 10 strong, 2 to 10
    moderate) and median rank in the post-training samples.

On current data the networks hold 40 to 50 modules each. Training moves 10 main-network modules in
BFR and 11 in HLRT; two change-network modules, the TCA cycle (M4, hub TUFM) and oxidative
phosphorylation (M10, hub DLAT), rise in both arms. No module differs between the arms. 38 of 40
pre-network modules are at least moderately preserved after training.
