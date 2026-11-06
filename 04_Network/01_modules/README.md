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

Inputs: `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds`,
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
  - `membership`: one row per network and protein: module (M0 is unassigned), its module at
    power 12, and kME with its own module;
  - `modules`: one row per network and module: size, hub protein (highest kME), WGCNA's
    connectivity hub (`chooseTopHubInEachModule()`), best Jaccard overlap with a module of the
    same network rebuilt at power 12, and top set and its BH q where one passes;
  - `enrichment`: one row per network, module and set with nominal p < 0.05: overlap and
    hypergeometric p and q (BH over every set tested for the module);
  - `eigengene_tests`: one row per main-network module and contrast: estimate, t, p, 95% interval
    and BH q;
  - `change_tests`: one row per change-network module and test (BFR training, HLRT training,
    BFR - HLRT): mean log2 change, t, p, 95% interval and BH q;
  - `preservation`: one row per pre-network module: size, Zsummary (above 10 strong, 2 to 10
    moderate) and median rank in the post-training samples.

On current data the networks hold 32 to 46 modules each. Training moves 7 main-network modules in
BFR and 10 in HLRT; three change-network modules rise in BFR and none in HLRT: the TCA cycle (M3,
hub IDH3B), oxidative phosphorylation (M9, hub DLAT) and M22 (no set, hub ATP5MD). No module
differs between the arms. All 32 pre-network modules are at least moderately preserved after
training.

At power 12, the FAQ's table value, the networks thin: 36% to 58% of proteins go unassigned and
15 to 23 modules form. The large main-network modules re-form (M1, M2, M3, M4, M5 and M13, Jaccard
0.62 to 0.73); most small ones do not (median best Jaccard 0.05 to 0.13 per network), so the small
modules of the fine split depend on the chosen power. The highest-kME hub matches the connectivity
hub in 61% to 79% of modules.
