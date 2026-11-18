# 01 · Modules

Groups proteins into signed WGCNA co-expression modules in four networks: all samples with
participant removed (main), each leg's training change (change), and the pre- and post-training
samples alone, all at soft power 8, `deepSplit = 2` and modules of at least 20 proteins. Names
each module by its BH top set and hub protein, tests main-network eigengenes
with the stage 02 contrasts and change-network modules for training, and checks whether
pre-training modules hold after training.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds`,
`02_Differential_Expression/01_Design/c_data/design.rds`, and `gene_sets.rds` and the
`protein_gene_map` sheet from `03_Pathway_Enrichment/00_build_gene_sets/`.

```sh
quarto render 04_Network/01_modules/a_script/01_modules.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/01_modules.html`: the walkthrough.
- `b_reports/01_modules.pdf`: scale-free fit by soft power per network (power 8 circled),
  then one protein dendrogram with module colours per network.
- `c_data/01_modules.xlsx`, read by `02_module_traits`, `07_Multivariate/02_splsda` and
  `08_Plasticity`:
  - `networks`: per network, samples, soft power, modules, share unassigned, median kME;
  - `soft_threshold`: per network and power, signed scale-free R² and mean connectivity;
  - `membership`: per network and protein, module (M0 unassigned) and own kME;
  - `modules`: per network and module, size, kME hub, connectivity hub, BH top set and its q;
  - `enrichment`: per network, module and set at nominal p < 0.05, overlap, p and BH q;
  - `eigengene_tests`: per main-network module and contrast, estimate, t, p, 95% CI, BH q;
  - `change_tests`: per change-network module and test (BFR training, HLRT training, BFR - HLRT),
    mean log2 change, t, p, 95% CI, BH q;
  - `preservation`: per pre-network module, size, Zsummary and median rank in post samples.

Training moves the same main-network modules in both arms, and no module differs between the
arms. Oxidative phosphorylation (M3, hub SUCLA2) rises with BFR (0.11, q 5 × 10⁻⁶) and HLRT
(0.09, q 1 × 10⁻⁴), and calcium ion transport (M4, hub MYLK2) falls with both (−0.09 and −0.13,
q 1 × 10⁻⁴ and 1 × 10⁻⁸); HLRT also raises ATP synthesis (M7, hub NDUFA10, q 0.022). No
between-arm contrast passes BH (smallest p 0.096, q ≥ 0.45). In the change network, organic acid
catabolism (M2, hub L2HGDH) and oxidative phosphorylation (M7, hub DLAT) rise with BFR (q 0.030
and 0.023); no module passes with HLRT (ATP synthesis M6, q 0.051) or between the arms (q ≥ 0.58).
All 8 pre-network modules are strongly preserved after training (Zsummary 13 to 55).

The networks hold 7 or 8 modules each and leave a large share unassigned (M0): 49% in main, 50%
in change, 39% in pre and 42% in post. Main-network modules are muscle cell differentiation (M1,
SYNPO2), cytoskeleton organisation (M2, AHNAK), M3, M4, haem metabolism (M5, SNCA), M6 (no set,
GSTP1) and M7. The kME hub matches the connectivity hub in 29% to 71% of modules. Power 8 follows
WGCNA guidance for signed networks; with these settings modules re-formed in 86% to 100% of
leave-one-participant-out refits (median Jaccard 0.62 to 0.76), a one-off check outside this step.
