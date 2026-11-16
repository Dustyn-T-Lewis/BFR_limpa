# 01 · Modules

Groups proteins into signed WGCNA co-expression modules in four networks: all samples with
participant removed (main), each leg's training change (change), and the pre- and post-training
samples alone. Names each module by its BH top set and hub protein, tests main-network eigengenes
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
- `b_reports/01_modules.pdf`: scale-free fit by soft power per network (chosen power circled),
  then one protein dendrogram with module colours per network.
- `c_data/01_modules.xlsx`, read by `02_module_traits`, `07_Multivariate/02_splsda` and
  `08_Plasticity`:
  - `networks`: per network, samples, soft power, modules, share unassigned, median kME;
  - `soft_threshold`: per network and power, signed scale-free R² and mean connectivity;
  - `membership`: per network and protein, module (M0 unassigned), module at power 12, own kME;
  - `modules`: per network and module, size, kME hub, connectivity hub, best Jaccard with a
    power-12 module, BH top set and its q;
  - `enrichment`: per network, module and set at nominal p < 0.05, overlap, p and BH q;
  - `eigengene_tests`: per main-network module and contrast, estimate, t, p, 95% CI, BH q;
  - `change_tests`: per change-network module and test (BFR training, HLRT training, BFR - HLRT),
    mean log2 change, t, p, 95% CI, BH q;
  - `preservation`: per pre-network module, size, Zsummary and median rank in post samples.

Training moves the same modules in both arms: 7 main-network modules with BFR and 10 with HLRT.
Aerobic respiration (M3, hub SUCLA2) and mitochondrial gene expression (M9, hub OGDH) rise in both,
and calcium ion transport (M8, hub TXLNB), M6 (hub CMBL) and M16 (hub PLEKHF1) fall in both; HLRT
also raises ATP synthesis (M4, hub NDUFA10). In the change network the TCA cycle (M3, hub IDH3B),
oxidative phosphorylation (M9, hub DLAT) and M22 (no set, hub ATP5MD) rise with BFR at BH; the
same three rise with HLRT by a little less (M3 0.15 against 0.19 log2) and miss BH (q 0.058 to
0.19). No module differs between the arms at BH; the one lead is main-network M22 (no set, hub
RBKS), which falls with HLRT (q 3.4 × 10⁻⁵) but not with BFR and sits higher in BFR after training
(p 0.0069, q 0.29). The networks hold 32 to 46 modules each, and all 32 pre-network modules are at
least moderately preserved after training.

Small modules depend on the chosen power: at power 12, 36% to 58% of proteins go unassigned, the
large main-network modules re-form (Jaccard 0.62 to 0.73) and most small ones do not (median best
Jaccard 0.05 to 0.13); the kME hub matches the connectivity hub in 61% to 79% of modules.
