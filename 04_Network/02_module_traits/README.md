# 02 · Module traits

Draws the module-trait heatmap of the YvO Figure 5 for each network from `01_modules`: each
module's contrasts and its associations with VL and RF cross-sectional area and echo intensity, on
one r-equivalent scale, sign(t) × √(t² / (t² + df)). Levels use arm as a covariate and participant
as a `duplicateCorrelation()` block; changes use the pooled and within-participant designs of
`03_Protein_Association`, which they reproduce for a single protein.

Inputs: the `membership`, `modules`, `eigengene_tests` and `change_tests` sheets of
`01_modules.xlsx`, `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds` and
`proteins.rds`, `02_Differential_Expression/02_Differential/c_data/fit.rds` and
`00_Input/phenotype.csv`.

```sh
quarto render 04_Network/02_module_traits/a_script/02_module_traits.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_module_traits.html`: the walkthrough.
- `b_reports/02_module_traits.pdf`: pages 1 to 4, one heatmap per network (main, change, pre,
  post), contrasts in purple and green, phenotype in red and blue, solid outline and stars at BH
  0.05, dashed at nominal 0.05, bold at |r| ≥ 0.40, module size bars labelled with top set and hub;
  page 5, main-network modules by training NES in each arm and by between-leg r before and after
  training.
- `c_data/02_module_traits.xlsx`:
  - `heatmap`: per cell, network, block, column, module, n, r-equivalent, t, df, p, BH q;
  - `module_nes`: per main-network module, size, and NES and BH q in each arm's training.

Three module levels track the phenotype at BH, all mitochondrial and all with RF muscle. Before
training, oxidative phosphorylation (main M3, hub SUCLA2) is higher where RF CSA is larger (r 0.35,
q 0.030). After training, ATP synthesis is higher where RF echo intensity is higher, in the main
network (M7, hub NDUFA10, r 0.37, q 0.019) and the post network (M6, hub NDUFS7, r 0.41, q 0.0056).
No module's change tracks the phenotype. Nominal hits run at or below chance (9 of 228 tests
against 11.4 expected), so these three are leads.
