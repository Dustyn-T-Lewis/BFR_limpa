# 02 · Module traits

The module-trait heatmap of the YvO Figure 5, one page per network from `01_modules`. Each module's
contrasts and its associations with VL and RF cross-sectional area and echo intensity share one
scale, the r-equivalent sign(t) × √(t² / (t² + df)).

Associations use both arms in one model per module and outcome. Levels (main network before and
after training, and the pre and post networks) regress the eigengene on the outcome with arm as a
covariate and participant as a block (`duplicateCorrelation()`). Changes regress each leg's
eigengene change on its outcome change, pooled (one intercept per arm) and within participant (one
intercept per participant); these match the slope t of `03_Protein_Association`'s designs.
Main-network eigengenes are recomputed on raw abundance so levels keep between-person differences.
Least-squares t, BH within each network, block and outcome across modules.

Inputs: the `membership`, `modules`, `eigengene_tests` and `change_tests` sheets of
`01_modules.xlsx`, `01_Preprocess/02_Quantification/c_data/proteins.rds`,
`02_Differential_Expression/02_Differential/c_data/fit.rds` and `00_Input/phenotype.csv`.

```sh
quarto render 04_Network/02_module_traits/a_script/02_module_traits.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_module_traits.html`: the walkthrough.
- `b_reports/02_module_traits.pdf`, A4 portrait:
  - pages 1 to 4, one per network (main, change, pre, post): the heatmap. Purple and green cells
    are contrasts, red and blue cells phenotype associations. A solid outline with stars marks BH
    p < 0.05, a dashed outline nominal p < 0.05, bold text |r| ≥ 0.40. Strips give n. The left bars
    give module size on a square-root scale, labelled with the BH top set ("no set" when none
    passes) and the hub protein.
  - page 5, main network: modules placed by fgsea NES in the two training contrasts (module
    proteins as sets, stage 02 moderated t as ranks), and by eigengene r-equivalent, BFR against
    HLRT, before and after training. Point size is module size.
- `c_data/02_module_traits.xlsx`:
  - `heatmap`: one row per cell: network, block, column, module, n, r-equivalent, t, df, p and BH
    q;
  - `module_nes`: one row per main-network module: size, and NES and BH q in each arm's training.

On current data, 7 associations pass BH, all with post-training echo intensity: in the main network
M5 (hub NDUFB9, a complex I subunit), M17 (TOLLIP) and M34 (SDHB) with RF echo; in the post network M24
(adipogenesis, hub AOC3) with VL echo, r = 0.49, M34 (glycogen metabolism) with VL echo, and M3
(complex I biogenesis) and M38 (oxidative phosphorylation) with RF echo. None involve CSA or change.
Nominal hits run at or below chance in the main and change networks and about twice chance in the
pre and post networks.
