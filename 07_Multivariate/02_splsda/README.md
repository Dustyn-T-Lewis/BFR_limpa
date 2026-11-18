# 02 · Multilevel sPLS-DA

Asks whether a weighted panel of proteins separates two groups in participants the model has not
seen: T1 against T2 in each leg (the positive control), and BFR against HLRT before training,
after training and on change, each with `multilevel = participant`. Each arm's training panel also
scores the other arm's legs, and the training panels are annotated by WGCNA module and
over-represented gene sets.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins_adjusted.rds`, the `membership` and
`modules` sheets of `04_Network/01_modules/c_data/01_modules.xlsx`, and
`03_Pathway_Enrichment/00_build_gene_sets/c_data/gene_sets.rds`.

```sh
quarto render 07_Multivariate/02_splsda/a_script/02_splsda.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_splsda.html`: the walkthrough.
- `b_reports/02_splsda.pdf`: one page per comparison (component scores, the ROC curve of each of
  the 50 cross-validation repeats, the 25 largest loadings per component shaded by stability), and
  a last page counting each training panel's proteins by main-network module.
- `c_data/02_splsda.xlsx` (`models` and `panels` are read by `08_Plasticity`):
  - `models`: per comparison and component, participants, samples, panel size, tuning BER,
    `perf()` BER under each distance, and cross-validated AUC with SD;
  - `panels`: per selected protein, comparison and component, gene, loading, stability across the
    250 refits, main-network module and its top set;
  - `transfer`: per training model, the arm it scores, transfer AUC and that arm's own CV AUC;
  - `panel_pathways`: per training panel and set at nominal p < 0.05, overlap, p and BH q;
  - `scores`: per sample and comparison, participant, group and the two component scores.

Training separates clearly: component 1 BER 0.006 and AUC 1.000 in the HLRT leg (20 proteins), 0.18
and 0.88 in the BFR leg (1,000 proteins). The HLRT panel scores BFR legs at AUC 0.91 and the BFR
panel HLRT legs at 0.98, though the scored legs belong to participants the panel was fitted on. The
BFR panel draws most from the oxidative phosphorylation module M3 (227 proteins) and the muscle
module M1 (128), with 383 in no module, and is enriched for oxidative phosphorylation,
mitochondrial protein degradation and the TCA cycle (BH q < 0.001); the HLRT panel draws 8
proteins each from M3 and the calcium transport module M4, and its top set, arginine and proline
metabolism, misses BH (q 0.058).

The arm comparisons sit near chance: BER 0.39 and AUC 0.67 before training (5 proteins), 0.44 and
0.60 after, 0.46 and 0.57 on change. Their panels pick the leads of stage 02 (PCYT1A, GPD2,
EIF2S1, FAHD2A and ARL8A on change). Tuning the panel size on the same data pulls BER below 0.5
even without signal, and the lowest BER falls before training, where stage 02 has no hit.
