# 02 · Test mito

Measures mitochondrial protein enrichment (MPE) per biopsy and tests whether it rose differently
between legs, then refits the mitochondrial proteins alone (mitochondrial normalisation) to ask
whether any protein, complex or MitoPathway changed beyond overall content. The stage README
defines both.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds` and
`quant_runs/proteins_slope_0.7.rds`, `design.rds` from `01_Design`, `fit.rds` from
`02_Differential`, the protein map from `00_build_gene_sets`, and the `membership` sheet of
`01_build_mito_sets`.

```sh
quarto render 05_Mitochondria/02_test_mito/a_script/02_test_mito.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_test_mito.html`: the walkthrough.
- `b_reports/02_test_mito.pdf`: MPE by arm and time with each leg's pair joined; complexes I to V
  by contrast for both fits, labelled with fry FDR; volcanoes of the normalised proteins (no
  baseline) ringed by up to twelve MitoPathways at fgsea FDR 0.05.
- `c_data/02_test_mito.xlsx`:
  - `mpe`: per biopsy, MPE;
  - `mpe_test`: per contrast, log2 MPE change, plain least-squares t, p and 95% CI;
  - `citrate_synthase`: per contrast, citrate synthase protein from the stage 02 fit;
  - `complex_tests`: per fit, contrast and complex, fry (read by `03_complex_subunits`);
  - `pathway_tests`: per method, contrast and MitoPathway at 20% coverage, fry and fgsea on the
    normalised fit;
  - `mito_proteins`: per MitoCarta protein and contrast, the normalised dpcDE fit.

| Contrast | log2 MPE change | 95% CI | p |
|---|---:|---|---:|
| BFR post − pre | 0.155 | 0.073 to 0.237 | 0.0003 |
| HLRT post − pre | 0.138 | 0.056 to 0.220 | 0.0011 |
| BFR − HLRT, baseline | 0.017 | −0.064 to 0.099 | 0.68 |
| BFR − HLRT, post | 0.034 | −0.048 to 0.116 | 0.41 |
| interaction | 0.017 | −0.099 to 0.133 | 0.77 |

MPE rises with training by a similar amount in both legs. In the standard fit complexes II to V
rise with BFR and I and III to V with HLRT (fry FDR 0.005 to 0.031); relative to total
mitochondrial content no complex moves at fry FDR 0.05. Within the mitochondrion both legs shift
the same way: mitophagy, autophagy, selenoproteins and ROS and glutathione metabolism fall relative
to content (fry, 9 MitoPathways with BFR and 16 with HLRT), and HLRT also raises the mitochondrial
ribosome, translation and the TCA cycle. After normalisation 26 proteins reach BH 0.05 with BFR and
46 with HLRT, ALDH6A1 and PRKN falling in both, none between legs or at baseline. Between legs,
fgsea puts complex I lower in BFR relative to content after training (CI subunits NES −2.85; fry
FDR 0.055) and in the interaction (NES −1.84). These are leads: fry confirms no between-leg
MitoPathway, and fgsea also calls 4 at baseline, CI subunits among them (NES −2.04).
