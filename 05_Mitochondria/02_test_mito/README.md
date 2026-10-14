# 02 · Test mito

Measures mitochondrial protein enrichment (MPE) per biopsy and tests whether it rose differently
between legs, then refits the mitochondrial proteins on their own (mitochondrial normalisation) to
ask whether any protein, complex or MitoPathway changed beyond overall content.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds` and
`quant_runs/proteins_slope_0.7.rds`, `design.rds` from `01_Design`, `fit.rds` from
`02_Differential`, the protein map from `00_build_gene_sets`, and `01_build_mito_sets`.

```sh
quarto render 05_Mitochondria/02_test_mito/a_script/02_test_mito.qmd --output-dir ../b_reports
```

Outputs:

- `b_reports/02_test_mito.html`: the walkthrough.
- `b_reports/02_test_mito.pdf`: MPE by arm and time; complexes I to V by contrast for both fits,
  coloured by fry; volcanoes of the normalised proteins with fgsea rings (no baseline volcano).
- `c_data/02_test_mito.xlsx`: MPE per biopsy, the MPE test, citrate synthase protein, complex and
  MitoPathway tests, and the mitochondria-normalised protein table.

| Contrast | log2 MPE change | 95% CI | p |
|---|---:|---|---:|
| BFR post − pre | 0.152 | 0.066 to 0.237 | 0.0007 |
| HLRT post − pre | 0.155 | 0.070 to 0.239 | 0.0005 |
| BFR − HLRT, baseline | 0.031 | −0.054 to 0.115 | 0.47 |
| BFR − HLRT, post | 0.028 | −0.058 to 0.113 | 0.52 |
| interaction | −0.003 | −0.123 to 0.117 | 0.96 |

After mitochondrial normalisation, 17 proteins reach BH 0.05 in BFR and 28 in HLRT, none between
legs or at baseline. No complex reaches fry FDR 0.05 on the normalised fit. fry calls 18
MitoPathways in BFR, 36 in HLRT and none between legs; fgsea's between-leg calls (8 after training,
4 in the interaction) are leads, since it treats correlated proteins as independent.
