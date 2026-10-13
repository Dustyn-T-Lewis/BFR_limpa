# 02_test_mito

This step measures mitochondrial content per biopsy and tests the mitochondrial proteome after
setting that content aside. It reads the unnormalised quantification
(`01_Preprocess/02_Quantification/c_data/quant_runs/proteins_slope_0.7.rds`, cut to the 3,042 rows
of `proteins.rds`), the design, the stage 02 fit and `01_build_mito_sets`. It writes
`c_data/02_test_mito.xlsx` and `b_reports/02_test_mito.pdf`.

```sh
Rscript 05_Mitochondria/02_test_mito/a_script/02_test_mito.R
```

MPE comes from raw intensities, as in both papers. Cyclic loess assumes the typical protein holds
still, and a mitochondrial expansion breaks that. The script checks that the raw matrix becomes
`proteins.rds` after loess and nothing else. log2 MPE goes through the saved design and contrasts
as one outcome row, so there is nothing for eBayes to borrow from and the test is the plain
least-squares t on 95 residual degrees of freedom.

| Contrast | log2 MPE change | 95% CI | p |
|---|---:|---|---:|
| BFR post − pre | 0.152 | 0.066 to 0.237 | 0.0007 |
| HLRT post − pre | 0.155 | 0.070 to 0.239 | 0.0005 |
| BFR − HLRT, baseline | 0.031 | −0.054 to 0.115 | 0.47 |
| BFR − HLRT, post | 0.028 | −0.058 to 0.113 | 0.52 |
| interaction | −0.003 | −0.123 to 0.117 | 0.96 |

Citrate synthase protein, read from the stage 02 fit for comparison with activity studies, rose
0.22 log2 in BFR (BH 0.076) and 0.18 in HLRT (BH 0.17).

Mitochondrial normalisation subsets the 660 MitoCarta proteins, normalises them alone with cyclic
loess and refits with `dpcDE`, `eBayes(robust = TRUE)` and the same contrasts. Protein hits at BH
0.05: 17 in BFR, 28 in HLRT, none between legs and none at baseline.

fry tests complexes I to V on both fits. On the standard fit, complexes III, IV and V rose in both
legs (FDR 0.0017 to 0.023) and complex I in HLRT (0.004); complex I in BFR and complex II in both
legs sat between 0.07 and 0.1. On the normalised fit none reaches FDR 0.05; the nearest is complex
I lower in BFR than HLRT after training (0.065).

MitoPathways are tested when at least 20% of their listed genes are measured, 141 of 149. There is
no size floor, since MitoCarta lists some functions that are one or two proteins. The 20% floor is
ours, not a published rule. fry is the confirmatory test: after normalisation it already reads
change relative to the rest of the mitochondrion, so a competitive test would repeat the question.
fry calls 18 pathways in BFR, 36 in HLRT and none between legs. fgsea draws the rings and calls 8
pathways after training and 4 in the interaction, led by complex I down and fatty-acid oxidation up
in BFR. fgsea treats proteins as independent, and complex I subunits correlate, so those are leads
that fry does not confirm.

The PDF has three pages: MPE by arm and time with each leg's pair joined, a complex × contrast grid
for both fits coloured by fry, and volcanoes of the normalised proteins with fgsea rings. The
baseline contrast is not drawn as a volcano.
