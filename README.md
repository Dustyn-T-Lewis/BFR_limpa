# BFR Proteomics — limpa pipeline

Label-free DIA-MS proteomics of human vastus lateralis from a unilateral resistance training trial. Each participant trained one leg with blood-flow restriction (`BFR`) and the contralateral leg with conventional high load (`HL`), with biopsies before (`T1`) and after (`T2`) the intervention.

131 runs from 33 participants; 32 contribute all four cells. Leg assignment is counterbalanced — 16 of the 33 analysed participants trained the right leg under BFR and 17 the left — so treatment is never confounded with limb. `Treatment` is derived by joining the leg letter to the protocol sheet, never inferred from it — `D` and `E` are Portuguese *direita* and *esquerda*. The derivation was done once, by hand, and the result lives in `00_Input/metadata.csv`.

A parallel proteoDA implementation lives in `../BFR_proteoDA`, and `../BFR_comparison` documents where the two disagree and why.

## Design

Both factors are within-subject. Every participant is their own control twice over. Comparing legs cancels everything systemic, because both legs sat in one body under the same diet, sleep and hormonal state. Comparing timepoints within a leg cancels whatever was specific to that limb at baseline.

The fit is a means model over the four `Treatment_Timepoint` cells with `duplicateCorrelation` blocking on `Subject_ID`, consensus rho 0.244. The blocking is what earns the result: without it `Training_BFR` returns 1 protein instead of 27. Five contrasts, each a linear combination of the four cell means:

| Contrast | Definition | Reads as |
|------------------------|------------------------|------------------------|
| `Training_BFR` | `BFR_T2 - BFR_T1` | what restricted training did |
| `Training_HL` | `HL_T2 - HL_T1` | what conventional training did |
| `Baseline_BFRvHL` | `BFR_T1 - HL_T1` | **negative control** — two legs, one person, one day |
| `Post_BFRvHL` | `BFR_T2 - HL_T2` | leg difference after the intervention |
| `Interaction` | `(BFR_T2 - BFR_T1) - (HL_T2 - HL_T1)` | **primary** — what restriction adds |

`Baseline_BFRvHL` is a control, not a result. Signal there means a technical asymmetry — biopsy order, limb dominance, run order — and invalidates everything downstream. Every stage reads it first and carries a `stopifnot` that fires if it is not empty. Never weaken one of those to make a render pass.

Selection is Benjamini-Hochberg within each contrast. The five are never pooled, because they share participants and overlap by construction. The Π score (`p^|log2FC|`, Xiao et al. 2014, PMID 22321699) is reported beside BH and never selects: it is a transformed raw p-value that controls no error rate, and BH cannot repair that. On this data it inverts, needing \|log2FC\| above 0.43 even at p = 0.001 against a median effect of 0.09.

## Structure

One preprocessing stage produces one matrix and one design. Everything downstream reads both. The number in a filename is always the run order; a `.R` file is slow and cached, a `.qmd` renders, and numbering restarts inside each sub-stage. Figures live with the stage that produced them.

| Stage | Sub-stages | Produces |
|------------------------|------------------------|------------------------|
| `00_Input/` | — | the study data, hand-verified: three CSVs, nothing computed |
| `01_Preprocess/` | `01_Filtering` `02_Quantification` `03_Design` | one matrix, one design |
| `02_Differential_Expression/` | `01_Differential` | 3,035 proteins × 5 contrasts |
| `03_Pathway_Enrichment/` | `01_Gene_Sets` `02_Set_Tests` `03_Enrichment` | 2,031 set scores × 131 samples |
| `04_Network/` | `01_Modules` `02_Mechanism` | 7 WGCNA modules; STRING null; INDRA paths |
| `05_Figures/` | — | the manuscript panels — **not written yet** |

Each sub-stage carries the same three directories, and each holds exactly one kind of thing:

|  |  |
|------------------------------------|------------------------------------|
| `a_script/` | the code |
| `b_reports/` | everything rendered — the HTML, the figures it drew, the proteoDA PDFs |
| `c_data/` | every table and matrix the step produces |

A figure belongs with the report that draws it, so `b_reports/` holds both, and a notebook shows the same figure inline.

Each sub-stage writes one canonical file per grain, not a scatter. Nothing is written twice: there are no workbook copies of tables that already exist as CSV.

**A shared-helpers directory is planned and does not exist yet.** Identifier parsing, the five contrasts and the figure theme are currently inlined in the stage that needs them. When enough is duplicated to justify it, it goes in `00_Setup/functions/`; until then, do not write `source(here("00_Setup", ...))` into a script.

## Quantification

No imputation. limpa's `dpcQuant()` fits an additive model per protein and returns an intensity, a standard error and an observation count for every value, so a poorly observed protein is down-weighted rather than filled in. `dpcDE()` carries all three into limma.

The detection-probability curve is estimated on 46,430 precursors: slope 0.498. Two protein groups are capped at 500 precursors (titin at 3,245, and `P20929`), because `dpcQuant` cost grows as roughly n\^2.5; the cap uses a systematic sample stratified on detection and intensity, matching `dpc()`'s own `subset` device.

Contaminant removal is **by accession, never by gene symbol**. A symbol filter deletes 39.4% of signal here, because cRAP's rabbit ALDOA collides with muscle ALDOA and bovine P00366 sits one digit from human GLUD1 P00367. Separately, 43 of 95 compound protein groups are one accession appearing with and without a `cRAP-` prefix; both limpa filters would delete those as ambiguous, taking about 59% of signal. They are resolved before either filter runs. **Do not reorder that.**

The blood threshold is 0.50, not 0.40. At 0.40 the filter takes PRKAA1 (AMPK) plus HINT2, MRPL21 and PPA2 — all mitochondrial — out of a study whose headline result is mitochondrial.

## Results as they stand

| Contrast          | Proteins (BH \< 0.05) | Pathway sets              | Modules |
|-------------------|-----------------------|---------------------------|---------|
| `Training_BFR`    | 27                    | 9 by over-representation  | 2       |
| `Training_HL`     | 140                   | 23 by over-representation | 2       |
| `Baseline_BFRvHL` | **0**                 | 0                         | 0       |
| `Post_BFRvHL`     | 0                     | 0                         | 0       |
| `Interaction`     | **0**                 | 0                         | 0       |

The negative control is clean at every level. The interaction is null: this design did not detect a proteome-wide difference between what restriction added and what conventional loading added. A power calculation carried out alongside this pipeline, but reproduced by no script in it, put the minimum detectable effect at 0.434 log2 (1.35×) at nominal alpha and 0.811 log2 (1.75×) at BH's effective bar. On that basis "BFR and high load are indistinguishable here" is supported and "identical" is not.

The `Training_BFR` enrichment is mitochondrial, and the drivers say how much so. Nine significant sets rest on roughly five enzymes — `IDH3B`, `IDH3G`, `SUCLA2`, `ME2`, `SHMT2` — appearing under TCA cycle, NAD binding, aerobic respiration and magnesium binding. Five enzymes wearing nine labels, not nine independent findings. The same happens in `Training_HL`, where seven T-cell sets are all the same five proteasome subunits: PSMA3, PSMA4, PSMA5, PSMB1 and PSMB2. Every enrichment row carries its driver genes for this reason.

## What is not tested

No competitive gene-set test runs. fgsea's permutation null assumes genes vary independently, which muscle mitochondrial and ribosomal sets do not, and it called 13 sets significant on the negative control. camera handles inter-gene correlation but accepts no `block` or `correlation` argument — limma drops both into `...` with a warning — so it would treat 131 samples as independent when 33 participants contribute four each. What runs is fry (self-contained, blocks correctly), singscore, and over-representation on the BH-selected proteins.

## Two results that are null, and stay reported

**Protein interaction enrichment.** An earlier pass reported 9 STRING edges against 2 expected, p = 0.00097, using STRING's whole-proteome background. Permuting against the proteins actually tested gives 9 observed against 8.41 expected, **p = 0.46**. The effect was the background. All four sets tested come back null, p 0.42 to 0.91.

**INDRA.** Every one of the 351 protein pairs submitted returns a literature path — 351 of 351. Without a random-pair control that is not a finding. Two direct edges rest on more than one publication: `IDH3B → IDH3G` (belief 0.975) and `SLC25A3 → FSCN1` (0.699). The first says two subunits of one enzyme are related. Kept as a documented failure.

## limpa against proteoDA

The two pipelines disagree on `Training_BFR` (27 here, 97 there). Four tests ruled out every tunable cause: the BH denominator, the filtering, the DPC slope (forcing 0.7 gives 28, forcing 0.9 gives 23), and sparsity. Per-protein correlation between them is 0.805, limpa's within-subject residual SD is 9% wider, and the proteoDA matrix is still 7.2% missing.

They disagree because they **measure** differently, not because either is mistuned. See `../BFR_comparison/comparison.html`.

## Running it

Open `BFR_limpa.Rproj` so the working directory is the project root and every `here::here()` path resolves. Then run the stages in order. A `.qmd` renders to its sibling `b_reports/`; a `.R` file is slow and cached, so run it before the `.qmd` that reads what it wrote.

```         
00_Input runs nothing.

quarto render 01_Preprocess/01_Filtering/a_script/01_filter.qmd --output-dir ../b_reports
Rscript      01_Preprocess/02_Quantification/a_script/run_dpcquant.R
quarto render 01_Preprocess/02_Quantification/a_script/02_quantify.qmd --output-dir ../b_reports
quarto render 01_Preprocess/03_Design/a_script/03_design.qmd --output-dir ../b_reports

quarto render 02_Differential_Expression/01_Differential/a_script/01_differential.qmd --output-dir ../b_reports

Rscript      03_Pathway_Enrichment/01_Gene_Sets/a_script/01_build_gene_sets.R
quarto render 03_Pathway_Enrichment/02_Set_Tests/a_script/02_set_tests.qmd --output-dir ../b_reports
quarto render 03_Pathway_Enrichment/03_Enrichment/a_script/03_enrichment.qmd --output-dir ../b_reports

quarto render 04_Network/01_Modules/a_script/01_modules.qmd --output-dir ../b_reports
Rscript      04_Network/02_Mechanism/a_script/01_run_ppi_null.R
Rscript      04_Network/02_Mechanism/a_script/02_run_indra.R
quarto render 04_Network/02_Mechanism/a_script/03_mechanism.qmd --output-dir ../b_reports
```

The `--output-dir` flag is the only thing putting reports in `b_reports/`; no `.qmd` sets it in its own YAML.

About 35 minutes end to end. The four cached `.R` steps are 23 of the 35: `dpcQuant` 9, the gene-set build 5, the STRING permutation 4, INDRA 5.

## Contributing

Branch before committing: `feature/…`, `fix/…`, `docs/…`. `main` requires a pull request.

`.gitignore` excludes `b_reports/` and `c_data/` — outputs regenerate from code, and four people committing their own render of the same HTML makes every pull a conflict in a file nobody reads. One tracked exception is deliberate: `03_Pathway_Enrichment/01_Gene_Sets/c_data/` pins set membership against an msigdbr release, and rebuilding it silently changes every enrichment number downstream. That exception is a literal path in `.gitignore` and must be updated if the directory ever moves again.

There is no `renv.lock`, by choice, so package versions are unpinned and msigdbr drift is a live risk the gene-set cache only documents.

## Licence

MIT. Copyright (c) 2026 Dustyn T. Lewis.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.