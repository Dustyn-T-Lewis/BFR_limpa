# 00 · Input

The study data, consolidated and verified by hand. Three files: one row per precursor, one
row per sample, one row per leg. No stage downstream writes back into this directory.

| File | Rows × cols | What it is |
|------------------|------------------|--------------------------------------|
| `precursors.csv` | 53,852 × 141 | 10 DIA-NN annotation columns, then the 131 sample columns |
| `metadata.csv` | 131 × 8 | one row per sample column of `precursors.csv`, same order |
| `phenotype.csv` | 70 × 14 | one row per leg — muscle size, echo intensity, training volume |

The annotation is not a separate file because it is not separate data: columns 1–10 describe
the precursor in the same row the intensities sit on. proteoDA's own example dataset,
`Lou_HF_DIANN_uni_prot_quan`, ships in exactly this shape — annotation columns first, then
sample columns — and `DAList()` performs the split in memory:

``` r
precursors <- read_csv(here("00_Input", "precursors.csv"))
metadata <- read_csv(here("00_Input", "metadata.csv"))

annotation <- precursors[, 1:10]
intensity <- as.matrix(precursors[, -(1:10)])
```

Splitting on disk would also force a `uniprot_id` key, which these rows cannot honestly
supply: they are precursors, not proteins. The protein-level split is available downstream
instead, where the key really is a protein — `dpcQuant()` in `01_Preprocess/02_Quantification` returns an
`EList` whose `$E` is the 3,758-row protein matrix, `$genes` its annotation and `$targets`
the sample sheet. (`dpcDE()` is one step further on and returns a limma `MArrayLM` fit, not
a matrix.)

## precursors.csv

Columns 1–10 describe the precursor: `Protein.Group`, `Protein.Ids`, `Protein.Names`,
`Genes`, `First.Protein.Description`, `Proteotypic`, `Stripped.Sequence`,
`Modified.Sequence`, `Precursor.Charge`, `Precursor.Id`. A precursor is one peptide at one
charge state. The DIA-NN names are kept as they came off the report so a field can be traced
back to it.

Columns 11–141 are the 131 samples, named `S02_L_T1` and so on, in the same order as the rows
of `metadata.csv`. Values are intensities on the **linear** scale, empty where the precursor
was not detected.

53,852 precursors over 3,758 protein groups and 37,475 distinct peptide sequences, charges
1–4. Median 6 precursors per protein group; the largest holds 3,248. **51.3% of cells are
empty** — 3,620,948 of 7,054,612 — which is ordinary for DIA and is what limpa's detection
probability curve models rather than imputes. There are no zeros: non-detection is always
empty. 1,983 precursors are detected in all 131 samples.

## metadata.csv

131 rows, one per sample column of `precursors.csv` and in the same order.

| Column | Values | Why it is here |
|------------------------|------------------------|------------------------|
| `sample_id` | `S02_L_T1` | the key — matches the sample columns of `precursors.csv` |
| `participant` | `S02` … `S37`, 33 total | the **block** |
| `leg` | `L` or `R` | the unit randomised to a treatment |
| `treatment` | `BFR` or `HLRT` | |
| `timepoint` | `T1` or `T2` | |
| `group` | `BFR_T1` `BFR_T2` `HLRT_T1` `HLRT_T2` | the **factor the model is fitted on** |
| `sex` | `F` or `M` | 14 F, 19 M, constant within participant |
| `run` | the original instrument filename | provenance; carries the reinjection suffix |

`group` is `treatment` and `timepoint` pasted together; all three were written and checked
against each other when the file was built.

## phenotype.csv

70 rows, one per participant leg — the unit randomised, and the grain the measurements were
taken at. Keyed by `participant` + `leg`, which is how it joins to `metadata.csv`.

| Column | Units | What it is |
|------------------------|------------------------|------------------------|
| `participant`, `leg` | | the join key |
| `treatment`, `sex` | | duplicated from `metadata.csv` so the file reads alone |
| `vl_csa_pre_cm2`, `vl_csa_post_cm2` | cm² | vastus lateralis cross-sectional area — muscle size |
| `vl_echo_pre_au`, `vl_echo_post_au` | A.U. | vastus lateralis echo intensity — muscle quality proxy |
| `rf_csa_pre_cm2`, `rf_csa_post_cm2` | cm² | rectus femoris cross-sectional area |
| `rf_echo_pre_au`, `rf_echo_post_au` | A.U. | rectus femoris echo intensity |
| `sessions_completed` | count, of 20 | sessions with a recorded volume load |
| `volume_load_total_kg` | kg | those sessions summed |

`pre` and `post` are `T1` and `T2`. Keeping this at leg grain rather than folding it into
`metadata.csv` is deliberate:

**It holds all 140 measurements.** 70 legs at two timepoints. A sample-grain copy fits only
131 and drops nine — S06's left leg at T2, whose sample was never acquired, and all eight for
S22 and S34, who trained and were scanned but contribute no MS data.

**It cannot be double-counted.** `volume_load_total_kg` is one number per leg. At sample grain
it appears twice per leg, and `sum()` over rows quietly returns twice the training volume.

**Higher echo intensity means worse muscle quality.** It rises with intramuscular fat and
fibrous tissue, so it moves opposite to CSA. Do not read the two as pointing the same way.

**`sessions_completed` is what makes the volume total readable.** S29's HLRT leg completed 11
sessions and S05's 19, so those two totals are not comparable with the other 68 at face value.

**Change scores are not stored.** Δ% and Δ-absolute both follow from the pre/post pair, and a
stored delta fixes a scale choice the analysis should be making itself.

**Pennation angle, fascicle length and per-session detail are not here.** They were in
`Auburn_data.xlsx`, which is no longer in this directory. `phenotype.csv` is now the record
of what was measured, and cannot be rebuilt from anything present.

**Participant names were dropped.** The workbook carried a full-name column; `participant` is
the key everything joins on.

## Design

Two factors, **both within-subject**:

-   **treatment** — each participant trained one leg with BFR and the other with HLRT, so BFR and HLRT are compared *inside* one body, under the same diet, sleep and hormonal state.
-   **timepoint** — each leg was sampled before (`T1`) and after (`T2`) training, so pre/post is compared *inside* one leg.

Every participant is therefore their own control twice over, and no between-subject variation enters any contrast. Assignment is counterbalanced: 17 participants trained BFR on the left, 16 on the right. `leg` is not a term in the model, and **treatment can never be read off the leg letter.**

## Model

Defined and asserted in `01_Preprocess/03_Design`, which is the authority. The short version,
because this directory supplies the labels it is built from.

The four `group` levels are the four cells of a 2 x 2, and the participant term is **fixed**:

``` r
design <- model.matrix(~ 0 + group + participant, data = targets)
```

$$y_{gj} = \mathbf{x}_j^{\top}\boldsymbol{\beta}_g + \varepsilon_{gj},
\qquad \varepsilon_{gj} \sim N(0,\ \sigma^2_g / w_{gj})$$

for protein $g$ in sample $j$, where $\boldsymbol{\beta}_g$ holds the four cell means
followed by 32 participant deviations, and $w_{gj}$ is the precision weight limpa derives from
each value's standard error. The matrix is 131 x 36, full rank, leaving 95 residual degrees of
freedom.

Every contrast this study asks for is within-participant, which is the paired case of the
limma User's Guide 9.4.2 — "the treatments are compared only within each block". A fixed block
term is what that section prescribes, so `dpcDE()` is called **without** `block =` and no
intra-block correlation is estimated. The price is that `sex` is collinear with the
participant term and cannot be tested by this model.

The group names are chosen to survive `makeContrasts()`, which parses its arguments as R code.
Boersma named them `2E-T1`, which R reads as a number minus a variable.

## What to know before using it

**33 participants, 66 legs, 131 samples.** 66 legs × 2 timepoints is 132.

**One sample is missing.** S06 left leg at T2 was never acquired, so `BFR_T2` has 32 samples where the other three groups have 33.

**Eight samples were reinjected**, flagged by a suffix on `run` — seven `…r` and one `…rR` (`CL_MDR_12D-T2rR`), so match case-insensitively. They fall 1/3/3/1 across the four groups, so no group effect is possible. Only the retained injection appears in `precursors.csv`.

**One injection was discarded.** S08 left leg T1 was run twice; the first recovered 2,858 precursors against 30,326 for the second, and the two correlate at r = 0.032 on log2. Those two figures come from the DIA-NN report, not from anything in this directory.

**Participant IDs are not contiguous.** S01, S19, S22 and S34 have no MS data.

**The batch data.** No acquisition date, instrument or plate is recorded anywhere in the inputs, so the prefix is kept in `run` and modelled nowhere.

**The phenotypes are in their own file, at their own grain.** `phenotype.csv` has 70 rows against 131 in `metadata.csv`; join on `participant` + `leg`, never `cbind`. See the section above for why it is not folded in.