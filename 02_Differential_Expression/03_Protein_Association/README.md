# 03 · Protein association

Tests whether each protein's training change tracks the same leg's outcome: change in VL and RF
cross-sectional area (CSA) and echo intensity, and total volume load. Each slope is a `dpcDE()`
coefficient with `eBayes(robust = TRUE)` and BH within each outcome and coefficient. The average
and difference t values rank proteins for `03_Pathway_Enrichment/05_set_association`.

```
arm:    ~ 0 + leg_id + treatment:post + treatment:post:x + blood
within: ~ 0 + leg_id + participant:post + bfr_post + post:x + blood
```

`post` is 1 at T2; `x` is the leg's outcome, centred. The arm model gives BFR, HLRT, difference
and average slopes from all 65 legs with both samples; the within slope comes from the 64 legs of
the 32 participants with both legs. Volume load has no within slope, because dose differs by arm by
design. The notebook explains each term.

Inputs: `01_Preprocess/02_Quantification/c_data/proteins.rds` (models), `proteins_adjusted.rds`
(figures), `00_Input/phenotype.csv`.

```sh
quarto render \
  02_Differential_Expression/03_Protein_Association/a_script/03_protein_association.qmd \
  --output-dir ../b_reports
```

Outputs:

- `b_reports/03_protein_association.html`: the report.
- `b_reports/03_association_<outcome>.pdf`: per outcome, the three smallest p per coefficient,
  protein change against outcome; panel headers give t, p and BH q.
- `c_data/association_t.rds`: moderated t, proteins by ranking (`average_<outcome>`,
  `difference_<outcome>`), for fgsea in `05_set_association`.
- `c_data/association_designs.rds`: per ranking, the arm design and contrast vector, for fry.
- `c_data/03_protein_association.xlsx`:
  - `outcomes`: per outcome and coefficient, legs used, formula, covariates, nominal and BH
    counts, best FDR and top five genes with direction.
  - `top_leads`: the 25 smallest p per outcome and coefficient.
  - `associations`: every protein, outcome and coefficient: slope, t, p and FDR.

Results: volume load is the only outcome with BH hits. Across BFR legs, CA14 (t -4.9, q 0.019) and
TOMM40 (t -4.6, q 0.022) fall more and HSPA12A (t 4.7, q 0.019) rises more in legs that lifted
more; the average slope keeps CA14 (q 0.004) and TOMM40 (q 0.018), and HSPA12A is the one arm
difference (t 4.8, q 0.034). No protein tracks a CSA or echo change after BH; the closest is AGMAT
on RF CSA (average, q 0.11). Nominal counts sit near the 152 expected by chance except VL echo in
HLRT legs (305) and its arm difference (361). The volume-load slopes compare people, so strength,
size and sex ride along with dose.
