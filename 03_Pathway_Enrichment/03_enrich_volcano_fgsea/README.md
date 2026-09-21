# 03_enrich_volcano_fgsea

Protein volcanoes with the surviving pathways ringed. Computes nothing.

| | |
|---|---|
| **Script** | `a_script/03_enrich_volcano_fgsea.R` |
| **Reads** | `fgsea.rds`, `design.rds`, `config.yml` |
| **Writes** | `b_reports/figures/`, six PNGs and six PDFs |

Needs `enrichVolcano`, which is not on CRAN.

## What a panel claims

Two separate claims share each panel. Point colours and the count badges read protein-level BH
FDR; the ring reads set-level fgsea FDR, and only sets marked `main` are ringed. Red is up, blue
is down, for points and arcs alike.

Every panel names its contrast in the coefficients `01_Design` fitted — `BFR_T2 - BFR_T1`, not
"post vs pre" — and the interaction carries the nested form, `(BFR_T2 - BFR_T1) - (HLRT_T2 -
HLRT_T1)`, because a difference of differences is what it is. A `stopifnot` evaluates each string
against the fitted contrast matrix, so a notation that drifts from the design stops the run.

## The control is not drawn

Four FDR panels — interaction, the two training responses, the post-training comparison — plus two
pi-ranked repeats. The pre-training control is fitted, tested and reported in `01_run_fgsea`, where
it returns six pathways on a comparison that cannot contain signal. That six is what calibrates
every ring here, and it is not visible on these figures.

## Two traps the script handles

`volcano_ring()` matches leading-edge genes against point labels, and a label carries its accession
when a symbol sits on more than one protein. The script translates `leadingEdge` into label space
first; skip that and the tick lines draw nothing, with no warning.

Ten overlapping databases mean two ringed sets can reduce to the same display name —
`GOBP_MUSCLE_CONTRACTION` and `REACTOME_MUSCLE_CONTRACTION` both render as "Muscle Contraction".
The script appends the database when that happens.
