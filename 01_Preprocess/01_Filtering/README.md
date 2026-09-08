# 01_Preprocess / 01_Filtering

Reads the DIA-NN report and removes rows. It drops no sample, and it filters nothing on
missingness.

| | |
|---|---|
| **Script** | `a_script/01_filter.qmd` |
| **Reads** | `00_Input/report.parquet`, `00_Input/metadata.csv` |
| **Writes** | `c_data/precursors_filtered.rds`, `c_data/filter_log.csv`, `c_data/contaminants_removed.csv` |
| **Read by** | `01_Preprocess/02_Quantification` |

## What it does

**Read.** `readDIANN()` with one change to the defaults: `Protein.Ids` is added, because the
repair below cannot run without it. It is added, not substituted, so `Protein.Names` still
comes through for later annotation. Everything else stays as shipped. `q.cutoffs = 0.01`
over `Q.Value`, `Lib.Q.Value` and `Lib.PG.Q.Value` is the set the vignette recommends for a
match-between-runs search. `intensity.column = "Precursor.Normalised"` is DIA-NN's
normalized column, which is why nothing here normalizes again. `log = TRUE` means `y$E`
arrives on the log2 scale, so never log it a second time.

**Map runs to samples.** Columns come back named by DIA-NN `Run`. One case-insensitive
match against `metadata$run` does three jobs: it renames columns to `sample_id`, it absorbs
a run the instrument and the sheet capitalised differently, and it drops the discarded S08
injection, which is in the report but not the sheet.

**Check what survived.** `metadata.csv` decides which runs enter, so a failed injection is
excluded by hand and the code has no way to notice. One chunk reports precursor count,
median intensity and missingness for each sample that made it, and fails if any of them sits
below half the median depth. It drops nothing. It checks that the exclusion list was right,
and it is what would catch a re-search that failed on a different run.

**Repair the annotation.** The search FASTA held every protein twice, once plain and once
with a `cRAP-` prefix, so myoglobin arrives as `P02144;cRAP-P02144` with `Proteotypic = 0`.
Drop every group that mentions `cRAP` and you delete most of the muscle. The rule instead:
drop a group only when *every* member is cRAP, then strip the prefix from the rest,
de-duplicate the members, and recompute `Proteotypic`. An assertion confirms no `cRAP`
string survives.

Dropping has to come first, and the contaminant list cannot do this job, because `Genes`
carries the prefix too. Six cRAP-only symbols appear in this dataset: `cRAP-MB`, `cRAP-ALB`,
`cRAP-LYZ`, `cRAP-AHSG`, `cRAP-CSN1S1`, and two rows with no symbol at all. Strip the prefix
first and `cRAP-MB` becomes `MB`, which `MUST_KEEP` requires to survive, so ten precursors
of a contaminant standard would merge into human myoglobin and nobody would see it.

**limpa's filters**, in the vignette's order: `filterNonProteotypicPeptides()` then
`filterCompoundProteins()`. We skip `filterSingletonPeptides()` on the vignette's own advice
to keep all proteins and hold on to as much information as possible.

**Contaminants.** A needle biopsy bleeds and passes through skin, so much of the signal is
blood, plasma, skin and antibody protein. These rows go before quantification, because
`dpc()` fits one curve across every precursor and `dpcQuant()` builds its prior from the
whole matrix. Leave contaminants in and they shape the missing-value model that every muscle
protein is then quantified through.

The list is a plain vector in the notebook: 133 symbols and three family patterns. It comes
from reading every gene symbol in this dataset against single-cell myonuclei expression
rather than whole-muscle tissue, because tissue measurements are themselves
blood-contaminated. Twelve blood-abundant symbols are kept on purpose, and they sit in the
same `MUST_KEEP` vector as the muscle sentinels rather than in a list of their own.

A group counts as contamination only when **every** one of its gene symbols does.
`filterCompoundProteins()` already removed groups holding several proteins, but one protein
can still carry several gene names, so testing the first symbol or any symbol would get it
wrong in both directions.

## The guard

One assertion, run after every filter: all 28 symbols in `MUST_KEEP` must survive. It exists
in order to fail. One line catches two separate faults, a repair that breaks and deletes
muscle through the `Proteotypic` flag, and a contaminant rule that reaches too far and
deletes muscle directly. Run it on the unrepaired annotation and it fails on thirteen
symbols.

## Not measured here

How much of each sample's signal was blood. Biopsy technique varies between legs and between
visits, and blood takes MS duty cycle away from muscle protein, so the amount could differ
between the groups and confound the comparison. Someone should check whether it does before
write-up. The haemoglobin rows are still in the matrix at the contaminant chunk, which is
where that measurement would have to be taken.
