# Half of a DIA matrix is empty, and not at random: faint precursors go missing more often than bright ones. `dpc()` fits that as a curve, and
#`dpcQuant()` uses it to give every protein a value in every sample without imputing anything. This stage is the limpa vignette's own sequence, unaltered.
library(here)
library(limpa)
library(limma)
library(dplyr)
library(tibble)
library(writexl)

precursors_file <- here("01_Preprocess", "01_Filtering", "c_data",
                        "precursors_filtered.rds")

y <- readRDS(precursors_file)
dim(y)

checkpoints <- here("01_Preprocess", "02_Quantification", "c_data", "proteins_quantified")
primary_file <- file.path(checkpoints, "proteins_slope_0.7.rds")
sensitivity_file <- file.path(checkpoints, "proteins_slope_fitted.rds")

# The missingness relationship: A downward trend means detection depends on intensity and the DPC has something to fit. 
plotAveVsMis(y)

# The curve
dpcest <- dpc(y$E)
plotDPC(dpcest)

dpc_parameters <- tibble(
  intercept = dpcest$dpc[1],
  slope = dpcest$dpc[2]
)
round(dpc_parameters, 3)

proteins <- readRDS(here("01_Preprocess", "02_Quantification", "c_data", "proteins_quantified.rds"))
dim(proteins)

# Detection Filter
group_sizes <- table(proteins$targets$group)
keep <- filterByDetection(proteins, n.samples = min(group_sizes))
proteins <- proteins[keep,]
# Removes protiens that do not add much to the information and those that have a high SE relative to others (0.66)
group_sizes
c(before = length(keep), kept = sum(keep), dropped = sum(!keep))

# Normlalization
# limpa places normalisation after dpcQuant() and before dpcDE() and performs none itself
#  because DIA-NN has usually normalised the precursors already. What matters is which kind is needed, and two numbers separate the cases.
# cyclic loess fixes: comparing a sample against the average of all samples gives a difference per protein, and fitting that against abundance gives a curve. 
# A flat curve means the sample disagrees by a constant. A curve that travels means the disagreement depends on abundance. Both numbers are reported before and after so the correction can be judged rather than assumed.
a <- rowMeans(proteins$E)

ma_width <- function(E) {
  a <- rowMeans(E)
  apply(E, 2, \(x) diff(range(lowess(a, x - a, f = 0.5)$y)))
}
spread <- function(E) diff(range(apply(E, 2, median)))

before <- c(ma_width = median(ma_width(proteins$E)), median_spread = spread(proteins$E))
proteins$E <- normalizeBetweenArrays(proteins$E, method = "cyclicloess")
after <- c(ma_width = median(ma_width(proteins$E)), median_spread = spread(proteins$E))

round(rbind(before, after), 3)
# The curvature is unrelated to treatment or timepoint, so cyclic loess is not absorbing the effect being tested

# Cyclic loess remaps $E non-linearly and leaves $other$standard.error untouched
# After normalisation the stored errors describe the unnormalised scale, and dpcDE() reads them as though they described the new one.

# What Survives
protein_quality <- tibble(
  protein = rownames(proteins),
  gene = proteins$genes$Genes,
  n_precursors = proteins$genes$NPrec,
  prop_observed = proteins$genes$PropObs,
  mean_log2 = rowMeans(proteins$E),
  mean_standard_error = rowMeans(proteins$other$standard.error),
  mean_n_observations = rowMeans(proteins$other$n.observations)
)
summary(protein_quality[c("mean_log2", "mean_standard_error", "prop_observed")])

# A protein with a high mean standard error and few observations is one whose estimate leans on the DPC rather than on data.
mds <- plotMDSUsingSEs(proteins,
  pch = 16,
  col = as.integer(factor(proteins$targets$participant))
)

r2 <- function(f) summary(lm(mds$x ~ factor(f)))$r.squared
round(c(
  participant = r2(proteins$targets$participant),
  group = r2(proteins$targets$group)
), 3)

# Ouputs
out <- here("01_Preprocess", "02_Quantification", "c_data")
saveRDS(proteins, file.path(out, "proteins.rds"))
write_xlsx(
  list(dpc_parameters = dpc_parameters, protein_quality = protein_quality),
  file.path(out, "02_quantify.xlsx")
)
