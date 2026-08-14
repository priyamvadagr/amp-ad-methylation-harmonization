# ============================================================================
# Run the vendored CorticalClock DNAm-age predictor
# (/home/ec2-user/CorticalClock/PredCorticalAge/CorticalClock.r, Shireby et
# al., https://github.com/gemmashireby/CorticalClock) on each cohort's
# PRE-FILTER noob beta matrix.
#
# CorticalClock() is used unmodified. Two things about it that shape this
# wrapper:
#   - it always writes to the CURRENT working directory under fixed names
#     (CorticalPred.csv, accuracy_statistics.csv, CorticalClockplot.pdf) --
#     so this script setwd()s into a per-cohort output dir before each call
#     (restored via on.exit even if the call errors) and renames the
#     fixed-name outputs to <cohort>_* afterward so cohorts don't collide.
#   - `dir` is pasted directly onto "CorticalClockCoefs.txt" with no
#     separator inserted, so it MUST end in "/".
#
# Reads <normalized_dir>/<cohort>_noob_beta.rds (pre probe-filter -- probe
# filtering may drop clock CpGs) plus each cohort's idat manifest +
# individual metadata (for sample IDs + ageDeath) and writes
# <plot_root>/<cohort>/<cohort>_CorticalPred.csv,
# <cohort>_accuracy_statistics.csv, <cohort>_CorticalClockplot.pdf,
# <cohort>_<tissue>_clock_accuracy.pdf.
#
# ageDeath is censored at "90+" in these cohorts' metadata; as.numeric()
# turns that into NA, which just drops those samples from the accuracy
# stats (correlation/RMSE/MAD) -- their predicted age is still written to
# <cohort>_CorticalPred.csv.
#
# Usage:
#   Rscript run_cortical_clock.R                  # all cohorts
#   Rscript run_cortical_clock.R MSBB ROSMAP       # named cohorts only
#   Rscript run_cortical_clock.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir:/path/normalized_dir"
#
# Cohort paths (manifest, individual_csv, normalized_dir) live in
# R/cohort_config.R, shared with the rest of the pipeline so stages can't
# drift apart.
# ============================================================================

source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cohort_config.R')
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cortical_clock.R')
source('/home/ec2-user/CorticalClock/PredCorticalAge/CorticalClock.r')

clock_dir <- '/home/ec2-user/CorticalClock/PredCorticalAge/'   # trailing slash required, see above
plot_root <- '/home/ec2-user/AMP-AD_methylation_harmonization/Results/Cortical_clock'

# Build one pheno row per sample (same "<Sentrix_ID>_<Sentrix_Row_Column>" ID
# minfi uses as colnames(beta)); mirrors pca_qc_per_cohort.R's build_pheno().
build_pheno <- function(cfg) {
  manifest <- read.csv(cfg$manifest, stringsAsFactors = FALSE)
  manifest$sample <- paste0(manifest$Sentrix_ID, "_", manifest$Sentrix_Row_Column)

  individual <- read.csv(cfg$individual_csv, stringsAsFactors = FALSE)
  individual$ageDeath <- suppressWarnings(as.numeric(individual$ageDeath))  # "90+" -> NA

  merge(manifest, individual, by = "individualID", all.x = TRUE)
}

# Run CorticalClock() for one cohort, isolated in its own output dir (its
# fixed output filenames would otherwise collide across cohorts) with the
# working directory always restored afterward.
run_clock <- function(cohort_name, cfg) {
  cat("\n==== ", cohort_name, " ====\n")

  beta  <- readRDS(file.path(cfg$normalized_dir, paste0(cohort_name, "_noob_beta.rds")))
  pheno <- build_pheno(cfg)

  out_dir <- file.path(plot_root, cohort_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)
  setwd(out_dir)

  CorticalClock(betas = beta, pheno = pheno,
               dir = clock_dir, IDcol = "sample", Agecol = "ageDeath")

  # CorticalClock() drops every pheno column except IDcol/Agecol before
  # writing CorticalPred.csv (ID, Age, brainpred only) -- join specimenID/
  # individualID back on via ID (== pheno$sample) since they're otherwise lost.
  if (file.exists("CorticalPred.csv")) {
    pred <- read.csv("CorticalPred.csv", stringsAsFactors = FALSE, row.names = 1)
    m <- match(pred$ID, pheno$sample)
    pred$individualID <- pheno$individualID[m]
    pred$specimenID   <- pheno$specimenID[m]
    pred$tissue       <- pheno$tissue[m]
    pred <- pred[c("ID", "individualID", "specimenID", "tissue", "Age", "brainpred")]
    write.csv(pred, "CorticalPred.csv", row.names = FALSE)

    # per-tissue accuracy plots (R/cortical_clock.R) -- a clock trained on
    # cortex can be systematically off in non-cortical tissue, which the
    # single pooled CorticalClockplot.pdf/accuracy_statistics.csv would hide.
    plot_clock_by_tissue(pred, tissue_col = "tissue", cohort = cohort_name,
                        out_prefix = cohort_name)
  }

  for (f in c("CorticalPred.csv", "accuracy_statistics.csv", "CorticalClockplot.pdf")) {
    if (file.exists(f)) file.rename(f, paste0(cohort_name, "_", f))
  }
}

cohorts <- select_cohorts(commandArgs(trailingOnly = TRUE))

for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  run_clock(cohort_name, cfg)
  gc()
}
