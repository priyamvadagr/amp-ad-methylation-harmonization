# ============================================================================
# Run noob normalization on each cohort's QC'd RGChannelSet (output of
# sample_level_QC.R's minfi phase).
#
# Reads <qc_dir>/sample_qc/<cohort>_qc.rds (res$rg from sample_qc()) and
# <qc_dir>/sample_qc/<cohort>_sample_qc.csv (qc_tab), drops samples flagged on
# drop_reasons, runs preprocessNoob(), maps to the genome, masks
# high-detection-p probe x sample values, and saves the result to
# <normalized_dir>/<cohort>_noob_gmset.rds and _noob_beta.rds.
#
# Usage:
#   Rscript normalize_data.R                  # all cohorts
#   Rscript normalize_data.R MSBB ROSMAP       # named cohorts only
#   Rscript normalize_data.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir:/path/normalized_dir"
#
# Cohort paths (qc_dir, normalized_dir) live in R/cohort_config.R, shared with
# compile_idat_files.R and sample_level_QC.R so the pipeline stages can't
# drift apart.
# ============================================================================

library(minfi)
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/noob_normalize.R')
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cohort_config.R')

cohorts <- select_cohorts(commandArgs(trailingOnly = TRUE))

for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  cat("\n==== ", cohort_name, " ====\n")

  qc_rds <- file.path(cfg$qc_dir, "sample_qc", paste0(cohort_name, "_qc.rds"))
  qc_csv <- file.path(cfg$qc_dir, "sample_qc", paste0(cohort_name, "_sample_qc.csv"))

  rg <- readRDS(qc_rds)
  qc <- read.csv(qc_csv, stringsAsFactors = FALSE)

  dir.create(cfg$normalized_dir, recursive = TRUE, showWarnings = FALSE)
  res <- normalize_noob(rg, qc_tab = qc,
                        out_prefix = file.path(cfg$normalized_dir, cohort_name),
                        save_gmset = TRUE, save_mvals = TRUE,
                        return_objects = FALSE)

  rm(rg, qc, res); gc()
}
