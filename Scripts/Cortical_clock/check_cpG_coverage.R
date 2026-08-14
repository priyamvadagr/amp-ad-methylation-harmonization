# ============================================================================
# Run subset_clock_cpgs() + plot_clock_missingness() (R/cortical_clock.R) on
# each cohort's PRE-FILTER noob beta matrix, to check coverage of the 347
# cortical-clock CpGs before running the clock itself.
#
# Reads <normalized_dir>/<cohort>_noob_beta.rds (output of normalize_data.R;
# deliberately the pre probe-filter beta -- filtering may drop clock CpGs)
# plus <qc_dir>/sample_qc/<cohort>_qc.rds (the QC'd RGChannelSet, res$rg from
# sample_qc() -- same file normalize_data.R reads) to recompute detection
# p-values, since normalize_data.R doesn't persist its detP matrix to disk.
# Writes <data_root>/<cohort>/Cortical_clock/<cohort>_clock_sample_coverage.csv,
# _clock_cpg_coverage.csv (data) and
# <plot_root>/<cohort>/<cohort>_clock_missingness.pdf (figure).
#
# Usage:
#   Rscript check_cpG_coverage.R                  # all cohorts
#   Rscript check_cpG_coverage.R MSBB ROSMAP       # named cohorts only
#   Rscript check_cpG_coverage.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir:/path/normalized_dir"
#
# Cohort paths (qc_dir, normalized_dir) live in R/cohort_config.R, shared
# with the rest of the pipeline so stages can't drift apart.
#
# Read-out: ideally 347/347 present and the per-sample histogram sits near 0%.
# Red bars in plot 2 are CpGs entirely absent (must be imputed with reference
# betas). The plot 3 detP heatmap and plot 4 UpSet co-occurrence plot show
# whether masking is scattered or concentrated in specific samples -> those
# samples' DNAmAge will lean on imputation, interpret with caution.
# ============================================================================

library(minfi)
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cortical_clock.R')
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cohort_config.R')

coef_file <- '/home/ec2-user/CorticalClock/PredCorticalAge/CorticalClockCoefs.txt'
data_root <- '/home/ec2-user/data/methyl_harmonization'
plot_root <- '/home/ec2-user/AMP-AD_methylation_harmonization/Results/Cortical_clock'

clock_cpgs <- clock_cpgs_from_coef(coef_file)
cohorts    <- select_cohorts(commandArgs(trailingOnly = TRUE))

for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  cat("\n==== ", cohort_name, " ====\n")

  # detP, restricted to clock CpGs right away (rg is the full, all-samples
  # RGChannelSet -- detectionP() is the expensive/memory-heavy step, same as
  # in noob_normalize.R, so rg is freed the instant detP is in hand).
  rg   <- readRDS(file.path(cfg$qc_dir, "sample_qc", paste0(cohort_name, "_qc.rds")))
  detP <- minfi::detectionP(rg)
  rm(rg); gc()
  detP <- detP[intersect(clock_cpgs, rownames(detP)), , drop = FALSE]

  beta <- readRDS(file.path(cfg$normalized_dir, paste0(cohort_name, "_noob_beta.rds")))
  detP <- detP[, colnames(beta), drop = FALSE]   # normalization may have dropped failed samples

  sub <- subset_clock_cpgs(beta, coef_file = coef_file)
  plot_clock_missingness(sub$mat, present = sub$present, detP = detP,
                         cohort = cohort_name,
                         out_prefix  = file.path(data_root, cohort_name, "Cortical_clock", cohort_name),
                         plot_prefix = file.path(plot_root, cohort_name, cohort_name))

  rm(beta, detP, sub); gc()
}
