# ============================================================================
# Re-render sample-QC plots (render_qc_plots(), Scripts/functions/sample_qc.R)
# from the QC'd RGChannelSet + qc_tab that sample_level_QC.R already computed
# and saved -- does NOT rerun sample_qc() itself. Rerunning sample_qc() here
# would recompute the minfi checks from scratch AND lose the ewastools merge
# (bisulfite/control flags), since this script has no access to the raw idats
# sample_level_QC.R's ewastools phase reads.
#
# Reads <qc_dir>/sample_qc/<cohort>_qc.rds and _sample_qc.csv (written by
# sample_level_QC.R's minfi phase) and writes the plots into that SAME
# directory, alongside the data files:
#   <qc_dir>/sample_qc/<cohort>_*.{pdf,png}             (render_qc_plots)
#   <qc_dir>/sample_qc/<cohort>_flagged_qcReport.pdf    (render_qc_plots)
#
# Usage:
#   Rscript render_qc_plots.R                  # all cohorts
#   Rscript render_qc_plots.R MSBB ROSMAP       # named cohorts only
#   Rscript render_qc_plots.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir"
#
# Cohort paths (qc_dir) live in Scripts/functions/cohort_config.R, shared with
# sample_level_QC.R so the two scripts can't drift apart.
# ============================================================================

source("/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/sample_qc.R")
source("/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cohort_config.R")

cohorts <- select_cohorts(commandArgs(trailingOnly = TRUE))

for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  cat("\n==== ", cohort_name, " ====\n")

  sample_qc_dir <- file.path(cfg$qc_dir, "sample_qc")
  qc_rds <- file.path(sample_qc_dir, paste0(cohort_name, "_qc.rds"))
  qc_csv <- file.path(sample_qc_dir, paste0(cohort_name, "_sample_qc.csv"))

  rg     <- readRDS(qc_rds)
  qc_tab <- read.csv(qc_csv, stringsAsFactors = FALSE)

  render_qc_plots(rg, cohort = cohort_name, out_dir = cfg$qc_dir,
                  subdir = "sample_qc", qc_tab = qc_tab)

  rm(rg, qc_tab)
  gc()
}
