# ============================================================================
# Run sample_qc() + render_qc_plots() (see R/sample_qc.R) for each cohort's
# saved RGChannelSet .rds, one cohort at a time. Each cohort's RGChannelSet
# is dropped from memory (rm + gc) before the next one is loaded.
#
# Per cohort, reads:
#   <data_root>/<cohort>/combined_idat/<cohort>.rds
# and writes:
#   <data_root>/<cohort>/qc/<cohort>_sample_qc.csv        (sample_qc)
#   <plot_root>/<cohort>/<cohort>_*.{pdf,png}             (render_qc_plots)
#   <plot_root>/<cohort>/<cohort>_flagged_qcReport.pdf    (render_qc_plots)
# ============================================================================

source("/home/ec2-user/AMP-AD_methylation_harmonization/R/sample_qc.R")

data_root <- "/home/ec2-user/data/methyl_harmonization"
plot_root <- "/home/ec2-user/AMP-AD_methylation_harmonization/Results/QC/raw_values"

cohorts <- c("ROSMAP", "ROSMAP_APOE4", "MSBB", "MOA-PAD")

for (cohort in cohorts) {

  rds_file <- file.path(data_root, cohort, "combined_idat", paste0(cohort, ".rds"))
  qc_dir   <- file.path(data_root, cohort, "qc")
  dir.create(qc_dir, showWarnings = FALSE, recursive = TRUE)

  rg  <- readRDS(rds_file)   # full object = original view
  res <- sample_qc(rg, out_prefix = file.path(qc_dir, cohort))
  render_qc_plots(res$rg, cohort = cohort, out_dir = plot_root, qc_tab = res$qc_tab)

  rm(rg, res)
  gc()
}
