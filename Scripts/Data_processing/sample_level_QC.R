# ============================================================================
# Run sample-level QC on each cohort's compiled RGChannelSet, in two phases so
# the raw idats (ewastools) and the RGChannelSet (minfi) are never resident in
# memory at the same time -- that co-residency is what was getting the job
# killed on the large cohorts.
#
#   phase=ewastools : read raw idats in memory-bounded chunks and write the
#                     control-probe QC (bisulfite + sample_failure) to
#                     <out_prefix>_ewastools.csv. Does NOT load the rds.
#   phase=minfi     : load the RGChannelSet, run detection-p / intensity / sex
#                     QC, MERGE the ewastools CSV from the previous phase, write
#                     <out_prefix>_sample_qc.csv, render figures, save the
#                     flagged RGChannelSet. Does NOT read raw idats.
#   phase=all       : run ewastools THEN minfi per cohort in one process
#                     (default). Safe on smaller cohorts because ewastools_qc()
#                     frees its idats before the rds is loaded; for the cohorts
#                     that get killed, run the two phases as SEPARATE Rscript
#                     calls so each process only holds one of the two.
#
# The ewastools checks need the raw idat files, so each cohort supplies its
# idat manifest and raw idat directory. Basenames are built from the manifest's
# grnFile column (same as compile_idat_files.R) -- this also covers ROSMAP,
# whose grnFile already includes the per-Sentrix-barcode subfolder prefix
# (e.g. "5815381027/5815381027_R06C01_Grn.idat").
#
# Reported sex is merged in from <cohort>_individual_metadata_processed.csv
# (not carried by the idat manifest), matched on individualID.
#
# Usage:
#   Rscript sample_level_QC.R                       # all cohorts, phase=all
#   Rscript sample_level_QC.R MSBB ROSMAP           # named cohorts, phase=all
#   Rscript sample_level_QC.R phase=ewastools MSBB  # only the ewastools phase
#   Rscript sample_level_QC.R phase=minfi     MSBB  # only the minfi phase
#   Rscript sample_level_QC.R chunk=50 phase=ewastools ROSMAP   # smaller chunks
#   Rscript sample_level_QC.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir"
#
# Recommended for large cohorts (two processes):
#   Rscript sample_level_QC.R phase=ewastools ROSMAP
#   Rscript sample_level_QC.R phase=minfi     ROSMAP
#
# Cohort paths (raw_dir, manifest, combined_dir, individual_csv, qc_dir) live
# in R/cohort_config.R, shared with compile_idat_files.R so the two stages
# can't drift apart.
# ============================================================================

library(minfi)
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/sample_qc.R')
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cohort_config.R')

# Attach reported sex onto pData(rg) from the cohort's processed individual
# metadata, matched on individualID.
attach_sex <- function(rg, individual_csv) {
  ind <- read.csv(individual_csv, stringsAsFactors = FALSE)
  pData(rg)$sex <- ind$sex[match(pData(rg)$individualID, ind$individualID)]
  rg
}

# Where render_qc_plots() writes its per-cohort figures (subfolder per cohort).
plot_root <- '/home/ec2-user/AMP-AD_methylation_harmonization/Results/QC/raw_values'

# --- parse args: phase=, chunk=, cohort names, ad-hoc specs -----------------
args      <- commandArgs(trailingOnly = TRUE)
phase     <- "all"
chunk_sz  <- 75
requested <- character(0)

for (arg in args) {
  if (grepl("^phase=", arg)) {
    phase <- sub("^phase=", "", arg)
  } else if (grepl("^chunk=", arg)) {
    chunk_sz <- as.integer(sub("^chunk=", "", arg))
  } else {
    requested <- c(requested, arg)
  }
}
if (!phase %in% c("all", "ewastools", "minfi"))
  stop("phase must be one of: all, ewastools, minfi (got '", phase, "')")

# --- select cohorts (shared with compile_idat_files.R; see R/cohort_config.R)
cohorts <- select_cohorts(requested)

# --- per-phase workers ------------------------------------------------------

# ewastools phase: read raw idats (chunked) -> <qc_dir>/sample_qc/<cohort>_ewastools.csv.
run_ewastools_phase <- function(cohort_name, cfg) {
  cat("\n==== ", cohort_name, " [ewastools] ====\n")
  manifest <- read.csv(cfg$manifest, stringsAsFactors = FALSE)
  ewastools_qc(targets   = manifest,
               raw_dir   = cfg$raw_dir,
               grn_col   = "grnFile",
               chunk_size = chunk_sz,
               out_csv   = file.path(cfg$qc_dir, "sample_qc", paste0(cohort_name, "_ewastools.csv")))
  rm(manifest); gc()
}

# minfi phase: load rg -> minfi QC + merge ewastools CSV -> figures + save.
run_minfi_phase <- function(cohort_name, cfg) {
  cat("\n==== ", cohort_name, " [minfi] ====\n")
  out_prefix <- file.path(cfg$qc_dir, "sample_qc", cohort_name)
  ew_csv     <- paste0(out_prefix, "_ewastools.csv")
  if (!file.exists(ew_csv))
    warning("[", cohort_name, "] no ewastools CSV (", ew_csv,
            "); bisulfite/control flags will be NA. Run phase=ewastools first.")

  rg <- readRDS(file.path(cfg$combined_dir, paste0(cohort_name, ".rds")))
  rg <- attach_sex(rg, cfg$individual_csv)

  dir.create(cfg$qc_dir, recursive = TRUE, showWarnings = FALSE)
  res <- sample_qc(rg,
                   sex_col       = "sex",
                   ewastools_csv = ew_csv,      # merged if present, else NA flags
                   out_prefix    = out_prefix)

  cat(cohort_name, "sample QC:", nrow(res$qc_tab), "samples,",
      sum(res$qc_tab$any_flag), "flagged (nothing dropped)\n")

  render_qc_plots(res$rg, cohort = cohort_name, out_dir = plot_root, qc_tab = res$qc_tab)
  saveRDS(res$rg, paste0(out_prefix, "_qc.rds"))

  rm(rg, res); gc()
}

# --- run --------------------------------------------------------------------
for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  dir.create(cfg$qc_dir, recursive = TRUE, showWarnings = FALSE)

  # ewastools first (frees idats) THEN minfi, so the two are never co-resident.
  if (phase %in% c("all", "ewastools")) run_ewastools_phase(cohort_name, cfg)
  if (phase %in% c("all", "minfi"))     run_minfi_phase(cohort_name, cfg)
}