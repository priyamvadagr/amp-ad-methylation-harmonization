# ============================================================================
# Compile each cohort's .idat files into a minfi RGChannelSet object.
#
#   RGChannelSet = the raw red/green intensity object that all downstream
#   minfi steps (QC, noob, filtering, beta/M-values) build on.
#
# Reads the per-cohort IDAT manifest written by Scripts/Misc/get_sample_info.Rmd
# (specimenID/individualID <-> Grn/Red IDAT file, already restricted to
# samples with both metadata and IDAT files present), builds one RGChannelSet
# per cohort, and checkpoints each to disk.
#
# Usage:
#   Rscript compile_idat_files.R          # process all cohorts
#   Rscript compile_idat_files.R MSBB ROSMAP
#                                          # process only the named cohort(s)
#   Rscript compile_idat_files.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir"
#                                          # run a cohort not in R/cohort_config.R
#                                          # (individual_csv/qc_dir are unused here
#                                          # but required for cohorts_all.R's shared
#                                          # config shape -- see sample_level_QC.R)
#
# Cohort paths (raw_dir, manifest, combined_dir) live in R/cohort_config.R,
# shared with sample_level_QC.R so the two stages can't drift apart.
# ============================================================================

library(minfi)
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/read_cohort_idats.R')
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cohort_config.R')

cohorts <- select_cohorts(commandArgs(trailingOnly = TRUE))

# Samples per read.metharray.exp() call -- bounds peak memory on large
# cohorts so the compile step doesn't get OOM-killed (see chunk_size doc in
# read_cohort_idats()).
chunk_size <- 75

for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  cat("\n==== ", cohort_name, " ====\n")

  manifest <- read.csv(cfg$manifest, stringsAsFactors = FALSE)

  # Drop any specimen the manifest couldn't resolve to both idat files
  unresolved <- is.na(manifest$grnFile) | is.na(manifest$redFile)
  if (any(unresolved)) {
    cat(sum(unresolved), "specimen(s) missing a Grn/Red file in the manifest; skipping\n")
    manifest <- manifest[!unresolved, ]
  }

  # Basename minfi needs: full path to each idat pair, minus the _Grn/_Red
  # suffix. Built from grnFile (not Sentrix_ID) so per-barcode subfolders
  # (e.g. ROSMAP) are preserved.
  manifest$Basename <- file.path(cfg$raw_dir, sub('_Grn\\.idat$', '', manifest$grnFile))

  rg <- read_cohort_idats(idat_dir = cfg$raw_dir, sample_sheet = manifest,
                          chunk_size = chunk_size)

  cat(cohort_name, "RGChannelSet:", ncol(rg), "samples,", nrow(rg), "probes\n")

  dir.create(cfg$combined_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(rg, file.path(cfg$combined_dir, paste0(cohort_name, '.rds')))
  rm(rg); gc()
}
