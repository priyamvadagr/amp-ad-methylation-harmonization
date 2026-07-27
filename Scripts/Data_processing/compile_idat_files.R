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
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/out_dir"
#                                          # run a cohort not hardcoded below
# ============================================================================

library(minfi)
source('/home/ec2-user/AMP-AD_methylation_harmonization/R/read_cohort_idats.R')

cohorts_all <- list(
  MSBB = list(
    raw_dir  = '/home/ec2-user/data/methyl_harmonization/methyl_data/MSBB/raw_data',
    manifest = '/home/ec2-user/data/methyl_harmonization/MSBB/metadata/MSBB_idat_manifest.csv',
    out_dir  = '/home/ec2-user/data/methyl_harmonization/MSBB/combined_idat/'
  ),
  ROSMAP = list(
    raw_dir  = '/home/ec2-user/data/methyl_harmonization/methyl_data/ROSMAP/raw_data',
    manifest = '/home/ec2-user/data/methyl_harmonization/ROSMAP/metadata/ROSMAP_idat_manifest.csv',
    out_dir  = '/home/ec2-user/data/methyl_harmonization/ROSMAP/combined_idat/'
  ),
  ROSMAP_APOE4 = list(
    raw_dir  = '/home/ec2-user/data/methyl_harmonization/methyl_data/ROSMAP_APOE4/raw_data',
    manifest = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/metadata/ROSMAP_APOE4_idat_manifest.csv',
    out_dir  = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/combined_idat/'
  ),
  `MOA-PAD` = list(
    raw_dir  = '/home/ec2-user/data/methyl_harmonization/methyl_data/MOA-PAD/raw_data',
    manifest = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/metadata/MOA-PAD_idat_manifest.csv',
    out_dir  = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/combined_idat/'
  )
)

# Select which cohort(s) to process from the command line; default to all
# built-in cohorts. Each arg is either:
#   - a name from cohorts_all (e.g. "MSBB"), or
#   - an ad-hoc spec "NAME:RAW_DIR:MANIFEST:OUT_DIR" for a cohort not
#     (yet) hardcoded above, so new cohorts can be run without editing this file.

requested <- commandArgs(trailingOnly = TRUE)

if (length(requested) == 0) {
  cohorts <- cohorts_all
} else {
  cohorts <- list()
  for (arg in requested) {
    if (grepl(":", arg, fixed = TRUE)) {
      parts <- strsplit(arg, ":", fixed = TRUE)[[1]]
      if (length(parts) != 4) {
        stop("Ad-hoc cohort spec must be 'NAME:RAW_DIR:MANIFEST:OUT_DIR', got: ", arg)
      }
      cohorts[[parts[1]]] <- list(raw_dir = parts[2], manifest = parts[3], out_dir = parts[4])
    } else if (arg %in% names(cohorts_all)) {
      cohorts[[arg]] <- cohorts_all[[arg]]
    } else {
      stop("Unknown cohort '", arg, "'. Valid built-in options: ",
           paste(names(cohorts_all), collapse = ", "),
           ". To run a cohort not listed above, pass 'NAME:RAW_DIR:MANIFEST:OUT_DIR'.")
    }
  }
}

rg_sets <- list()

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

  rg <- read_cohort_idats(idat_dir = cfg$raw_dir, sample_sheet = manifest)

  cat(cohort_name, "RGChannelSet:", ncol(rg), "samples,", nrow(rg), "probes\n")

  dir.create(cfg$out_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(rg, file.path(cfg$out_dir, paste0(cohort_name, '.rds')))
  rg_sets[[cohort_name]] <- rg
}
