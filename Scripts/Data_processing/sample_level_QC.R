# ============================================================================
# Run sample-level QC (detection p-values, intensity, predicted-vs-reported
# sex, bisulfite conversion, control probes) on each cohort's compiled
# RGChannelSet, as written by compile_idat_files.R to
# <cohort>/combined_idat/<cohort>.rds. Reported sex is merged in from
# <cohort>/metadata/<cohort>_individual_metadata_processed.csv (not carried by
# the idat manifest used to compile the RGChannelSet), matched on
# individualID.
#
# The bisulfite-conversion and control-probe checks (ewastools) need the raw
# idat files, not just the compiled RGChannelSet, so each cohort also supplies
# its idat manifest (<cohort>/metadata/<cohort>_idat_manifest.csv, the same
# file compile_idat_files.R used) and the raw idat directory
# (methyl_data/<cohort>/raw_data). Basenames are built from the manifest's
# grnFile column, same as compile_idat_files.R -- this also covers ROSMAP,
# whose idats sit in per-Sentrix-barcode subdirectories under raw_data/,
# because grnFile there already includes that subfolder prefix (e.g.
# "5815381027/5815381027_R06C01_Grn.idat").
#
# Usage:
#   Rscript sample_level_QC.R          # process all cohorts
#   Rscript sample_level_QC.R MSBB ROSMAP
#                                       # process only the named cohort(s)
#   Rscript sample_level_QC.R \
#     "NewCohort:/path/to/rg.rds:/path/individual_metadata.csv:/path/idat_manifest.csv:/path/raw_data:/path/out_prefix"
#                                       # run a cohort not hardcoded below
# ============================================================================

library(minfi)
source('/home/ec2-user/AMP-AD_methylation_harmonization/R/sample_qc.R')

# Attach reported sex onto pData(rg) from the cohort's processed individual
# metadata, matched on individualID.
attach_sex <- function(rg, individual_csv) {
  ind <- read.csv(individual_csv, stringsAsFactors = FALSE)
  pData(rg)$sex <- ind$sex[match(pData(rg)$individualID, ind$individualID)]
  rg
}

# Where render_qc_plots() writes its per-cohort figures (in a subfolder named
# after the cohort).
plot_root <- '/home/ec2-user/AMP-AD_methylation_harmonization/Results/QC/raw_values'

cohorts_all <- list(
  MSBB = list(
    rds_in         = '/home/ec2-user/data/methyl_harmonization/MSBB/combined_idat/MSBB.rds',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/MSBB/metadata/MSBB_individual_metadata_processed.csv',
    idat_manifest  = '/home/ec2-user/data/methyl_harmonization/MSBB/metadata/MSBB_idat_manifest.csv',
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/MSBB/raw_data',
    out_prefix     = '/home/ec2-user/data/methyl_harmonization/MSBB/qc/MSBB'
  ),
  ROSMAP = list(
    rds_in         = '/home/ec2-user/data/methyl_harmonization/ROSMAP/combined_idat/ROSMAP.rds',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/ROSMAP/metadata/ROSMAP_individual_metadata_processed.csv',
    idat_manifest  = '/home/ec2-user/data/methyl_harmonization/ROSMAP/metadata/ROSMAP_idat_manifest.csv',
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/ROSMAP/raw_data',
    out_prefix     = '/home/ec2-user/data/methyl_harmonization/ROSMAP/qc/ROSMAP'
  ),
  ROSMAP_APOE4 = list(
    rds_in         = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/combined_idat/ROSMAP_APOE4.rds',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/metadata/ROSMAP_APOE4_individual_metadata_processed.csv',
    idat_manifest  = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/metadata/ROSMAP_APOE4_idat_manifest.csv',
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/ROSMAP_APOE4/raw_data',
    out_prefix     = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/qc/ROSMAP_APOE4'
  ),
  `MOA-PAD` = list(
    rds_in         = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/combined_idat/MOA-PAD.rds',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/metadata/MOA-PAD_individual_metadata_processed.csv',
    idat_manifest  = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/metadata/MOA-PAD_idat_manifest.csv',
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/MOA-PAD/raw_data',
    out_prefix     = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/qc/MOA-PAD'
  )
)

# Select which cohort(s) to process from the command line; default to all
# built-in cohorts. Each arg is either:
#   - a name from cohorts_all (e.g. "MSBB"), or
#   - an ad-hoc spec "NAME:RDS_IN:INDIVIDUAL_CSV:IDAT_MANIFEST:RAW_DIR:OUT_PREFIX"
#     for a cohort not (yet) hardcoded above, so new cohorts can be run
#     without editing this file.

requested <- commandArgs(trailingOnly = TRUE)

if (length(requested) == 0) {
  cohorts <- cohorts_all
} else {
  cohorts <- list()
  for (arg in requested) {
    if (grepl(":", arg, fixed = TRUE)) {
      parts <- strsplit(arg, ":", fixed = TRUE)[[1]]
      if (length(parts) != 6) {
        stop("Ad-hoc cohort spec must be ",
             "'NAME:RDS_IN:INDIVIDUAL_CSV:IDAT_MANIFEST:RAW_DIR:OUT_PREFIX', got: ", arg)
      }
      cohorts[[parts[1]]] <- list(rds_in = parts[2], individual_csv = parts[3],
                                   idat_manifest = parts[4], raw_dir = parts[5],
                                   out_prefix = parts[6])
    } else if (arg %in% names(cohorts_all)) {
      cohorts[[arg]] <- cohorts_all[[arg]]
    } else {
      stop("Unknown cohort '", arg, "'. Valid built-in options: ",
           paste(names(cohorts_all), collapse = ", "),
           ". To run a cohort not listed above, pass ",
           "'NAME:RDS_IN:INDIVIDUAL_CSV:IDAT_MANIFEST:RAW_DIR:OUT_PREFIX'.")
    }
  }
}

qc_results <- list()

for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  cat("\n==== ", cohort_name, " ====\n")

  rg <- readRDS(cfg$rds_in)
  rg <- attach_sex(rg, cfg$individual_csv)

  idat_manifest <- read.csv(cfg$idat_manifest, stringsAsFactors = FALSE)

  dir.create(dirname(cfg$out_prefix), recursive = TRUE, showWarnings = FALSE)
  res <- sample_qc(rg, targets = idat_manifest, raw_dir = cfg$raw_dir, grn_col = "grnFile",
                   sex_col = "sex", out_prefix = cfg$out_prefix)

  cat(cohort_name, "sample QC:", nrow(res$qc_tab), "samples,",
      sum(res$qc_tab$any_flag), "flagged (nothing dropped)\n")

  render_qc_plots(res$rg, cohort = cohort_name, out_dir = plot_root, qc_tab = res$qc_tab)

  saveRDS(res$rg, paste0(cfg$out_prefix, "_qc.rds"))
  qc_results[[cohort_name]] <- res

  rm(rg, res, idat_manifest)
  gc()
}
