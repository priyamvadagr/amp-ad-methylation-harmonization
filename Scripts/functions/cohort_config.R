# ============================================================================
# Canonical per-cohort file locations, shared by compile_idat_files.R,
# sample_level_QC.R and normalize_data.R so raw_dir / idat-manifest paths
# can't drift between pipeline stages (each used to hardcode its own copy).
# Add a new cohort here once; all scripts pick it up.
#
#   raw_dir        directory of raw *_Grn.idat / *_Red.idat files
#   manifest       idat manifest CSV (specimenID/individualID <-> Grn/Red idat),
#                  written by Scripts/Misc/get_sample_info.Rmd
#   combined_dir   where compile_idat_files.R saves <cohort>.rds
#   individual_csv processed individual metadata CSV (used by sample_level_QC.R
#                  to attach reported sex)
#   qc_dir         where sample_level_QC.R writes <cohort>_ewastools.csv,
#                  <cohort>_sample_qc.csv, <cohort>_qc.rds (under qc_dir/sample_qc/)
#                  and pca_qc_per_cohort.R writes <cohort>_pca.rds and friends
#                  (under qc_dir/pca_qc/)
#   normalized_dir where normalize_data.R writes <cohort>_noob_gmset.rds,
#                  <cohort>_noob_beta.rds
# ============================================================================

cohorts_all <- list(
  MSBB = list(
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/MSBB/raw_data',
    manifest       = '/home/ec2-user/data/methyl_harmonization/MSBB/metadata/MSBB_idat_manifest.csv',
    combined_dir   = '/home/ec2-user/data/methyl_harmonization/MSBB/combined_idat',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/MSBB/metadata/MSBB_individual_metadata_processed.csv',
    qc_dir         = '/home/ec2-user/data/methyl_harmonization/MSBB/qc',
    normalized_dir = '/home/ec2-user/data/methyl_harmonization/MSBB/normalized'
  ),
  ROSMAP = list(
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/ROSMAP/raw_data',
    manifest       = '/home/ec2-user/data/methyl_harmonization/ROSMAP/metadata/ROSMAP_idat_manifest.csv',
    combined_dir   = '/home/ec2-user/data/methyl_harmonization/ROSMAP/combined_idat',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/ROSMAP/metadata/ROSMAP_individual_metadata_processed.csv',
    qc_dir         = '/home/ec2-user/data/methyl_harmonization/ROSMAP/qc',
    normalized_dir = '/home/ec2-user/data/methyl_harmonization/ROSMAP/normalized'
  ),
  ROSMAP_APOE4 = list(
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/ROSMAP_APOE4/raw_data',
    manifest       = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/metadata/ROSMAP_APOE4_idat_manifest.csv',
    combined_dir   = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/combined_idat',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/metadata/ROSMAP_APOE4_individual_metadata_processed.csv',
    qc_dir         = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/qc',
    normalized_dir = '/home/ec2-user/data/methyl_harmonization/ROSMAP_APOE4/normalized'
  ),
  `MOA-PAD` = list(
    raw_dir        = '/home/ec2-user/data/methyl_harmonization/methyl_data/MOA-PAD/raw_data',
    manifest       = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/metadata/MOA-PAD_idat_manifest.csv',
    combined_dir   = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/combined_idat',
    individual_csv = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/metadata/MOA-PAD_individual_metadata_processed.csv',
    qc_dir         = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/qc',
    normalized_dir = '/home/ec2-user/data/methyl_harmonization/MOA-PAD/normalized'
  )
)

# ---------------------------------------------------------------------------
# select_cohorts(): pick cohort configs from command-line args, shared by
# compile_idat_files.R and sample_level_QC.R.
#
#   requested  character vector of args (cohort names and/or ad-hoc specs)
#   all_cfg    named list of built-in cohort configs (cohorts_all above)
#
# Each element of `requested` is either a name in `all_cfg`, or an ad-hoc
# "NAME:RAW_DIR:MANIFEST:COMBINED_DIR:INDIVIDUAL_CSV:QC_DIR:NORMALIZED_DIR"
# spec for a cohort not (yet) hardcoded in cohorts_all, so new cohorts can be
# run without editing this file. Empty `requested` selects all of `all_cfg`.
# ---------------------------------------------------------------------------
select_cohorts <- function(requested, all_cfg = cohorts_all) {
  if (length(requested) == 0) return(all_cfg)

  fields <- c("raw_dir", "manifest", "combined_dir", "individual_csv", "qc_dir",
              "normalized_dir")
  cohorts <- list()
  for (arg in requested) {
    if (grepl(":", arg, fixed = TRUE)) {
      parts <- strsplit(arg, ":", fixed = TRUE)[[1]]
      if (length(parts) != length(fields) + 1) {
        stop("Ad-hoc cohort spec must be 'NAME:", paste(toupper(fields), collapse = ":"),
             "', got: ", arg)
      }
      cfg <- as.list(parts[-1])
      names(cfg) <- fields
      cohorts[[parts[1]]] <- cfg
    } else if (arg %in% names(all_cfg)) {
      cohorts[[arg]] <- all_cfg[[arg]]
    } else {
      stop("Unknown cohort '", arg, "'. Valid built-in options: ",
           paste(names(all_cfg), collapse = ", "),
           ". To run a cohort not listed above, pass ",
           "'NAME:", paste(toupper(fields), collapse = ":"), "'.")
    }
  }
  cohorts
}
