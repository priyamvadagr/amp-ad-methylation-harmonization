# ============================================================================
# Run Zhou-mask probe filtering (R/filter_probes.R) on each cohort's noob
# output (normalize_data.R), per Lundin et al. Artherosclerosis(2025) and
# Methylation(2024).
#
# Reads <normalized_dir>/<cohort>_noob_gmset.rds and _noob_beta.rds, applies
# filter_probes() (low call rate -> cross-reactive -> SNP artefact ->
# non-autosomal), and writes <normalized_dir>/<cohort>_filtered_gmset.rds,
# _filtered_beta.rds, _probe_filter_summary.csv.
#
# The Zhou mask manifests depend on the array platform -- HM450 for ROSMAP,
# EPIC for everything else (see PLATFORM below). Two SEPARATE manifests per
# platform are used (never merged -- see R/filter_probes.R's file header):
#   mapping mask  base ".rds" GRanges manifest -> drop_crossreactive()
#   snp mask      ".pop" tsv manifest          -> drop_snp()
# Downloaded from http://zwdzwd.github.io/InfiniumAnnotation (hg19, to match
# ilmn12.hg19 / ilm10b4.hg19) and cached in mask_dir; each platform's pair of
# masks is loaded at most once even if several cohorts share it.
#
# Usage:
#   Rscript filter_probes_per_cohort.R                  # all cohorts
#   Rscript filter_probes_per_cohort.R MSBB ROSMAP       # named cohorts only
#   Rscript filter_probes_per_cohort.R callrate=0.95 MOA-PAD
#   Rscript filter_probes_per_cohort.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir:/path/normalized_dir"
#                                          # ad-hoc cohort not in cohort_config.R;
#                                          # platform defaults to EPIC -- see
#                                          # PLATFORM below to add it explicitly.
#
# Cohort paths (normalized_dir) live in R/cohort_config.R, shared with
# compile_idat_files.R / sample_level_QC.R / normalize_data.R so the pipeline
# stages can't drift apart.
# ============================================================================

library(minfi)
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/filter_probes.R')
source('/home/ec2-user/AMP-AD_methylation_harmonization/Scripts/functions/cohort_config.R')

mask_dir <- '/home/ec2-user/data/methyl_harmonization/Zhou_mask_tables'

# cohort -> array platform, so the matching Zhou manifest gets loaded.
# Cohorts not listed here (e.g. a new ad-hoc cohort) default to EPIC.
PLATFORM <- list(
  ROSMAP       = "450k",
  MSBB         = "epic",
  ROSMAP_APOE4 = "epic",
  `MOA-PAD`    = "epic"
)

mapping_mask_path <- list(
  "450k" = file.path(mask_dir, "HM450.hg19.manifest.rds"),
  "epic" = file.path(mask_dir, "EPIC.hg19.manifest.rds")
)
snp_mask_path <- list(
  "450k" = file.path(mask_dir, "HM450.hg19.manifest.pop.tsv.gz"),
  "epic" = file.path(mask_dir, "EPIC.hg19.manifest.pop.tsv.gz")
)

# --- parse args: callrate=, cohort names, ad-hoc specs ----------------------
args         <- commandArgs(trailingOnly = TRUE)
min_callrate <- 0.99
requested    <- character(0)

for (arg in args) {
  if (grepl("^callrate=", arg)) {
    min_callrate <- as.numeric(sub("^callrate=", "", arg))
  } else {
    requested <- c(requested, arg)
  }
}

cohorts <- select_cohorts(requested)

# --- load each platform's masks at most once --------------------------------
mapping_masks <- list()
get_mapping_mask <- function(platform) {
  if (is.null(mapping_masks[[platform]])) {
    message("[mask] loading ", platform, " mapping (Zhou base) manifest ...")
    mapping_masks[[platform]] <<- readRDS(mapping_mask_path[[platform]])
  }
  mapping_masks[[platform]]
}

snp_masks <- list()
get_snp_mask <- function(platform) {
  if (is.null(snp_masks[[platform]])) {
    message("[mask] loading ", platform, " SNP (Zhou .pop) manifest ...")
    snp_masks[[platform]] <<- load_zhou_mask(snp_mask_path[[platform]])
  }
  snp_masks[[platform]]
}

# --- run ---------------------------------------------------------------------
for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  cat("\n==== ", cohort_name, " ====\n")

  platform <- PLATFORM[[cohort_name]]
  if (is.null(platform)) {
    platform <- "epic"
    message("[", cohort_name, "] no PLATFORM entry; defaulting to EPIC mask.")
  }
  mapping_mask <- get_mapping_mask(platform)
  snp_mask     <- get_snp_mask(platform)

  gmset <- readRDS(file.path(cfg$normalized_dir, paste0(cohort_name, "_noob_gmset.rds")))
  beta  <- readRDS(file.path(cfg$normalized_dir, paste0(cohort_name, "_noob_beta.rds")))

  res <- filter_probes(gmset, beta = beta,
                        mapping_mask_tbl = mapping_mask, snp_mask_tbl = snp_mask,
                        min_callrate = min_callrate, autosomal_only = TRUE,
                        out_prefix = file.path(cfg$normalized_dir, cohort_name))

  print(res$summary)
  rm(gmset, beta, res); gc()
}
