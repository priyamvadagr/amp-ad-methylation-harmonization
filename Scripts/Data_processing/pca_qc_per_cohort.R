# ============================================================================
# Run pca_qc() (R/pca_qc.R) on each cohort's filtered beta matrix (output of
# filter_probes_per_cohort.R), testing PCs against technical and biological
# covariates assembled from that cohort's idat manifest + individual metadata.
#
# Reads <normalized_dir>/<cohort>_filtered_beta.rds and builds one pheno row
# per sample (matching colnames(beta), i.e. "<Sentrix_ID>_<Sentrix_Row_Column>")
# by:
#   1. manifest.csv       -> Sentrix_ID, array_row/array_col (split out of
#                            Sentrix_Row_Column, e.g. "R06C01" -> "R06"/"C01")
#   2. individual_csv      -> sex, race, ageDeath, Braak, apoe4Status, tissue
#                             (matched on individualID; columns common to all
#                             cohorts -- see COVARIATES below)
# and writes <qc_dir>/<cohort>_pca.rds, _pca_scores.csv, _pca_varexplained.csv,
# _pca_assoc.csv, _pca_qc.pdf.
#
# Usage:
#   Rscript pca_qc_per_cohort.R                    # all cohorts
#   Rscript pca_qc_per_cohort.R MSBB ROSMAP        # named cohorts only
#   Rscript pca_qc_per_cohort.R ntop=10000 npcs=15 MOA-PAD
#   Rscript pca_qc_per_cohort.R \
#     "NewCohort:/raw/dir:/path/manifest.csv:/path/combined_dir:/path/individual.csv:/path/qc_dir:/path/normalized_dir"
#
# Cohort paths (manifest, individual_csv, qc_dir, normalized_dir) live in
# R/cohort_config.R, shared with the rest of the pipeline so stages can't
# drift apart.
# ============================================================================

source('/home/ec2-user/AMP-AD_methylation_harmonization/R/pca_qc.R')
source('/home/ec2-user/AMP-AD_methylation_harmonization/R/cohort_config.R')

# Where pca_qc() ALSO copies each cohort's _pca_qc.pdf (data files stay in
# qc_dir); mirrors render_qc_plots()'s Results/QC/<subfolder>/<cohort>/ layout.
plot_root <- '/home/ec2-user/AMP-AD_methylation_harmonization/Results/QC/PCA'

# technical (from the idat manifest) + biological (from individual metadata)
# covariates common to every cohort's individual_csv -- see cohort_config.R.
TECH_COVARIATES <- c("Sentrix_ID", "array_row", "array_col")
BIO_COVARIATES  <- c("sex", "race", "ageDeath", "Braak", "apoe4Status", "tissue")
COLOR_BY        <- c("Sentrix_ID", "sex", "tissue", "race", "array_row", "array_col")

# Build one pheno row per sample (keyed on the same
# "<Sentrix_ID>_<Sentrix_Row_Column>" ID minfi uses as colnames(beta)).
build_pheno <- function(cfg) {
  manifest <- read.csv(cfg$manifest, stringsAsFactors = FALSE)
  manifest$sample     <- paste0(manifest$Sentrix_ID, "_", manifest$Sentrix_Row_Column)
  manifest$array_row  <- sub("C[0-9]+$", "", manifest$Sentrix_Row_Column)
  manifest$array_col  <- sub("^R[0-9]+", "", manifest$Sentrix_Row_Column)

  individual <- read.csv(cfg$individual_csv, stringsAsFactors = FALSE)

  pheno <- merge(manifest, individual, by = "individualID", all.x = TRUE)
  rownames(pheno) <- pheno$sample
  pheno
}

# --- parse args: ntop=, npcs=, r2flag=, cohort names, ad-hoc specs ----------
args      <- commandArgs(trailingOnly = TRUE)
n_top     <- 20000
n_pcs     <- 10
r2_flag   <- 0.05
requested <- character(0)

for (arg in args) {
  if (grepl("^ntop=", arg)) {
    n_top <- as.integer(sub("^ntop=", "", arg))
  } else if (grepl("^npcs=", arg)) {
    n_pcs <- as.integer(sub("^npcs=", "", arg))
  } else if (grepl("^r2flag=", arg)) {
    r2_flag <- as.numeric(sub("^r2flag=", "", arg))
  } else {
    requested <- c(requested, arg)
  }
}

cohorts <- select_cohorts(requested)

# --- run ---------------------------------------------------------------------
for (cohort_name in names(cohorts)) {
  cfg <- cohorts[[cohort_name]]
  cat("\n==== ", cohort_name, " ====\n")

  beta  <- readRDS(file.path(cfg$normalized_dir, paste0(cohort_name, "_filtered_beta.rds")))
  pheno <- build_pheno(cfg)

  covariates <- intersect(c(TECH_COVARIATES, BIO_COVARIATES), names(pheno))
  missing    <- setdiff(c(TECH_COVARIATES, BIO_COVARIATES), names(pheno))
  if (length(missing))
    message("[", cohort_name, "] covariate(s) not found, skipping: ",
            paste(missing, collapse = ", "))

  # skip covariates that are "missing or unknown" (the metadata's NA
  # placeholder) for every sample -- e.g. Braak/apoe4Status in MOA-PAD --
  # lm() can't associate a PC with a column that has no variation anyway.
  all_unknown <- vapply(covariates, function(cv) {
    x <- pheno[[cv]]
    all(is.na(x) | x == "missing or unknown")
  }, logical(1))
  if (any(all_unknown))
    message("[", cohort_name, "] covariate(s) all 'missing or unknown', skipping: ",
            paste(covariates[all_unknown], collapse = ", "))
  covariates <- covariates[!all_unknown]

  color_by <- intersect(COLOR_BY, covariates)

  dir.create(cfg$qc_dir, recursive = TRUE, showWarnings = FALSE)
  res <- pca_qc(beta, pheno,
               sample_col   = "sample",
               covariates   = covariates,
               color_by     = color_by,
               n_top        = n_top,
               n_pcs        = n_pcs,
               r2_flag      = r2_flag,
               out_prefix   = file.path(cfg$qc_dir, cohort_name),
               plot_dir     = plot_root)

  print(res$assoc_r2)
  rm(beta, pheno, res); gc()
}
