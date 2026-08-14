# ============================================================================
# noob normalization for one cohort (minfi). Run AFTER sample_qc().
#
# noob = Normal-exponential out-of-band background correction + dye-bias
# correction. It is a WITHIN-SAMPLE method (each array normalized on its own).
#
# This version is MEMORY-CONSCIOUS: on large cohorts (≈1000 EPIC samples) each
# probe x sample matrix is several GB, so the function frees every large object
# the instant it is no longer needed (rm + gc between steps), writes each
# matrix to disk one at a time, and by default returns only file paths + counts
# rather than the matrices (so the caller doesn't retain them across cohorts).
# Peak memory is ~2 large matrices instead of ~5.
#
# normalize_noob():
#   1. drops samples failing the chosen QC criteria (default: detp, intensity,
#      sex, bisulfite -- control-only samples KEPT), using qc_tab flags;
#   2. computes detP (while rg is alive) for masking, if needed;
#   3. preprocessNoob() -> MethylSet, then frees rg;
#   4. mapToGenome() -> GenomicMethylSet, then frees the MethylSet;
#   5. saves the GenomicMethylSet, then derives + masks + saves the beta matrix,
#      freeing the GenomicMethylSet and detP as soon as possible;
#   6. optionally derives M-values from the (masked) beta and saves them.
#
#   rg            QC'd RGChannelSet / RGChannelSetExtended
#   qc_tab        data.frame from sample_qc() (or its _sample_qc.csv). If NULL,
#                 no samples are dropped.
#   drop_reasons  flags that force a drop; subset of
#                 c("detp","intensity","sex","bisulfite","control","bead")
#   detP          optional precomputed detection-p matrix (res$detP from
#                 sample_qc()); recomputed from rg if NULL and mask_detp = TRUE
#   detp_cut      per-probe detection-p threshold for masking
#   mask_detp     TRUE = set beta (and M) to NA where detP > detp_cut
#   dye_method    preprocessNoob dyeMethod: "single" (per-sample, default) or
#                 "reference"
#   out_prefix    output path prefix
#   save_gmset    TRUE = write <out_prefix>_noob_gmset.rds (the biggest file;
#                 set FALSE if you only need the beta/M matrices)
#   save_mvals    TRUE = also write <out_prefix>_noob_mvals.rds. M is derived
#                 from the masked beta as log2(beta/(1-beta)) to avoid holding
#                 a second GenomicMethylSet in memory.
#   return_objects TRUE = also return the beta/gmset objects (old behavior;
#                 memory-heavy). Default FALSE returns only paths + counts.
#
# Returns (invisibly) a small list: paths, n_samples, n_probes, n_masked
# [+ gmset/beta/mvals if return_objects = TRUE].
# ============================================================================

library(minfi)

normalize_noob <- function(rg,
                           qc_tab         = NULL,
                           drop_reasons   = c("detp", "intensity", "sex", "bisulfite"),
                           detP           = NULL,
                           detp_cut       = 0.01,
                           mask_detp      = TRUE,
                           dye_method     = "single",
                           out_prefix     = "cohort",
                           save_gmset     = TRUE,
                           save_mvals     = FALSE,
                           return_objects = FALSE) {

  stopifnot(is(rg, "RGChannelSet"))
  tag <- basename(out_prefix)

  ## --- 1. drop failed samples (control-only kept) ------------------------
  if (!is.null(qc_tab)) {
    flag_col <- c(detp = "fail_detp", intensity = "fail_intensity",
                  sex = "sex_mismatch", bisulfite = "fail_bisulfite",
                  control = "fail_control", bead = "fail_bead")
    use <- flag_col[drop_reasons]
    if (anyNA(use)) stop("drop_reasons must be in: ", paste(names(flag_col), collapse = ", "))
    drop <- Reduce(`|`, lapply(use, function(cc) qc_tab[[cc]] %in% TRUE))
    keep_ids <- qc_tab$sample[!drop]
    rg <- rg[, colnames(rg) %in% keep_ids]
    message(sprintf("[%s] dropped %d sample(s) on {%s}; %d retained.",
                    tag, sum(drop), paste(drop_reasons, collapse = ","), ncol(rg)))
  }
  n_samples <- ncol(rg)

  ## --- 2. detP now, while rg is alive (only if we'll mask) ---------------
  if (mask_detp && is.null(detP)) {
    message(sprintf("[%s] computing detection p-values ...", tag))
    detP <- minfi::detectionP(rg)
  }

  ## --- 3. noob normalization, then free rg -------------------------------
  message(sprintf("[%s] preprocessNoob (dyeMethod = %s) on %d samples ...",
                  tag, dye_method, n_samples))
  mset <- preprocessNoob(rg, dyeCorr = TRUE, dyeMethod = dye_method)
  rm(rg); gc()

  ## --- 4. map to genome, then free the MethylSet -------------------------
  gmset <- mapToGenome(mset)
  rm(mset); gc()

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  paths <- list()

  ## --- 5a. save the GenomicMethylSet (largest object) --------------------
  if (save_gmset) {
    p <- paste0(out_prefix, "_noob_gmset.rds")
    saveRDS(gmset, p); paths$gmset <- p
    message(sprintf("[%s] wrote %s", tag, basename(p)))
  }

  ## --- 5b. derive beta, free gmset (unless returning it), mask, save -----
  beta <- getBeta(gmset)
  n_probes <- nrow(beta)
  if (!return_objects) { rm(gmset); gc() }

  n_masked <- 0L
  if (mask_detp) {
    detP <- detP[rownames(beta), colnames(beta), drop = FALSE]   # align to beta
    mask <- detP > detp_cut
    mask[is.na(mask)] <- FALSE
    rm(detP); gc()
    beta[mask] <- NA
    n_masked <- sum(mask)
    message(sprintf("[%s] masked %d probe x sample values (detP > %.3f); %.3f%% of matrix.",
                    tag, n_masked, detp_cut, 100 * n_masked / length(beta)))
    rm(mask); gc()
  }

  p <- paste0(out_prefix, "_noob_beta.rds")
  saveRDS(beta, p); paths$beta <- p
  message(sprintf("[%s] wrote %s", tag, basename(p)))

  ## --- 6. optional M-values, derived from the masked beta ----------------
  if (save_mvals) {
    mvals <- log2(beta / (1 - beta))          # NAs in beta propagate automatically
    p <- paste0(out_prefix, "_noob_mvals.rds")
    saveRDS(mvals, p); paths$mvals <- p
    message(sprintf("[%s] wrote %s", tag, basename(p)))
    if (!return_objects) rm(mvals)
  }

  message(sprintf("[%s] done: %d samples x %d probes.", tag, n_samples, n_probes))

  out <- list(paths = paths, n_samples = n_samples,
              n_probes = n_probes, n_masked = n_masked)
  if (return_objects) {
    out$beta  <- beta
    if (exists("gmset")) out$gmset <- gmset
    if (save_mvals)      out$mvals <- log2(beta / (1 - beta))
  } else {
    rm(beta); gc()
  }
  invisible(out)
}

