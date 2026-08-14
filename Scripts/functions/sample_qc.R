# ============================================================================
# Sample-level QC for one cohort's RGChannelSet (minfi + ewastools).
# Run BEFORE normalization. This file defines three functions:
#
#   ewastools_qc()     Reads raw idats (in memory-bounded chunks) and computes
#                       the ewastools control-probe QC: bisulfite conversion +
#                       the 17-metric sample_failure(). Returns / writes a small
#                       per-sample table. Run this SEPARATELY (ideally in its
#                       own R process) from sample_qc() so the raw idats and the
#                       RGChannelSet are never both resident in memory -- this
#                       is what was getting the job killed on large cohorts.
#
#   sample_qc()        Computes the minfi per-sample QC flags (detection p,
#                       intensity, sex mismatch) from the RGChannelSet, MERGES
#                       in the ewastools flags from ewastools_qc() (if given),
#                       writes them to a CSV, and copies them into pData(rg).
#                       FLAGS every problem; drops NOTHING.
#
#   render_qc_plots()   Renders the diagnostic figures and, given qc_tab, a
#                       control-probe report restricted to flagged samples.
#
# Typical flow (two steps, ideally two Rscript calls):
#   ew  <- ewastools_qc(targets = manifest, raw_dir = cfg$raw_dir,
#                       out_csv = paste0(cfg$out_prefix, "_ewastools.csv"))
#   res <- sample_qc(rg, sex_col = "sex", ewastools_tab = ew,
#                    out_prefix = cfg$out_prefix)
#   render_qc_plots(res$rg, cohort = "MSBB", out_dir = "…", qc_tab = res$qc_tab)
#
# sample_qc() flags (all flag-only; none force removal):
#   fail_detp      - too many probes with signal indistinguishable from noise
#   fail_intensity - overall array signal too low
#   sex_mismatch   - predicted sex disagrees with recorded sex
#   fail_bisulfite - bisulfite conversion control below threshold   (ewastools)
#   fail_control   - ewastools sample_failure() over the 17 controls (ewastools)
#   fail_beadcount - too many probes measured by < min_beadcount beads
#                    (requires rg to be an RGChannelSetExtended, i.e. compiled
#                    with extended = TRUE in read_cohort_idats(); NA otherwise)
#
# `flag_reason` records every triggered check per sample ("none" if clean).
# ============================================================================

library(minfi)
library(ewastools)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ---------------------------------------------------------------------------
# idat_has_valid_magic(): TRUE if `path` starts with the IDAT binary magic
# header. Mirrors the helper of the same name in read_cohort_idats.R.
# ---------------------------------------------------------------------------
idat_has_valid_magic <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- tryCatch(readBin(con, "raw", n = 4), error = function(e) raw(0))
  length(magic) == 4 && rawToChar(magic) == "IDAT"
}

# ---------------------------------------------------------------------------
# resolve_basenames(): build idat path prefixes (no _Grn/_Red.idat) from either
# an explicit `basenames` vector or a manifest `targets` (+ raw_dir/grn_col).
# ---------------------------------------------------------------------------
resolve_basenames <- function(basenames = NULL, targets = NULL,
                              raw_dir = NULL, grn_col = "grnFile") {
  if (!is.null(basenames)) return(basenames)
  if (is.null(targets)) return(NULL)
  tgm <- as.data.frame(targets, stringsAsFactors = FALSE)
  if ("Basename" %in% names(tgm)) return(tgm$Basename)
  if (!is.null(raw_dir) && grn_col %in% names(tgm)) {
    tgm <- tgm[!is.na(tgm[[grn_col]]) & nzchar(tgm[[grn_col]]), ]
    return(file.path(raw_dir, sub("_Grn\\.idat$", "", tgm[[grn_col]])))
  }
  NULL
}

# ---------------------------------------------------------------------------
# safe_sample_failure(): wraps ewastools::sample_failure(), which errors on a
# single-sample chunk -- sapply() over the 17 threshold-metrics collapses to a
# plain vector instead of a 1-row matrix when each metric evaluates to length
# 1, so its internal apply(failed, 1, any, na.rm = TRUE) has no dim to iterate
# over ("dim(X) must have a positive length"). Reimplements the same logic
# with an explicit reshape so it works regardless of chunk size.
# ---------------------------------------------------------------------------
safe_sample_failure <- function(metrics) {
  failed <- sapply(metrics, function(metric) metric < attr(metric, "threshold"))
  if (is.null(dim(failed)))
    failed <- matrix(failed, nrow = length(metrics[[1]]), dimnames = list(NULL, names(metrics)))
  apply(failed, 1, any, na.rm = TRUE)
}

# ---------------------------------------------------------------------------
# ewastools_qc(): control-probe QC (bisulfite conversion + sample_failure)
# from raw idats, processed in memory-bounded chunks.
#
#   basenames      idat path prefixes (no _Grn/_Red.idat). If NULL, derived
#                  from `targets` (+ raw_dir/grn_col).
#   targets        idat manifest data.frame (Basename col, or grn_col+raw_dir)
#   raw_dir        directory of raw idats (used with grn_col)
#   grn_col        manifest column with bare "*_Grn.idat" filenames
#   bisulfite_min  flag sample if any bisulfite-conversion metric < this
#   chunk_size     number of samples to read_idats() at once (bounds memory);
#                  NULL / Inf = all at once. control_metrics/sample_failure are
#                  per-sample, so chunking does not change the result.
#   check_magic    TRUE = skip idats with a bad magic header (flags -> NA)
#   out_csv        optional path to write the per-sample table
#
# Returns a data.frame: sample, bisulfite_min, fail_bisulfite, fail_control,
# corrupt_idat_file (sample = basename() of each idat prefix, matching
# colnames(rg)). Samples with a corrupt idat (bad magic header) are still
# included as a row -- corrupt_idat_file = TRUE and all other metrics NA,
# since they can't be read -- rather than silently dropped from the table.
# ---------------------------------------------------------------------------
ewastools_qc <- function(basenames     = NULL,
                         targets       = NULL,
                         raw_dir       = NULL,
                         grn_col       = "grnFile",
                         bisulfite_min = 1,
                         chunk_size    = 75,
                         check_magic   = TRUE,
                         out_csv       = NULL) {

  basenames <- resolve_basenames(basenames, targets, raw_dir, grn_col)
  if (is.null(basenames) || !length(basenames))
    stop("ewastools_qc(): no basenames; supply `basenames` or `targets` (+ raw_dir/grn_col).")

  corrupt_samples <- character(0)
  if (isTRUE(check_magic)) {
    ok <- vapply(paste0(basenames, "_Grn.idat"), idat_has_valid_magic, logical(1)) &
          vapply(paste0(basenames, "_Red.idat"), idat_has_valid_magic, logical(1))
    if (any(!ok)) {
      corrupt_samples <- basename(basenames[!ok])
      warning(sum(!ok), " sample(s) have a corrupt IDAT file (bad magic header); ",
              "flagged as corrupt_idat_file (other metrics NA):\n  ",
              paste(head(corrupt_samples, 10), collapse = "\n  "))
      basenames <- basenames[ok]
    }
  }
  if (!length(basenames)) stop("ewastools_qc(): no valid idats left after magic-header check.")

  if (is.null(chunk_size) || !is.finite(chunk_size) || chunk_size < 1)
    chunk_size <- length(basenames)
  groups <- split(seq_along(basenames),
                  ceiling(seq_along(basenames) / chunk_size))

  message(sprintf("[ewastools_qc] %d samples in %d chunk(s) of up to %d",
                  length(basenames), length(groups), chunk_size))

  parts <- vector("list", length(groups))
  for (i in seq_along(groups)) {
    bn   <- basenames[groups[[i]]]
    meth <- read_idats(bn)
    ctrl <- control_metrics(meth)
    failed <- safe_sample_failure(ctrl)

    bis_names <- grep("bisulfite", names(ctrl), ignore.case = TRUE, value = TRUE)
    if (length(bis_names)) {
      bis_mat <- do.call(cbind, lapply(bis_names, function(n) ctrl[[n]]))
      bis_min <- apply(bis_mat, 1, min, na.rm = TRUE)
    } else {
      if (i == 1L) warning("No bisulfite-conversion metric found in control_metrics().")
      bis_min <- rep(NA_real_, length(failed))
    }
    metrics_df <- as.data.frame(lapply(ctrl, round, 3))   # 17 columns, one row per sample
    parts[[i]] <- cbind(data.frame(sample = basename(bn),
             bisulfite_min     = round(bis_min, 2),
             fail_bisulfite    = bis_min < bisulfite_min,
             fail_control      = failed,
             corrupt_idat_file = FALSE,
             stringsAsFactors = FALSE),
            metrics_df)
    message(sprintf("  chunk %d/%d done (%d samples)", i, length(groups), length(bn)))
    rm(meth, ctrl); gc()
  }

  out <- do.call(rbind, parts)

  if (length(corrupt_samples)) {
    na_metrics <- as.data.frame(matrix(NA_real_,
                                       nrow = length(corrupt_samples),
                                       ncol = ncol(metrics_df),
                                       dimnames = list(NULL, names(metrics_df))))
    corrupt_rows <- cbind(data.frame(sample            = corrupt_samples,
             bisulfite_min     = NA_real_,
             fail_bisulfite    = NA,
             fail_control      = NA,
             corrupt_idat_file = TRUE,
             stringsAsFactors = FALSE),
            na_metrics)
    out <- rbind(out, corrupt_rows)
  }

  if (!is.null(out_csv)) {
    dir.create(dirname(out_csv), recursive = TRUE, showWarnings = FALSE)
    write.csv(out, out_csv, row.names = FALSE)
    message("[ewastools_qc] wrote ", out_csv)
  }
  out
}

# ---------------------------------------------------------------------------
# beadcount_from_rgset(): per-CpG functional bead counts from a compiled
# RGChannelSetExtended (i.e. one read with extended = TRUE in
# read_cohort_idats()). Type II probes use their single address; Type I
# probes use the minimum bead count across AddressA and AddressB (both bead
# types are required). Returns a CpG x sample matrix. Runs directly off the
# already-compiled .rds -- no separate re-read of the raw idats needed.
#
#   Example:
#     bc            <- beadcount_from_rgset(rg)          # CpG x samples
#     frac_low_bead <- colMeans(bc < 3, na.rm = TRUE)     # per-sample metric
# ---------------------------------------------------------------------------
beadcount_from_rgset <- function(rg) {
  stopifnot(is(rg, "RGChannelSetExtended"))
  nb  <- getNBeads(rg)                                 # addresses x samples
  rn  <- rownames(nb)
  tI  <- getProbeInfo(rg, type = "I")
  tII <- getProbeInfo(rg, type = "II")
  bcII <- nb[match(as.character(tII$AddressA), rn), , drop = FALSE]
  rownames(bcII) <- tII$Name
  bcI <- pmin(nb[match(as.character(tI$AddressA), rn), , drop = FALSE],
              nb[match(as.character(tI$AddressB), rn), , drop = FALSE])
  rownames(bcI) <- tI$Name
  rbind(bcI, bcII)
}

# ---------------------------------------------------------------------------
# sample_qc(): minfi per-sample QC flags + merge of ewastools flags.
#
# The RGChannelSet is obtained: `rg` if supplied > `rds` if it exists.
# sample_qc() no longer reads raw idats itself -- run ewastools_qc() separately
# and pass its result as `ewastools_tab` (or a CSV path via `ewastools_csv`).
#
#   rg             RGChannelSet (annotation set); if supplied, `rds` is ignored
#   rds            path to an existing combined-idat RGChannelSet .rds
#   ewastools_tab  data.frame from ewastools_qc() (sample, bisulfite_min,
#                  fail_bisulfite, fail_control). If NULL, those flags = NA.
#   ewastools_csv  alternative to ewastools_tab: path to the CSV it wrote
#   sex_col        pData column with reported sex ("female"/"male"/...)
#   detp           per-probe detection-p threshold (probe "fails" above this)
#   max_fail_frac  flag sample if > this fraction of probes fail detection
#   min_intensity  flag sample if median log2 (M+U)/2 intensity below this
#   min_beadcount        per-probe bead-count threshold (probe "fails" below this)
#   max_beadcount_frac   flag sample if > this fraction of probes fail bead-count.
#                        Requires rg to be an RGChannelSetExtended (extended =
#                        TRUE in read_cohort_idats()); fail_beadcount = NA otherwise.
#   out_prefix     prefix for the output CSV (plots are in render_qc_plots())
#
# Example:
#   ew  <- ewastools_qc(targets = manifest, raw_dir = cfg$raw_dir,
#                       out_csv = paste0(cfg$out_prefix, "_ewastools.csv"))
#   res <- sample_qc(rg, sex_col = "sex", ewastools_tab = ew,
#                    out_prefix = cfg$out_prefix)
#   res$qc_tab[res$qc_tab$any_flag, c("sample","flag_reason")]  # who & why
# ---------------------------------------------------------------------------
sample_qc <- function(rg                  = NULL,
                      rds                 = NULL,
                      ewastools_tab       = NULL,
                      ewastools_csv       = NULL,
                      sex_col             = "sex",
                      detp                = 0.01,
                      max_fail_frac       = 0.01,
                      min_intensity       = 10.5,
                      min_beadcount       = 3,
                      max_beadcount_frac  = 0.01,
                      out_prefix          = "cohort") {

  ## --- obtain the RGChannelSet: object > existing rds ---------------------
  if (is.null(rg)) {
    if (!is.null(rds) && file.exists(rds)) {
      message("[sample_qc] using existing combined-idat rds: ", rds)
      rg <- readRDS(rds)
    } else {
      stop("Provide `rg` or an existing `rds` (as written by compile_idat_files.R).")
    }
  }
  stopifnot(is(rg, "RGChannelSet"))
  pd <- as.data.frame(pData(rg))

  ## optional: load ewastools table from CSV -------------------------------
  if (is.null(ewastools_tab) && !is.null(ewastools_csv) && file.exists(ewastools_csv))
    ewastools_tab <- read.csv(ewastools_csv, stringsAsFactors = FALSE)

  ## detection p-values ----------------------------------------------------
  ## minfi:: qualified -- ewastools also exports detectionP() and masks minfi's.
  detP        <- minfi::detectionP(rg)
  frac_failed <- colMeans(detP > detp)
  mean_detP   <- colMeans(detP)

  ## intensity QC ----------------------------------------------------------
  mset <- preprocessRaw(rg)
  qc   <- getQC(mset)
  qc_intensity <- (qc$mMed + qc$uMed) / 2

  ## predicted sex vs reported ---------------------------------------------
  gmset         <- mapToGenome(mset)
  predicted_sex <- getSex(gmset)$predictedSex
  reported_sex  <- NA_character_
  if (sex_col %in% names(pd)) {
    r <- tolower(as.character(pd[[sex_col]]))
    reported_sex <- ifelse(r == "female", "F",
                    ifelse(r == "male",   "M", NA_character_))
  }
  sex_mismatch <- !is.na(reported_sex) & reported_sex != predicted_sex

  ## bead-count QC (requires RGChannelSetExtended; NA if compiled without
  ## extended = TRUE in read_cohort_idats()) ------------------------------
  if (is(rg, "RGChannelSetExtended")) {
    bc            <- beadcount_from_rgset(rg)
    frac_low_bead <- colMeans(bc < min_beadcount, na.rm = TRUE)[colnames(rg)]
  } else {
    message("[sample_qc] rg is not an RGChannelSetExtended; bead-count QC ",
            "skipped (fail_beadcount = NA). Re-compile with extended = TRUE ",
            "in read_cohort_idats() to enable it.")
    frac_low_bead <- rep(NA_real_, ncol(rg))
  }

  ## ewastools flags: merged from ewastools_qc() output, matched by sample --
  fail_bisulfite    <- rep(NA, ncol(rg))
  fail_control      <- rep(NA, ncol(rg))
  bisulfite_min_val <- rep(NA_real_, ncol(rg))
  if (!is.null(ewastools_tab)) {
    m <- match(colnames(rg), ewastools_tab$sample)
    if (anyNA(m))
      warning(sum(is.na(m)), " sample(s) in rg not found in `ewastools_tab`; ",
              "their ewastools flags stay NA.")
    bisulfite_min_val <- ewastools_tab$bisulfite_min[m]
    fail_bisulfite    <- ewastools_tab$fail_bisulfite[m]
    fail_control      <- ewastools_tab$fail_control[m]
  }

  ## Assemble QC table: ALL checks are flags -------------------------------
  qc_tab <- data.frame(
    sample         = colnames(rg),
    frac_failed    = round(frac_failed, 4),
    mean_detP      = signif(mean_detP, 3),
    qc_intensity   = round(qc_intensity, 2),
    predicted_sex  = predicted_sex,
    reported_sex   = reported_sex,
    bisulfite_min  = round(bisulfite_min_val, 2),
    frac_low_bead  = round(frac_low_bead, 4),
    fail_detp      = frac_failed  > max_fail_frac,
    fail_intensity = qc_intensity < min_intensity,
    sex_mismatch   = sex_mismatch,
    fail_bisulfite = fail_bisulfite,
    fail_control   = fail_control,
    fail_beadcount = frac_low_bead > max_beadcount_frac,
    stringsAsFactors = FALSE
  )

  ## flag_reason: every triggered check, semicolon-joined ("none" if clean)
  flag_mat <- cbind(
    detP      = qc_tab$fail_detp      %in% TRUE,
    intensity = qc_tab$fail_intensity %in% TRUE,
    sex       = qc_tab$sex_mismatch   %in% TRUE,
    bisulfite = qc_tab$fail_bisulfite %in% TRUE,
    control   = qc_tab$fail_control   %in% TRUE,
    beadcount = qc_tab$fail_beadcount %in% TRUE
  )
  qc_tab$flag_reason <- apply(flag_mat, 1, function(x) {
    r <- colnames(flag_mat)[x]; if (length(r)) paste(r, collapse = ";") else "none"
  })
  qc_tab$any_flag <- qc_tab$flag_reason != "none"

  ## copy flags into pData(rg) so they travel with the object --------------
  pData(rg)$qc_flag_reason <- qc_tab$flag_reason
  pData(rg)$qc_any_flag    <- qc_tab$any_flag

  ## report + return (NOTHING dropped) ------------------------------------
  message(sprintf(
    "[%s] %d samples | flagged: %d (detP %d, intensity %d, sex %d, bisulfite %d, control %d, beadcount %d)",
    out_prefix, nrow(qc_tab), sum(qc_tab$any_flag),
    sum(qc_tab$fail_detp %in% TRUE), sum(qc_tab$fail_intensity %in% TRUE),
    sum(qc_tab$sex_mismatch %in% TRUE), sum(qc_tab$fail_bisulfite %in% TRUE),
    sum(qc_tab$fail_control %in% TRUE), sum(qc_tab$fail_beadcount %in% TRUE)))

  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
  write.csv(qc_tab, paste0(out_prefix, "_sample_qc.csv"), row.names = FALSE)

  list(qc_tab = qc_tab,
       rg = rg,          # full object (flags in pData); nothing removed
       detP = detP,
       qc = qc)
}

# ============================================================================
# render_qc_plots(): render sample-QC figures from a saved RGChannelSet, with
# customizable titles.
#
# Produces, in <out_dir>/<cohort>/, one PDF and one PNG per panel:
#   <file_prefix>_qc_intensity.{pdf,png}         intensity QC scatter
#   <file_prefix>_beta_density.{pdf,png}         beta-value density
#   <file_prefix>_failed_fraction_hist.{pdf,png} failed-probe histogram
#   <file_prefix>_sex_prediction.{pdf,png}       predicted vs reported sex
#   <file_prefix>_flagged_qcReport.pdf           minfi control-probe plots for
#                                                 ONLY the flagged samples
#                                                 (requires qc_tab)
#
# Pass qc_tab (from sample_qc()) so the flagged set is known; the control-probe
# report is then restricted to those samples instead of the whole cohort.
# Example usage:
# render_qc_plots(res$rg, cohort = "MSBB",
#                 out_dir = "/…/Results/QC/raw_values",
#                 qc_tab  = res$qc_tab)         # flagged_qcReport uses any_flag
#
# # restrict the control-probe report to a specific criterion instead:
# render_qc_plots(res$rg, cohort = "MSBB", out_dir = "…",
#                 qc_tab = res$qc_tab, flag_col = "sex_mismatch")
# ============================================================================

render_qc_plots <- function(rg,
                            cohort             = "Cohort",
                            out_dir            = ".",
                            file_prefix        = cohort,
                            qc_tab             = NULL,      # from sample_qc(); needed for flagged report
                            flag_col           = "any_flag", # which flag defines "to be dropped"
                            detp               = 0.01,
                            max_fail_frac      = 0.01,
                            min_intensity      = 10.5,
                            min_beadcount      = 3,
                            max_beadcount_frac = 0.01,
                            sex_col       = "sex",
                            titles = list(
                              qc      = paste0(cohort, ": Median Methylated vs Unmethylated Intensity"),
                              density = paste0(cohort, ": Beta-Value Density (raw)"),
                              hist    = paste0(cohort, ": Per-Sample Fraction of Failed Probes"),
                              sex     = paste0(cohort, ": Predicted Sex (X vs Y intensity)"))) {

  pal_ok  <- "#2E7D6F"   # teal
  pal_bad <- "#C0392B"   # red
  pal_m   <- "#2C7FB8"   # blue (male)
  pd      <- as.data.frame(pData(rg))

  ## --- shared computations -----------------------------------------------
  mset        <- preprocessRaw(rg)
  qc          <- getQC(mset)
  qc_int      <- (qc$mMed + qc$uMed) / 2
  detP        <- minfi::detectionP(rg)   # ewastools also exports detectionP(); qualify to avoid masking
  frac_failed <- colMeans(detP > detp)

  ## bead-count QC (requires RGChannelSetExtended; skipped otherwise) ------
  if (is(rg, "RGChannelSetExtended")) {
    bc             <- beadcount_from_rgset(rg)
    frac_low_bead  <- colMeans(bc < min_beadcount, na.rm = TRUE)[colnames(rg)]
    fail_beadcount <- frac_low_bead > max_beadcount_frac
  } else {
    fail_beadcount <- rep(NA, ncol(rg))
  }

  bad <- qc_int < min_intensity | frac_failed > max_fail_frac | fail_beadcount %in% TRUE

  gmset <- mapToGenome(mset)
  sx    <- getSex(gmset)                       # xMed, yMed, predictedSex
  pred  <- sx$predictedSex
  rep_sex <- if (sex_col %in% names(pd)) {
    r <- tolower(as.character(pd[[sex_col]]))
    ifelse(r == "female", "F", ifelse(r == "male", "M", NA_character_))
  } else rep(NA_character_, ncol(rg))
  mism <- !is.na(rep_sex) & rep_sex != pred

  cohort_dir <- file.path(out_dir, cohort)
  dir.create(cohort_dir, showWarnings = FALSE, recursive = TRUE)
  ## --- Panel 1: intensity QC scatter --------------------------------------
  plot_qc_intensity <- function() {
    par(mar = c(5, 5, 4, 2), cex.main = 1.2, cex.lab = 1.05, font.main = 1)
    plot(qc$mMed, qc$uMed,
         xlim = range(c(qc$mMed, qc$uMed, 9)), ylim = range(c(qc$mMed, qc$uMed, 9)),
         pch = 21, cex = 1.1, lwd = 1.1, bg = ifelse(bad, pal_bad, pal_ok), col = "grey30",
         xlab = "Median methylated intensity (log2)",
         ylab = "Median unmethylated intensity (log2)", main = titles$qc)
    abline(a = 2 * min_intensity, b = -1, lty = 2, col = "grey50")
    legend("topleft", bty = "n", pch = 21, pt.bg = c(pal_ok, pal_bad), col = "grey30",
           legend = c("pass", "low intensity / detP / bead-count fail"))
    mtext(sprintf("n = %d samples", ncol(rg)), side = 3, line = 0.2, cex = 0.85, col = "grey40")
  }

  ## --- Panel 2: beta density (color by reported sex if present) ----------
  plot_beta_density <- function() {
    par(mar = c(5, 5, 4, 2), cex.main = 1.2, cex.lab = 1.05, font.main = 1)
    grp <- if (sex_col %in% names(pd)) pd[[sex_col]] else NULL
    densityPlot(rg, sampGroups = grp, main = titles$density,
                xlab = "Beta value", pal = c("#2E7D6F", "#8E44AD", "grey60"))
  }

  ## --- Panel 3: failed-fraction histogram ---------------------------------
  plot_failed_fraction <- function() {
    par(mar = c(5, 5, 4, 2), cex.main = 1.2, cex.lab = 1.05, font.main = 1)
    h <- hist(frac_failed, breaks = 50, plot = FALSE)
    barplot(h$counts, space = 0, col = "#BFD8D2", border = "white",
            names.arg = round(h$mids, 3), las = 2, cex.names = 0.6,
            xlab = "Fraction of probes with detP > threshold", ylab = "Number of samples",
            main = titles$hist)
    abline(v = which.min(abs(h$mids - max_fail_frac)), col = pal_bad, lty = 2, lwd = 2)
    legend("topright", bty = "n", lty = 2, lwd = 2, col = pal_bad,
           legend = sprintf("threshold = %.3f", max_fail_frac))
  }

  ## --- Panel 4: sex prediction (X vs Y median intensity) ------------------
  plot_sex_prediction <- function() {
    par(mar = c(5, 5, 4, 2), cex.main = 1.2, cex.lab = 1.05, font.main = 1)
    plot(sx$xMed, sx$yMed, pch = 19, cex = 1.1,
         col = ifelse(pred == "M", pal_m, pal_bad),
         xlab = "X chr, median total intensity (log2)",
         ylab = "Y chr, median total intensity (log2)", main = titles$sex)
    if (any(mism))                                  # ring the mismatches
      points(sx$xMed[mism], sx$yMed[mism], pch = 1, cex = 2.6, lwd = 2, col = "black")
    legend("topright", bty = "n", pch = c(19, 19, 1), pt.cex = c(1.1, 1.1, 2.2),
           col = c(pal_m, pal_bad, "black"),
           legend = c("predicted M", "predicted F", "reported != predicted"))
    mtext(sprintf("%d predicted-vs-reported mismatch(es)", sum(mism)),
          side = 3, line = 0.2, cex = 0.85, col = "grey40")
  }

  ## --- write each panel as its own PDF and PNG ----------------------------
  save_panel <- function(plot_fn, name, width = 8.5, height = 6.5, png_res = 300) {
    pdf(paste0(cohort_dir, "/", name, ".pdf"), width = width, height = height)
    plot_fn()
    dev.off()

    png(paste0(cohort_dir, "/", name, ".png"), width = width, height = height, units = "in", res = png_res)
    plot_fn()
    dev.off()

    message("wrote ", name, ".pdf and ", name, ".png")
  }

  save_panel(plot_qc_intensity,    paste0(file_prefix, "_qc_intensity"))
  save_panel(plot_beta_density,    paste0(file_prefix, "_beta_density"))
  save_panel(plot_failed_fraction, paste0(file_prefix, "_failed_fraction_hist"))
  save_panel(plot_sex_prediction,  paste0(file_prefix, "_sex_prediction"))

  ## --- control-probe report for FLAGGED samples only ---------------------
  if (!is.null(qc_tab)) {
    if (!flag_col %in% names(qc_tab)) stop("flag_col '", flag_col, "' not in qc_tab")
    flagged_ids <- qc_tab$sample[qc_tab[[flag_col]] %in% TRUE]
    sel <- colnames(rg) %in% flagged_ids
    if (any(sel)) {
      rep_path <- file.path(cohort_dir, paste0(file_prefix, "_flagged_qcReport.pdf"))
      qcReport(rg[, sel], sampNames = colnames(rg)[sel], pdf = rep_path)
      message("wrote ", rep_path, " (", sum(sel), " flagged samples)")
    } else {
      message("no flagged samples (", flag_col, "); skipping control-probe report")
    }
  } else {
    message("qc_tab not supplied; skipping flagged control-probe report")
  }

  invisible(list(qc = qc, frac_failed = frac_failed, bad = bad,
                 predicted_sex = pred, sex_mismatch = mism))
}