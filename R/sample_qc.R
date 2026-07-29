# ============================================================================
# Sample-level QC for one cohort's RGChannelSet (minfi + ewastools).
# Run BEFORE normalization. This file defines two functions, meant to be used
# together:
#
#   sample_qc()        Computes per-sample QC flags (detection p, intensity,
#                       sex mismatch, bisulfite conversion, ewastools control
#                       probes), writes them to a CSV, and copies them into
#                       pData(rg). FLAGS every problem; drops NOTHING -- you
#                       decide what to remove downstream using the flags.
#
#   render_qc_plots()   Renders the diagnostic figures (intensity scatter,
#                       beta density, failed-probe histogram, sex prediction)
#                       and, given the qc_tab from sample_qc(), a
#                       control-probe report restricted to flagged samples.
#
# Typical flow: res <- sample_qc(rg, ...); render_qc_plots(res$rg, qc_tab = res$qc_tab, ...)
#
# sample_qc() flags (all flag-only; none force removal):
#   fail_detp      - too many probes with signal indistinguishable from noise
#   fail_intensity - overall array signal too low
#   sex_mismatch   - predicted sex disagrees with recorded sex
#   fail_bisulfite - bisulfite conversion control below threshold
#   fail_control   - ewastools sample_failure() over the 17 Illumina controls
#
# `flag_reason` records every triggered check per sample ("none" if clean).
# The flags are also copied into pData(rg) so they travel with the object.
# ============================================================================

library(minfi)
library(ewastools)
library(IlluminaHumanMethylationEPICmanifest)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19)

# ---------------------------------------------------------------------------
# sample_qc(): compute per-sample QC flags for one cohort's RGChannelSet.
#
# The RGChannelSet is obtained in this priority order: `rg` if supplied >
# `rds` if it exists. sample_qc() never builds a RGChannelSet from raw idats
# itself -- that's compile_idat_files.R's job; pass its output (directly, or
# via the .rds it wrote) here.
#
#   rg             RGChannelSet (annotation already set); if supplied, `rds`
#                  is ignored
#   rds            path to an existing combined-idat RGChannelSet .rds, as
#                  written by compile_idat_files.R; read if `rg` is not supplied
#   targets        idat manifest data.frame/list (e.g. the one compile_idat_files.R
#                  used); needs a Basename column, or `raw_dir` + `grn_col` to
#                  build one -- used only to derive `basenames` for ewastools
#   raw_dir        directory of raw idats, used with `grn_col` to build
#                  Basename paths when `targets` only has bare filenames
#   grn_col        `targets` column with bare "*_Grn.idat" filenames (paired
#                  with `raw_dir`) when no Basename column is present
#   basenames      idat path prefixes (no _Grn/_Red.idat) for ewastools;
#                  matched to rg by basename(). If NULL (and not derivable
#                  from `targets`/`raw_dir`), ewastools flags = NA.
#   sex_col        pData column with reported sex ("female"/"male"/...)
#   detp           per-probe detection-p threshold (probe "fails" above this)
#   max_fail_frac  flag sample if > this fraction of probes fail detection
#   min_intensity  flag sample if median log2 (M+U)/2 intensity below this
#   bisulfite_min  flag sample if any bisulfite-conversion metric < this
#   out_prefix     prefix for the output CSV (plots are in render_qc_plots())
#
# Example:
#   manifest <- read.csv(cfg$manifest, stringsAsFactors = FALSE)
#   manifest$Basename <- file.path(cfg$raw_dir, sub("_Grn\\.idat$", "", manifest$grnFile))
#
#   res <- sample_qc(rg_msbb, sex_col = "sex",
#                    basenames = manifest$Basename, out_prefix = "MSBB")
#
#   res$qc_tab[res$qc_tab$any_flag, c("sample","flag_reason")]  # who & why
#
#   # YOU decide what to remove, e.g. drop only technical failures:
#   keep <- !(res$qc_tab$fail_detp | res$qc_tab$fail_intensity | res$qc_tab$sex_mismatch)
#   rg_clean <- res$rg[, keep]
#
# Example (loading a compiled rds instead of passing `rg` directly):
#   res <- sample_qc(rds = "/…/MSBB/combined_idat/MSBB.rds",
#                    targets = manifest, raw_dir = cfg$raw_dir,
#                    sex_col = "sex", out_prefix = "MSBB")
# ---------------------------------------------------------------------------
sample_qc <- function(rg            = NULL,
                      rds           = NULL,
                      targets       = NULL,
                      raw_dir       = NULL,
                      grn_col       = "grnFile",
                      basenames     = NULL,
                      sex_col       = "sex",
                      detp          = 0.01,
                      max_fail_frac = 0.01,
                      min_intensity = 10.5,
                      bisulfite_min = 1,
                      out_prefix    = "cohort") {

  ## --- obtain the RGChannelSet: object > existing rds ---------------------
  if (is.null(rg)) {
    if (!is.null(rds) && file.exists(rds)) {
      message("[sample_qc] using existing combined-idat rds: ", rds)
      rg <- readRDS(rds)
    } else {
      stop("Provide `rg` or an existing `rds` (as written by compile_idat_files.R); ",
           "sample_qc() does not build RGChannelSets from raw idats.")
    }
  }
  ## if we loaded an rds (or were handed rg) and only a manifest is available,
  ## build ewastools Basenames from raw_dir + grn_col
  if (is.null(basenames) && !is.null(targets) && !is.null(raw_dir)) {
    tgm <- as.data.frame(targets, stringsAsFactors = FALSE)
    if (grn_col %in% names(tgm)) {
      tgm <- tgm[!is.na(tgm[[grn_col]]) & nzchar(tgm[[grn_col]]), ]
      basenames <- file.path(raw_dir, sub("_Grn\\.idat$", "", tgm[[grn_col]]))
    } else if ("Basename" %in% names(tgm)) {
      basenames <- tgm$Basename
    }
  }

  stopifnot(is(rg, "RGChannelSet"))
  pd <- as.data.frame(pData(rg))

  ## detection p-values ------------------------------------------------------
  ## minfi:: qualified -- ewastools also exports detectionP() (for its own
  ## read_idats() objects), and it masks minfi's version since ewastools is
  ## library()'d after minfi above.
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

  ## ewastools control metrics ---------------------------------------------
  fail_bisulfite    <- rep(NA, ncol(rg))
  fail_control      <- rep(NA, ncol(rg))
  bisulfite_min_val <- rep(NA_real_, ncol(rg))

  if (!is.null(basenames)) {
    ord <- match(colnames(rg), basename(basenames))
    if (anyNA(ord))
      warning(sum(is.na(ord)), " sample(s) in rg not found in `basenames`; ",
              "their ewastools flags will be NA.")
    meth <- read_idats(basenames)
    ctrl <- control_metrics(meth)
    ctrl_failed <- sample_failure(ctrl)

    bis_names <- grep("bisulfite", names(ctrl), ignore.case = TRUE, value = TRUE)
    if (length(bis_names)) {
      bis_mat <- do.call(cbind, lapply(bis_names, function(n) ctrl[[n]]))
      bis_min <- apply(bis_mat, 1, min, na.rm = TRUE)
    } else {
      warning("No bisulfite-conversion metric found in control_metrics().")
      bis_min <- rep(NA_real_, length(ctrl_failed))
    }
    bisulfite_min_val <- bis_min[ord]
    fail_bisulfite    <- bis_min[ord] < bisulfite_min
    fail_control      <- ctrl_failed[ord]
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
    fail_detp      = frac_failed  > max_fail_frac,
    fail_intensity = qc_intensity < min_intensity,
    sex_mismatch   = sex_mismatch,
    fail_bisulfite = fail_bisulfite,
    fail_control   = fail_control,
    stringsAsFactors = FALSE
  )

  ## flag_reason: every triggered check, semicolon-joined ("none" if clean)
  flag_mat <- cbind(
    detP      = qc_tab$fail_detp      %in% TRUE,
    intensity = qc_tab$fail_intensity %in% TRUE,
    sex       = qc_tab$sex_mismatch   %in% TRUE,
    bisulfite = qc_tab$fail_bisulfite %in% TRUE,
    control   = qc_tab$fail_control   %in% TRUE
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
    "[%s] %d samples | flagged: %d (detP %d, intensity %d, sex %d, bisulfite %d, control %d)",
    out_prefix, nrow(qc_tab), sum(qc_tab$any_flag),
    sum(qc_tab$fail_detp %in% TRUE), sum(qc_tab$fail_intensity %in% TRUE),
    sum(qc_tab$sex_mismatch %in% TRUE), sum(qc_tab$fail_bisulfite %in% TRUE),
    sum(qc_tab$fail_control %in% TRUE)))

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
# res <- sample_qc(rg, sex_col = "sex", basenames = manifest$Basename, out_prefix = "MSBB")
# render_qc_plots(res$rg, cohort = "MSBB",
#                 out_dir = "/…/Results/QC/raw_values",
#                 qc_tab  = res$qc_tab)         # flagged_qcReport uses any_flag
#
# # restrict the control-probe report to a specific criterion instead:
# render_qc_plots(res$rg, cohort = "MSBB", out_dir = "…",
#                 qc_tab = res$qc_tab, flag_col = "sex_mismatch")
# ============================================================================

render_qc_plots <- function(rg,
                            cohort        = "Cohort",
                            out_dir       = ".",
                            file_prefix   = cohort,
                            qc_tab        = NULL,           # from sample_qc(); needed for flagged report
                            flag_col      = "any_flag",     # which flag defines "to be dropped"
                            detp          = 0.01,
                            max_fail_frac = 0.01,
                            min_intensity = 10.5,
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
  bad         <- qc_int < min_intensity | frac_failed > max_fail_frac

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
           legend = c("pass", sprintf("low intensity / detP fail")))
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

# ---------------------------------------------------------------------------




