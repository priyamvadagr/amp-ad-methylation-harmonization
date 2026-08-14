# ============================================================================
# Subset a beta matrix to the cortical-clock CpGs and inspect their coverage
# (missingness) across samples BEFORE running the clock -- so you know how many
# of the 347 sites are present and how much has to be imputed.
#
# The clock CpG list is read from the authors' coefficient file
# (CorticalClockCoefs.csv, https://github.com/gemmashireby/CorticalClock),
# not hard-coded.
#
# Two kinds of "missing" are tracked separately:
#   - ABSENT : the CpG is not on your platform / was removed pre-clock
#              (whole row NA)  -> a coverage problem
#   - MASKED : the CpG is present but NA in some samples (detP masking)
#              (scattered NA)  -> per-sample data loss
# ============================================================================

# --- read the clock CpG ids from the coefficient file -----------------------
# CorticalClockCoefs.txt (as shipped at
# https://github.com/gemmashireby/CorticalClock) is whitespace-delimited,
# not comma-delimited -- read.table (default sep = any whitespace) parses it
# correctly where read.csv would collapse "probe coef" into one column.
clock_cpgs_from_coef <- function(coef_file, probe_col = NULL,
                                 intercept_label = "(Intercept)") {
  cf <- read.table(coef_file, header = TRUE, stringsAsFactors = FALSE)
  if (is.null(probe_col)) probe_col <- names(cf)[1]
  ids <- as.character(cf[[probe_col]])
  setdiff(ids, intercept_label)
}

# --- subset any probe x sample matrix to a fixed CpG id list (NA rows for
# ids absent from the matrix) -- shared by subset_clock_cpgs() and
# plot_clock_missingness()'s detP heatmap.
subset_to_clock <- function(mat, ids) {
  present <- ids %in% rownames(mat)
  out <- matrix(NA_real_, length(ids), ncol(mat),
               dimnames = list(ids, colnames(mat)))
  out[ids[present], ] <- mat[ids[present], , drop = FALSE]
  out
}

# --- subset beta to the 347 clock CpGs (NA rows for absent sites) -----------
subset_clock_cpgs <- function(beta, coef_file, probe_col = NULL,
                              intercept_label = "(Intercept)") {
  clock_cpgs <- clock_cpgs_from_coef(coef_file, probe_col, intercept_label)
  present    <- clock_cpgs %in% rownames(beta)
  mat <- subset_to_clock(beta, clock_cpgs)
  message(sprintf("[clock] %d/%d clock CpGs present on this matrix (%d absent)",
                  sum(present), length(clock_cpgs), sum(!present)))
  list(mat = mat, clock_cpgs = clock_cpgs, present = present)
}

# --- coverage summary + plots ----------------------------------------------
#   clock_mat   347 x samples beta (from subset_clock_cpgs)
#   detP        full detection-p matrix (probes x samples); subset here to the
#               clock CpGs and shown as a heatmap (median detP per sample too).
#   out_prefix  where the two coverage CSVs are written (data)
#   plot_prefix where the missingness PDF is written (figures); defaults to
#               out_prefix so existing single-prefix callers still work
plot_clock_missingness <- function(clock_mat, present = NULL, detP = NULL,
                                   cohort = "Cohort", out_prefix = NULL,
                                   plot_prefix = out_prefix,
                                   upset_top = 20) {

  n_cpg  <- nrow(clock_mat); n_samp <- ncol(clock_mat)
  na_mat <- is.na(clock_mat)
  per_sample <- colMeans(na_mat)
  per_cpg    <- rowMeans(na_mat)
  if (is.null(present)) present <- per_cpg < 1

  ## clock-CpG detP + per-sample median detP
  detP_clock <- if (!is.null(detP)) subset_to_clock(detP, rownames(clock_mat)) else NULL
  med_detp   <- if (!is.null(detP_clock)) apply(detP_clock, 2, median, na.rm = TRUE) else NA

  message(sprintf("[%s] clock coverage: %d CpGs x %d samples | absent %d | present-masked %d",
                  cohort, n_cpg, n_samp, sum(!present), sum(present & per_cpg > 0)))
  message(sprintf("  per-sample missing: median %.4f, max %.4f",
                  median(per_sample), max(per_sample)))
  if (!is.null(detP_clock))
    message(sprintf("  per-sample median detP at clock CpGs: median %.2e, max %.2e",
                    median(med_detp, na.rm = TRUE), max(med_detp, na.rm = TRUE)))

  if (!is.null(out_prefix)) {
    dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
    write.csv(data.frame(sample = colnames(clock_mat),
                         missing_frac = round(per_sample, 5),
                         n_missing = colSums(na_mat),
                         median_detP = signif(med_detp, 3)),
              paste0(out_prefix, "_clock_sample_coverage.csv"), row.names = FALSE)
    write.csv(data.frame(cpg = rownames(clock_mat), present = present,
                         missing_frac_samples = round(per_cpg, 5)),
              paste0(out_prefix, "_clock_cpg_coverage.csv"), row.names = FALSE)
  }

  if (!is.null(plot_prefix)) {
    dir.create(dirname(plot_prefix), recursive = TRUE, showWarnings = FALSE)
    pdf(paste0(plot_prefix, "_clock_missingness.pdf"), width = 9, height = 7)
  }
  on.exit(if (!is.null(plot_prefix) && grDevices::dev.cur() > 1) dev.off(), add = TRUE)

  ## 1. per-sample missingness histogram
  par(mar = c(5, 5, 4, 2))
  hist(100 * per_sample, breaks = 40, col = "#2E7D6F", border = "white",
       xlab = "% of 347 clock CpGs missing (per sample)", ylab = "samples",
       main = paste0(cohort, ": per-sample clock CpG missingness"))
  abline(v = 100 * median(per_sample), col = "#C0392B", lty = 2, lwd = 2)

  ## 2. per-CpG missingness, absent sites flagged
  o <- order(per_cpg, decreasing = TRUE)
  par(mar = c(5, 5, 4, 2))
  barplot(100 * per_cpg[o], col = ifelse(!present[o], "#C0392B", "#2C7FB8"),
          border = NA, names.arg = FALSE,
          xlab = "clock CpGs (sorted by missingness)", ylab = "% of samples missing",
          main = paste0(cohort, ": per-CpG missingness"))
  legend("topright", bty = "n", fill = c("#C0392B", "#2C7FB8"),
         legend = c("absent on platform", "present (some masked)"))

  ## 3. detP heatmap at clock CpGs (replaces the binary raster) ------------
  if (!is.null(detP_clock)) {
    ord  <- order(med_detp, na.last = TRUE)                 # samples worst->best
    logp <- log10(pmax(detP_clock[, ord, drop = FALSE], 1e-16))
    par(mar = c(5, 5, 4, 5))
    image(x = seq_len(n_samp), y = seq_len(n_cpg),
          z = t(logp[n_cpg:1, , drop = FALSE]),
          col = colorRampPalette(c("#08306B", "#F7F7F7", "#C0392B"))(50),
          zlim = c(-16, 0), axes = TRUE,
          xlab = "samples (ordered by median detP)", ylab = "clock CpGs",
          main = paste0(cohort, ": detection-p at clock CpGs  log10(detP)"))
    mtext("blue = strong signal (low detP), red = weak (detP -> 1); white rows = absent",
          side = 3, line = 0.2, cex = 0.75, col = "grey40")
  } else {
    message("[", cohort, "] no detP supplied; skipping detP heatmap.")
  }

  ## 4. UpSet plot of masking co-occurrence -------------------------------
  ## sets = clock CpGs with variable masking; elements = samples; membership =
  ## CpG is NA (masked) in that sample. Reveals whether the same samples lose
  ## the same CpGs (structured) vs scattered dropouts.
  masked_cpgs <- rownames(clock_mat)[present & per_cpg > 0]
  if (length(masked_cpgs) >= 2 && requireNamespace("UpSetR", quietly = TRUE)) {
    bin  <- as.data.frame(t(is.na(clock_mat[masked_cpgs, , drop = FALSE])) * 1)
    keep <- names(sort(colMeans(bin), decreasing = TRUE))[seq_len(min(upset_top, ncol(bin)))]
    bin  <- bin[, keep, drop = FALSE]
    bin  <- bin[rowSums(bin) > 0, , drop = FALSE]           # samples with any masking
    if (ncol(bin) >= 2 && nrow(bin) >= 1) {
      print(UpSetR::upset(bin, sets = rev(keep), nsets = length(keep),
                          order.by = "freq", keep.order = TRUE,
                          mainbar.y.label = "samples sharing this missing-CpG set",
                          sets.x.label   = "samples missing the CpG"))
      grid::grid.text(paste0(cohort, ": co-occurrence of masked clock CpGs"),
                      x = 0.65, y = 0.97, gp = grid::gpar(fontsize = 12))
    }
  } else if (length(masked_cpgs) < 2) {
    message("[", cohort, "] <2 masked clock CpGs; skipping UpSet (coverage is clean).")
  } else {
    message("[", cohort, "] package 'UpSetR' not installed; skipping UpSet plot ",
            "(install.packages('UpSetR')).")
  }

  invisible(list(per_sample = per_sample, per_cpg = per_cpg,
                 present = present, median_detP = med_detp))
}

# --- accuracy plot (predicted vs chronological age), one PDF per tissue -----
# Same geom_abline/geom_point/theme as CorticalClock.r's plotAge(), just run
# once per distinct value of tissue_col instead of over all samples pooled --
# a clock trained on cortex can be systematically off in non-cortical tissue
# (e.g. cerebellum), which a single pooled plot/correlation would hide.
#   data       data.frame with Age, brainpred (as written to
#              <cohort>_CorticalPred.csv by run_cortical_clock.R) plus a
#              tissue column
#   tissue_col name of the tissue column in `data`
#   out_prefix one PDF per tissue written to
#              <out_prefix>_<tissue>_clock_accuracy.pdf (tissue name
#              sanitized for use in a filename)
plot_clock_by_tissue <- function(data, tissue_col = "tissue", cohort = "Cohort",
                                 out_prefix = NULL) {
  if (!requireNamespace("ggplot2", quietly = TRUE))
    stop("package 'ggplot2' is required for plot_clock_by_tissue()")
  stopifnot(all(c("Age", "brainpred", tissue_col) %in% names(data)))

  clock_theme <- ggplot2::theme_minimal() +
    ggplot2::theme(panel.border = ggplot2::element_blank(),
                  panel.grid.major = ggplot2::element_blank(),
                  panel.grid.minor = ggplot2::element_blank(),
                  axis.line = ggplot2::element_line(colour = "black"),
                  axis.title = ggplot2::element_text(size = 12, face = "bold"),
                  legend.text = ggplot2::element_text(size = 12, face = "bold"),
                  plot.title = ggplot2::element_text(size = 12, face = "bold.italic"),
                  plot.subtitle = ggplot2::element_text(size = 11),
                  axis.text = ggplot2::element_text(size = 14, face = "bold"),
                  strip.text.x = ggplot2::element_text(face = "bold", size = 12))

  tissues <- sort(unique(data[[tissue_col]][!is.na(data[[tissue_col]])]))
  if (!length(tissues)) stop("no non-NA values in tissue_col '", tissue_col, "'")

  plots <- list()
  for (ts in tissues) {
    sub <- data[!is.na(data[[tissue_col]]) & data[[tissue_col]] == ts, ]
    sub <- sub[complete.cases(sub[c("Age", "brainpred")]), ]
    if (!nrow(sub)) {
      message("[", cohort, "] tissue '", ts, "': no samples with complete Age + brainpred, skipping")
      next
    }

    corr <- if (nrow(sub) >= 2) cor(sub$Age, sub$brainpred) else NA_real_
    rmse <- sqrt(mean((sub$Age - sub$brainpred)^2))
    mad  <- median(abs(sub$Age - sub$brainpred))

    p <- ggplot2::ggplot(data = sub, ggplot2::aes(x = Age, y = brainpred)) +
      ggplot2::geom_abline(intercept = 0, slope = 1) +
      ggplot2::geom_point() +
      ggplot2::xlab("Chronological Age") +
      ggplot2::ylab("Predicted Age") +
      ggplot2::ggtitle(paste0(cohort, ": Cortical Clock -- ", ts),
                       subtitle = sprintf("n=%d, r=%.2f, RMSE=%.1f yrs, MAD=%.1f yrs",
                                         nrow(sub), corr, rmse, mad)) +
      clock_theme

    plots[[ts]] <- p
    message(sprintf("[%s] %s: n=%d, r=%.2f, RMSE=%.1f, MAD=%.1f", cohort, ts, nrow(sub), corr, rmse, mad))

    if (!is.null(out_prefix)) {
      dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
      safe_ts <- gsub("[^A-Za-z0-9]+", "_", tolower(ts))
      f <- paste0(out_prefix, "_", safe_ts, "_clock_accuracy.pdf")
      ggplot2::ggsave(f, plot = p, height = 5, width = 10)
      message("[", cohort, "] wrote ", basename(f))
    }
  }

  invisible(plots)
}

# ---------------------------------------------------------------------------
# Usage:
#   source("cortical_clock.R")
#   beta  <- readRDS("/…/ROSMAP_noob_beta.rds")           # pre probe-filter
#   detP  <- readRDS("/…/ROSMAP_qc.rds")$detP              # optional, probes x samples
#   coefs <- "/…/CorticalClock/CorticalClockCoefs.txt"
#
#   sub <- subset_clock_cpgs(beta, coef_file = coefs)      # 347 x samples matrix
#   plot_clock_missingness(sub$mat, present = sub$present, detP = detP,
#                          cohort = "ROSMAP",
#                          out_prefix  = "/…/data/…/ROSMAP/Cortical_clock/ROSMAP",   # CSVs
#                          plot_prefix = "/…/Results/Cortical_clock/ROSMAP/ROSMAP") # PDF
#
# Read-out: ideally 347/347 present and the per-sample histogram sits near 0%.
# Red bars in plot 2 are CpGs entirely absent (must be imputed with reference
# betas). The plot 3 detP heatmap and plot 4 UpSet co-occurrence plot (both
# skipped if detP is NULL / UpSetR isn't installed) show whether masking is
# scattered or concentrated in specific samples -> those samples' DNAmAge
# will lean on imputation, interpret with caution.
# ---------------------------------------------------------------------------