# ============================================================================
# pca_qc(): PCA on a cohort's filtered methylation matrix, then test each
# principal component for association with technical and biological factors --
# to detect residual/unwanted clustering (batch, chip, plate, scan date) and
# confirm expected biology (sex, region, diagnosis) BEFORE modelling.
#
# Run AFTER probe filtering (filter_probes_zhou). PCA is on M-values of the
# most-variable, complete-case probes.
#
#   beta        filtered (masked) beta matrix: probes x samples
#   pheno       data.frame of sample metadata (covariates to test). Rows must
#               correspond to samples in `beta`.
#   sample_col  column in `pheno` holding the sample IDs that match
#               colnames(beta); if NULL, rownames(pheno) are used
#   covariates  character vector of `pheno` columns to test against the PCs
#               (mix of categorical + numeric; type auto-detected)
#   color_by    subset of `covariates` to color the PC1-PC2 scatter by
#               (default: first up to 4 covariates)
#   n_top       number of most-variable probes to use (default 20000)
#   n_pcs       number of PCs to compute/associate (default 10)
#   scale_probes scale each probe to unit variance before PCA (default TRUE)
#   r2_flag     R^2 threshold for flagging a PC~covariate association as a
#               concern (default 0.05 = the covariate explains >=5% of the PC)
#   out_prefix  writes <out_prefix>_pca.rds, _pca_scores.csv,
#               _pca_varexplained.csv, _pca_assoc.csv, plus figures both
#               combined (_pca_qc.pdf, multi-page) and one-per-file
#               (_pca_elbow.pdf, _pca_scatter_<covariate>.pdf, _pca_heatmap.pdf)
#   plot_dir    if given, ALSO copies the PDF to
#               <plot_dir>/<basename(out_prefix)>/<basename(out_prefix)>_pca_qc.pdf
#               (mirrors render_qc_plots()'s <plot_root>/<cohort>/ layout, so
#               figures collect under Results/QC/ next to the other QC plots
#               while the data files stay with the rest of that cohort's QC
#               outputs at out_prefix)
#
# Association: for each PC ~ covariate we fit lm() and report R^2 (variance of
# that PC explained by the covariate) as the PRIMARY, decision-relevant metric,
# with the ANOVA F-test p-value kept as a secondary column. R^2 is preferred
# for this batch-detection step because it is an effect size -- unlike the
# p-value it is not inflated by large sample size, so a factor that explains a
# trivial fraction of variance won't be flagged just because n is large.
#
# Returns list(pca, var_explained, assoc_p, assoc_r2, scores).
# ============================================================================

pca_qc <- function(beta,
                   pheno,
                   sample_col   = NULL,
                   covariates,
                   color_by     = NULL,
                   n_top        = 20000,
                   n_pcs        = 10,
                   scale_probes = TRUE,
                   r2_flag      = 0.05,
                   out_prefix   = "cohort",
                   plot_dir     = NULL) {

  ## --- align pheno rows to beta columns ----------------------------------
  ids <- colnames(beta)
  key <- if (is.null(sample_col)) rownames(pheno) else as.character(pheno[[sample_col]])
  ph  <- pheno[match(ids, key), , drop = FALSE]
  if (anyNA(match(ids, key)))
    warning(sum(is.na(match(ids, key))), " sample(s) in beta not found in pheno.")
  covariates <- intersect(covariates, names(ph))
  if (!length(covariates)) stop("None of `covariates` are columns of `pheno`.")

  ## --- most-variable, complete-case probes -> M-values -------------------
  cc <- rowSums(is.na(beta)) == 0
  beta_cc <- beta[cc, , drop = FALSE]
  message(sprintf("[pca_qc] %d complete-case probes of %d", nrow(beta_cc), nrow(beta)))

  v <- if (requireNamespace("matrixStats", quietly = TRUE))
         matrixStats::rowVars(beta_cc) else apply(beta_cc, 1, var)
  top <- order(v, decreasing = TRUE)[seq_len(min(n_top, nrow(beta_cc)))]
  M   <- log2(beta_cc[top, , drop = FALSE] / (1 - beta_cc[top, , drop = FALSE]))
  M[!is.finite(M)] <- NA                       # guard against beta at 0/1
  M   <- M[rowSums(is.na(M)) == 0, , drop = FALSE]

  ## --- PCA (samples as rows) ---------------------------------------------
  n_pcs <- min(n_pcs, ncol(M) - 1L)
  pca   <- prcomp(t(M), center = TRUE, scale. = scale_probes)
  ve    <- (pca$sdev^2) / sum(pca$sdev^2)
  scores <- pca$x[, seq_len(n_pcs), drop = FALSE]

  ## --- associate each PC with each covariate -----------------------------
  assoc_p  <- matrix(NA_real_, n_pcs, length(covariates),
                     dimnames = list(paste0("PC", seq_len(n_pcs)), covariates))
  assoc_r2 <- assoc_p
  for (cv in covariates) {
    x <- ph[[cv]]
    if (is.character(x) || is.logical(x)) x <- factor(x)
    for (i in seq_len(n_pcs)) {
      d <- data.frame(pc = scores[, i], x = x)
      d <- d[stats::complete.cases(d), ]
      if (nrow(d) < 3 || length(unique(d$x)) < 2) next
      fit <- try(stats::lm(pc ~ x, data = d), silent = TRUE)
      if (inherits(fit, "try-error")) next
      aov_p <- tryCatch(stats::anova(fit)[["Pr(>F)"]][1], error = function(e) NA_real_)
      assoc_p[i, cv]  <- aov_p
      assoc_r2[i, cv] <- summary(fit)$r.squared
    }
  }

  ## --- outputs -----------------------------------------------------------
  if (is.null(color_by)) color_by <- head(covariates, 4)
  dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)

  saveRDS(pca, paste0(out_prefix, "_pca.rds"))
  write.csv(data.frame(sample = ids, scores, ph[, covariates, drop = FALSE],
                       check.names = FALSE),
            paste0(out_prefix, "_pca_scores.csv"), row.names = FALSE)
  write.csv(data.frame(PC = paste0("PC", seq_along(ve)),
                       var_explained = round(ve, 5)),
            paste0(out_prefix, "_pca_varexplained.csv"), row.names = FALSE)
  ## long format, R^2-primary (p demoted to a secondary column), sorted by R^2
  assoc_long <- data.frame(
    PC          = rep(rownames(assoc_r2), times = ncol(assoc_r2)),
    covariate   = rep(colnames(assoc_r2), each = nrow(assoc_r2)),
    r2          = round(as.vector(assoc_r2), 4),
    pct_var     = round(100 * as.vector(assoc_r2), 2),
    p_value     = signif(as.vector(assoc_p), 3),
    stringsAsFactors = FALSE
  )
  assoc_long$flagged <- assoc_long$r2 >= r2_flag
  assoc_long <- assoc_long[order(-assoc_long$r2), ]
  write.csv(assoc_long, paste0(out_prefix, "_pca_assoc.csv"), row.names = FALSE)

  ## --- figures -------------------------------------------------------------
  ## Each figure is its own draw_*() closure so it can be rendered into the
  ## combined multi-page PDF AND into its own single-plot PDF without
  ## duplicating the plotting code.
  tag <- basename(out_prefix)

  draw_elbow <- function() {
    par(mar = c(5, 5, 4, 2))
    plot(seq_len(n_pcs), 100 * ve[seq_len(n_pcs)], type = "b", pch = 19, lwd = 2,
         col = "#2E7D6F", xaxt = "n",
         xlab = "PC", ylab = "% variance explained",
         main = paste0(tag, ": elbow plot"))
    axis(1, at = seq_len(n_pcs), labels = paste0("PC", seq_len(n_pcs)), las = 2)
  }

  draw_scatter <- function(cv) {
    x <- ph[[cv]]
    par(mar = c(5, 5, 4, 8), xpd = NA)
    if (is.numeric(x)) {
      pal <- colorRampPalette(c("#2C7FB8", "#EDD9A3", "#C0392B"))(100)
      col <- pal[cut(x, 100, labels = FALSE)]
      plot(scores[, 1], scores[, 2], pch = 19, col = col,
           xlab = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
           ylab = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
           main = paste0(tag, ": PC1 vs PC2 by ", cv))
    } else {
      f   <- factor(x)
      pal <- grDevices::hcl.colors(max(2, nlevels(f)), "Dark 3")
      plot(scores[, 1], scores[, 2], pch = 19, col = pal[as.integer(f)],
           xlab = sprintf("PC1 (%.1f%%)", 100 * ve[1]),
           ylab = sprintf("PC2 (%.1f%%)", 100 * ve[2]),
           main = paste0(tag, ": PC1 vs PC2 by ", cv))
      legend(par("usr")[2], par("usr")[4], legend = levels(f), col = pal,
             pch = 19, bty = "n", title = cv, cex = 0.8)
    }
  }

  draw_heatmap <- function() {
    r2m <- assoc_r2
    r2m[!is.finite(r2m)] <- 0
    par(mar = c(10, 5, 4, 5), xpd = FALSE)
    image(x = seq_len(ncol(r2m)), y = seq_len(nrow(r2m)), z = t(r2m[nrow(r2m):1, , drop = FALSE]),
          col = colorRampPalette(c("white", "#FEE08B", "#C0392B"))(50),
          zlim = c(0, max(r2m, na.rm = TRUE)),
          axes = FALSE, xlab = "", ylab = "",
          main = paste0(tag, ": PC ~ covariate  R^2 (variance explained)"))
    axis(1, at = seq_len(ncol(r2m)), labels = colnames(r2m), las = 2, cex.axis = 0.8)
    axis(2, at = seq_len(nrow(r2m)), labels = rev(rownames(r2m)), las = 2, cex.axis = 0.8)
    ## annotate cells at/above the R^2 flag with their % variance explained
    for (i in seq_len(nrow(assoc_r2))) for (j in seq_len(ncol(assoc_r2))) {
      r2 <- assoc_r2[i, j]; if (is.na(r2) || r2 < r2_flag) next
      text(j, nrow(assoc_r2) - i + 1, sprintf("%.0f%%", 100 * r2), cex = 0.75, font = 2)
    }
    mtext(sprintf("cells labelled where R^2 >= %.0f%%", 100 * r2_flag),
          side = 1, line = 8.5, cex = 0.8, col = "grey40")
  }

  ## open `path`, run draw_fn(), close -- on.exit is scoped to this call so the
  ## device is closed even if draw_fn() errors, without leaking into later plots.
  emit_pdf <- function(path, draw_fn) {
    pdf(path, width = 9, height = 7)
    on.exit(dev.off())
    draw_fn()
  }

  ## combined multi-page PDF (unchanged output file)
  emit_pdf(paste0(out_prefix, "_pca_qc.pdf"), function() {
    draw_elbow()
    for (cv in color_by) draw_scatter(cv)
    draw_heatmap()
  })

  ## same figures, each also saved as its own single-plot PDF
  emit_pdf(paste0(out_prefix, "_pca_elbow.pdf"), draw_elbow)
  for (cv in color_by) {
    cv_safe <- gsub("[^A-Za-z0-9]+", "_", cv)
    emit_pdf(paste0(out_prefix, "_pca_scatter_", cv_safe, ".pdf"), function() draw_scatter(cv))
  }
  emit_pdf(paste0(out_prefix, "_pca_heatmap.pdf"), draw_heatmap)

  message(sprintf("[pca_qc] wrote %s_pca_qc.pdf (+ separate per-plot PDFs) and CSVs", tag))

  if (!is.null(plot_dir)) {
    cohort_dir <- file.path(plot_dir, basename(out_prefix))
    dir.create(cohort_dir, recursive = TRUE, showWarnings = FALSE)
    file.copy(paste0(out_prefix, "_pca_qc.pdf"),
              file.path(cohort_dir, paste0(basename(out_prefix), "_pca_qc.pdf")),
              overwrite = TRUE)
    message(sprintf("[pca_qc] copied plot to %s", cohort_dir))
  }

  invisible(list(pca = pca, var_explained = ve,
                 assoc_p = assoc_p, assoc_r2 = assoc_r2, scores = scores))
}

# ---------------------------------------------------------------------------
# Usage (per cohort, after filter_probes_zhou):
#
#   source("pca_qc.R")
#   beta  <- readRDS("/…/ROSMAP/normalized/ROSMAP_filtered_beta.rds")
#   pheno <- read.csv("/…/ROSMAP_individual_metadata_processed.csv")  # + technical cols
#   # pheno needs one row per SAMPLE (matching colnames(beta)); merge in Sentrix
#   # chip/plate/row/scan-date from the idat manifest if not already present.
#
#   res <- pca_qc(beta, pheno,
#                 sample_col = "specimenID",           # or NULL if rownames match
#                 covariates = c("Sentrix_ID","plate","array_row","scan_date",
#                                "sex","region","ageDeath","Braak","diagnosis"),
#                 color_by   = c("Sentrix_ID","sex","region"),
#                 out_prefix = "/…/ROSMAP/qc/ROSMAP")
#
#   res$assoc_r2          # PC x covariate R^2 (variance explained); primary
#   res$assoc_p           # p-values (secondary)
#
# Read-out: judge on R^2 (the heatmap and _pca_assoc.csv are R^2-driven). A
# technical factor (chip/plate/row/scan_date) explaining a large share (>= r2_flag)
# of an early, high-variance PC is a batch effect -> add it as a covariate /
# random effect (or ComBat) in the EWAS. Sex/region loading a PC is expected
# biology. Use p_value only as supporting evidence, not the decision.
# ---------------------------------------------------------------------------
