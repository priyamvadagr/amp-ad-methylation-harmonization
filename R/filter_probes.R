# ============================================================================
# Probe filtering per Lundin et al. Artherosclerosis(2025) and Methylation(2024)
#
# Each filtering step is its OWN function, with a common interface --
#   drop_*(gmset, beta, ...) -> list(gmset, beta, n_dropped, dropped)
# -- so they compose and can be called independently, in any order:
#
#   drop_low_callrate()   remove probes detected (non-NA) in < min_callrate of
#                          samples. Uses the MASKED beta from normalize_noob().
#   drop_crossreactive()  remove probes that co-hybridize to alternate genomic
#                          sequences               (Zhou MASK.mapping).
#   drop_snp()            remove probes with SNP-introduced artefact; by default
#                          the population-specific masks for White (EUR) and
#                          Black/African American (AFR)  (MASK.snp5.EUR/AFR).
#   drop_nonautosomal()   restrict to chr1..22 (drop chrX / chrY / chrM).
#
#   filter_probes()         just chains the four and returns a per-step summary.
#
# Works for 450K (ROSMAP) and EPIC (others): load the matching platform's Zhou
# ".pop" manifest (see load_zhou_mask()).  
# ============================================================================

library(minfi)

# --- load a Zhou mask table from a downloaded InfiniumAnnotation TSV --------
# http://zwdzwd.github.io/InfiniumAnnotation  -- use the ".pop" manifest (has
# the per-population SNP masks), hg19 to match ilmn12.hg19 / ilm10b4.hg19:
#   HM450.hg19.manifest.pop.tsv.gz   (450K -> ROSMAP)
#   EPIC.hg19.manifest.pop.tsv.gz    (EPIC -> MOA-PAD / ROSMAP_APOE4 / MSBB)
# Columns include: probeID, MASK.mapping, MASK.snp5.EUR, MASK.snp5.AFR,
# MASK.snp5.GMAF1p, MASK.general, ...
load_zhou_mask <- function(tsv_path) {
  stopifnot(file.exists(tsv_path))
  read.delim(tsv_path, header = TRUE, stringsAsFactors = FALSE,
             na.strings = c("NA", ""), check.names = FALSE)
}

# --- internal helpers -------------------------------------------------------
# align gmset & beta to the same probes/order (beta from gmset if NULL)
.align_gb <- function(gmset, beta) {
  stopifnot(is(gmset, "GenomicMethylSet") || is(gmset, "GenomicRatioSet"))
  if (is.null(beta)) beta <- getBeta(gmset)
  common <- intersect(rownames(gmset), rownames(beta))
  list(gmset = gmset[common, ], beta = beta[common, , drop = FALSE])
}

# TRUE for probe_ids flagged in ANY of the given Zhou MASK columns
.zhou_mask_hits <- function(probe_ids, mask_tbl, cols, id_col = "probeID") {
  if (is(mask_tbl, "GRanges")) {
    mdf <- as.data.frame(S4Vectors::mcols(mask_tbl), stringsAsFactors = FALSE)
    mdf[[id_col]] <- names(mask_tbl)
  } else mdf <- as.data.frame(mask_tbl, stringsAsFactors = FALSE)
  if (!id_col %in% names(mdf))
    stop("mask_tbl has no id column '", id_col, "'. Columns: ",
         paste(head(names(mdf), 20), collapse = ", "))
  miss <- setdiff(cols, names(mdf))
  if (length(miss))
    stop("mask cols not found: ", paste(miss, collapse = ", "),
         "\nAvailable MASK columns: ",
         paste(grep("^MASK", names(mdf), value = TRUE), collapse = ", "))
  m   <- match(probe_ids, mdf[[id_col]])
  hit <- rep(FALSE, length(probe_ids))
  for (cc in cols) { v <- as.logical(mdf[[cc]][m]); v[is.na(v)] <- FALSE; hit <- hit | v }
  hit
}

# subset both objects by a keep logical, package the drop info
.finish <- function(gmset, beta, keep, tag, detail) {
  dropped <- rownames(gmset)[!keep]
  message(sprintf("[%s] removed %d (%s); %d retained", tag, sum(!keep), detail, sum(keep)))
  list(gmset = gmset[keep, ], beta = beta[keep, , drop = FALSE],
       n_dropped = sum(!keep), dropped = dropped)
}

# ---------------------------------------------------------------------------
# Step 1: low call-rate probes (uses the MASKED beta) -----------------------
# ---------------------------------------------------------------------------
drop_low_callrate <- function(gmset, beta = NULL, min_callrate = 0.99) {
  a <- .align_gb(gmset, beta)
  keep <- rowMeans(!is.na(a$beta)) >= min_callrate
  .finish(a$gmset, a$beta, keep, "drop_low_callrate",
          sprintf("call rate < %.2f", min_callrate))
}

# ---------------------------------------------------------------------------
# Step 2: cross-hybridizing probes (co-hybridize to alternate sequences) ----
# ---------------------------------------------------------------------------
drop_crossreactive <- function(gmset, beta = NULL, mask_tbl,
                               cols = "MASK.mapping", id_col = "probeID") {
  a   <- .align_gb(gmset, beta)
  hit <- .zhou_mask_hits(rownames(a$gmset), mask_tbl, cols, id_col)
  .finish(a$gmset, a$beta, !hit, "drop_crossreactive", paste(cols, collapse = "+"))
}

# ---------------------------------------------------------------------------
# Step 3: SNP-introduced artefact (population-specific by default) -----------
# ---------------------------------------------------------------------------
drop_snp <- function(gmset, beta = NULL, mask_tbl,
                     cols = c("MASK.snp5.EUR", "MASK.snp5.AFR"),
                     id_col = "probeID") {
  a   <- .align_gb(gmset, beta)
  hit <- .zhou_mask_hits(rownames(a$gmset), mask_tbl, cols, id_col)
  .finish(a$gmset, a$beta, !hit, "drop_snp", paste(cols, collapse = "+"))
}

# ---------------------------------------------------------------------------
# Step 4: restrict to autosomes (uses the object's own hg19 coordinates) ----
# ---------------------------------------------------------------------------
drop_nonautosomal <- function(gmset, beta = NULL) {
  a    <- .align_gb(gmset, beta)
  chr  <- as.character(GenomeInfoDb::seqnames(granges(a$gmset)))
  keep <- chr %in% paste0("chr", 1:22)
  .finish(a$gmset, a$beta, keep, "drop_nonautosomal", "chrX/Y/other")
}

# ---------------------------------------------------------------------------
# Wrapper: chain the four steps, return a per-step summary -------------------
# (steps are applied in sequence, so each n_removed is *additional* removal on
#  the already-reduced set; the retained set is order-independent.)
# ---------------------------------------------------------------------------
filter_probes <- function(gmset,
                               beta               = NULL,
                               mask_tbl,
                               crossreactive_cols = "MASK.mapping",
                               snp_cols           = c("MASK.snp5.EUR", "MASK.snp5.AFR"),
                               id_col             = "probeID",
                               min_callrate       = 0.99,
                               autosomal_only     = TRUE,
                               out_prefix         = NULL) {

  a <- .align_gb(gmset, beta); gmset <- a$gmset; beta <- a$beta
  n0 <- nrow(beta)

  r <- drop_low_callrate(gmset, beta, min_callrate);            gmset <- r$gmset; beta <- r$beta; n1 <- nrow(beta)
  r <- drop_crossreactive(gmset, beta, mask_tbl, crossreactive_cols, id_col); gmset <- r$gmset; beta <- r$beta; n2 <- nrow(beta)
  r <- drop_snp(gmset, beta, mask_tbl, snp_cols, id_col);       gmset <- r$gmset; beta <- r$beta; n3 <- nrow(beta)
  n4 <- n3
  if (autosomal_only) { r <- drop_nonautosomal(gmset, beta);    gmset <- r$gmset; beta <- r$beta; n4 <- nrow(beta) }

  summary <- data.frame(
    step        = c("input", "low_call_rate", "cross_reactive", "snp_artefact", "non_autosomal"),
    n_removed   = c(0L, n0 - n1, n1 - n2, n2 - n3, n3 - n4),
    n_remaining = c(n0, n1, n2, n3, n4),
    stringsAsFactors = FALSE
  )
  message(sprintf("[filter_probes_zhou] %d -> %d probes retained", n0, n4))

  if (!is.null(out_prefix)) {
    dir.create(dirname(out_prefix), recursive = TRUE, showWarnings = FALSE)
    saveRDS(gmset, paste0(out_prefix, "_filtered_gmset.rds"))
    saveRDS(beta,  paste0(out_prefix, "_filtered_beta.rds"))
    write.csv(summary, paste0(out_prefix, "_probe_filter_summary.csv"), row.names = FALSE)
    message("[filter_probes_zhou] wrote _filtered_gmset.rds, _filtered_beta.rds, _probe_filter_summary.csv")
  }

  list(gmset = gmset, beta = beta, summary = summary)
}

