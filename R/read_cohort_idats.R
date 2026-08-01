library(minfi)

# ---------------------------------------------------------------------------
# idat_has_valid_magic(): TRUE if `path` starts with the IDAT binary magic
# header. Catches files that are corrupt/truncated at the source (e.g. a bad
# Synapse upload) — same size as a real IDAT but zeroed or garbled content —
# which otherwise abort the whole read.metharray.exp() bplapply batch with a
# cryptic "Unknown magic" BiocParallel error.
# ---------------------------------------------------------------------------
idat_has_valid_magic <- function(path) {
  con <- file(path, "rb")
  on.exit(close(con))
  magic <- tryCatch(readBin(con, "raw", n = 4), error = function(e) raw(0))
  length(magic) == 4 && rawToChar(magic) == "IDAT"
}

# ---------------------------------------------------------------------------
# read_cohort_idats(): read all idat pairs for one cohort into one object.
#
#   idat_dir     folder containing the *_Grn.idat / *_Red.idat files
#   sample_sheet (optional) data.frame of samples to read. If supplied, only
#                these samples are read (recommended — lets you restrict to
#                metadata-matched samples and attach phenotype up front).
#   sentrix_col  name of the column holding the full "<slide>_<array>" id
#                (e.g. "5822062002_R03C02"). Used to build the file prefix.
#   recursive    search idat_dir subfolders (TRUE if idats are in per-chip dirs)
#   force        TRUE lets minfi bind arrays whose manifests differ slightly;
#                safe within one platform and avoids spurious read failures.
#   extended     TRUE (default) reads into an RGChannelSetExtended, which also
#                carries the per-probe bead counts (NBeads) and bead SDs. This
#                lets bead-count QC (probes measured by < 3 functional beads)
#                run off the compiled .rds directly, with no separate re-read
#                of the raw idats. It is a subclass of RGChannelSet, so all
#                downstream minfi steps (preprocessRaw, detectionP, getSex,
#                preprocessNoob, ...) work unchanged. Trade-off: the extended
#                object / .rds is larger (adds the NBeads + SD assays), so set
#                extended = FALSE if you don't need bead QC and want to save
#                memory/disk.
#   chunk_size   read.metharray.exp() loads every listed sample's raw idats
#                into memory at once before assembling the RGChannelSet -- on
#                large cohorts that's what gets the compile step killed by the
#                OOM killer. If set, samples are instead read chunk_size at a
#                time and the chunks combined with cbind() as they complete,
#                bounding peak memory to one chunk's raw idats. NULL (default)
#                = read everything in a single call (original behavior).
# ---------------------------------------------------------------------------
read_cohort_idats <- function(idat_dir,
                              sample_sheet = NULL,
                              sentrix_col  = "Sentrix_ID",
                              recursive    = TRUE,
                              force        = TRUE,
                              extended     = TRUE,
                              chunk_size   = NULL) {

  if (is.null(sample_sheet)) {
    if (is.null(chunk_size)) {
      # --- read every idat pair found under idat_dir, in one call ---
      return(read.metharray.exp(base = idat_dir,
                                recursive = recursive,
                                extended = extended,
                                force = force,
                                verbose = TRUE))
    }

    # Chunking needs an explicit sample list -- build one from the idat pairs
    # found under idat_dir instead of handing base= to read.metharray.exp().
    grn_files <- list.files(idat_dir, pattern = "_Grn\\.idat$",
                            recursive = recursive, full.names = TRUE)
    ss <- data.frame(Basename = sub("_Grn\\.idat$", "", grn_files),
                     stringsAsFactors = FALSE)
  } else {
    # --- read only the samples listed in sample_sheet ---
    ss <- as.data.frame(sample_sheet, stringsAsFactors = FALSE)

    # Build the Basename column minfi needs: path to each idat WITHOUT the
    # trailing _Grn.idat / _Red.idat. If a Basename column already exists, use it.
    if (!"Basename" %in% names(ss)) {
      stopifnot(sentrix_col %in% names(ss))
      ss$Basename <- file.path(idat_dir, ss[[sentrix_col]])
    }

    # Warn about any listed samples whose idat files are missing on disk
    missing <- !file.exists(paste0(ss$Basename, "_Grn.idat")) |
               !file.exists(paste0(ss$Basename, "_Red.idat"))
    if (any(missing)) {
      warning(sum(missing), " sample(s) have no _Grn.idat/_Red.idat and will be skipped:\n  ",
              paste(head(ss$Basename[missing], 10), collapse = "\n  "))
      ss <- ss[!missing, , drop = FALSE]
    }

    # Drop samples whose Grn/Red idat exists but is corrupt (bad magic header)
    # before handing off to read.metharray.exp(), so one bad file doesn't
    # abort the whole cohort's parallel read.
    corrupt <- !vapply(paste0(ss$Basename, "_Grn.idat"), idat_has_valid_magic, logical(1)) |
               !vapply(paste0(ss$Basename, "_Red.idat"), idat_has_valid_magic, logical(1))
    if (any(corrupt)) {
      warning(sum(corrupt), " sample(s) have a corrupt IDAT file (bad magic header) and will be skipped:\n  ",
              paste(head(ss$Basename[corrupt], 10), collapse = "\n  "))
      ss <- ss[!corrupt, , drop = FALSE]
    }
  }

  if (is.null(chunk_size) || !is.finite(chunk_size) || chunk_size < 1 || chunk_size >= nrow(ss)) {
    return(read.metharray.exp(targets = ss, extended = extended,
                              force = force, verbose = TRUE))
  }

  # --- chunked read: read + cbind chunk_size samples at a time, so peak
  # memory is bounded to one chunk's raw idats instead of the whole cohort's ---
  groups <- split(seq_len(nrow(ss)), ceiling(seq_len(nrow(ss)) / chunk_size))
  message(sprintf("[read_cohort_idats] %d samples in %d chunk(s) of up to %d",
                  nrow(ss), length(groups), chunk_size))

  rg <- NULL
  for (i in seq_along(groups)) {
    chunk <- ss[groups[[i]], , drop = FALSE]
    message(sprintf("  chunk %d/%d (%d samples)", i, length(groups), nrow(chunk)))
    rg_chunk <- read.metharray.exp(targets = chunk, extended = extended,
                                   force = force, verbose = TRUE)
    rg <- if (is.null(rg)) rg_chunk else cbind(rg, rg_chunk)
    rm(rg_chunk); gc()
  }

  rg
}

# NOTE: bead-count QC (beadcount_from_rgset(), which turns the NBeads assay
# enabled by extended = TRUE above into per-sample fail_beadcount flags) now
# lives in R/sample_qc.R, alongside the other QC-flag logic (ewastools_qc(),
# sample_qc()).