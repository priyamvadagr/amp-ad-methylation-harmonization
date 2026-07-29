
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
# ---------------------------------------------------------------------------
read_cohort_idats <- function(idat_dir,
                              sample_sheet = NULL,
                              sentrix_col  = "Sentrix_ID",
                              recursive    = TRUE,
                              force        = TRUE) {

  if (is.null(sample_sheet)) {
    # --- read every idat pair found under idat_dir ---
    rg <- read.metharray.exp(base = idat_dir,
                             recursive = recursive,
                             force = force,
                             verbose = TRUE)
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

    rg <- read.metharray.exp(targets = ss, force = force, verbose = TRUE)
  }

  rg
}