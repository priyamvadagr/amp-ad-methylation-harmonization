# ============================================================================
# build_idat_manifest() : map each specimen to its green/red IDAT file.
#
# One row per specimenID/individualID with the matching Grn/Red IDAT file
# path (relative to that cohort's raw_data directory). grn_files/red_files
# should come from list.files(..., recursive = TRUE) so subfolder structure
# (e.g. ROSMAP's per-barcode folders) is preserved.
# ============================================================================

build_idat_manifest <- function(meth_df, sentrix_id_col, sentrix_rc_col, grn_files, red_files) {
  sentrix_id <- meth_df[[sentrix_id_col]]
  sentrix_rc <- meth_df[[sentrix_rc_col]]
  sample_name <- paste0(sentrix_id, '_', sentrix_rc)

  grn_lookup <- setNames(grn_files, sub('_Grn\\.idat$', '', basename(grn_files)))
  red_lookup <- setNames(red_files, sub('_Red\\.idat$', '', basename(red_files)))

  data.frame(
    specimenID         = meth_df$specimenID,
    individualID       = meth_df$individualID,
    Sentrix_ID         = sentrix_id,
    Sentrix_Row_Column = sentrix_rc,
    grnFile            = unname(grn_lookup[sample_name]),
    redFile            = unname(red_lookup[sample_name]),
    stringsAsFactors = FALSE
  )
}
