# ---------------------------------------------------------------------------
# Usage (per cohort, after normalize_noob):
#
#   source("filter_probes_zhou.R")
#   mask  <- load_zhou_mask("/…/annotation/HM450.hg19.manifest.pop.tsv.gz")  # or EPIC
#   gmset <- readRDS("/…/ROSMAP_noob_gmset.rds")
#   beta  <- readRDS("/…/ROSMAP_noob_beta.rds")                 # MASKED beta
#
#   # all steps at once:
#   res <- filter_probes_zhou(gmset, beta = beta, mask_tbl = mask,
#                             min_callrate = 0.99, autosomal_only = TRUE,
#                             out_prefix = "/…/ROSMAP/normalized/ROSMAP")
#   res$summary
#
#   # or step by step / in your own order (each returns $gmset and $beta):
#   s <- drop_low_callrate(gmset, beta, 0.99)
#   s <- drop_crossreactive(s$gmset, s$beta, mask)
#   s <- drop_snp(s$gmset, s$beta, mask, cols = c("MASK.snp5.EUR","MASK.snp5.AFR"))
#   s <- drop_nonautosomal(s$gmset, s$beta)
#   final_beta <- s$beta
#
# If the protocol also uses a population-agnostic SNP column, add it, e.g.:
#   drop_snp(..., cols = c("MASK.snp5.GMAF1p","MASK.snp5.EUR","MASK.snp5.AFR"))
# ---------------------------------------------------------------------------