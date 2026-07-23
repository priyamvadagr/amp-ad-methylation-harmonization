# ============================================================================
# harmonize_adkp() : ONE generic harmonizer for ADKP individual metadata.
# You harmonize each study by handing it a MAPPING TABLE (one row per output
# column). Target = ADKP harmonized data dictionary (syn73713784).
#
# Quick start:
#   source("harmonize_adkp.R")
#   moa_pad_h <- harmonize_adkp(moa_pad_ind_df, mapping_moa_pad,
#                               keep_extra = keep_moa_pad)
#   validate_harmonized(moa_pad_h, "MOA-PAD")
#
# To onboard a new study, copy a mapping table below and edit the `source`
# column (and `map`/`const`) to match that study's raw columns/encodings.
# ============================================================================

suppressPackageStartupMessages({ library(dplyr); library(tibble) })
MISS <- "missing or unknown"

# ---------------------------------------------------------------------------
# 1. Canonical dictionary: output column -> allowed values (NULL = numeric/id)
# ---------------------------------------------------------------------------
ADKP_DICT <- list(
  individualID          = NULL,
  dataContributionGroup = c("Allen Institute","Banner Sun Health Research Institute",
                            "Columbia University","Emory University","Mayo Clinic",
                            "Mount Sinai School of Medicine","NIMH Human Brain Collection Core",
                            "Rush Alzheimer's Disease Center","Stanley Medical Research Institute"),
  cohort                = c("Banner Sun Health Research Institute","Biggs Institute Brain Bank",
                            "CLINCOR","Columbia ADRC","EFIGA","Emory","FBS","HBCC","LATC","MAP",
                            "MARS","Mayo Clinic Brain Bank","Mount Sinai Brain Bank","ROS",
                            "SEA-AD","SMRI","UFL","University of Kentucky","UPenn","WHICAP"),
  species               = c("Human"),
  sex                   = c("female","male", MISS),
  race                  = c("American Indian or Alaska Native","Asian",
                            "Black or African American",
                            "Native Hawaiian or Other Pacific Islander","White","Other", MISS),
  isHispanic            = c("True","False", MISS),
  ageDeath              = NULL,
  PMI                   = NULL,
  apoeGenotype          = c("22","23","24","33","34","44", MISS),
  apoe4Status           = c("yes","no", MISS),
  amyCerad              = c("None/No AD/C0","Sparse/Possible/C1","Moderate/Probable/C2",
                            "Frequent/Definite/C3", MISS),
  amyAny                = c("no","yes", MISS),
  amyThal               = c("None","Phase 1","Phase 2","Phase 3","Phase 4","Phase 5", MISS),
  amyA                  = c("None","Thal Phase 1 or 2","Thal Phase 3","Thal Phase 4 or 5", MISS),
  Braak                 = c("None","Stage I","Stage II","Stage III","Stage IV","Stage V",
                            "Stage VI", MISS),
  bScore                = c("None","Low (Stage I-II)","Moderate (Stage III-IV)",
                            "High (Stage V-VI)", MISS)
)
DICT_COLS <- names(ADKP_DICT)

# ---------------------------------------------------------------------------
# 2. Value-map registry: referenced by name from the `map` column of a table.
#    Left = raw source value (as character); right = dictionary value.
#    Covers ROSMAP numeric codes, MSBB letter codes, and string spellings.
# ---------------------------------------------------------------------------
MAPS <- list(
  sex = c(female="female", male="male", Female="female", Male="male",
          F="female", M="male", "0"="female", "1"="male"),          # ROSMAP msex 0/1
  race = c("White"="White","Asian"="Asian",
           "Black"="Black or African American","African American"="Black or African American",
           "Black or African American"="Black or African American",
           "American Indian or Alaska Native"="American Indian or Alaska Native",
           "Native Hawaiian or Other Pacific Islander"="Native Hawaiian or Other Pacific Islander",
           "Other"="Other",
           "1"="White","2"="Black or African American","3"="American Indian or Alaska Native",
           "4"="Native Hawaiian or Other Pacific Islander","5"="Asian","6"="Other","7"=MISS,
           "W"="White","B"="Black or African American","A"="Asian","H"="Other","O"="Other","U"=MISS),
  hispanic = c("Hispanic or Latino"="True","Hispanic"="True","Yes"="True","True"="True",
               "Not Hispanic or Latino"="False","Not Hispanic"="False","No"="False","False"="False",
               "1"="True","2"="False"),                              # ROSMAP spanish 1/2
  cerad = c("0"="None/No AD/C0","1"="Sparse/Possible/C1","2"="Moderate/Probable/C2",
            "3"="Frequent/Definite/C3",
            "None/No AD/C0"="None/No AD/C0","Sparse/Possible/C1"="Sparse/Possible/C1",
            "Moderate/Probable/C2"="Moderate/Probable/C2","Frequent/Definite/C3"="Frequent/Definite/C3"),
  braak = c("0"="None","1"="Stage I","2"="Stage II","3"="Stage III","4"="Stage IV",
            "5"="Stage V","6"="Stage VI",
            "I"="Stage I","II"="Stage II","III"="Stage III","IV"="Stage IV","V"="Stage V","VI"="Stage VI",
            "None"="None","Stage I"="Stage I","Stage II"="Stage II","Stage III"="Stage III",
            "Stage IV"="Stage IV","Stage V"="Stage V","Stage VI"="Stage VI"),
  thal = c("0"="None","1"="Phase 1","2"="Phase 2","3"="Phase 3","4"="Phase 4","5"="Phase 5",
           "None"="None","Phase 1"="Phase 1","Phase 2"="Phase 2","Phase 3"="Phase 3",
           "Phase 4"="Phase 4","Phase 5"="Phase 5")
)

# ---------------------------------------------------------------------------
# 3. Type handlers + derivers (deriver reads an already-harmonized column)
# ---------------------------------------------------------------------------
.as_miss <- function(x){ x<-as.character(x); x[is.na(x)|trimws(x)==""]<-MISS; x }
.map_vals <- function(x, lookup, col){
  xc<-trimws(as.character(x)); out<-unname(lookup[xc]); out[is.na(x)|xc==""]<-MISS
  bad<-is.na(out)&!(is.na(x)|xc==""); if(any(bad)){
    warning(sprintf("[%s] unmapped: %s -> '%s'. Add to MAPS.",col,
                    paste(unique(xc[bad]),collapse=", "),MISS)); out[bad]<-MISS}; out }
.censor_age <- function(x){ a<-suppressWarnings(as.numeric(as.character(x)))
  o<-as.character(a); o[!is.na(a)&a>=90]<-"90+"; o[is.na(a)]<-NA_character_; o }
.clean_apoe <- function(x){ g<-trimws(as.character(x)); o<-rep(MISS,length(g))
  ok<-g%in%c("22","23","24","33","34","44"); o[ok]<-g[ok]; o }

DERIVERS <- list(
  apoe4 = function(g){ g<-trimws(as.character(g)); o<-rep(MISS,length(g))
    ok<-g%in%c("22","23","24","33","34","44"); o[ok]<-ifelse(grepl("4",g[ok]),"yes","no"); o },
  amyAny = function(c){ o<-rep(MISS,length(c)); o[c=="None/No AD/C0"]<-"no"
    o[c%in%c("Sparse/Possible/C1","Moderate/Probable/C2","Frequent/Definite/C3")]<-"yes"; o },
  amyA = function(t){ o<-rep(MISS,length(t)); o[t=="None"]<-"None"
    o[t%in%c("Phase 1","Phase 2")]<-"Thal Phase 1 or 2"; o[t=="Phase 3"]<-"Thal Phase 3"
    o[t%in%c("Phase 4","Phase 5")]<-"Thal Phase 4 or 5"; o },
  bScore = function(b){ o<-rep(MISS,length(b)); o[b=="None"]<-"None"
    o[b%in%c("Stage I","Stage II")]<-"Low (Stage I-II)"
    o[b%in%c("Stage III","Stage IV")]<-"Moderate (Stage III-IV)"
    o[b%in%c("Stage V","Stage VI")]<-"High (Stage V-VI)"; o }
)

# ---------------------------------------------------------------------------
# 4. THE generic function
#
#   mapping: a data.frame / tibble with columns (extra columns ignored):
#     target   output (dictionary) column name              [required]
#     source   source column in df, OR harmonized column to read for derive
#     type     one of: id | string | map | numeric | age | apoe | const | derive
#     map      name of a MAPS entry            (used when type == "map")
#     const    fixed value                     (used when type == "const")
#     fn       name of a DERIVERS entry        (used when type == "derive")
#     default  fallback for blank/absent source (used when type == "string")
#   Any dictionary column not listed is auto-filled: MISS (string) or NA (numeric).
#
#   keep_extra: character vector of non-dictionary source columns to carry through.
# ---------------------------------------------------------------------------
harmonize_adkp <- function(df, mapping, keep_extra = character()) {
  mapping <- as.data.frame(mapping, stringsAsFactors = FALSE)
  for (c in c("source","type","map","const","fn","default"))
    if (!c %in% names(mapping)) mapping[[c]] <- NA_character_
  bad_t <- setdiff(mapping$target, DICT_COLS)
  if (length(bad_t)) warning("ignoring non-dictionary target(s): ",
                             paste(bad_t, collapse=", "))
  mapping <- mapping[mapping$target %in% DICT_COLS, , drop = FALSE]

  n <- nrow(df); out <- vector("list", length(DICT_COLS)); names(out) <- DICT_COLS
  # defaults for any unmapped column
  for (dc in DICT_COLS)
    out[[dc]] <- if (dc %in% c("ageDeath","PMI")) rep(NA_real_, n) else rep(MISS, n)

  get_src <- function(col) if (!is.na(col) && col %in% names(df)) df[[col]] else rep(NA, n)

  # pass 1: non-derive rows
  for (i in which(mapping$type != "derive" | is.na(mapping$type))) {
    r <- mapping[i, ]; dc <- r$target; src <- get_src(r$source)
    out[[dc]] <- switch(r$type,
      id      = as.character(src),
      const   = rep(r$const, n),
      string  = { v <- .as_miss(src)
                  if (!is.na(r$default)) v[is.na(src) | trimws(as.character(src))==""] <- r$default; v },
      map     = .map_vals(src, MAPS[[r$map]], dc),
      numeric = suppressWarnings(as.numeric(src)),
      age     = .censor_age(src),
      apoe    = .clean_apoe(src),
      .as_miss(src)
    )
  }
  # pass 2: derive rows (read already-harmonized `source` column)
  for (i in which(mapping$type == "derive")) {
    r <- mapping[i, ]; fn <- DERIVERS[[r$fn]]
    if (is.null(fn)) stop("Unknown deriver: ", r$fn)
    out[[r$target]] <- fn(out[[r$source]])
  }

  res <- as_tibble(out)[DICT_COLS]
  for (ec in keep_extra) res[[ec]] <- if (ec %in% names(df)) df[[ec]] else NA
  res
}

# ---------------------------------------------------------------------------
# 5. Mapping table (edit / copy per study)
# ---------------------------------------------------------------------------
# ---- Identity table: for files already in dictionary column names -----------
# (e.g. harmonized ROSMAP / MSBB). Re-validates + re-derives grouped columns.
mapping_identity <- tribble(
  ~target,                 ~source,                 ~type,     ~map,        ~const, ~fn,      ~default,
  "individualID",          "individualID",          "id",      NA,          NA,     NA,       NA,
  "dataContributionGroup", "dataContributionGroup", "string",  NA,          NA,     NA,       NA,
  "cohort",                "cohort",                "string",  NA,          NA,     NA,       NA,
  "species",               "species",               "string",  NA,          NA,     NA,       "Human",
  "sex",                   "sex",                   "map",     "sex",       NA,     NA,       NA,
  "race",                  "race",                  "map",     "race",      NA,     NA,       NA,
  "isHispanic",            "isHispanic",            "map",     "hispanic",  NA,     NA,       NA,
  "ageDeath",              "ageDeath",              "age",     NA,          NA,     NA,       NA,
  "PMI",                   "PMI",                   "numeric", NA,          NA,     NA,       NA,
  "apoeGenotype",          "apoeGenotype",          "apoe",    NA,          NA,     NA,       NA,
  "apoe4Status",           "apoeGenotype",          "derive",  NA,          NA,     "apoe4",  NA,
  "amyCerad",              "amyCerad",              "map",     "cerad",     NA,     NA,       NA,
  "amyAny",                "amyCerad",              "derive",  NA,          NA,     "amyAny", NA,
  "amyThal",               "amyThal",               "map",     "thal",      NA,     NA,       NA,
  "amyA",                  "amyThal",               "derive",  NA,          NA,     "amyA",   NA,
  "Braak",                 "Braak",                 "map",     "braak",     NA,     NA,       NA,
  "bScore",                "Braak",                 "derive",  NA,          NA,     "bScore", NA
)

# ---------------------------------------------------------------------------
# 6. Validate + combine
# ---------------------------------------------------------------------------
validate_harmonized <- function(df, study = "") {
  ok <- TRUE
  for (col in names(ADKP_DICT)) {
    allowed <- ADKP_DICT[[col]]
    if (is.null(allowed) || !col %in% names(df)) next
    bad <- setdiff(unique(as.character(df[[col]])), c(allowed, NA))
    if (length(bad)) { ok <- FALSE
      message(sprintf("[%s] INVALID %s: %s", study, col, paste(bad, collapse=", "))) }
  }
  if ("ageDeath" %in% names(df)) {
    a <- df$ageDeath; bad <- a[!is.na(a) & a!="90+" & is.na(suppressWarnings(as.numeric(a)))]
    if (length(bad)) { ok <- FALSE
      message(sprintf("[%s] INVALID ageDeath: %s", study, paste(unique(bad), collapse=", "))) }
  }
  if (ok) message(sprintf("[%s] all dictionary columns conform.", study))
  invisible(ok)
}

combine_studies <- function(...) {
  dfs <- list(...); all_cols <- Reduce(union, lapply(dfs, names))
  fill <- function(d){ for (c in setdiff(all_cols, names(d))) d[[c]] <- NA; d[all_cols] }
  bind_rows(lapply(dfs, fill))
}

# ---------------------------------------------------------------------------
# 7. Example
# ---------------------------------------------------------------------------
# moa_pad_h <- harmonize_adkp(moa_pad_ind_df, mapping_moa_pad, keep_extra = keep_moa_pad)
# rosmap_h  <- harmonize_adkp(rosmap_ind_df, mapping_identity,
#                             keep_extra = setdiff(names(rosmap_ind_df), DICT_COLS))
# msbb_h    <- harmonize_adkp(msbb_ind_df,   mapping_identity,
#                             keep_extra = setdiff(names(msbb_ind_df), DICT_COLS))
# validate_harmonized(moa_pad_h,"MOA-PAD"); validate_harmonized(rosmap_h,"ROSMAP")
# validate_harmonized(msbb_h,"MSBB")
# combined <- combine_studies(rosmap_h, msbb_h, moa_pad_h)