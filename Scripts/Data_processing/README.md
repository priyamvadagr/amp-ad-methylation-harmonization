# Data_processing

Main pipeline: raw IDATs -> compiled RGChannelSet -> sample QC -> noob normalization -> probe filtering -> PCA QC. Shared helper functions live in [`../functions`](../functions); cohort paths are centralized in `functions/cohort_config.R`.

- **compile_idat_files.R** — Reads each cohort's IDAT manifest (built by `Scripts/Misc/get_sample_info.Rmd`) and assembles the raw `.idat` files into one minfi `RGChannelSet` per cohort, read in memory-bounded chunks. Saves `<combined_dir>/<cohort>.rds`.

- **sample_level_QC.R** — Runs sample-level QC on each cohort's compiled `RGChannelSet`, in two phases (`ewastools` on raw idats, `minfi` on the RGChannelSet) to avoid memory usage for large cohorts. Writes `<cohort>_ewastools.csv`, `<cohort>_sample_qc.csv`, the flagged `<cohort>_qc.rds`, and QC figures, all under `<qc_dir>/sample_qc/`. Nothing is dropped — every sample is flagged, not removed.

- **render_qc_plots.R** — Re-renders the sample-QC plots from the `<qc_dir>/sample_qc/<cohort>_qc.rds` + `_sample_qc.csv` already written by `sample_level_QC.R`, without recomputing QC. Requires `sample_level_QC.R` to have been run first.

- **normalize_data.R** — Runs noob normalization (`functions/noob_normalize.R`) on each cohort's QC'd `RGChannelSet`, dropping samples flagged in the QC step, masking high detection-p values, and writing `<cohort>_noob_gmset.rds` / `<cohort>_noob_beta.rds`.

- **filter_probes_per_cohort.R** — Applies Zhou-mask probe filtering (`functions/filter_probes.R`: low call rate -> cross-reactive -> SNP artefact -> non-autosomal) to each cohort's noob output, using the platform-appropriate (450K or EPIC) Zhou manifests. Writes `<cohort>_filtered_gmset.rds`, `_filtered_beta.rds`, `_probe_filter_summary.csv`.

- **pca_qc_per_cohort.R** — Runs PCA (`functions/pca_qc.R`) on each cohort's filtered beta matrix and tests PCs against technical (Sentrix ID, array row/col) and biological (sex, race, age, Braak, APOE4, tissue) covariates, to catch batch effects before modelling. Writes PCA scores, variance-explained, PC~covariate association tables, and figures under `<qc_dir>/pca_qc/`.

- **check-qc_flags.Rmd** — Diagnostic report investigating which ewastools control-probe metrics drive "control-only" sample QC flags in each cohort, to distinguish systematic control-probe artifacts from genuine sample failure. Rendered output: `check-qc_flags.pdf`.


